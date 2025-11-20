# Copyright (C) 2025  Technische Universitaet Berlin
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA

import traceback
from flask import jsonify
from xisearch_common.config import Crosslinker, Modification, ModificationConfig, Loss, \
    FragmentationConfig, Config
from xisearch_common.mock_context import MockContext
from xisearch_common.spectra_reader import Spectrum
from xisearch_common.fragment_peptides import fragment_crosslinked_peptide_pair, \
    fragment_linear_peptide, fragment_noncovalent_peptide_pair
from xisearch_common.fragmentation import spread_charges, include_losses
from xisearch_common.filters import IsotopeDetector
from xisearch_common import const
import numpy as np
import re
import os


def annotate_request(json_request):
    """
    Annotate the json request.

    :param json_request: JSON annotation request
    :return: JSON annotation response
    """
    try:
        # create Config object
        if 'config' in json_request['annotation'].keys():
            # add default losses if there were no losses defined
            # ToDo: maybe this should happen on the frontend request generation?
            if 'losses' not in json_request['annotation']['config']['fragmentation'].keys():
                json_request['annotation']['config']['fragmentation']['losses'] = [
                    {"name": 'H2O',
                     "mass": 18.01056027,
                     "specificity": ['S', 'T', 'D', 'E', 'cterm']},
                    {"name": 'NH3',
                     "mass": 17.02654493,
                     "specificity": ['R', 'K', 'N', 'Q', 'nterm']}
                ]
            config = Config(**json_request['annotation']['config'])
        else:
            # create from xi1 json style format
            config = create_config_from_json_format(json_request['annotation'])

        # set return mod syntax
        return_mod_syntax = json_request['annotation'].get('returnModSyntax',
                                                           config.mod_peptide_syntax)
        # set crosslinker
        is_crosslinked = False
        crosslinker = None
        if len(config.crosslinker) > 0:
            is_crosslinked = True
            try:
                crosslinker_idx = json_request['annotation']['crosslinkerID']
            except KeyError:
                if len(config.crosslinker) == 1:
                    crosslinker_idx = 0
                else:
                    raise ValueError(
                        "More than 1 crosslinker in config without defined crosslinkerID!")
            crosslinker = config.crosslinker[crosslinker_idx]

        # create peptide database and set up Context
        ctx = MockContext(config)
        # xi2 style annotation
        if 'base_sequence' in json_request['Peptides'][0].keys():
            base_seqs = [p['base_sequence'].encode('ascii') for p in json_request['Peptides']]
            mod_ids = [p['modification_ids'] for p in json_request['Peptides']]
            mod_pos = [p['modification_positions'] for p in json_request['Peptides']]

            ctx.setup_peptide_db_xi2annotator(base_seqs, mod_ids, mod_pos)

            # get the peptide indices of input peptides in the peptide db
            peps_db = ctx.peptide_db.unmodified_sequences[ctx.peptide_db.peptides['sequence_index']]
            pep_idx = [np.where(peps_db == p)[0][0] for p in base_seqs]
        # xi1 style annotation
        else:
            peptides = []
            for peptide in json_request['Peptides']:
                pep_seq = ''.join([aa['Modification'] + aa['aminoAcid']
                                   for aa in peptide['sequence']])
                peptides.append(pep_seq.encode('ascii'))
            peptides = np.array(peptides)
            ctx.setup_peptide_db(peptides)

            # get the peptide indices of input peptides in the peptide db
            peps_db = ctx.peptide_db.mod_pep_sequence(range(ctx.peptide_db.peptides.size))
            pep_idx = [np.where(peps_db == p)[0][0] for p in peptides]

        # create Spectrum object
        precursor = {
            'mz': json_request['annotation'].get('precursorMZ', None),
            'charge': json_request['annotation']['precursorCharge'],
            'intensity': json_request['annotation'].get('precursorIntensity', -1)
        }
        mz_array = np.array([p['mz'] for p in json_request['peaks']])
        int_array = np.array([p['intensity'] for p in json_request['peaks']])
        spectrum = Spectrum(precursor, mz_array, int_array, -1)

        # Process Spectrum for annotation:
        # detect and reduce isotope clusters to monoisotopic peaks
        detector = IsotopeDetector(ctx)
        full_match_spectrum = detector.process(spectrum)

        # create unique clusters, cluster indices and cluster ids
        unique_clusters, cluster_indices, cluster_ids = np.unique(
            full_match_spectrum.isotope_cluster_peaks['cluster_id'], return_index=True,
            return_inverse=True)

        # create response peaks block with clusterIds mz-ordered
        json_request['peaks'] = [
            {'mz': float(m), 'intensity': float(i), 'clusterIds':
                [int(cid) for cid in full_match_spectrum.isotope_cluster_peaks['cluster_id'][
                    np.where(full_match_spectrum.isotope_cluster_peaks['peak_id'] == peakid)]]}
            for peakid, (m, i) in enumerate(zip(
                full_match_spectrum.mz_values, full_match_spectrum.int_values))
        ]

        # create clusters
        json_request['clusters'] = [
            {'charge': int(c), 'firstPeakId':
                int(full_match_spectrum.isotope_cluster_peaks['peak_id'][cluster_indices[i]])}
            for i, c in enumerate(full_match_spectrum.isotope_cluster_charge_values)
        ]

        # create fragments
        n_peptides = len(ctx.peptide_db.peptides)
        # Linear
        if n_peptides == 1:
            fragments = fragment_linear_peptide(
                0, ctx, add_precursor=config.fragmentation.add_precursor)
            fragments = include_losses(fragments, [0], ctx)
            fragments = spread_charges(fragments, ctx, precursor['charge'])
            # overwrite LinkSite with empty list for linears
            json_request['LinkSite'] = []
        elif n_peptides == 2:
            link1_pos = json_request['LinkSite'][0]['linkSite']
            link2_pos = json_request['LinkSite'][1]['linkSite']
            # noncovalently associated peptides NAPs
            if link1_pos == -1 and link2_pos == -1:
                fragments = fragment_noncovalent_peptide_pair(
                    pep_idx[0], pep_idx[1], ctx,
                    add_precursor=config.fragmentation.add_precursor)
            # Crosslinked peptide
            elif is_crosslinked:
                fragments = fragment_crosslinked_peptide_pair(
                    pep_idx[0], pep_idx[1], link1_pos, link2_pos,
                    crosslinker, ctx,
                    add_precursor=config.fragmentation.add_precursor)
            else:
                raise ValueError("2 peptides with crosslink positions defined but no crosslinker.")
            fragments = include_losses(fragments, [pep_idx[0], pep_idx[1]], ctx)
            fragments = spread_charges(fragments, ctx, precursor['charge'])
        else:
            raise ValueError("Unsupported number of peptides given!")

        # annotate the spectrum with fragments
        annotations = full_match_spectrum.annotate_spectrum(fragments, ctx)

        if len(annotations) == 0:
            json_request['fragments'] = []
        else:
            # generate peptides in return_mod_syntax as list of modified amino acids
            # for assembly of fragment sequences
            if return_mod_syntax == 'modX':
                pep_to_aa_re = const.PEPTIDE_TO_AMINO_ACID
            elif return_mod_syntax == 'Xmod':
                pep_to_aa_re = re.compile(b'([A-Z][^A-Z\\-]*)')

            return_pep1 = ctx.peptide_db.mod_pep_sequence([pep_idx[0]],
                                                          mod_peptide_syntax=return_mod_syntax)
            pep_mod_arr1 = pep_to_aa_re.findall(return_pep1)
            return_pep2 = ctx.peptide_db.mod_pep_sequence([pep_idx[1]],
                                                          mod_peptide_syntax=return_mod_syntax) \
                if n_peptides == 2 else None
            pep_mod_arr2 = pep_to_aa_re.findall(return_pep2) if n_peptides == 2 else []
            pep_mod_arr = [pep_mod_arr1, pep_mod_arr2]

            # create unique fragments
            fragment_cols = ['ion_type', 'idx', 'pep_id', 'LN', 'loss', 'nlosses', 'stub',
                             'ranges']
            annotations.sort(order=fragment_cols)
            if len(annotations) == 1:
                fragments = annotations[fragment_cols]
                frag_indices = [0]
                frag_counts = [1]
            else:
                fragments, frag_indices, frag_counts = np.unique(
                    annotations[fragment_cols], return_index=True, return_counts=True)

            json_request['fragments'] = []
            for i, f in enumerate(fragments):
                name = f['ion_type'].decode('ascii')
                if f['ion_type'] != b'P':
                    name += str(f['idx'])
                if not f['LN']:
                    name += '+P'
                name += f['stub'].decode('ascii')
                if f['nlosses'] > 0:
                    name += '_' + f['loss'].decode('ascii')
                # get all annotations for this fragment
                f_annotations = annotations[frag_indices[i]: frag_indices[i] + frag_counts[i]]

                # loop over fragment annotations for clusterIds and clusterInfo
                cluster_ids = []
                cluster_info = []
                for annotation in f_annotations:
                    # get cluster id and convert to int (int16 not JSON serializable)
                    cluster_id = int(annotation['cluster_id'])
                    # append to cluster ids
                    cluster_ids.append(cluster_id)
                    cluster_info.append({
                        'Clusterid': cluster_id,
                        'calcMZ': annotation['frag_mz'],
                        'error': annotation['rel_error'] / 1e-6,
                        'errorUnit': 'ppm',
                        'matchedMissingMonoIsotopic': int(annotation['missing_monoisotopic_peak']),
                        'matchedCharge': int(annotation['frag_charge'])
                    })

                # reformat ranges according to annotator format
                f_pep_id = f['pep_id'] - 1  # change pep_id from 1-based to 0-based

                if f['LN']:
                    ranges = [{
                        'peptideId': int(f_pep_id),
                        'from': int(f['ranges'][f_pep_id][0]),
                        'to': int(f['ranges'][f_pep_id][1]) - 1
                    }]
                else:
                    ranges = [{'peptideId': pep_id, 'from': int(r[0]), 'to': int(r[1]) - 1}
                              for pep_id, r in enumerate(f['ranges'])]

                # assemble fragment sequence
                frag_sequence = [pep_mod_arr[i][r[0]:r[1]]
                                 for i, r in enumerate(f['ranges']) if r.sum() > 0]
                frag_sequence_str = ' + '.join(
                    [''.join([aa.decode() for aa in seq]) for seq in frag_sequence])
                fragment = {
                    'name': name,
                    'ionNumber': int(f['idx']),
                    'peptideId': int(f_pep_id),
                    'range': ranges,
                    'type': f['ion_type'].decode('ascii'),
                    # ToDo: change class to primary True, False? or just have field nlosses
                    'class': 'lossy' if f['nlosses'] > 0 else 'non-lossy',
                    'nlosses': int(f['nlosses']),
                    'stub': f['stub'].decode('ascii'),
                    'clusterIds': cluster_ids,
                    'clusterInfo': cluster_info,
                    'sequence': frag_sequence_str
                }
                json_request['fragments'].append(fragment)

        # theoretical calculated precursor m/z and error
        peptides_mass = ctx.peptide_db.peptide_mass(np.arange(n_peptides)).sum()
        if is_crosslinked:
            cl_mass = crosslinker.mass
            peptides_mass += cl_mass
        calc_mz = (peptides_mass / precursor['charge']) + const.PROTON_MASS
        json_request['annotation']['calculatedMZ'] = calc_mz
        if precursor['mz'] is None:
            json_request['annotation']['precursorError'] = ''
        else:
            json_request['annotation']['precursorError'] = {
                'tolerance': (precursor['mz'] - calc_mz) / calc_mz * 1e6,
                'unit': 'ppm'
            }

        # write out crosslinker modMass for xispec
        if is_crosslinked:
            try:
                xi2_stubs = config.crosslinker[crosslinker_idx].cleavage_stubs
                stubs = []
                for s in xi2_stubs:
                    stubs.append(f"{s.name}:{s.mass}:{''.join(s.pairs_with)}")
            except (KeyError, TypeError):
                stubs = []

            json_request['annotation']['crosslinker'] = {
                'name': crosslinker.name,
                'modMass': crosslinker.mass,
                'specificity': crosslinker.specificity,
            }
            if len(stubs) > 0:
                json_request['annotation']['crosslinker'].update({
                    'stubs1': stubs,
                    'stubs2': stubs,
                    'cleavage_stubs': [s.to_dict() for s in xi2_stubs]
                })

        # write out modifications with masses (it's possible to just give composition and not mass)
        json_request['annotation']['modifications'] = []
        for mod in config.modification.modifications:
            json_request['annotation']['modifications'].append({
                'id': mod.name,
                'aminoAcids': mod.specificity,
                'mass': mod.mass
            })

        # write out losses in annotation block for backwards compatibility
        json_request['annotation']['losses'] = []
        for loss in config.fragmentation.losses:
            specificity = loss.specificity.copy()
            if return_mod_syntax != config.mod_peptide_syntax:
                if return_mod_syntax == 'modX':
                    specificity = [s[1:] + s[0] for s in specificity]
                elif return_mod_syntax == 'Xmod':
                    specificity = [s[-1] + s[:-1] for s in specificity]
                else:
                    raise Exception
            old_loss = {
                "id": loss.name,
                "specificity": specificity,
                "mass": loss.mass
            }
            if loss.cterm:
                old_loss['specificity'].append('CTerm')
            if loss.nterm:
                old_loss['specificity'].append('NTerm')
            json_request['annotation']['losses'].append(old_loss)

        # ToDo: the version should come from a central place
        json_request['annotation']['xiVersion'] = const.VERSION

        return jsonify(json_request)
    except Exception as e:
        debug_value = os.environ.get("XI2ANNOTATOR_DEBUG", "false")
        if debug_value.lower() != "false" and debug_value != "0":
            return jsonify({'error': str(e), "stacktrace": traceback.format_exc()}), 400
        raise e


def create_config_from_json_format(annotation_json):
    """
    Create a Config from the xi1 json format.

    :param annotation_json: annotation block of the json request
    :return: xi2 config
    :rtype: Config
    """
    # Crosslinker
    if 'crosslinker' in annotation_json.keys():
        a_stubs1 = annotation_json['crosslinker'].get('stubs1', [])
        a_stubs2 = annotation_json['crosslinker'].get('stubs2', [])
        # merge the stubs removing duplicates
        stubs = a_stubs1.copy()
        for s2 in a_stubs2:
            if s2 not in a_stubs1:
                stubs.append(s2)
        # convert stubs from "[name]:[mass]:[pairs with]" to the xi2 cleavage stub format"
        xi2_stubs = []
        for s in stubs:
            parts = s.split(':')
            if len(parts) != 3:
                raise ValueError(f"Crosslinker stub {s} not in expected format"
                                 " [name]:[mass]:[pairs with]")
            xi2_stubs.append({
                'name': parts[0],
                'mass': float(parts[1]),
                'pairs_with': list(parts[2])
            })
        if (len(xi2_stubs) > 0):
            annotation_json['crosslinker']['cleavage_stubs'] = xi2_stubs
        # create crosslinker with specificity X (any amino acid)
        crosslinker = [Crosslinker(mass=annotation_json['crosslinker']['modMass'],
                                   specificity=['X'], name='MockCL', cleavage_stubs=xi2_stubs)]
    else:
        crosslinker = []

    # Modifications
    modifications = []
    for mod in annotation_json['modifications']:
        modifications.append(
            Modification(name=mod['id'], mass=mod['mass'], specificity=mod['aminoAcids'],
                         type='variable'))
    mod_config = ModificationConfig(modifications=modifications)

    # MS2 fragment tolerance
    frag_tol = annotation_json['fragmentTolerance']
    ms2_tol = f"{frag_tol['tolerance']} {frag_tol['unit']}"

    # Fragmentation config
    ion_types = [i['type'] for i in annotation_json['ions']]
    nterm_ions = [i[0].lower() for i in ion_types if i in ('AIon', 'BIon', 'CIon')]
    cterm_ions = [i[0].lower() for i in ion_types if i in ('XIon', 'YIon', 'ZIon')]
    add_precursor = 'PeptideIon' in ion_types

    # set standard losses if no losses defined
    if 'losses' not in annotation_json.keys():
        losses = [Loss.H2O, Loss.NH3]
    else:
        losses = []
        for loss in annotation_json['losses']:
            # swap AA to end of string and lowercase n/cterm
            specificity = [s[1:] + s[0] if s not in ('CTerm', 'NTerm') else
                           s.lower() for s in loss['specificity']]
            losses.append(Loss(name=loss['id'], specificity=specificity, mass=loss['mass']))

    frag_config = FragmentationConfig(nterm_ions=nterm_ions, cterm_ions=cterm_ions,
                                      add_precursor=add_precursor, losses=losses)
    # create Config
    config = Config(crosslinker=crosslinker, ms2_tol=ms2_tol,
                    fragmentation=frag_config, modification=mod_config)

    return config
