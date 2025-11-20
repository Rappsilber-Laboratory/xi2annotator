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

from xi2annotator import create_app
import pytest
from flask import url_for
import os
import json
import copy
from recursive_diff import recursive_eq


@pytest.fixture
def app():
    app = create_app()
    return app


def recursive_ppm_error_adjust(d):
    """
    Go through the json dictionary and "transform" it back into a Da error.

    The reason for that is, that ppm error should be compared with a different tolerance then
    e.g. m/z values. The transformation is done by simply dividing the error by 100000. This is
    based on ppm error = (mz1-mz2)/mz2*1000000. As the mz1 and 2 are more complicated to get I
    just used 100000 instead of 1000000.

    :param d: dictionary representation of the JSON
    """
    if isinstance(d, dict):
        for key in d:
            if isinstance(d[key], dict):
                recursive_ppm_error_adjust(d[key])
            elif isinstance(d[key], set) or isinstance(d[key], list):
                for v in d[key]:
                    recursive_ppm_error_adjust(v)
            elif key == "error":
                if "errorUnit" in d and d["errorUnit"] == "ppm":
                    d[key] /= 100000

    elif isinstance(d, set) or isinstance(d, list):
        for v in d:
            recursive_ppm_error_adjust(v)

    return d


def test_general(client, config):
    url = url_for('xi2annotator.annotate')
    # GET not allowed
    assert client.get(url)._status_code == 405
    # check cors headers
    assert config['CORS_HEADERS'] == 'Content-Type'


def check_result(res_json, exp):
    # check peak clusterIds
    if 'peak_cluster_ids' in exp.keys():
        assert [p['clusterIds'] for p in res_json['peaks']] == exp['peak_cluster_ids']

    # check calculated precursor and error
    if 'precursor' in exp.keys():
        assert res_json['annotation']['precursorError']['tolerance'] == exp['precursor']['error']
        assert res_json['annotation']['precursorError']['unit'] == exp['precursor']['errorUnit']
        assert res_json['annotation']['calculatedMZ'] == exp['precursor']['calcMZ']

    # check losses
    if 'losses' in exp.keys():
        assert res_json['annotation']['losses'] == exp['losses']

    # check clusters
    if 'cluster' in exp.keys():
        assert [p['charge'] for p in res_json['clusters']] == exp['cluster']['charges']
        assert [p['firstPeakId'] for p in res_json['clusters']] == exp['cluster']['firstPeakIds']

    # check fragments
    if 'fragments' in exp.keys():
        recursive_eq(recursive_ppm_error_adjust(copy.deepcopy(exp['fragments'])),
                     recursive_ppm_error_adjust(copy.deepcopy(res_json['fragments'])),
                     abs_tol=1e-8, rel_tol=1e-6)


def check_unchanged_values(res, req, req_format):
    # check that Peptides and LinkSite are unchanged
    assert res['Peptides'] == req['Peptides']
    if 'LinkSite' in req.keys():
        assert res['LinkSite'] == req['LinkSite']
    else:
        assert res['LinkSite'] == []
    # check unchanged peaks
    assert [p['intensity'] for p in res['peaks']] == [p['intensity'] for p in req['peaks']]
    assert [p['mz'] for p in res['peaks']] == [p['mz'] for p in req['peaks']]

    if req_format == 'xi1':
        # check unchanged annotation block values
        assert res['annotation']['crosslinker']['modMass'] == \
               req['annotation']['crosslinker']['modMass']
        assert res['annotation']['fragmentTolerance'] == req['annotation']['fragmentTolerance']
        # check if same ions are defined (independent of order)
        for i in res['annotation']['ions']:
            assert i in req['annotation']['ions']
        assert res['annotation']['precursorCharge'] == req['annotation']['precursorCharge']


# the default losses in backwards compatible json format
default_losses = [
    {
        "id": "H2O",
        "specificity": ["S", "T", "D", "E", "CTerm"],
        "mass": 18.01056027
    },
    {
        "id": "NH3",
        "specificity": ["R", "K", "N", "Q", "NTerm"],
        "mass": 17.02654493
    }
]

# expected values for AKT-KMR_1-0:3 (BS3)
exp_simple_synthetic = {
    # synthetic spectrum with 3 isotope peaks -> 4 peaks per cluster
    'peak_cluster_ids': [
        [0, 1], [1], [1], [0, 1],
        [2, 3], [3], [2, 3], [3],
        [4, 5], [5], [5], [4, 5],
        [6, 7], [7], [7], [6, 7],
        [8, 9], [9], [8, 9], [9],
        [10], [10], [10], [10],
        [11, 12], [12], [11, 12], [12],
        [13, 14], [14], [14], [13, 14],
        [15], [15], [15], [15],
        [16, 17], [17], [16, 17], [17],
        [18], [18], [18], [18],
        [19, 20], [20], [20], [19, 20],
        [21, 22], [22], [22], [21, 22],
        [23, 24], [24], [24], [23, 24],
        [25, 26], [26], [26], [25, 26],
        [27, 28], [28], [27, 28], [28],
        [29], [29], [29], [29],
        [30, 31], [31], [30, 31], [31],
        [32, 33], [33], [32, 33], [33],
        [34, 35], [35], [34, 35], [35],
        [36], [36], [36], [36],
        [37], [37], [37], [37],
        [38], [38], [38], [38],
        [39], [39], [39], [39]
    ],
    'precursor': {
        'calcMZ': 297.50911753191565,
        'error': 0.0015733445603505064,
        'errorUnit': 'ppm',
    },
    'cluster': {
        'charges': [1, 3, 1, 2, 1, 3, 1, 3, 1, 2, 1, 1, 2, 1, 3, 1, 1, 2, 1, 1, 3, 1, 3, 1,
                    3, 1, 3, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1],
        'firstPeakIds': [0, 0, 4, 4, 8, 8, 12, 12, 16, 16, 20, 24, 24, 28, 28, 32, 36, 36, 40,
                         44, 44, 48, 48, 52, 52, 56, 56, 60, 60, 64, 68, 68, 72, 72, 76, 76, 80,
                         84, 88, 92]
    },
    # no losses defined in request so we expect the default losses
    'losses': default_losses,
    # manually verified matching fragments on spectrumviewer.org
    'fragments': [
        {
            "class": "non-lossy",
            "clusterIds": [
                1,
                3,
                10
            ],
            "clusterInfo": [
                {
                    "Clusterid": 1,
                    "calcMZ": 24.68631439511567,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 3,
                    "calcMZ": 36.52583335923401,
                    "error": -1.9453156038135122e-10,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 10,
                    "calcMZ": 72.04439025158901,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 0
                }
            ],
            "sequence": "A",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                20,
                28,
                36
            ],
            "clusterInfo": [
                {
                    "Clusterid": 20,
                    "calcMZ": 195.79173078848567,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 28,
                    "calcMZ": 293.183957949289,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 36,
                    "calcMZ": 585.3606394316989,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b1+P",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 2
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 0
                }
            ],
            "sequence": "AKT + K",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                24,
                33,
                38
            ],
            "clusterInfo": [
                {
                    "Clusterid": 24,
                    "calcMZ": 257.82303661121233,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 33,
                    "calcMZ": 386.23091668337895,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 38,
                    "calcMZ": 771.4545568998789,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b2+P",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 1
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 2
                }
            ],
            "sequence": "AK + KMR",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                22,
                31,
                37
            ],
            "clusterInfo": [
                {
                    "Clusterid": 22,
                    "calcMZ": 239.471892426149,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 31,
                    "calcMZ": 358.7042004057839,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 37,
                    "calcMZ": 716.4011243446889,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b2+P",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 2
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 1
                }
            ],
            "sequence": "AKT + KM",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                5,
                9,
                15
            ],
            "clusterInfo": [
                {
                    "Clusterid": 5,
                    "calcMZ": 40.693357517582335,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 9,
                    "calcMZ": 60.536398042934,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 15,
                    "calcMZ": 120.065519618989,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 2,
                    "peptideId": 0,
                    "to": 2
                }
            ],
            "sequence": "T",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                7,
                12,
                18
            ],
            "clusterInfo": [
                {
                    "Clusterid": 7,
                    "calcMZ": 59.044501702645675,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 12,
                    "calcMZ": 88.06311432052901,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 18,
                    "calcMZ": 175.11895217417901,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 2,
                    "peptideId": 1,
                    "to": 2
                }
            ],
            "sequence": "R",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                26,
                35,
                39
            ],
            "clusterInfo": [
                {
                    "Clusterid": 26,
                    "calcMZ": 273.830079733679,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 35,
                    "calcMZ": 410.24148136707896,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 39,
                    "calcMZ": 819.4756862672789,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y2+P",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 1,
                    "peptideId": 0,
                    "to": 2
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 2
                }
            ],
            "sequence": "KT + KMR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                14,
                17,
                29
            ],
            "clusterInfo": [
                {
                    "Clusterid": 14,
                    "calcMZ": 102.724663340309,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 17,
                    "calcMZ": 153.58335677702397,
                    "error": 1.8505722251967272e-10,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 29,
                    "calcMZ": 306.15943708716895,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 1,
                    "peptideId": 1,
                    "to": 2
                }
            ],
            "sequence": "MR",
            "stub": "",
            "type": "y"
        }
    ]
}


def test_annotate_simple_synthetic_xi1_format(client):
    """
    Test the annotation of the synthetic spectrum of AKT-KMR_1-0:3 (BS3) using the xi1 json style.

    no modifications defined
    no losses in the request -> use standard losses
    """
    url = url_for('xi2annotator.annotate')

    # Test annotation with old config format
    current_dir = os.path.dirname(__file__)
    # Load request with synthetic spectrum for AKT-KMR_1-0:3 BS3 with 3 extra isotope peaks
    json_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                             'xi1_format_AKT-KMR_1-0_z3_BS3.json')
    with open(json_file) as f:
        request = json.load(f)

    # request annotation
    res = client.post(url, json=request)
    # status ok
    assert res._status_code == 200
    assert res.mimetype == 'application/json'

    check_unchanged_values(res.json, request, 'xi1')
    check_result(res.json, exp_simple_synthetic)
    # check that the modifications have been written out
    exp_mods = []  # no modifications in this case
    assert res.json['annotation']['modifications'] == exp_mods


def sort_json_recursively(obj):
    """
    Recursively sort JSON object by keys, and sort lists of dictionaries by their values.

    Args:
        obj: The JSON object to sort (dict, list, or primitive)

    Returns:
        Sorted version of the object
    """
    if isinstance(obj, dict):
        # Sort dictionary by keys
        sorted_dict = {}
        for key in sorted(obj.keys()):
            sorted_dict[key] = sort_json_recursively(obj[key])
        return sorted_dict
    elif isinstance(obj, list):
        # Sort each item in the list recursively
        sorted_list = [sort_json_recursively(item) for item in obj]

        # If list contains dictionaries, sort by their string representation
        # This ensures consistent ordering of complex objects
        if sorted_list and isinstance(sorted_list[0], dict):
            try:
                sorted_list.sort(key=lambda x: json.dumps(x, sort_keys=True))
            except TypeError:
                # If sorting fails, keep original order
                pass
        elif sorted_list and not isinstance(sorted_list[0], (dict, list)):
            # Sort primitive values
            try:
                sorted_list.sort()
            except TypeError:
                # If sorting fails (mixed types), keep original order
                pass

        return sorted_list
    else:
        # Return primitive values as-is
        return obj


def test_annotate_simple_synthetic_xi1_stubs(client):
    """
    Test the annotation of the synthetic spectrum of AKT-KMR_1-0:3 (BS3) using the xi1 json style.

    no modifications defined
    no losses in the request -> use standard losses
    """
    url = url_for('xi2annotator.annotate')

    # Test annotation with old config format
    current_dir = os.path.dirname(__file__)
    # Load request with synthetic spectrum for AKT-KMR_1-0:3 BS3 with 3 extra isotope peaks
    json_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                             'x1_format_KKK-KKR_1-1_z1_BS3_stubs.json')
    with open(json_file) as f:
        request = json.load(f)

    # request annotation
    res = client.post(url, json=request)
    # read expected response
    json_exp_resp_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                                      'x1_format_KKK-KKR_1-1_z1_BS3_stubs_expected_response.json')

    with open(json_exp_resp_file) as f:
        exp_res = json.load(f)

    exp_res = sort_json_recursively(exp_res)
    res_json = sort_json_recursively(res.json)
    # status ok
    assert res._status_code == 200
    assert res.mimetype == 'application/json'

    check_unchanged_values(res_json, request, 'xi1')

    assert set(res_json['annotation']['crosslinker']['stubs1']) == set(["b:138.06808:o", "o:0.0:b"])
    assert set(res_json['annotation']['crosslinker']['stubs2']) == set(["b:138.06808:o", "o:0.0:b"])
    check_result(res_json, exp_res)


def test_annotate_simple_synthetic_xi2_format(client):
    """
    Test the annotation of the synthetic spectrum of AKT-KMR_1-0:3 (BS3) using a xi2 json config.

    no modifications defined
    no losses in the request -> use standard losses
    """
    url = url_for('xi2annotator.annotate')

    # Test annotation with old config format
    current_dir = os.path.dirname(__file__)
    # Load request with synthetic spectrum for AKT-KMR_1-0:3 BS3 with 3 extra isotope peaks
    json_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                             'xi2_format_AKT-KMR_1-0_z3_BS3.json')
    with open(json_file) as f:
        request = json.load(f)

    # request annotation
    res = client.post(url, json=request)
    # status ok
    assert res._status_code == 200
    assert res.mimetype == 'application/json'

    check_unchanged_values(res.json, request, 'xi2')
    check_result(res.json, exp_simple_synthetic)
    # check that the modifications have been written out
    exp_mods = []  # no modifications in this case
    assert res.json['annotation']['modifications'] == exp_mods


exp_loss_and_mod = {
    # we don't check this for 'peak_cluster_ids', 'precursor' and 'cluster' values atm
    # 'peak_cluster_ids': [],
    # 'precursor': {
    #     'calcMZ': ,
    #     'error': ,
    #     'errorUnit': 'ppm',
    # },
    # 'cluster': {
    #     'charges': ,
    #     'firstPeakIds': ,
    # },
    # loss as defined in request
    'losses': [
        {
            "id": "H2O",
            "specificity": ["S", "T", "D", "E", "CTerm"],
            "mass": 18.01056027
        }
    ],
    # manually verified matching fragments on spectrumviewer.org
    'fragments': [
        {
            "class": "non-lossy",
            "clusterIds": [
                30
            ],
            "clusterInfo": [
                {
                    "Clusterid": 30,
                    "calcMZ": 243.108781413299,
                    "error": -2.185907460485464,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 1
                }
            ],
            "sequence": "QN",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                51
            ],
            "clusterInfo": [
                {
                    "Clusterid": 51,
                    "calcMZ": 403.13946619800896,
                    "error": -4.3811091621034794,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b3",
            "ionNumber": 3,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 2
                }
            ],
            "sequence": "QNCcm",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                66
            ],
            "clusterInfo": [
                {
                    "Clusterid": 66,
                    "calcMZ": 532.1820592859791,
                    "error": 2.87629767701505,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b4",
            "ionNumber": 4,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 3
                }
            ],
            "sequence": "QNCcmE",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                71
            ],
            "clusterInfo": [
                {
                    "Clusterid": 71,
                    "calcMZ": 645.2661232631091,
                    "error": -2.6706238665404034,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b5",
            "ionNumber": 5,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 4
                }
            ],
            "sequence": "QNCcmEL",
            "stub": "",
            "type": "b"
        },
        {
            "class": "lossy",
            "clusterIds": [
                70
            ],
            "clusterInfo": [
                {
                    "Clusterid": 70,
                    "calcMZ": 627.255562993109,
                    "error": 0.28219262028307385,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b5_H2Ox1",
            "ionNumber": 5,
            "nlosses": 1,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 4
                }
            ],
            "sequence": "QNCcmEL",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                86
            ],
            "clusterInfo": [
                {
                    "Clusterid": 86,
                    "calcMZ": 792.3345371760992,
                    "error": -1.624536150762437,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b6",
            "ionNumber": 6,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 5
                }
            ],
            "sequence": "QNCcmELF",
            "stub": "",
            "type": "b"
        },
        {
            "class": "lossy",
            "clusterIds": [
                144
            ],
            "clusterInfo": [
                {
                    "Clusterid": 144,
                    "calcMZ": 1680.3764885898188,
                    "error": -17.608193772598554,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "b6+P_H2Ox1",
            "ionNumber": 6,
            "nlosses": 1,
            "peptideId": 1,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 5
                }
            ],
            "sequence": "QNCcmELFEQLGEYKFQNALLVR + KQTALV",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                57
            ],
            "clusterInfo": [
                {
                    "Clusterid": 57,
                    "calcMZ": 461.19220336547403,
                    "error": -1.9370784404474806,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "b7",
            "ionNumber": 7,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 6
                }
            ],
            "sequence": "QNCcmELFE",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                104
            ],
            "clusterInfo": [
                {
                    "Clusterid": 104,
                    "calcMZ": 1049.435707769349,
                    "error": -5.248315173714634,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b8",
            "ionNumber": 8,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 7
                }
            ],
            "sequence": "QNCcmELFEQ",
            "stub": "",
            "type": "b"
        },
        {
            "class": "lossy",
            "clusterIds": [
                103
            ],
            "clusterInfo": [
                {
                    "Clusterid": 103,
                    "calcMZ": 1031.425147499349,
                    "error": 1.6991059945413365,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b8_H2Ox1",
            "ionNumber": 8,
            "nlosses": 1,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 7
                }
            ],
            "sequence": "QNCcmELFEQ",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                117
            ],
            "clusterInfo": [
                {
                    "Clusterid": 117,
                    "calcMZ": 1219.5412354670489,
                    "error": -2.817015898184973,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b10",
            "ionNumber": 10,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 9
                }
            ],
            "sequence": "QNCcmELFEQLG",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                129
            ],
            "clusterInfo": [
                {
                    "Clusterid": 129,
                    "calcMZ": 1348.583828555019,
                    "error": -9.809219670915793,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b11",
            "ionNumber": 11,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 10
                }
            ],
            "sequence": "QNCcmELFEQLGE",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                138
            ],
            "clusterInfo": [
                {
                    "Clusterid": 138,
                    "calcMZ": 1511.647157087569,
                    "error": 2.0791309772884,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "b12",
            "ionNumber": 12,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 11
                }
            ],
            "sequence": "QNCcmELFEQLGEY",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                17
            ],
            "clusterInfo": [
                {
                    "Clusterid": 17,
                    "calcMZ": 175.11895217417901,
                    "error": -1.0973922389307718,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 20,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "R",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                11
            ],
            "clusterInfo": [
                {
                    "Clusterid": 11,
                    "calcMZ": 147.11280416457902,
                    "error": -3.8349114627194485,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 9,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "K",
            "stub": "",
            "type": "y"
        },
        {
            "class": "lossy",
            "clusterIds": [
                6
            ],
            "clusterInfo": [
                {
                    "Clusterid": 6,
                    "calcMZ": 129.10224389457903,
                    "error": -4.600187890701663,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y1_H2Ox1",
            "ionNumber": 1,
            "nlosses": 1,
            "peptideId": 1,
            "range": [
                {
                    "from": 9,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "K",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                35
            ],
            "clusterInfo": [
                {
                    "Clusterid": 35,
                    "calcMZ": 274.187366087169,
                    "error": -1.298700133697604,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 19,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "VR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                31
            ],
            "clusterInfo": [
                {
                    "Clusterid": 31,
                    "calcMZ": 246.181218077569,
                    "error": -2.3887994933091172,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 8,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "VK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                49
            ],
            "clusterInfo": [
                {
                    "Clusterid": 49,
                    "calcMZ": 387.271430064299,
                    "error": -2.065900649824338,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y3",
            "ionNumber": 3,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 18,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "LVR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                47
            ],
            "clusterInfo": [
                {
                    "Clusterid": 47,
                    "calcMZ": 359.265282054699,
                    "error": -2.1768167927146775,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y3",
            "ionNumber": 3,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 7,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "LVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                63
            ],
            "clusterInfo": [
                {
                    "Clusterid": 63,
                    "calcMZ": 500.355494041429,
                    "error": -1.1472671646986707,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y4",
            "ionNumber": 4,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 17,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "LLVR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                61
            ],
            "clusterInfo": [
                {
                    "Clusterid": 61,
                    "calcMZ": 488.307875142669,
                    "error": -3.7172094930026893,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y4",
            "ionNumber": 4,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 6,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "ELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "lossy",
            "clusterIds": [
                58
            ],
            "clusterInfo": [
                {
                    "Clusterid": 58,
                    "calcMZ": 470.297314872669,
                    "error": -6.453093762214033,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y4_H2Ox1",
            "ionNumber": 4,
            "nlosses": 1,
            "peptideId": 1,
            "range": [
                {
                    "from": 6,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "ELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                67
            ],
            "clusterInfo": [
                {
                    "Clusterid": 67,
                    "calcMZ": 571.3926078261389,
                    "error": -1.2212726054408456,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y5",
            "ionNumber": 5,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 16,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "ALLVR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                69
            ],
            "clusterInfo": [
                {
                    "Clusterid": 69,
                    "calcMZ": 587.376289055659,
                    "error": -1.5646795352608827,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y5",
            "ionNumber": 5,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 5,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "VELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                75
            ],
            "clusterInfo": [
                {
                    "Clusterid": 75,
                    "calcMZ": 685.435535267279,
                    "error": -0.4307732292114382,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y6",
            "ionNumber": 6,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 15,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "NALLVR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                78
            ],
            "clusterInfo": [
                {
                    "Clusterid": 78,
                    "calcMZ": 700.460353032789,
                    "error": 2.5796852072673273,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y6",
            "ionNumber": 6,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 4,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "LVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "lossy",
            "clusterIds": [
                24
            ],
            "clusterInfo": [
                {
                    "Clusterid": 24,
                    "calcMZ": 222.151261808849,
                    "error": -15.538101476000366,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                }
            ],
            "name": "y6_H2Ox2",
            "ionNumber": 6,
            "nlosses": 2,
            "peptideId": 1,
            "range": [
                {
                    "from": 4,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "LVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                87
            ],
            "clusterInfo": [
                {
                    "Clusterid": 87,
                    "calcMZ": 813.494112772559,
                    "error": -0.7163820240096549,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y7",
            "ionNumber": 7,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 14,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "QNALLVR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                82
            ],
            "clusterInfo": [
                {
                    "Clusterid": 82,
                    "calcMZ": 771.497466817499,
                    "error": -3.119669995760977,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y7",
            "ionNumber": 7,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 3,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "ALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                60,
                95
            ],
            "clusterInfo": [
                {
                    "Clusterid": 60,
                    "calcMZ": 480.784901576214,
                    "error": -0.23207096072373026,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                },
                {
                    "Clusterid": 95,
                    "calcMZ": 960.562526685549,
                    "error": -0.2776347624579975,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y8",
            "ionNumber": 8,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 13,
                    "peptideId": 0,
                    "to": 20
                }
            ],
            "sequence": "FQNALLVR",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                90
            ],
            "clusterInfo": [
                {
                    "Clusterid": 90,
                    "calcMZ": 872.545145285909,
                    "error": 1.7703543470826593,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y8",
            "ionNumber": 8,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 2,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "TALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                85,
                116
            ],
            "clusterInfo": [
                {
                    "Clusterid": 85,
                    "calcMZ": 785.4771773238724,
                    "error": 1.3912640102859326,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 116,
                    "calcMZ": 1177.712127752369,
                    "error": -1.8065173307117677,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "y9+P",
            "ionNumber": 9,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 12,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "KFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "lossy",
            "clusterIds": [
                84
            ],
            "clusterInfo": [
                {
                    "Clusterid": 84,
                    "calcMZ": 779.4736572338724,
                    "error": -17.495272288036325,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                }
            ],
            "name": "y9+P_H2Ox1",
            "ionNumber": 9,
            "nlosses": 1,
            "peptideId": 0,
            "range": [
                {
                    "from": 12,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "KFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                100
            ],
            "clusterInfo": [
                {
                    "Clusterid": 100,
                    "calcMZ": 1000.603722791189,
                    "error": 2.375774501699233,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1
                }
            ],
            "name": "y9",
            "ionNumber": 9,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 1,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "QTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "lossy",
            "clusterIds": [
                62
            ],
            "clusterInfo": [
                {
                    "Clusterid": 62,
                    "calcMZ": 491.800219494034,
                    "error": 0.5296987590542648,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "y9_H2Ox1",
            "ionNumber": 9,
            "nlosses": 1,
            "peptideId": 1,
            "range": [
                {
                    "from": 1,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "QTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                88,
                123
            ],
            "clusterInfo": [
                {
                    "Clusterid": 88,
                    "calcMZ": 839.8316201680557,
                    "error": -0.8216385752317645,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 123,
                    "calcMZ": 1259.243792018644,
                    "error": -1.6611744745002,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "y10+P",
            "ionNumber": 10,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 11,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "YKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                127
            ],
            "clusterInfo": [
                {
                    "Clusterid": 127,
                    "calcMZ": 1323.765088562629,
                    "error": -4.976991527585901,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "y11+P",
            "ionNumber": 11,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 10,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "EYKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                92,
                132
            ],
            "clusterInfo": [
                {
                    "Clusterid": 92,
                    "calcMZ": 901.852972437569,
                    "error": 1.2060640308799626,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 132,
                    "calcMZ": 1352.275820422914,
                    "error": -1.9376431011335047,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "y12+P",
            "ionNumber": 12,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 9,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "GEYKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "lossy",
            "clusterIds": [
                91,
                128
            ],
            "clusterInfo": [
                {
                    "Clusterid": 91,
                    "calcMZ": 895.8494523475689,
                    "error": -1.3655875358367577,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 1,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 128,
                    "calcMZ": 1343.270540287914,
                    "error": -10.431450160644264,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 1,
                    "matchedCharge": 2
                }
            ],
            "name": "y12+P_H2Ox1",
            "ionNumber": 12,
            "nlosses": 1,
            "peptideId": 0,
            "range": [
                {
                    "from": 9,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "GEYKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                93,
                134
            ],
            "clusterInfo": [
                {
                    "Clusterid": 93,
                    "calcMZ": 939.5476604299458,
                    "error": 3.7324894364820502,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 1,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 134,
                    "calcMZ": 1408.8178524114792,
                    "error": 1.1696214088294483,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "y13+P",
            "ionNumber": 13,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 8,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "LGEYKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                99,
                137
            ],
            "clusterInfo": [
                {
                    "Clusterid": 99,
                    "calcMZ": 982.2338529317058,
                    "error": 6.763497524220899,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 1,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 137,
                    "calcMZ": 1472.847141164119,
                    "error": -1.1141476080134431,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2
                }
            ],
            "name": "y14+P",
            "ionNumber": 14,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 7,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "QLGEYKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "lossy",
            "clusterIds": [
                98
            ],
            "clusterInfo": [
                {
                    "Clusterid": 98,
                    "calcMZ": 976.2303328417057,
                    "error": -11.780387630739282,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 1,
                    "matchedCharge": 3
                }
            ],
            "name": "y14+P_H2Ox1",
            "ionNumber": 14,
            "nlosses": 1,
            "peptideId": 0,
            "range": [
                {
                    "from": 7,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "QLGEYKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                106,
                142
            ],
            "clusterInfo": [
                {
                    "Clusterid": 106,
                    "calcMZ": 1074.2708552653592,
                    "error": -7.498235031411179,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3
                },
                {
                    "Clusterid": 142,
                    "calcMZ": 1610.902644664599,
                    "error": 11.342982113670729,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 1,
                    "matchedCharge": 2
                }
            ],
            "name": "y16+P",
            "ionNumber": 16,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 5,
                    "peptideId": 0,
                    "to": 20
                },
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 9
                }
            ],
            "sequence": "FEQLGEYKFQNALLVR + KQTALVELVK",
            "stub": "",
            "type": "y"
        }
    ]
}


def test_annotate_loss_and_mod_xi1_format(client):
    """
    Test the annotation of the synthetic spectrum of QNCcmELFEQLGEYKFQNALLVR-KQTALVELVK_12-0:4 (BS3)
    using the xi1 json style.

    cm modification defined in request
    H2O loss defined in request
    """
    url = url_for('xi2annotator.annotate')

    # Test annotation with old config format
    current_dir = os.path.dirname(__file__)
    # Load request with synthetic spectrum for AKT-KMR_1-0:3 BS3 with 3 extra isotope peaks
    json_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                             'xi1_format_QNCcmELFEQLGEYKFQNALLVR-KQTALVELVK_12-0_z4_BS3.json')
    with open(json_file) as f:
        request = json.load(f)

    # request annotation
    res = client.post(url, json=request)
    assert res._status_code == 200
    check_unchanged_values(res.json, request, 'xi1')
    check_result(res.json, exp_loss_and_mod)
    # check that the modifications have been written out
    exp_mods = [{'aminoAcids': ['C'], 'id': 'cm', 'mass': 57.0215}]
    assert res.json['annotation']['modifications'] == exp_mods


def test_annotate_loss_and_mod_xi2_format(client):
    """
    Test the annotation of the synthetic spectrum of QNCcmELFEQLGEYKFQNALLVR-KQTALVELVK_12-0:4 (BS3)
    using the xi2 json style.

    cm modification defined in request
    H2O loss defined in request
    """
    url = url_for('xi2annotator.annotate')

    # Test annotation with old config format
    current_dir = os.path.dirname(__file__)
    # Load request with synthetic spectrum for AKT-KMR_1-0:3 BS3 with 3 extra isotope peaks
    json_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                             'xi2_format_QNCcmELFEQLGEYKFQNALLVR-KQTALVELVK_12-0_z4_BS3.json')
    with open(json_file) as f:
        request = json.load(f)

    # request annotation
    res = client.post(url, json=request)
    assert res._status_code == 200
    check_unchanged_values(res.json, request, 'xi2')
    check_result(res.json, exp_loss_and_mod)
    # check that the modifications have been written out
    exp_mods = [{'aminoAcids': ['C'], 'id': 'cm', 'mass': 57.0215}]
    assert res.json['annotation']['modifications'] == exp_mods


exp_noncov_synthetic = {
    # synthetic spectrum without isotope peaks so we should have a single peak per cluster
    'peak_cluster_ids': [[x] for x in range(30)],
    # we don't check this for 'precursor' values atm
    # 'precursor': {
    #     'calcMZ': ,
    #     'error': '',
    #     'errorUnit': 'ppm',
    # },
    # check clusters - only single peaks so all clusters should have charge 0
    'cluster': {
        'charges': [0] * 30,
        'firstPeakIds': list(range(30)),
    },
    # no losses defined -> default losses
    'losses': default_losses,
    # manually verified matching fragments on spectrumviewer.org
    'fragments': [
        {
            "class": "non-lossy",
            "clusterIds": [
                18,
                24,
                29
            ],
            "clusterInfo": [
                {
                    "Clusterid": 18,
                    "calcMZ": 138.43013323339233,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 24,
                    "calcMZ": 207.14156161664897,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 29,
                    "calcMZ": 413.27584676641897,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "P",
            "ionNumber": 3,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 2
                }
            ],
            "sequence": "LAKsda",
            "stub": "",
            "type": "P"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                16,
                21,
                28
            ],
            "clusterInfo": [
                {
                    "Clusterid": 16,
                    "calcMZ": 121.73773732687233,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 21,
                    "calcMZ": 182.102967756869,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 28,
                    "calcMZ": 363.198659046859,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "P",
            "ionNumber": 3,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 2
                }
            ],
            "sequence": "TSR",
            "stub": "",
            "type": "P"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                1,
                3,
                14
            ],
            "clusterInfo": [
                {
                    "Clusterid": 1,
                    "calcMZ": 38.70196445925566,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 3,
                    "calcMZ": 57.549308455444,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 14,
                    "calcMZ": 114.09134044400899,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "b1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 0
                }
            ],
            "sequence": "L",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                0,
                2,
                13
            ],
            "clusterInfo": [
                {
                    "Clusterid": 0,
                    "calcMZ": 34.689835956349,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 2,
                    "calcMZ": 51.531115701084005,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 13,
                    "calcMZ": 102.054954935289,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "b1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 0
                }
            ],
            "sequence": "T",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                5,
                10,
                22
            ],
            "clusterInfo": [
                {
                    "Clusterid": 5,
                    "calcMZ": 62.38100238749234,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 10,
                    "calcMZ": 93.067865347799,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 22,
                    "calcMZ": 185.128454228719,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "b2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 0,
                    "peptideId": 0,
                    "to": 1
                }
            ],
            "sequence": "LA",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                6,
                11,
                23
            ],
            "clusterInfo": [
                {
                    "Clusterid": 6,
                    "calcMZ": 63.70051209110567,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 11,
                    "calcMZ": 95.047129903219,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 23,
                    "calcMZ": 189.086983339559,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "b2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 0,
                    "peptideId": 1,
                    "to": 1
                }
            ],
            "sequence": "TS",
            "stub": "",
            "type": "b"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                7,
                15,
                25
            ],
            "clusterInfo": [
                {
                    "Clusterid": 7,
                    "calcMZ": 77.056407312779,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 15,
                    "calcMZ": 115.08097273572899,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 25,
                    "calcMZ": 229.15466900457898,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "y1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 2,
                    "peptideId": 0,
                    "to": 2
                }
            ],
            "sequence": "Ksda",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                4,
                9,
                20
            ],
            "clusterInfo": [
                {
                    "Clusterid": 4,
                    "calcMZ": 59.04450170264567,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 9,
                    "calcMZ": 88.063114320529,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 20,
                    "calcMZ": 175.118952174179,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "y1",
            "ionNumber": 1,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 2,
                    "peptideId": 1,
                    "to": 2
                }
            ],
            "sequence": "R",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                12,
                19,
                27
            ],
            "clusterInfo": [
                {
                    "Clusterid": 12,
                    "calcMZ": 100.73544524101567,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 19,
                    "calcMZ": 150.59952962808399,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 27,
                    "calcMZ": 300.191782789289,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "y2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 0,
            "range": [
                {
                    "from": 1,
                    "peptideId": 0,
                    "to": 2
                }
            ],
            "sequence": "AKsda",
            "stub": "",
            "type": "y"
        },
        {
            "class": "non-lossy",
            "clusterIds": [
                8,
                17,
                26
            ],
            "clusterInfo": [
                {
                    "Clusterid": 8,
                    "calcMZ": 88.05517783740233,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 3,
                },
                {
                    "Clusterid": 17,
                    "calcMZ": 131.57912852266398,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 2,
                },
                {
                    "Clusterid": 26,
                    "calcMZ": 262.150980578449,
                    "error": 0,
                    "errorUnit": "ppm",
                    "matchedMissingMonoIsotopic": 0,
                    "matchedCharge": 1,
                }
            ],
            "name": "y2",
            "ionNumber": 2,
            "nlosses": 0,
            "peptideId": 1,
            "range": [
                {
                    "from": 1,
                    "peptideId": 1,
                    "to": 2
                }
            ],
            "sequence": "SR",
            "stub": "",
            "type": "y"
        }
    ]
}


def test_annotate_synthetic_noncovalent_xi1_format(client):
    """
    Test the annotation of the synthetic spectrum of noncovalently associated peptides
    LAsdaK & TSR using the xi1 json style.
    """
    url = url_for('xi2annotator.annotate')

    # Test annotation with old config format
    current_dir = os.path.dirname(__file__)
    # Load request with synthetic spectrum for LAsdaK TSR as noncovalently associated peptides
    json_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                             'xi1_format_LAsdaK-TSR_z3_NAP.json')
    with open(json_file) as f:
        request = json.load(f)

    # request annotation
    res = client.post(url, json=request)
    assert res._status_code == 200

    check_unchanged_values(res.json, request, 'xi1')
    check_result(res.json, exp_noncov_synthetic)
    # check that the modifications have been written out
    exp_mods = [{'aminoAcids': ['K'], 'id': 'sda', 'mass': 82.04186484}]
    assert res.json['annotation']['modifications'] == exp_mods


def test_annotate_synthetic_noncovalent_xi2_format(client):
    """
    Test the annotation of the synthetic spectrum of noncovalently associated peptides
    LAsdaK & TSR using the xi2 json style.
    """
    url = url_for('xi2annotator.annotate')

    # Test annotation with old config format
    current_dir = os.path.dirname(__file__)
    # Load request with synthetic spectrum for LAsdaK TSR as noncovalently associated peptides
    json_file = os.path.join(current_dir, '../fixtures', 'annotation_requests',
                             'xi2_format_LAsdaK-TSR_z3_NAP.json')
    with open(json_file) as f:
        request = json.load(f)

    # request annotation
    res = client.post(url, json=request)
    assert res._status_code == 200

    check_unchanged_values(res.json, request, 'xi2')
    check_result(res.json, exp_noncov_synthetic)
    # check that the modifications have been written out
    exp_mods = [{'aminoAcids': ['K'], 'id': 'sda', 'mass': 82.04186484}]
    assert res.json['annotation']['modifications'] == exp_mods


def test_annotate_with_composition_modification(client):
    """
    Test that the annotator reports the modification mass when chemical composition was defined.
    Also checks that the annotator doesn't fail when there are no annotated fragments
    """
    url = url_for('xi2annotator.annotate')

    # basic request with a modification that has a composition defined and no mass
    request = {
        "Peptides": [
            {
                "sequence": [
                    {"aminoAcid": "A", "Modification": ""},
                    {"aminoAcid": "M", "Modification": "ox"},
                ]
            },
        ],
        "peaks": [{"mz": 104.31995, "intensity": 1537.41}],
        "annotation": {
            "precursorCharge": 2,
            "config": {
                "ms2_tol": "20 ppm",
                "fragmentation": {"nterm_ions": ["b"], "cterm_ions": ["y"]},
                "modification": {
                    "modifications": [
                        {
                            "specificity": ["M"],
                            "name": "ox",
                            "composition": "O",
                            "type": "variable"
                        }
                    ]
                }
            }
        }
    }

    # request annotation
    res = client.post(url, json=request)
    assert res._status_code == 200

    check_unchanged_values(res.json, request, 'xi2')
    # check that the modifications has been written out with a mass
    exp_mods = [{'aminoAcids': ['M'], 'id': 'ox', 'mass': 15.99491461956}]
    assert res.json['annotation']['modifications'] == exp_mods
