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

"""
The conftest.py file serves as a means of providing fixtures for an entire directory.
Fixtures defined in a conftest.py can be used by any test in that package without needing to
import them (pytest will automatically discover them).
"""

import pytest

from xicommon.config import Config, Crosslinker, Enzyme, DigestionConfig


@pytest.fixture()
def search_config():
    # Basic Search config to use for tests
    return Config(
        reporting_requirements={'report_top_ranking_only': False},
        digestion=DigestionConfig(
            enzymes=[Enzyme.trypsin], missed_cleavages=0, min_peptide_length=3
        ),
        ms1_tol="5ppm",
        ms2_tol="15ppm",
        fragmentation={"add_precursor": False},
        top_n_alpha_scores=10,
        crosslinker=[Crosslinker.BS3],
    )
