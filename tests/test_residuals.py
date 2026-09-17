import os
import sys
from math import isclose
import pytest

build_src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "src"))
if build_src not in sys.path:
    sys.path.insert(0, build_src)

from dynare_preprocessor import DynareModel

declared_steady_state = [
    "model_KR2000_IRF.mod",
    "model_KR2000_STAT.mod",
    "Occbin_example.mod",
    "RBC.mod",
]


@pytest.mark.parametrize("filename", declared_steady_state)
def test_steady_state(filename):
    path = os.path.join(os.path.dirname(__file__), "modfiles", filename)
    with open(path, encoding="utf-8") as f:
        mod_string = f.read()
    model = DynareModel(mod_string)
    endo = [model.context[x] for x in model.endogenous]
    exo = [model.context[x] for x in model.exogenous]
    exo_det = [model.context[x] for x in model.exogenous_det]
    params = [model.context[x] for x in model.parameters]
    for x in model.residuals(endo, endo, endo, exo, exo_det, params):
        assert isclose(x, 0, abs_tol=1e-10)

