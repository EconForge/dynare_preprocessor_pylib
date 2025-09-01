from dynare_preprocessor import DynareModel
import pytest

declared_steady_state = [
    "model_KR2000_IRF.mod",
    "model_KR2000_STAT.mod",
    "Occbin_example.mod",
    "RBC.mod"
]

@pytest.mark.parametrize("filename", declared_steady_state)
def test_steady_state(filename):
    from math import isclose
    f = filename
    filename = "tests/modfiles/" + f
    mod_string = open(filename).read()
    model = DynareModel(mod_string)
    endo = [model.context[x] for x in model.endogenous]
    exo = [model.context[x] for x in model.exogenous]
    exo_det = [model.context[x] for x in model.exogenous_det]
    params = [model.context[x] for x in model.parameters]
    for x in model.residuals(endo,endo,endo,exo,exo_det,params):
        assert(isclose(x,0, abs_tol=1e-10))