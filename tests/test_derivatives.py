import pytest
from dynare_preprocessor import DynareModel
import numpy as np

files = [
    "fs2000.mod",
    "model_KR2000_IRF.mod",
    "model_KR2000_STAT.mod",
    "RBC.mod",
    "fs2000_nonstationary.mod",
    "example2.mod",
]

def jacobian(func, initial, delta: float = 1e-10):
    f = func
    f0 = f(initial)
    nrow = len(f0)
    ncol = len(initial)
    output = np.zeros((nrow, ncol))
    for j in range(ncol):
        ej = np.zeros(ncol)
        ej[j] = 1
        x = (initial + delta * ej).reshape(ncol)
        dj = (f(x) - f0) / (delta)
        output[:, j] = dj

    return output

@pytest.mark.parametrize("filename", files)
def test_jacobians(filename):
    f = filename
    filename = "tests/modfiles/" + f
    mod_string = open(filename).read()
    model = DynareModel(mod_string, 1, 1)
    endo = [model.context[x] for x in model.endogenous]
    exo = [model.context[x] for x in model.exogenous]
    exo_det = [model.context[x] for x in model.exogenous_det]
    params = [model.context[x] for x in model.parameters]
    lengths = [len(model.endogenous)] * 3 + [len(model.exogenous), len(model.exogenous_det), len(model.parameters)]
    n = len(model.equations)
    preprocessor_jacobians = [np.zeros((n, l)) for l in lengths]
    sparse_jacobs = model.jacobians(endo,endo,endo,exo,exo_det,params)
    for i, sparse_jacob in enumerate(sparse_jacobs):
        for (k,l), v in sparse_jacob.items():
            preprocessor_jacobians[i][k,l] = v
    
    y = np.array(endo)
    e = np.array(exo)
    ed = np.array(exo_det)
    p = np.array(params)
    dyn = lambda u,v,w,x,y,z: np.array(model.dynamic_function(u,v,w,x,y,z))

    fin_diff_jacobians = [
        jacobian(func, init) for (func,init) in [
            (lambda u: dyn(u, y, y, e, ed, p), y),
            (lambda u: dyn(y, u, y, e, ed, p), y),
            (lambda u: dyn(y, y, u, e, ed, p), y),
            (lambda u: dyn(y, y, y, u, ed, p), e),
            (lambda u: dyn(y, y, y, e, u, p), ed),
            (lambda u: dyn(y, y, y, e, ed, u), p)
        ]
    ]
    assert(len(preprocessor_jacobians) == 6)
    for i in range(6):
        if not np.allclose(preprocessor_jacobians[i], fin_diff_jacobians[i], atol=1e-3):
            print(preprocessor_jacobians[i])
            print(fin_diff_jacobians[i])
            assert(False)