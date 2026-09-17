import os
import sys
import numpy as np
import pytest

build_src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "src"))
if build_src not in sys.path:
    sys.path.insert(0, build_src)

from dynare_preprocessor import DynareModel

files = [
    "fs2000.mod",
    "model_KR2000_IRF.mod",
    "model_KR2000_STAT.mod",
    "RBC.mod",
    "fs2000_nonstationary.mod",
    "example2.mod",
]


def finite_diff_jac(func, initial, delta: float = 1e-7):
    f = func
    f0 = f(initial)
    nrow = len(f0)
    ncol = len(initial)
    output = np.zeros((nrow, ncol))
    for j in range(ncol):
        ej = np.zeros(ncol)
        ej[j] = 1.0
        x_plus = initial + delta * ej
        x_minus = initial - delta * ej
        output[:, j] = (f(x_plus) - f(x_minus)) / (2.0 * delta)
    return output


@pytest.mark.parametrize("filename", files)
def test_jacobians(filename):
    path = os.path.join(os.path.dirname(__file__), "modfiles", filename)
    with open(path, encoding="utf-8") as f:
        mod_string = f.read()
    model = DynareModel(mod_string, 1, 1)
    endo = [model.context[x] for x in model.endogenous]
    exo = [model.context[x] for x in model.exogenous]
    exo_det = [model.context[x] for x in model.exogenous_det]
    params = [model.context[x] for x in model.parameters]

    # Test both jacobian_blocks and legacy jacobians
    jb = model.jacobian_blocks(endo, endo, endo, exo, exo_det, params)
    blocks = [jb.lead, jb.curr, jb.lag, jb.exo, jb.exo_det, jb.params]

    y = np.array(endo)
    e = np.array(exo)
    ed = np.array(exo_det)
    p = np.array(params)
    dyn = lambda u, v, w, x, y_v, z: np.array(model.residuals(u, v, w, x, y_v, z))

    fin_diff_jacobians = [
        finite_diff_jac(func, init)
        for (func, init) in [
            (lambda u: dyn(u, y, y, e, ed, p), y),
            (lambda u: dyn(y, u, y, e, ed, p), y),
            (lambda u: dyn(y, y, u, e, ed, p), y),
            (lambda u: dyn(y, y, y, u, ed, p), e),
            (lambda u: dyn(y, y, y, e, u, p), ed),
            (lambda u: dyn(y, y, y, e, ed, u), p),
        ]
    ]

    assert len(blocks) == 6
    for i in range(6):
        if blocks[i].shape[1] > 0:
            assert np.allclose(blocks[i], fin_diff_jacobians[i], atol=1e-4)


@pytest.mark.parametrize("filename", files)
def test_derivatives(filename):
    deriv_order = 3
    path = os.path.join(os.path.dirname(__file__), "modfiles", filename)
    with open(path, encoding="utf-8") as f:
        mod_string = f.read()
    model = DynareModel(mod_string, deriv_order, 1)
    endo = [model.context[x] for x in model.endogenous]
    exo = [model.context[x] for x in model.exogenous]
    exo_det = [model.context[x] for x in model.exogenous_det]
    params = [model.context[x] for x in model.parameters]
    derivatives = model.derivatives(endo, endo, endo, exo, exo_det, params)
    for i in range(deriv_order):
        for lst, val in derivatives[i]:
            assert len(lst) == i + 2

