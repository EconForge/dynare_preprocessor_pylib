"""Arbitrary-order perturbation with truncated power series (TaylorPerturbation)."""

import itertools
import os
import sys

import numpy as np
import pytest

build_src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "src"))
if build_src not in sys.path:
    sys.path.insert(0, build_src)

from dynare_preprocessor import (
    DynareModel,
    EvaluationException,
    ModelSemanticException,
    TaylorPerturbation,
)

MODFILES = os.path.join(os.path.dirname(__file__), "modfiles")


def load(name):
    return DynareModel(os.path.join(MODFILES, name), stochastic=True)


@pytest.fixture(scope="module")
def rbc():
    return load("RBC.mod")


@pytest.fixture(scope="module")
def tp3(rbc):
    return TaylorPerturbation(rbc, 3).solve()


def params(model):
    return [model.context[p] for p in model.parameters]


def euler_residual(model, tp, w, order, nq=8):
    """E f at w = (xhat, u, sigma), by Gauss-Hermite quadrature over eps'."""
    nx, m = tp.nx, tp.m
    x, wq = np.polynomial.hermite_e.hermegauss(nq)
    wq = wq / wq.sum()
    L = np.linalg.cholesky(tp.Sigma)
    sigma = w[nx + m]
    y0 = tp.g(w, order)
    y_lag = tp.ybar.copy()
    y_lag[tp.states] += w[:nx]
    p = params(model)
    res = np.zeros(tp.n)
    for nodes in itertools.product(range(nq), repeat=m):
        eps = L @ x[list(nodes)]
        weight = np.prod(wq[list(nodes)])
        w1 = np.concatenate([(y0 - tp.ybar)[tp.states], sigma * eps, [sigma]])
        y1 = tp.g(w1, order)
        res += weight * model.residuals(y1, y0, y_lag, w[nx : nx + m], [], p)
    return res


def test_dimensions(rbc, tp3):
    assert tp3.n == len(rbc.endogenous)
    assert tp3.m == len(rbc.exogenous)
    assert [rbc.endogenous[s] for s in tp3.states] == ["k", "a"]
    assert tp3.n_w == tp3.nx + tp3.m + 1
    # monomials of degree <= 3 in n_w = 4 variables
    assert tp3.exponents.shape == (35, 4)
    assert tp3.coefficients.shape == (tp3.n, 35)
    assert tp3.exponents.sum(axis=1).max() == 3
    assert np.allclose(tp3.Sigma, [[0.009**2]])


def test_order_k_residuals_vanish(tp3):
    assert tp3.solved
    assert tp3.max_residual < 1e-12
    for k in range(1, 4):
        assert np.abs(tp3.residual(k)).max() < 1e-12


def test_first_order(rbc, tp3):
    # g_x solves the quadratic matrix equation, g_u = -(f_0 + f_p g_x)^{-1} f_e
    X = np.zeros((tp3.n, tp3.n))
    X[:, tp3.states] = tp3.coefficients[:, 1 : 1 + tp3.nx]
    assert np.abs(tp3.f_m + tp3.f_0 @ X + tp3.f_p @ X @ X).max() < 1e-12
    g_u = -np.linalg.solve(tp3.f_0 + tp3.f_p @ X, tp3.f_e)
    assert np.allclose(tp3.coef(u=[0]), g_u[:, 0], atol=1e-12)
    # technology: log(a) = rho log(a(-1)) + epsilon
    a = rbc.endogenous.index("a")
    assert np.isclose(tp3.coef(u=[0])[a], 1.0)
    assert np.isclose(tp3.coef(u=[0, 0])[a], 0.5)  # exp(u) = 1 + u + u^2/2 + ...


def test_steady_state_and_decision_rule(rbc, tp3):
    assert np.abs(rbc.residuals(tp3.ybar, tp3.ybar, tp3.ybar, [0.0], [], params(rbc))).max() < 1e-12
    assert np.allclose(tp3.g(np.zeros(tp3.n_w)), tp3.ybar)
    # coefficients and exponents describe the same polynomial as g
    rng = np.random.default_rng(0)
    w = rng.normal(size=(5, tp3.n_w)) * 0.1
    mon = np.prod(w[:, None, :] ** tp3.exponents[None, :, :], axis=2)
    assert np.allclose(tp3.g(w), tp3.ybar + mon @ tp3.coefficients.T, atol=1e-14)


@pytest.mark.parametrize("states_only", [False, True])
@pytest.mark.parametrize("sylvester", ["schur", "block"])
def test_solver_variants_agree(rbc, tp3, states_only, sylvester):
    tp = TaylorPerturbation(rbc, 3, states_only=states_only, sylvester=sylvester).solve()
    assert tp.max_residual < 1e-12
    pos = {s: p for p, s in enumerate(tp.states)}
    rng = np.random.default_rng(0)
    for _ in range(5):
        x = rng.normal(size=tp3.nx) * 0.05
        u = rng.normal(size=tp3.m) * 0.01
        x_full = np.zeros(tp.nx)
        for p, s in enumerate(tp3.states):
            x_full[pos[s]] = x[p]
        w_ref = np.concatenate([x, u, [0.7]])
        w = np.concatenate([x_full, u, [0.7]])
        assert np.allclose(tp.g(w), tp3.g(w_ref), atol=1e-13)


def test_convergence_rates(rbc, tp3):
    """The Euler residual of the order-k truncation along a ray w = t d is O(t^(k+1))."""
    rng = np.random.default_rng(0)
    d = np.zeros(tp3.n_w)
    scale = np.maximum(np.abs(tp3.ybar[tp3.states]), 1.0) * 0.05
    d[: tp3.nx] = rng.normal(size=tp3.nx) * scale
    d[tp3.nx : tp3.nx + tp3.m] = rng.normal(size=tp3.m) * np.sqrt(np.diag(tp3.Sigma))
    d[-1] = 1.0
    ts = [0.8, 0.4, 0.2]
    table = np.array(
        [[np.abs(euler_residual(rbc, tp3, t * d, k)).max() for k in (1, 2, 3)] for t in ts]
    )
    rates = np.log2(table[:-1] / table[1:])
    for k in range(3):
        assert np.allclose(rates[:, k], k + 2, atol=0.1)


def test_steady_state_operator():
    # Occbin_example uses STEADY_STATE(Invest) inside the model block
    model = load("Occbin_example.mod")
    tp = TaylorPerturbation(model, 3).solve()
    assert tp.max_residual < 1e-12


def test_unit_root():
    # Gali (2015) has an eigenvalue exactly on the unit circle
    model = load("Gali_2015.mod")
    tp = TaylorPerturbation(model, 2).solve()
    assert tp.max_residual < 1e-12


def test_explicit_inputs(rbc, tp3):
    tp = TaylorPerturbation(rbc, 2, ybar=tp3.ybar, params=params(rbc), Sigma=4 * tp3.Sigma).solve()
    # the risk correction (sigma^2 term) scales with the variance
    assert np.allclose(tp.coef(s=2), 4 * tp3.coef(s=2), rtol=1e-10, atol=1e-16)
    assert np.abs(tp.coef(s=2)).max() > 0


def test_errors(rbc):
    with pytest.raises(ModelSemanticException):
        TaylorPerturbation(rbc, 0)
    with pytest.raises(ModelSemanticException):
        TaylorPerturbation(rbc, 2, sylvester="nope")
    with pytest.raises(ModelSemanticException):
        TaylorPerturbation(rbc, 2, ybar=[1.0, 2.0])
    tp = TaylorPerturbation(rbc, 2)
    with pytest.raises(ValueError):
        tp.coef(x=[0, 0, 1])  # degree 3 > order 2
    deterministic = DynareModel(
        """
        var x;
        parameters rho;
        rho = 0.5;
        model;
        x = rho * x(-1);
        end;
        """
    )
    with pytest.raises(ModelSemanticException):
        TaylorPerturbation(deterministic, 2)
    explosive = DynareModel(
        """
        var x;
        varexo e;
        parameters rho;
        rho = 1.5;
        model;
        x = rho * x(-1) + e;
        end;
        shocks;
        var e; stderr 0.1;
        end;
        """,
        stochastic=True,
    )
    with pytest.raises(EvaluationException, match="Blanchard-Kahn"):
        TaylorPerturbation(explosive, 2).solve()
