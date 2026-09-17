import os
import sys
import numpy as np
import pytest

# Ensure build/src is in sys.path if not installed
build_src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "src"))
if build_src not in sys.path:
    sys.path.insert(0, build_src)

import dynare_preprocessor as dp


def finite_difference_jacobian(func, x0, delta=1e-7):
    """Compute numerical Jacobian using central differences."""
    f0 = func(x0)
    m = len(f0)
    n = len(x0)
    jac = np.zeros((m, n))
    for j in range(n):
        dx = np.zeros(n)
        dx[j] = delta
        f_plus = func(x0 + dx)
        f_minus = func(x0 - dx)
        jac[:, j] = (f_plus - f_minus) / (2.0 * delta)
    return jac


def test_model_loading_and_properties():
    """Test loading a model from file and inspecting properties."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "RBC.mod")
    model = dp.DynareModel(modfile)

    assert "c" in model.endogenous
    assert "k" in model.endogenous
    assert "y" in model.endogenous
    assert "epsilon" in model.exogenous
    assert len(model.exogenous_det) == 0
    assert "alpha" in model.parameters
    assert "beta" in model.parameters

    assert len(model.equations) == len(model.endogenous)
    assert np.array(model.lead_lag_incidence).shape == (len(model.endogenous), 3)
    assert model.max_endo_lag >= 0
    assert model.max_endo_lead >= 0

    # Test context dictionary
    assert "alpha" in model.context
    assert model.context["alpha"] == 0.33


def test_in_memory_string_and_macro():
    """Test loading a model from in-memory string with macroprocessor directives."""
    mod_text = """
    @#define HAS_SHOCK = 1
    var y c;
    @#if HAS_SHOCK
    varexo e;
    @#endif
    parameters alpha beta;
    alpha = 0.35;
    beta = 0.99;
    model;
    c = y;
    y = alpha * c(-1) + e;
    end;
    initval;
    y = 0;
    c = 0;
    e = 0;
    end;
    """
    model = dp.DynareModel(mod_text)
    assert model.endogenous == ["y", "c"]
    assert model.exogenous == ["e"]
    assert model.parameters == ["alpha", "beta"]
    assert len(model.equations) == 2


def test_residuals_at_steady_state():
    """Test residual evaluation at steady state."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "RBC.mod")
    model = dp.DynareModel(modfile)

    y_ss = np.array([model.context[v] for v in model.endogenous])
    e_ss = np.array([model.context[v] for v in model.exogenous])
    ed_ss = np.array([model.context[v] for v in model.exogenous_det])
    p_ss = np.array([model.context[v] for v in model.parameters])

    # Dynamic residuals with y_{t+1} = y_t = y_{t-1} = y_ss
    res = model.residuals(y_ss, y_ss, y_ss, e_ss, ed_ss, p_ss)
    assert isinstance(res, np.ndarray)
    assert res.shape == (len(model.endogenous),)
    assert np.allclose(res, 0.0, atol=1e-10)

    # Static residuals
    sres = model.static_residuals(y_ss, e_ss, ed_ss, p_ss)
    assert isinstance(sres, np.ndarray)
    assert sres.shape == (len(model.endogenous),)
    assert np.allclose(sres, 0.0, atol=1e-10)


def test_jacobian_blocks_and_unpacking():
    """Test jacobian_blocks access patterns and unpacking."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "RBC.mod")
    model = dp.DynareModel(modfile, derivs_order=1, params_derivs_order=1)

    y = np.array([model.context[v] for v in model.endogenous])
    e = np.array([model.context[v] for v in model.exogenous])
    ed = np.array([model.context[v] for v in model.exogenous_det])
    p = np.array([model.context[v] for v in model.parameters])

    jb = model.jacobian_blocks(y, y, y, e, ed, p)
    assert isinstance(jb, dp.JacobianBlocks)
    assert len(jb) == 6

    # Test named property access
    assert jb.lead.shape == (len(model.endogenous), len(model.endogenous))
    assert jb.curr.shape == (len(model.endogenous), len(model.endogenous))
    assert jb.lag.shape == (len(model.endogenous), len(model.endogenous))
    assert jb.exo.shape == (len(model.endogenous), len(model.exogenous))
    assert jb.exo_det.shape == (len(model.endogenous), len(model.exogenous_det))
    assert jb.params.shape == (len(model.endogenous), len(model.parameters))

    # Test indexing
    assert np.array_equal(jb[0], jb.lead)
    assert np.array_equal(jb[1], jb.curr)
    assert np.array_equal(jb[2], jb.lag)
    assert np.array_equal(jb[3], jb.exo)
    assert np.array_equal(jb[4], jb.exo_det)
    assert np.array_equal(jb[5], jb.params)

    # Test sequence unpacking
    lead, curr, lag, exo, exo_det, params = jb
    assert np.array_equal(lead, jb.lead)
    assert np.array_equal(curr, jb.curr)
    assert np.array_equal(lag, jb.lag)
    assert np.array_equal(exo, jb.exo)
    assert np.array_equal(exo_det, jb.exo_det)
    assert np.array_equal(params, jb.params)


def test_jacobian_blocks_vs_finite_differences():
    """Verify analytic Jacobian blocks against finite-difference derivatives."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "RBC.mod")
    model = dp.DynareModel(modfile, derivs_order=1, params_derivs_order=1)

    y = np.array([model.context[v] for v in model.endogenous])
    e = np.array([model.context[v] for v in model.exogenous])
    ed = np.array([model.context[v] for v in model.exogenous_det])
    p = np.array([model.context[v] for v in model.parameters])

    jb = model.jacobian_blocks(y, y, y, e, ed, p)

    # Partial derivative wrt future endogenous y_{t+1}
    fd_lead = finite_difference_jacobian(lambda u: model.residuals(u, y, y, e, ed, p), y)
    assert np.allclose(jb.lead, fd_lead, atol=1e-5)

    # Partial derivative wrt present endogenous y_t
    fd_curr = finite_difference_jacobian(lambda u: model.residuals(y, u, y, e, ed, p), y)
    assert np.allclose(jb.curr, fd_curr, atol=1e-5)

    # Partial derivative wrt past endogenous y_{t-1}
    fd_lag = finite_difference_jacobian(lambda u: model.residuals(y, y, u, e, ed, p), y)
    assert np.allclose(jb.lag, fd_lag, atol=1e-5)

    # Partial derivative wrt shocks e_t
    fd_exo = finite_difference_jacobian(lambda u: model.residuals(y, y, y, u, ed, p), e)
    assert np.allclose(jb.exo, fd_exo, atol=1e-5)

    # Partial derivative wrt parameters
    fd_params = finite_difference_jacobian(lambda u: model.residuals(y, y, y, e, ed, u), p)
    assert np.allclose(jb.params, fd_params, atol=1e-5)


def test_dynamic_jacobian_matrix():
    """Verify combined dynamic Jacobian matrix shape and non-zero entries."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "RBC.mod")
    model = dp.DynareModel(modfile)

    y = np.array([model.context[v] for v in model.endogenous])
    e = np.array([model.context[v] for v in model.exogenous])
    ed = np.array([model.context[v] for v in model.exogenous_det])
    p = np.array([model.context[v] for v in model.parameters])

    jac = model.jacobian(y, y, y, e, ed, p)
    assert isinstance(jac, np.ndarray)
    # Total dynamic Jacobian columns = non-zero entries in lead_lag_incidence + exo + exo_det
    assert jac.shape[0] == len(model.endogenous)
    assert jac.shape[1] > 0


def test_static_jacobian_vs_finite_differences():
    """Verify analytic static Jacobian against finite differences."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "RBC.mod")
    model = dp.DynareModel(modfile)

    y = np.array([model.context[v] for v in model.endogenous])
    e = np.array([model.context[v] for v in model.exogenous])
    ed = np.array([model.context[v] for v in model.exogenous_det])
    p = np.array([model.context[v] for v in model.parameters])

    sjac = model.static_jacobian(y, e, ed, p)
    assert isinstance(sjac, np.ndarray)
    assert sjac.shape == (len(model.endogenous), len(model.endogenous))

    fd_sjac = finite_difference_jacobian(lambda u: model.static_residuals(u, e, ed, p), y)
    assert np.allclose(sjac, fd_sjac, atol=1e-5)


def test_higher_order_derivatives():
    """Verify computing pass at higher orders and COO format structure."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "RBC.mod")
    order = 3
    model = dp.DynareModel(modfile, derivs_order=order, params_derivs_order=0)

    y = np.array([model.context[v] for v in model.endogenous])
    e = np.array([model.context[v] for v in model.exogenous])
    ed = np.array([model.context[v] for v in model.exogenous_det])
    p = np.array([model.context[v] for v in model.parameters])

    derivs = model.derivatives(y, y, y, e, ed, p)
    assert len(derivs) == order
    for ord_idx in range(order):
        # Order ord_idx+1 has tuples with ord_idx+2 indices (eq_idx, var_1, ..., var_{ord_idx+1})
        for coords, val in derivs[ord_idx]:
            assert len(coords) == ord_idx + 2
            assert isinstance(val, float)


def test_exceptions():
    """Test structured exception throwing and Python hierarchy."""
    # Syntax error -> ParserException
    with pytest.raises(dp.ParserException) as exc_info:
        dp.DynareModel("var x; model; x = ; end;")
    assert "syntax error" in str(exc_info.value).lower()
    assert issubclass(dp.ParserException, dp.DynareException)
    assert issubclass(dp.ParserException, dp.SourceFileException)

    # Equation count mismatch -> ModelSemanticException
    with pytest.raises(dp.ModelSemanticException) as exc_info:
        dp.DynareModel("var x; model; x = 1; x = 2; end;")
    assert "2 equations but 1 endogenous" in str(exc_info.value)
    assert issubclass(dp.ModelSemanticException, dp.DynareException)

    # Mixed stochastic and perfect foresight -> StatementException
    with pytest.raises(dp.StatementException) as exc_info:
        dp.DynareModel("""
        var y;
        varexo e;
        model;
        y = e;
        end;
        perfect_foresight_setup;
        perfect_foresight_solver;
        stoch_simul;
        """)
    assert "cannot mix perfect foresight" in str(exc_info.value)
    assert issubclass(dp.StatementException, dp.DynareException)
