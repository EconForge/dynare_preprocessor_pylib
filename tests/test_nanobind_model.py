import json
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


def test_source_location_class():
    """Test SourceLocation properties, equality, and formatting."""
    loc = dp.SourceLocation("test.mod", 3, 5, 3, 10)
    assert loc.filename == "test.mod"
    assert loc.begin_line == 3
    assert loc.begin_column == 5
    assert loc.end_line == 3
    assert loc.end_column == 10
    assert loc.line == 3
    assert loc.column == 5
    assert loc.col == 5
    assert "test.mod: line 3, cols 5-9" in str(loc)
    assert repr(loc) == "<SourceLocation test.mod: line 3, cols 5-9>"

    loc2 = dp.SourceLocation("test.mod", 3, 5, 3, 10)
    assert loc == loc2
    assert loc != "not_a_loc"
    assert loc != None


def test_exceptions():
    """Test structured exception throwing, Python hierarchy, and location info."""
    # Syntax error -> ParserException with location
    with pytest.raises(dp.ParserException) as exc_info:
        dp.DynareModel("var x; model; x = ; end;")
    exc = exc_info.value
    assert "syntax error" in str(exc).lower()
    assert issubclass(dp.ParserException, dp.DynareException)
    assert issubclass(dp.ParserException, dp.SourceFileException)

    # Check location metadata on exception
    assert exc.location is not None
    assert isinstance(exc.location, dp.SourceLocation)
    assert exc.line == 1
    assert exc.column == 19 or exc.column == 20
    assert exc.filename == "in_memory.mod"
    assert "syntax error" in exc.message.lower()

    # Macro error -> MacroException with location
    with pytest.raises(dp.MacroException) as exc_info:
        dp.DynareModel("@#error \"explicit macro failure\"\nvar y;\nmodel;\ny=0;\nend;")
    m_exc = exc_info.value
    assert issubclass(dp.MacroException, dp.SourceFileException)
    assert issubclass(dp.MacroException, dp.DynareException)
    assert m_exc.location is not None
    assert m_exc.line == 1
    assert "explicit macro failure" in m_exc.message
    assert len(m_exc.backtrace) > 0

    # Equation count mismatch -> ModelSemanticException
    with pytest.raises(dp.ModelSemanticException) as exc_info:
        dp.DynareModel("var x; model; x = 1; x = 2; end;")
    sem_exc = exc_info.value
    assert "2 equations but 1 endogenous" in str(sem_exc)
    assert issubclass(dp.ModelSemanticException, dp.DynareException)
    assert issubclass(dp.EquationException, dp.ModelSemanticException)
    assert sem_exc.message is not None

    # Undeclared variable -> ParserException with location and undeclared_variables
    with pytest.raises(dp.ParserException) as exc_info:
        dp.DynareModel("var x;\nmodel;\nx = undefined_var;\nend;", strict=True)
    undec_exc = exc_info.value
    assert issubclass(dp.ParserException, dp.SourceFileException)
    assert undec_exc.location is not None
    assert undec_exc.line == 3
    assert undec_exc.column == 5
    assert len(undec_exc.undeclared_variables) == 1
    assert undec_exc.undeclared_variables[0][0] == "undefined_var"

    # Semantic equation error with equation number, line, and tag -> ModelSemanticException
    with pytest.raises(dp.ModelSemanticException) as exc_info:
        dp.DynareModel("""
        var y;
        varexo e;
        model;
        [name='eq_var']
        y = y(+1) + e;
        end;
        var_model(model_name=my_var, eqtags=['eq_var']);
        """)
    var_exc = exc_info.value
    assert var_exc.equation_number == 1
    assert var_exc.equation_tag == "eq_var"
    assert var_exc.equation_lineno == 5
    assert var_exc.line == 5
    assert var_exc.location is not None
    assert var_exc.location.line == 5
    assert "leaded endogenous variables on the RHS" in var_exc.message

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
    stmt_exc = exc_info.value
    assert "cannot mix perfect foresight" in str(stmt_exc)
    assert issubclass(dp.StatementException, dp.DynareException)
    assert stmt_exc.statement_name == "model"


def test_matched_irfs_json_serialization():
    """Test JSON serialization of matched_irfs and matched_irfs_weights statements."""
    mod_text = """
    var ghat;
    varexo eps_g;
    parameters a;
    a = 0.5;
    model;
    ghat = a * ghat(-1) + eps_g;
    end;
    matched_irfs;
    var ghat; varexo eps_g; periods 2:5; values 1.0; weights 1.0;
    end;
    matched_irfs_weights;
    ghat(2), eps_g, ghat(2), eps_g, 1;
    end;
    """
    model = dp.DynareModel(mod_text)
    data = json.loads(model.json_string)
    assert isinstance(data, dict)
    assert "transformed_modfile" in data
    statements = data["transformed_modfile"]["statements"]
    stmt_names = [stmt.get("statementName") for stmt in statements if isinstance(stmt, dict)]
    assert "matched_irfs" in stmt_names
    assert "matched_irfs_weights" in stmt_names

    matched_irfs_stmt = next(s for s in statements if s.get("statementName") == "matched_irfs")
    assert matched_irfs_stmt["contents"][0]["var"] == "ghat"
    assert matched_irfs_stmt["contents"][0]["varexo"] == "eps_g"
    pvw = matched_irfs_stmt["contents"][0]["periods_values_weights"][0]
    assert pvw["period1"] == 2
    assert pvw["period2"] == 5

    matched_irfs_weights_stmt = next(s for s in statements if s.get("statementName") == "matched_irfs_weights")
    item = matched_irfs_weights_stmt["contents"][0]
    assert item["endo1"] == "ghat"
    assert item["exo1"] == "eps_g"
    assert item["endo2"] == "ghat"
    assert item["exo2"] == "eps_g"


def test_perfect_foresight_models():
    """Test loading deterministic models with perfect_foresight_setup, perfect_foresight_solver, and simul."""
    base_txt = """
    var c;
    varexo x;
    parameters a;
    a = 0.5;
    model;
    c = a * c(-1) + x;
    end;
    initval;
    c = 0;
    x = 0;
    end;
    """
    m1 = dp.DynareModel(base_txt + "perfect_foresight_setup(periods=10);")
    assert len(m1.equations) == 1

    m2 = dp.DynareModel(base_txt + "perfect_foresight_setup(periods=10); perfect_foresight_solver;")
    assert len(m2.equations) == 1

    m3 = dp.DynareModel(base_txt + "simul(periods=10);")
    assert len(m3.equations) == 1


def test_deterministic_shocks_trajectories():
    """Test that deterministic shocks in shocks blocks populate model.trajectories."""
    mod_text = """
    var c;
    varexo x;
    parameters a;
    a = 0.5;
    model;
    c = a * c(-1) + x;
    end;
    initval;
    c = 0;
    x = 0;
    end;
    shocks;
    var x;
    periods 1, 3:5;
    values 1.2, 2.5;
    end;
    """
    model = dp.DynareModel(mod_text)
    assert "x" in model.trajectories
    assert model.trajectories["x"] == [(1, 1, 1.2), (3, 5, 2.5)]


def test_ramst_deterministic_model():
    """Test loading ramst.mod (deterministic setup with shocks and solver)."""
    modfile = os.path.join(os.path.dirname(__file__), "modfiles", "ramst.mod")
    model = dp.DynareModel(modfile)
    assert "c" in model.endogenous
    assert "k" in model.endogenous
    assert "x" in model.exogenous
    assert "x" in model.trajectories
    assert model.trajectories["x"] == [(1, 1, 1.2)]


def test_stochastic_parameter_override():
    """Test explicit stochastic constructor parameter."""
    txt = """
    var c;
    varexo x;
    parameters a;
    a = 0.5;
    model;
    c = a * c(-1) + x;
    end;
    initval;
    c = 0;
    x = 0;
    end;
    perfect_foresight_setup(periods=10);
    perfect_foresight_solver;
    """
    # Auto-detected as deterministic
    m = dp.DynareModel(txt)
    assert len(m.equations) == 1

    # Explicit stochastic=False succeeds
    m_det = dp.DynareModel(txt, stochastic=False)
    assert len(m_det.equations) == 1

    # Explicit stochastic=True fails with StatementException because mixing stochastic context with solver
    with pytest.raises(dp.StatementException):
        dp.DynareModel(txt, stochastic=True)



