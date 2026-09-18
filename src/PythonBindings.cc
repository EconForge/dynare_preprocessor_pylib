/*
 * Copyright © 2026 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#ifdef COMPILER
# undef COMPILER
#endif

#include "DynLib.hh"
#include "Exceptions.hh"

namespace nb = nanobind;

namespace
{

struct InputSpan
{
  const double* data {nullptr};
  size_t size {0};
  std::vector<double> storage;

  InputSpan(nb::handle h)
  {
    if (nb::isinstance<nb::ndarray<>>(h))
      {
        nb::ndarray<const double, nb::c_contig, nb::device::cpu> arr
            = nb::cast<nb::ndarray<const double, nb::c_contig, nb::device::cpu>>(h);
        data = arr.data();
        size = arr.size();
      }
    else if (nb::isinstance<nb::sequence>(h))
      {
        nb::sequence seq = nb::cast<nb::sequence>(h);
        size_t n = nb::len(seq);
        storage.reserve(n);
        for (size_t i = 0; i < n; i++)
          storage.push_back(nb::cast<double>(seq[i]));
        data = storage.data();
        size = storage.size();
      }
    else
      throw nb::type_error("Expected a numpy array or a sequence of floats");
  }

  [[nodiscard]] std::span<const double>
  span() const
  {
    return {data, size};
  }
};

nb::object
to_numpy_1d(const std::vector<double>& vec)
{
  size_t n = vec.size();
  auto* ptr = new double[n];
  std::copy(vec.begin(), vec.end(), ptr);
  nb::capsule owner(ptr, [](void* p) noexcept { delete[] static_cast<double*>(p); });
  size_t shape[1] = {n};
  return nb::ndarray<nb::numpy, double, nb::ndim<1>>(ptr, 1, shape, owner).cast();
}

nb::object
to_numpy_2d(const std::vector<double>& vec, size_t rows, size_t cols)
{
  size_t total = rows * cols;
  auto* ptr = new double[total];
  std::copy(vec.begin(), vec.end(), ptr);
  nb::capsule owner(ptr, [](void* p) noexcept { delete[] static_cast<double*>(p); });
  size_t shape[2] = {rows, cols};
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(ptr, 2, shape, owner).cast();
}

struct ExceptionTypes
{
  PyObject* dynare {nullptr};
  PyObject* source {nullptr};
  PyObject* parser {nullptr};
  PyObject* macro {nullptr};
  PyObject* semantic {nullptr};
  PyObject* equation {nullptr};
  PyObject* statement {nullptr};
  PyObject* fileio {nullptr};
  PyObject* eval {nullptr};
  PyObject* internal {nullptr};
};

static ExceptionTypes g_exc_types;

void
set_source_location_attrs(nb::object& exc_obj,
                          const std::optional<SourceLocation>& loc,
                          const std::string& msg)
{
  if (loc)
    {
      nb::setattr(exc_obj, "location", nb::cast(*loc));
      nb::setattr(exc_obj, "filename", nb::cast(loc->filename));
      nb::setattr(exc_obj, "line", nb::cast(loc->begin_line));
      nb::setattr(exc_obj, "column", nb::cast(loc->begin_column));
      nb::setattr(exc_obj, "begin_line", nb::cast(loc->begin_line));
      nb::setattr(exc_obj, "begin_column", nb::cast(loc->begin_column));
      nb::setattr(exc_obj, "end_line", nb::cast(loc->end_line));
      nb::setattr(exc_obj, "end_column", nb::cast(loc->end_column));
    }
  else
    {
      nb::setattr(exc_obj, "location", nb::none());
      nb::setattr(exc_obj, "filename", nb::none());
      nb::setattr(exc_obj, "line", nb::none());
      nb::setattr(exc_obj, "column", nb::none());
      nb::setattr(exc_obj, "begin_line", nb::none());
      nb::setattr(exc_obj, "begin_column", nb::none());
      nb::setattr(exc_obj, "end_line", nb::none());
      nb::setattr(exc_obj, "end_column", nb::none());
    }
  nb::setattr(exc_obj, "message", nb::cast(msg));
}

void
set_semantic_attrs(nb::object& exc_obj, const ModelSemanticException& e)
{
  if (e.getEquationNumber())
    nb::setattr(exc_obj, "equation_number", nb::cast(*e.getEquationNumber()));
  else
    nb::setattr(exc_obj, "equation_number", nb::none());

  if (e.getEquationLineno())
    {
      int l = *e.getEquationLineno();
      nb::setattr(exc_obj, "equation_lineno", nb::cast(l));
      nb::setattr(exc_obj, "line", nb::cast(l));
      nb::setattr(exc_obj, "begin_line", nb::cast(l));
      nb::setattr(exc_obj, "end_line", nb::cast(l));
    }
  else
    {
      nb::setattr(exc_obj, "equation_lineno", nb::none());
      nb::setattr(exc_obj, "line", nb::none());
      nb::setattr(exc_obj, "begin_line", nb::none());
      nb::setattr(exc_obj, "end_line", nb::none());
    }

  if (e.getLocation())
    nb::setattr(exc_obj, "location", nb::cast(*e.getLocation()));
  else
    nb::setattr(exc_obj, "location", nb::none());

  nb::setattr(exc_obj, "column", nb::none());
  nb::setattr(exc_obj, "begin_column", nb::none());
  nb::setattr(exc_obj, "end_column", nb::none());
  nb::setattr(exc_obj, "filename", nb::none());

  if (e.getEquationTag())
    nb::setattr(exc_obj, "equation_tag", nb::cast(*e.getEquationTag()));
  else
    nb::setattr(exc_obj, "equation_tag", nb::none());

  if (e.getSymbolName())
    nb::setattr(exc_obj, "symbol_name", nb::cast(*e.getSymbolName()));
  else
    nb::setattr(exc_obj, "symbol_name", nb::none());

  nb::setattr(exc_obj, "message", nb::cast(e.getMessage()));
}

void
set_statement_attrs(nb::object& exc_obj, const StatementException& e)
{
  nb::setattr(exc_obj, "statement_name", nb::cast(e.getStatementName()));
  if (e.getOptionName())
    nb::setattr(exc_obj, "option_name", nb::cast(*e.getOptionName()));
  else
    nb::setattr(exc_obj, "option_name", nb::none());
  nb::setattr(exc_obj, "message", nb::cast(e.getMessage()));
  nb::setattr(exc_obj, "location", nb::none());
  nb::setattr(exc_obj, "line", nb::none());
  nb::setattr(exc_obj, "column", nb::none());
  nb::setattr(exc_obj, "filename", nb::none());
}

void
set_fileio_attrs(nb::object& exc_obj, const FileIOException& e)
{
  nb::setattr(exc_obj, "path", nb::cast(e.getPath().string()));
  nb::setattr(exc_obj, "action", nb::cast(e.getAction()));
  nb::setattr(exc_obj, "message", nb::cast(std::string(e.what())));
  nb::setattr(exc_obj, "location", nb::none());
  nb::setattr(exc_obj, "line", nb::none());
  nb::setattr(exc_obj, "column", nb::none());
  nb::setattr(exc_obj, "filename", nb::none());
}

void
set_dynare_attrs(nb::object& exc_obj, const DynareException& e)
{
  nb::setattr(exc_obj, "message", nb::cast(std::string(e.what())));
  nb::setattr(exc_obj, "location", nb::none());
  nb::setattr(exc_obj, "line", nb::none());
  nb::setattr(exc_obj, "column", nb::none());
  nb::setattr(exc_obj, "filename", nb::none());
}

template <typename ExcType, typename SetAttrsFunc>
void
raise_py_exception(PyObject* py_type, const ExcType& e, SetAttrsFunc set_attrs)
{
  nb::str msg = nb::str(e.what());
  PyObject* exc = PyObject_CallOneArg(py_type, msg.ptr());
  if (!exc)
    return;
  nb::object exc_obj = nb::steal(exc);
  set_attrs(exc_obj, e);
#if PY_VERSION_HEX >= 0x030C0000
  PyErr_SetRaisedException(exc_obj.release().ptr());
#else
  PyErr_SetObject(py_type, exc_obj.ptr());
#endif
}

} // namespace

NB_MODULE(_dynare_preprocessor, m)
{
  m.doc() = "Dynare Preprocessor Model API and Python Wrapper";

  // Exception hierarchy
  nb::exception<DynareException> ex_dynare(m, "DynareException");
  [[maybe_unused]] nb::exception<SourceFileException> ex_source(m, "SourceFileException", ex_dynare);
  [[maybe_unused]] nb::exception<ParserException> ex_parser(m, "ParserException", ex_source);
  [[maybe_unused]] nb::exception<MacroException> ex_macro(m, "MacroException", ex_source);
  nb::exception<ModelSemanticException> ex_semantic(m, "ModelSemanticException", ex_dynare);
  [[maybe_unused]] nb::exception<EquationException> ex_equation(m, "EquationException", ex_semantic);
  [[maybe_unused]] nb::exception<StatementException> ex_statement(m, "StatementException", ex_dynare);
  [[maybe_unused]] nb::exception<FileIOException> ex_fileio(m, "FileIOException", ex_dynare);
  nb::exception<EvaluationException> ex_eval(m, "EvaluationException", ex_dynare);
  [[maybe_unused]] nb::exception<InternalCompilerException> ex_internal(m, "InternalCompilerException", ex_dynare);

  // Aliases for compatibility with dynare_lite
  m.attr("PreprocessorException") = ex_dynare;
  m.attr("UnsupportedFeatureException") = ex_semantic;
  m.attr("EvalException") = ex_eval;

  // SourceLocation class
  nb::class_<SourceLocation>(m, "SourceLocation")
      .def(nb::init<>())
      .def(nb::init<string, int, int, int, int>(),
           nb::arg("filename") = "",
           nb::arg("begin_line") = 1,
           nb::arg("begin_column") = 1,
           nb::arg("end_line") = 1,
           nb::arg("end_column") = 1)
      .def_rw("filename", &SourceLocation::filename)
      .def_rw("begin_line", &SourceLocation::begin_line)
      .def_rw("begin_column", &SourceLocation::begin_column)
      .def_rw("end_line", &SourceLocation::end_line)
      .def_rw("end_column", &SourceLocation::end_column)
      .def_prop_ro("line", [](const SourceLocation& sl) { return sl.begin_line; })
      .def_prop_ro("column", [](const SourceLocation& sl) { return sl.begin_column; })
      .def_prop_ro("col", [](const SourceLocation& sl) { return sl.begin_column; })
      .def("__repr__",
           [](const SourceLocation& sl) {
             return "<SourceLocation " + sl.to_string() + ">";
           })
      .def("__str__", &SourceLocation::to_string)
      .def("__eq__",
           [](const SourceLocation& a, nb::handle other) {
             if (!other.is_valid() || other.is_none() || !nb::isinstance<SourceLocation>(other))
               return false;
             const auto& b = nb::cast<const SourceLocation&>(other);
             return a.filename == b.filename && a.begin_line == b.begin_line
                    && a.begin_column == b.begin_column && a.end_line == b.end_line
                    && a.end_column == b.end_column;
           },
           nb::arg("other").none());

  // Register custom translator for detailed exception attributes
  g_exc_types.dynare = ex_dynare.ptr();
  g_exc_types.source = ex_source.ptr();
  g_exc_types.parser = ex_parser.ptr();
  g_exc_types.macro = ex_macro.ptr();
  g_exc_types.semantic = ex_semantic.ptr();
  g_exc_types.equation = ex_equation.ptr();
  g_exc_types.statement = ex_statement.ptr();
  g_exc_types.fileio = ex_fileio.ptr();
  g_exc_types.eval = ex_eval.ptr();
  g_exc_types.internal = ex_internal.ptr();

  nb::register_exception_translator(
      [](const std::exception_ptr& p, void* payload) {
        auto* t = static_cast<ExceptionTypes*>(payload);
        try
          {
            std::rethrow_exception(p);
          }
        catch (const ParserException& e)
          {
            raise_py_exception(t->parser, e, [](nb::object& o, const ParserException& ex) {
              set_source_location_attrs(o, ex.getLocation(), ex.getMessage());
              nb::setattr(o, "undeclared_variables", nb::cast(ex.getUndeclaredVariables()));
            });
          }
        catch (const MacroException& e)
          {
            raise_py_exception(t->macro, e, [](nb::object& o, const MacroException& ex) {
              set_source_location_attrs(o, ex.getLocation(), ex.getMessage());
              nb::setattr(o, "backtrace", nb::cast(ex.getBacktrace()));
            });
          }
        catch (const SourceFileException& e)
          {
            raise_py_exception(t->source, e, [](nb::object& o, const SourceFileException& ex) {
              set_source_location_attrs(o, ex.getLocation(), ex.getMessage());
            });
          }
        catch (const EquationException& e)
          {
            raise_py_exception(t->equation, e, [](nb::object& o, const EquationException& ex) {
              set_semantic_attrs(o, ex);
            });
          }
        catch (const ModelSemanticException& e)
          {
            raise_py_exception(t->semantic, e, [](nb::object& o, const ModelSemanticException& ex) {
              set_semantic_attrs(o, ex);
            });
          }
        catch (const StatementException& e)
          {
            raise_py_exception(t->statement, e, [](nb::object& o, const StatementException& ex) {
              set_statement_attrs(o, ex);
            });
          }
        catch (const FileIOException& e)
          {
            raise_py_exception(t->fileio, e, [](nb::object& o, const FileIOException& ex) {
              set_fileio_attrs(o, ex);
            });
          }
        catch (const EvaluationException& e)
          {
            raise_py_exception(t->eval, e, [](nb::object& o, const EvaluationException& ex) {
              set_dynare_attrs(o, ex);
            });
          }
        catch (const InternalCompilerException& e)
          {
            raise_py_exception(t->internal, e, [](nb::object& o, const InternalCompilerException& ex) {
              set_dynare_attrs(o, ex);
            });
          }
        catch (const DynareException& e)
          {
            raise_py_exception(t->dynare, e, [](nb::object& o, const DynareException& ex) {
              set_dynare_attrs(o, ex);
            });
          }
      },
      &g_exc_types);

  // SymbolType enum
  nb::enum_<SymbolType>(m, "SymbolType")
      .value("endogenous", SymbolType::endogenous)
      .value("exogenous", SymbolType::exogenous)
      .value("exogenous_det", SymbolType::exogenousDet)
      .value("parameter", SymbolType::parameter)
      .export_values();

  // JacobianBlocks class
  nb::class_<JacobianBlocks>(m, "JacobianBlocks")
      .def_prop_ro("lead",
                   [](const JacobianBlocks& b) {
                     return to_numpy_2d(b.lead, b.n_eq, b.n_endo);
                   })
      .def_prop_ro("curr",
                   [](const JacobianBlocks& b) {
                     return to_numpy_2d(b.curr, b.n_eq, b.n_endo);
                   })
      .def_prop_ro("lag",
                   [](const JacobianBlocks& b) {
                     return to_numpy_2d(b.lag, b.n_eq, b.n_endo);
                   })
      .def_prop_ro("exo",
                   [](const JacobianBlocks& b) {
                     return to_numpy_2d(b.exo, b.n_eq, b.n_exo);
                   })
      .def_prop_ro("exo_det",
                   [](const JacobianBlocks& b) {
                     return to_numpy_2d(b.exo_det, b.n_eq, b.n_exo_det);
                   })
      .def_prop_ro("params",
                   [](const JacobianBlocks& b) {
                     return to_numpy_2d(b.params, b.n_eq, b.n_params);
                   })
      .def("__len__", [](const JacobianBlocks&) { return 6; })
      .def("__getitem__",
           [](const JacobianBlocks& b, int idx) -> nb::object {
             if (idx < 0)
               idx += 6;
             switch (idx)
               {
               case 0:
                 return to_numpy_2d(b.lead, b.n_eq, b.n_endo);
               case 1:
                 return to_numpy_2d(b.curr, b.n_eq, b.n_endo);
               case 2:
                 return to_numpy_2d(b.lag, b.n_eq, b.n_endo);
               case 3:
                 return to_numpy_2d(b.exo, b.n_eq, b.n_exo);
               case 4:
                 return to_numpy_2d(b.exo_det, b.n_eq, b.n_exo_det);
               case 5:
                 return to_numpy_2d(b.params, b.n_eq, b.n_params);
               default:
                 throw nb::index_error("Index out of range for JacobianBlocks (must be 0-5)");
               }
           })
      .def("__iter__", [](const JacobianBlocks& b) {
        nb::list lst;
        lst.append(to_numpy_2d(b.lead, b.n_eq, b.n_endo));
        lst.append(to_numpy_2d(b.curr, b.n_eq, b.n_endo));
        lst.append(to_numpy_2d(b.lag, b.n_eq, b.n_endo));
        lst.append(to_numpy_2d(b.exo, b.n_eq, b.n_exo));
        lst.append(to_numpy_2d(b.exo_det, b.n_eq, b.n_exo_det));
        lst.append(to_numpy_2d(b.params, b.n_eq, b.n_params));
        return nb::iter(lst);
      });

  // DynareModel class
  nb::class_<DynareModel>(m, "DynareModel")
      .def(nb::init<const std::string&, int, int, bool>(),
           nb::arg("modfile_content_or_path"),
           nb::arg("derivs_order") = 1,
           nb::arg("params_derivs_order") = 0,
           nb::arg("strict") = false)
      .def_ro("endogenous", &DynareModel::endogenous)
      .def_ro("exogenous", &DynareModel::exogenous)
      .def_ro("exogenous_det", &DynareModel::exogenous_det)
      .def_ro("parameters", &DynareModel::parameters)
      .def_ro("equations", &DynareModel::equations)
      .def_ro("context", &DynareModel::context)
      .def_ro("covariances", &DynareModel::covariances)
      .def_ro("trajectories", &DynareModel::trajectories)
      .def_ro("symbol_info", &DynareModel::symbol_info)
      .def_ro("json_string", &DynareModel::json_string)
      .def_ro("lead_lag_incidence", &DynareModel::lead_lag_incidence)
      .def_ro("max_endo_lag", &DynareModel::max_endo_lag)
      .def_ro("max_endo_lead", &DynareModel::max_endo_lead)
      .def("residuals",
           [](const DynareModel& model,
              nb::handle endo_future,
              nb::handle endo_present,
              nb::handle endo_past,
              nb::handle exo,
              nb::handle exo_det,
              nb::handle params) {
             InputSpan s_fut(endo_future);
             InputSpan s_pres(endo_present);
             InputSpan s_past(endo_past);
             InputSpan s_exo(exo);
             InputSpan s_exo_det(exo_det);
             InputSpan s_params(params);
             return to_numpy_1d(model.residuals(
                 s_fut.span(), s_pres.span(), s_past.span(),
                 s_exo.span(), s_exo_det.span(), s_params.span()));
           },
           nb::arg("endo_future"), nb::arg("endo_present"), nb::arg("endo_past"),
           nb::arg("exo"), nb::arg("exo_det"), nb::arg("params"),
           "Evaluate dynamic residuals F(y_{t+1}, y_t, y_{t-1}, e_t, params) as a 1D NumPy array")
      .def("jacobians",
           [](const DynareModel& model,
              nb::handle endo_future,
              nb::handle endo_present,
              nb::handle endo_past,
              nb::handle exo,
              nb::handle exo_det,
              nb::handle params) {
             InputSpan s_fut(endo_future);
             InputSpan s_pres(endo_present);
             InputSpan s_past(endo_past);
             InputSpan s_exo(exo);
             InputSpan s_exo_det(exo_det);
             InputSpan s_params(params);
             return model.jacobians(
                 s_fut.span(), s_pres.span(), s_past.span(),
                 s_exo.span(), s_exo_det.span(), s_params.span());
           },
           nb::arg("endo_future"), nb::arg("endo_present"), nb::arg("endo_past"),
           nb::arg("exo"), nb::arg("exo_det"), nb::arg("params"),
           "Return sparse maps for jacobian blocks (compatibility with dynare_lite)")
      .def("jacobian_blocks",
           [](const DynareModel& model,
              nb::handle endo_future,
              nb::handle endo_present,
              nb::handle endo_past,
              nb::handle exo,
              nb::handle exo_det,
              nb::handle params) {
             InputSpan s_fut(endo_future);
             InputSpan s_pres(endo_present);
             InputSpan s_past(endo_past);
             InputSpan s_exo(exo);
             InputSpan s_exo_det(exo_det);
             InputSpan s_params(params);
             return model.jacobian_blocks(
                 s_fut.span(), s_pres.span(), s_past.span(),
                 s_exo.span(), s_exo_det.span(), s_params.span());
           },
           nb::arg("endo_future"), nb::arg("endo_present"), nb::arg("endo_past"),
           nb::arg("exo"), nb::arg("exo_det"), nb::arg("params"),
           "Return structured JacobianBlocks with 2D NumPy arrays: lead, curr, lag, exo, exo_det, params")
      .def("jacobian",
           [](const DynareModel& model,
              nb::handle endo_future,
              nb::handle endo_present,
              nb::handle endo_past,
              nb::handle exo,
              nb::handle exo_det,
              nb::handle params) {
             InputSpan s_fut(endo_future);
             InputSpan s_pres(endo_present);
             InputSpan s_past(endo_past);
             InputSpan s_exo(exo);
             InputSpan s_exo_det(exo_det);
             InputSpan s_params(params);
             auto [mat, n_eq, n_cols] = model.jacobian(
                 s_fut.span(), s_pres.span(), s_past.span(),
                 s_exo.span(), s_exo_det.span(), s_params.span());
             return to_numpy_2d(mat, n_eq, n_cols);
           },
           nb::arg("endo_future"), nb::arg("endo_present"), nb::arg("endo_past"),
           nb::arg("exo"), nb::arg("exo_det"), nb::arg("params"),
           "Evaluate the combined dynamic Jacobian as a 2D NumPy array")
      .def("static_residuals",
           [](const DynareModel& model,
              nb::handle endo,
              nb::handle exo,
              nb::handle exo_det,
              nb::handle params) {
             InputSpan s_endo(endo);
             InputSpan s_exo(exo);
             InputSpan s_exo_det(exo_det);
             InputSpan s_params(params);
             return to_numpy_1d(model.static_residuals(
                 s_endo.span(), s_exo.span(), s_exo_det.span(), s_params.span()));
           },
           nb::arg("endo"), nb::arg("exo"), nb::arg("exo_det"), nb::arg("params"),
           "Evaluate static (steady-state) residuals as a 1D NumPy array")
      .def("static_jacobian",
           [](const DynareModel& model,
              nb::handle endo,
              nb::handle exo,
              nb::handle exo_det,
              nb::handle params) {
             InputSpan s_endo(endo);
             InputSpan s_exo(exo);
             InputSpan s_exo_det(exo_det);
             InputSpan s_params(params);
             auto [mat, n_eq, n_endo] = model.static_jacobian(
                 s_endo.span(), s_exo.span(), s_exo_det.span(), s_params.span());
             return to_numpy_2d(mat, n_eq, n_endo);
           },
           nb::arg("endo"), nb::arg("exo"), nb::arg("exo_det"), nb::arg("params"),
           "Evaluate static (steady-state) Jacobian as a 2D NumPy array")
      .def("derivatives",
           [](const DynareModel& model,
              nb::handle endo_future,
              nb::handle endo_present,
              nb::handle endo_past,
              nb::handle exo,
              nb::handle exo_det,
              nb::handle params) {
             InputSpan s_fut(endo_future);
             InputSpan s_pres(endo_present);
             InputSpan s_past(endo_past);
             InputSpan s_exo(exo);
             InputSpan s_exo_det(exo_det);
             InputSpan s_params(params);
             return model.derivatives(
                 s_fut.span(), s_pres.span(), s_past.span(),
                 s_exo.span(), s_exo_det.span(), s_params.span());
           },
           nb::arg("endo_future"), nb::arg("endo_present"), nb::arg("endo_past"),
           nb::arg("exo"), nb::arg("exo_det"), nb::arg("params"),
           "Evaluate higher-order derivatives in COO format");
}

