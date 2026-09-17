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

} // namespace

NB_MODULE(dynare_preprocessor, m)
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
      .def(nb::init<const std::string&, int, int>(),
           nb::arg("modfile_content_or_path"),
           nb::arg("derivs_order") = 1,
           nb::arg("params_derivs_order") = 0)
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

