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

#ifndef DYNLIB_HH
#define DYNLIB_HH

#ifdef COMPILER
# undef COMPILER
#endif

#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "CommonEnums.hh"
#include "Exceptions.hh"
#include "ExprNode.hh"
#include "ModFile.hh"
#include "ParsingDriver.hh"
#include "WarningConsolidation.hh"

using EvaluationException = MathEvaluationException;

enum class ExprNodeType
{
  NumConstNode,
  VariableNode,
  UnaryOpNode,
  BinaryOpNode,
  TrinaryOpNode
};

ExprNodeType expression_type(expr_t expression);

// (equation id, type specific id) -> symbolic partial derivative
using symb_jacobian_t = std::map<std::pair<int, int>, expr_t>;

// (equation id, type specific id) -> evaluated partial derivative
using jacobian_t = std::map<std::pair<int, int>, double>;

//! Matrix in COO format: vector of ((eq, var1, var2, ...), value)
using coo_matrix_t = std::vector<std::pair<std::vector<int>, double>>;

//! Vector of COO matrices per derivation order
using derivatives_t = std::vector<coo_matrix_t>;

//! Structured dynamic Jacobian blocks
struct JacobianBlocks
{
  int n_eq {0};
  int n_endo {0};
  int n_exo {0};
  int n_exo_det {0};
  int n_params {0};

  // Dense row-major storage: (eq * n_cols + col)
  std::vector<double> lead;     // n_eq x n_endo (df / dy_{t+1})
  std::vector<double> curr;     // n_eq x n_endo (df / dy_t)
  std::vector<double> lag;      // n_eq x n_endo (df / dy_{t-1})
  std::vector<double> exo;      // n_eq x n_exo (df / d_eps)
  std::vector<double> exo_det;  // n_eq x n_exo_det
  std::vector<double> params;   // n_eq x n_params
};

class DynareModel
{
public:
  DynareModel(const std::string& modfile_content_or_path,
              int derivs_order = 1,
              int params_derivs_order = 0);

  // Introspection / symbols
  std::vector<std::string> endogenous;
  std::vector<std::string> exogenous;
  std::vector<std::string> exogenous_det;
  std::vector<std::string> parameters;
  std::vector<std::string> equations;
  std::map<std::string, double> context;
  std::map<std::pair<std::string, std::string>, double> covariances;
  std::map<std::string, std::vector<std::tuple<int, int, double>>> trajectories;
  std::vector<std::tuple<SymbolType, int, int>> symbol_info;
  std::string json_string;

  // Structural incidence
  int max_endo_lag {0};
  int max_endo_lead {0};
  std::vector<std::vector<int>> lead_lag_incidence;

  // Evaluators
  [[nodiscard]] std::vector<double> residuals(
      std::span<const double> endo_future,
      std::span<const double> endo_present,
      std::span<const double> endo_past,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

  [[nodiscard]] std::vector<jacobian_t> jacobians(
      std::span<const double> endo_future,
      std::span<const double> endo_present,
      std::span<const double> endo_past,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

  [[nodiscard]] JacobianBlocks jacobian_blocks(
      std::span<const double> endo_future,
      std::span<const double> endo_present,
      std::span<const double> endo_past,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

  // Returns (dense_data, n_eq, n_cols)
  [[nodiscard]] std::tuple<std::vector<double>, int, int> jacobian(
      std::span<const double> endo_future,
      std::span<const double> endo_present,
      std::span<const double> endo_past,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

  [[nodiscard]] std::vector<double> static_residuals(
      std::span<const double> endo,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

  [[nodiscard]] std::tuple<std::vector<double>, int, int> static_jacobian(
      std::span<const double> endo,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

  [[nodiscard]] derivatives_t derivatives(
      std::span<const double> endo_future,
      std::span<const double> endo_present,
      std::span<const double> endo_past,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

private:
  std::unique_ptr<WarningConsolidation> warnings;
  std::unique_ptr<ParsingDriver> driver;
  std::unique_ptr<ModFile> mod_file;

  void set_mod_file(const std::string& modfile_content_or_path, int derivs_order, int params_derivs_order);
  void set_json_string();
  void set_symbols();
  void set_equations();
  void set_context();
  void set_exogenous();
  void set_lead_lag_incidence();
  void set_symbolic_derivatives();

  double evaluate_with_lags(
      expr_t expression,
      const std::vector<std::span<const double>>& endo,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;

  std::vector<symb_jacobian_t> symb_jacob_endo; // 0: past (lag -1), 1: present (lag 0), 2: future (lag +1)
  symb_jacobian_t symb_jacob_exo;
  symb_jacobian_t symb_jacob_exo_det;
  symb_jacobian_t symb_jacob_params;

  std::vector<std::map<std::vector<int>, expr_t>> symb_derivatives;

  void eval_symb_jacob(
      const symb_jacobian_t& symb_jacob,
      jacobian_t& out,
      const std::vector<std::span<const double>>& endo,
      std::span<const double> exo,
      std::span<const double> exo_det,
      std::span<const double> params) const;
};

#endif

