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

#include "DynLib.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <typeinfo>

#include <boost/algorithm/string.hpp>

#include "Shocks.hh"

// Prototype for macro-expansion function from MacroExpandModFile.cc
stringstream macroExpandModFile(const filesystem::path& filename, const istream& modfile,
                                bool debug, bool save_macro, filesystem::path save_macro_file,
                                bool line_macro, const vector<pair<string, string>>& defines,
                                vector<filesystem::path> paths);

ExprNodeType
expression_type(expr_t expression)
{
  const std::type_info& type = typeid(*expression);
  if (type == typeid(NumConstNode))
    return ExprNodeType::NumConstNode;
  else if (type == typeid(VariableNode))
    return ExprNodeType::VariableNode;
  else if (type == typeid(UnaryOpNode))
    return ExprNodeType::UnaryOpNode;
  else if (type == typeid(BinaryOpNode))
    return ExprNodeType::BinaryOpNode;
  else if (type == typeid(TrinaryOpNode))
    return ExprNodeType::TrinaryOpNode;
  else
    throw EvaluationException("Unknown expression node type during evaluation");
  __builtin_unreachable();
}

DynareModel::DynareModel(const std::string& modfile_content_or_path,
                         int derivs_order,
                         int params_derivs_order)
{
  set_mod_file(modfile_content_or_path, derivs_order, params_derivs_order);
  set_json_string();
  set_symbols();
  set_equations();
  set_context();
  set_exogenous();
  set_lead_lag_incidence();
  set_symbolic_derivatives();
}

void
DynareModel::set_mod_file(const std::string& modfile_content_or_path,
                          int derivs_order,
                          int params_derivs_order)
{
  filesystem::path filepath;
  string content;
  vector<filesystem::path> paths;

  // Check if input is an existing file path
  std::error_code ec;
  if (filesystem::exists(modfile_content_or_path, ec))
    {
      filepath = filesystem::absolute(modfile_content_or_path);
      ifstream f(filepath);
      if (!f.is_open())
        throw SourceFileException("Cannot open file: " + modfile_content_or_path);
      stringstream ss;
      ss << f.rdbuf();
      content = ss.str();
      paths.push_back(filepath.parent_path());
    }
  else
    {
      filepath = "in_memory.mod";
      content = modfile_content_or_path;
    }

  // 1. Macro-expand MOD file
  istringstream modfile_stream(content);
  stringstream macro_output
      = macroExpandModFile(filepath, modfile_stream, false, false, "", false, {}, move(paths));

  // 2. Parse into AST
  const bool nostrict = true;
  const bool nowarn = true;
  warnings = make_unique<WarningConsolidation>(nowarn);
  driver = make_unique<ParsingDriver>(*warnings, nostrict);
  const bool debug = false;
  mod_file = driver->parse(macro_output, debug);

  // 3. Checking pass
  const bool stochastic = true;
  mod_file->checkPass(nostrict, stochastic);

  // 4. Transformation pass
  const bool compute_xrefs = false;
  const bool transform_unary_ops = false;
  string exclude_eqs = "";
  string include_eqs = "";
  mod_file->transformPass(nostrict, stochastic, compute_xrefs,
                          transform_unary_ops, exclude_eqs, include_eqs);

  // 5. Evaluate parameters initialization, initval, and steady states
  const bool warn_uninit = false;
  mod_file->evalAllExpressions(warn_uninit);

  // 6. Computing pass for derivatives
  const bool no_tmp_terms = true;
  mod_file->mod_file_struct.order_option = derivs_order;

  mod_file->static_model = static_cast<StaticModel>(mod_file->dynamic_model);
  mod_file->static_model.computingPass(
      derivs_order,
      params_derivs_order,
      mod_file->global_eval_context,
      no_tmp_terms,
      false,
      false
  );
  mod_file->dynamic_model.computingPass(
      derivs_order,
      params_derivs_order,
      mod_file->global_eval_context,
      no_tmp_terms,
      false,
      false
  );
}

void
DynareModel::set_json_string()
{
  const string basename = "model";
  auto outputpoint = JsonOutputPointType::computingpass;
  auto json_output_mode = JsonFileOutputType::standardout;
  bool onlyjson = false;
  bool jsonderivsimple = true;

  // RAII stream buffer redirection
  stringstream buffer;
  streambuf* old_buf = cout.rdbuf(buffer.rdbuf());
  struct StreamGuard
  {
    streambuf* old;
    ~StreamGuard() { cout.rdbuf(old); }
  } guard {old_buf};

  mod_file->writeJsonOutput(basename, outputpoint, json_output_mode, onlyjson, jsonderivsimple);

  json_string = buffer.str();
  boost::replace_all(json_string, ", ,", ",");
  string prefix = "//-- BEGIN JSON --// \n";
  string suffix = "\n//-- END JSON --// \nJSON written after Computing step.\n";
  if (json_string.ends_with(suffix))
    json_string.erase(json_string.length() - suffix.length(), suffix.length());
  if (json_string.starts_with(prefix))
    json_string.erase(0, prefix.length());
}

void
DynareModel::set_symbols()
{
  const SymbolTable& table = mod_file->symbol_table;
  endogenous.clear();
  for (int id : table.endo_ids)
    endogenous.push_back(table.getName(id));

  exogenous.clear();
  for (int id : table.exo_ids)
    exogenous.push_back(table.getName(id));

  exogenous_det.clear();
  for (int id : table.exo_det_ids)
    exogenous_det.push_back(table.getName(id));

  parameters.clear();
  for (int id : table.param_ids)
    parameters.push_back(table.getName(id));

  symbol_info.clear();
  const DynamicModel& dm = mod_file->dynamic_model;
  for (const auto& [id, lag] : dm.inv_deriv_id_table)
    {
      SymbolType type = table.getType(id);
      int tsid = table.getTypeSpecificID(id);
      symbol_info.emplace_back(type, tsid, lag);
    }
}

void
DynareModel::set_equations()
{
  const DynamicModel& dm = mod_file->dynamic_model;
  equations.clear();
  for (expr_t eq : dm.equations)
    equations.push_back(eq->toString());
}

void
DynareModel::set_context()
{
  const SymbolTable& table = mod_file->symbol_table;
  eval_context_t new_context = mod_file->global_eval_context;

  // Add steady state to eval context for uninitialized variables
  for (const auto& [vect, expr] : mod_file->steady_state_model.def_table)
    {
      try
        {
          double val = expr->eval(new_context);
          for (int id : vect)
            if (new_context[id] == 0)
              new_context[id] = val;
        }
      catch (const ExprNode::EvalExternalFunctionException&)
        {
          throw ModelSemanticException("External functions are not supported yet in steady state evaluation.");
        }
      catch (const ExprNode::EvalException&)
        {
          throw EvaluationException("Evaluation error in steady state");
        }
    }

  context.clear();
  for (const auto& [id, val] : new_context)
    context[table.getName(id)] = val;
}

void
DynareModel::set_exogenous()
{
  const SymbolTable& table = mod_file->symbol_table;
  const eval_context_t& eval_ctx = mod_file->global_eval_context;
  covariances.clear();
  trajectories.clear();

  try
    {
      for (const auto& statement : mod_file->statements)
        {
          const auto& type = typeid(*statement);
          if (type == typeid(ShocksStatement))
            {
              auto* shock = static_cast<ShocksStatement*>(statement.get());
              for (const auto& [id, expr] : shock->var_shocks)
                {
                  string s = table.getName(id);
                  covariances[{s, s}] = expr->eval(eval_ctx);
                }
              for (const auto& [id, expr] : shock->std_shocks)
                {
                  string s = table.getName(id);
                  covariances[{s, s}] = std::pow(expr->eval(eval_ctx), 2);
                }
              for (const auto& [key, expr] : shock->covar_shocks)
                {
                  const auto& [id1, id2] = key;
                  string s1 = table.getName(id1);
                  string s2 = table.getName(id2);
                  double covar = expr->eval(eval_ctx);
                  covariances[{s1, s2}] = covar;
                }
              for (const auto& [key, expr] : shock->corr_shocks)
                {
                  const auto& [id1, id2] = key;
                  string s1 = table.getName(id1);
                  string s2 = table.getName(id2);
                  double corr = expr->eval(eval_ctx);
                  double std_1 = std::sqrt(covariances[{s1, s1}]);
                  double std_2 = std::sqrt(covariances[{s2, s2}]);
                  covariances[{s1, s2}] = corr * std_1 * std_2;
                }
            }
          else if (type == typeid(ShocksSurpriseStatement))
            {
              auto* shock = static_cast<ShocksSurpriseStatement*>(statement.get());
              if (shock->overwrite)
                covariances.clear();
              for (const auto& [id, trajectory] : shock->surprise_shocks)
                {
                  string s = table.getName(id);
                  vector<tuple<int, int, double>> traj;
                  for (const auto& [period_range, expr] : trajectory)
                    {
                      if (holds_alternative<pair<int, int>>(period_range))
                        {
                          auto [p1, p2] = get<pair<int, int>>(period_range);
                          traj.emplace_back(p1, p2, expr->eval(eval_ctx));
                        }
                    }
                  trajectories[s] = traj;
                }
            }
        }
    }
  catch (const ExprNode::EvalException&)
    {
      throw EvaluationException("Evaluation error in shocks block");
    }
}

void
DynareModel::set_lead_lag_incidence()
{
  const DynamicModel& dm = mod_file->dynamic_model;
  const SymbolTable& table = mod_file->symbol_table;
  max_endo_lag = dm.max_endo_lag;
  max_endo_lead = dm.max_endo_lead;

  lead_lag_incidence.clear();
  int num_endo = table.endo_nbr();
  lead_lag_incidence.resize(num_endo);

  for (int endoID = 0; endoID < num_endo; endoID++)
    {
      lead_lag_incidence[endoID].resize(max_endo_lag + max_endo_lead + 1, 0);
      for (int lag = -max_endo_lag; lag <= max_endo_lead; lag++)
        {
          try
            {
              int varID = dm.getDerivID(table.getID(SymbolType::endogenous, endoID), lag);
              lead_lag_incidence[endoID][lag + max_endo_lag] = dm.getLegacyJacobianCol(varID) + 1;
            }
          catch (...)
            {
              lead_lag_incidence[endoID][lag + max_endo_lag] = 0;
            }
        }
    }
}

void
DynareModel::set_symbolic_derivatives()
{
  const DynamicModel& dm = mod_file->dynamic_model;
  symb_jacob_endo = vector<symb_jacobian_t>(3);
  symb_jacob_exo.clear();
  symb_jacob_exo_det.clear();

  if (dm.derivatives.size() > 1)
    {
      for (const auto& [vect, expr] : dm.derivatives[1])
        {
          int eq = vect[0];
          for (size_t i = 1; i < vect.size(); i++)
            {
              int derivID = vect[i];
              SymbolType st = dm.getTypeByDerivID(derivID);
              int tsid = dm.getTypeSpecificIDByDerivID(derivID);
              switch (st)
                {
                case SymbolType::endogenous:
                  {
                    int lag = dm.getLagByDerivID(derivID);
                    if (lag > 1 || lag < -1)
                      throw ModelSemanticException("Unsupported lag value: only -1, 0, and 1 are supported");
                    symb_jacob_endo[lag + 1][{eq, tsid}] = expr;
                    break;
                  }
                case SymbolType::exogenous:
                  symb_jacob_exo[{eq, tsid}] = expr;
                  break;
                case SymbolType::exogenousDet:
                  symb_jacob_exo_det[{eq, tsid}] = expr;
                  break;
                default:
                  break;
                }
            }
        }
    }

  symb_jacob_params.clear();
  if (dm.params_derivatives.contains({0, 1}))
    {
      for (const auto& [indices, expr] : dm.params_derivatives.at({0, 1}))
        {
          int eq = indices[0];
          int param = dm.getTypeSpecificIDByDerivID(indices[1]);
          symb_jacob_params[{eq, param}] = expr;
        }
    }

  symb_derivatives.clear();
  for (size_t i = 1; i < dm.derivatives.size(); i++)
    symb_derivatives.push_back(dm.derivatives[i]);
}

double
DynareModel::evaluate_with_lags(
    expr_t expression,
    const vector<span<const double>>& endo,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  try
    {
      switch (expression_type(expression))
        {
        case ExprNodeType::NumConstNode:
          {
            auto* expr = static_cast<NumConstNode*>(expression);
            return mod_file->num_constants.getDouble(expr->id);
          }
        case ExprNodeType::VariableNode:
          {
            auto* expr = static_cast<VariableNode*>(expression);
            int lag = expr->lag;
            if (lag > 1 || lag < -1)
              throw ModelSemanticException("Unsupported lag value: only -1, 0, and 1 are supported");
            int id = expr->symb_id;
            SymbolType type = mod_file->symbol_table.getType(id);
            int sid = mod_file->symbol_table.getTypeSpecificID(id);
            switch (type)
              {
              case SymbolType::endogenous:
                return endo[lag + 1][sid];
              case SymbolType::exogenous:
                return exo[sid];
              case SymbolType::exogenousDet:
                return exo_det[sid];
              case SymbolType::parameter:
                return params[sid];
              default:
                throw ModelSemanticException("Unsupported variable type in evaluation");
              }
          }
        case ExprNodeType::UnaryOpNode:
          {
            auto* expr = static_cast<UnaryOpNode*>(expression);
            double arg = evaluate_with_lags(expr->arg, endo, exo, exo_det, params);
            return expr->eval_opcode(expr->op_code, arg);
          }
        case ExprNodeType::BinaryOpNode:
          {
            auto* expr = static_cast<BinaryOpNode*>(expression);
            double arg1 = evaluate_with_lags(expr->arg1, endo, exo, exo_det, params);
            double arg2 = evaluate_with_lags(expr->arg2, endo, exo, exo_det, params);
            BinaryOpcode opcode = expr->op_code;
            if (opcode == BinaryOpcode::equal)
              opcode = BinaryOpcode::minus;
            return expr->eval_opcode(arg1, opcode, arg2, expr->powerDerivOrder);
          }
        case ExprNodeType::TrinaryOpNode:
          {
            auto* expr = static_cast<TrinaryOpNode*>(expression);
            double arg1 = evaluate_with_lags(expr->arg1, endo, exo, exo_det, params);
            double arg2 = evaluate_with_lags(expr->arg2, endo, exo, exo_det, params);
            double arg3 = evaluate_with_lags(expr->arg3, endo, exo, exo_det, params);
            return expr->eval_opcode(arg1, expr->op_code, arg2, arg3);
          }
        default:
          throw EvaluationException("Unknown expression node type");
        }
      __builtin_unreachable();
    }
  catch (const ExprNode::EvalExternalFunctionException&)
    {
      throw ModelSemanticException("External functions are not supported yet in evaluation.");
    }
  catch (const ExprNode::EvalException&)
    {
      throw EvaluationException("Evaluation error in equation");
    }
}

vector<double>
DynareModel::residuals(
    span<const double> endo_future,
    span<const double> endo_present,
    span<const double> endo_past,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  vector<span<const double>> endo {endo_past, endo_present, endo_future};
  vector<double> res;
  res.reserve(mod_file->dynamic_model.equations.size());
  for (expr_t eq : mod_file->dynamic_model.equations)
    res.push_back(evaluate_with_lags(eq, endo, exo, exo_det, params));
  return res;
}

void
DynareModel::eval_symb_jacob(
    const symb_jacobian_t& symb_jacob,
    jacobian_t& out,
    const vector<span<const double>>& endo,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  for (const auto& [key, expr] : symb_jacob)
    out[key] = evaluate_with_lags(expr, endo, exo, exo_det, params);
}

vector<jacobian_t>
DynareModel::jacobians(
    span<const double> endo_future,
    span<const double> endo_present,
    span<const double> endo_past,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  vector<span<const double>> endo {endo_past, endo_present, endo_future};
  vector<jacobian_t> res(6);
  eval_symb_jacob(symb_jacob_endo[2], res[0], endo, exo, exo_det, params);
  eval_symb_jacob(symb_jacob_endo[1], res[1], endo, exo, exo_det, params);
  eval_symb_jacob(symb_jacob_endo[0], res[2], endo, exo, exo_det, params);
  eval_symb_jacob(symb_jacob_exo, res[3], endo, exo, exo_det, params);
  eval_symb_jacob(symb_jacob_exo_det, res[4], endo, exo, exo_det, params);
  eval_symb_jacob(symb_jacob_params, res[5], endo, exo, exo_det, params);
  return res;
}

JacobianBlocks
DynareModel::jacobian_blocks(
    span<const double> endo_future,
    span<const double> endo_present,
    span<const double> endo_past,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  JacobianBlocks blocks;
  blocks.n_eq = static_cast<int>(mod_file->dynamic_model.equations.size());
  blocks.n_endo = static_cast<int>(endogenous.size());
  blocks.n_exo = static_cast<int>(exogenous.size());
  blocks.n_exo_det = static_cast<int>(exogenous_det.size());
  blocks.n_params = static_cast<int>(parameters.size());

  blocks.lead.assign(blocks.n_eq * blocks.n_endo, 0.0);
  blocks.curr.assign(blocks.n_eq * blocks.n_endo, 0.0);
  blocks.lag.assign(blocks.n_eq * blocks.n_endo, 0.0);
  blocks.exo.assign(blocks.n_eq * blocks.n_exo, 0.0);
  blocks.exo_det.assign(blocks.n_eq * blocks.n_exo_det, 0.0);
  blocks.params.assign(blocks.n_eq * blocks.n_params, 0.0);

  vector<span<const double>> endo {endo_past, endo_present, endo_future};

  // Future / Lead (symb_jacob_endo[2])
  for (const auto& [key, expr] : symb_jacob_endo[2])
    blocks.lead[key.first * blocks.n_endo + key.second]
        = evaluate_with_lags(expr, endo, exo, exo_det, params);

  // Present / Current (symb_jacob_endo[1])
  for (const auto& [key, expr] : symb_jacob_endo[1])
    blocks.curr[key.first * blocks.n_endo + key.second]
        = evaluate_with_lags(expr, endo, exo, exo_det, params);

  // Past / Lag (symb_jacob_endo[0])
  for (const auto& [key, expr] : symb_jacob_endo[0])
    blocks.lag[key.first * blocks.n_endo + key.second]
        = evaluate_with_lags(expr, endo, exo, exo_det, params);

  // Exogenous
  for (const auto& [key, expr] : symb_jacob_exo)
    blocks.exo[key.first * blocks.n_exo + key.second]
        = evaluate_with_lags(expr, endo, exo, exo_det, params);

  // Exogenous deterministic
  for (const auto& [key, expr] : symb_jacob_exo_det)
    blocks.exo_det[key.first * blocks.n_exo_det + key.second]
        = evaluate_with_lags(expr, endo, exo, exo_det, params);

  // Parameters
  for (const auto& [key, expr] : symb_jacob_params)
    blocks.params[key.first * blocks.n_params + key.second]
        = evaluate_with_lags(expr, endo, exo, exo_det, params);

  return blocks;
}

tuple<vector<double>, int, int>
DynareModel::jacobian(
    span<const double> endo_future,
    span<const double> endo_present,
    span<const double> endo_past,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  const DynamicModel& dm = mod_file->dynamic_model;
  int n_eq = static_cast<int>(dm.equations.size());
  int n_cols = dm.getJacobianColsNbr();
  vector<double> mat(n_eq * n_cols, 0.0);

  vector<span<const double>> endo {endo_past, endo_present, endo_future};

  if (dm.derivatives.size() > 1)
    {
      for (const auto& [vect, expr] : dm.derivatives[1])
        {
          int eq = vect[0];
          int derivID = vect[1];
          try
            {
              int col = dm.getJacobianCol(derivID);
              mat[eq * n_cols + col] = evaluate_with_lags(expr, endo, exo, exo_det, params);
            }
          catch (const DataTree::UnknownDerivIDException&)
            {
              // Ignore derivatives not belonging to endogenous/exogenous
            }
        }
    }

  return {mat, n_eq, n_cols};
}

vector<double>
DynareModel::static_residuals(
    span<const double> endo,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  // In static model, lag is 0 for all endogenous variables
  vector<span<const double>> endo_static {endo, endo, endo};
  vector<double> res;
  res.reserve(mod_file->static_model.equations.size());
  for (expr_t eq : mod_file->static_model.equations)
    res.push_back(evaluate_with_lags(eq, endo_static, exo, exo_det, params));
  return res;
}

tuple<vector<double>, int, int>
DynareModel::static_jacobian(
    span<const double> endo,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  const StaticModel& sm = mod_file->static_model;
  int n_eq = static_cast<int>(sm.equations.size());
  int n_endo = static_cast<int>(endogenous.size());
  vector<double> mat(n_eq * n_endo, 0.0);
  vector<span<const double>> endo_static {endo, endo, endo};

  if (sm.derivatives.size() > 1)
    {
      for (const auto& [vect, expr] : sm.derivatives[1])
        {
          int eq = vect[0];
          int derivID = vect[1];
          if (sm.getTypeByDerivID(derivID) == SymbolType::endogenous)
            {
              int tsid = sm.getTypeSpecificIDByDerivID(derivID);
              mat[eq * n_endo + tsid] = evaluate_with_lags(expr, endo_static, exo, exo_det, params);
            }
        }
    }

  return {mat, n_eq, n_endo};
}

derivatives_t
DynareModel::derivatives(
    span<const double> endo_future,
    span<const double> endo_present,
    span<const double> endo_past,
    span<const double> exo,
    span<const double> exo_det,
    span<const double> params) const
{
  vector<span<const double>> endo {endo_past, endo_present, endo_future};
  derivatives_t res;
  for (const auto& fixed_order_derivs : symb_derivatives)
    {
      coo_matrix_t coords;
      for (const auto& [vect, expr] : fixed_order_derivs)
        coords.emplace_back(vect, evaluate_with_lags(expr, endo, exo, exo_det, params));
      res.push_back(coords);
    }
  return res;
}

