/*
 * Copyright © 2003-2026 Dynare Team
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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "DynamicModel.hh"
#include "StaticModel.hh"

StaticModel::StaticModel(SymbolTable& symbol_table_arg, NumericalConstants& num_constants_arg,
                         ExternalFunctionsTable& external_functions_table_arg,
                         HeterogeneityTable& heterogeneity_table_arg,
                         DatabaseTable& database_table_arg) :
    ModelTree {symbol_table_arg, num_constants_arg, external_functions_table_arg,
               heterogeneity_table_arg, database_table_arg}
{
}

void
StaticModel::copyHelper(const StaticModel& m)
{
  auto f = [this](const ExprNode* e) { return e->clone(*this); };

  for (const auto& it : m.ramsey_multipliers_derivatives_temporary_terms)
    ramsey_multipliers_derivatives_temporary_terms.insert(f(it));
  for (const auto& it : m.ramsey_multipliers_derivatives_temporary_terms_idxs)
    ramsey_multipliers_derivatives_temporary_terms_idxs.emplace(f(it.first), it.second);
}

StaticModel::StaticModel(const StaticModel& m) :
    ModelTree {m},
    ramsey_multipliers_derivatives {m.ramsey_multipliers_derivatives},
    ramsey_multipliers_derivatives_sparse_colptr {m.ramsey_multipliers_derivatives_sparse_colptr},
    static_mfs {m.static_mfs}
{
  copyHelper(m);
}

StaticModel&
StaticModel::operator=(const StaticModel& m)
{
  ModelTree::operator=(m);

  ramsey_multipliers_derivatives = m.ramsey_multipliers_derivatives;
  ramsey_multipliers_derivatives_sparse_colptr = m.ramsey_multipliers_derivatives_sparse_colptr;

  static_mfs = m.static_mfs;

  copyHelper(m);

  return *this;
}

StaticModel::StaticModel(const DynamicModel& m) :
    ModelTree {m.symbol_table, m.num_constants, m.external_functions_table, m.heterogeneity_table,
               m.database_table}
{
  // Convert model local variables (need to be done first)
  for (int it : m.local_variables_vector)
    AddLocalVariable(it, m.local_variables_table.at(it)->toStatic(*this));

  // Convert equations
  int static_only_index = 0;
  set<int> dynamic_equations = m.equation_tags.getDynamicEqns();
  for (int i = 0; i < static_cast<int>(m.equations.size()); i++)
    try
      {
        // If equation is dynamic, replace it by an equation marked [static]
        if (dynamic_equations.contains(i))
          {
            auto [static_only_equations, static_only_equations_lineno,
                  static_only_complementarity_conditions, static_only_equations_equation_tags]
                = m.getStaticOnlyEquationsInfo();

            addEquation(static_only_equations[static_only_index]->toStatic(*this),
                        static_only_equations_lineno[static_only_index],
                        static_only_complementarity_conditions[static_only_index],
                        static_only_equations_equation_tags.getTagsByEqn(static_only_index));
            static_only_index++;
          }
        else
          addEquation(m.equations[i]->toStatic(*this), m.equations_lineno[i],
                      m.complementarity_conditions[i], m.equation_tags.getTagsByEqn(i));
      }
    catch (const DataTree::DivisionByZeroException& e)
      {
        cerr << "...division by zero error encountered when converting equation " << i;
        if (!e.message.empty())
          cerr << " (" << e.message << ")";
        cerr << '\n';
        exit(EXIT_FAILURE);
      }

  // Convert auxiliary equations
  for (auto aux_eq : m.aux_equations)
    addAuxEquation(aux_eq->toStatic(*this));

  cutoff = m.cutoff;

  user_set_add_flags = m.user_set_add_flags;
  user_set_subst_flags = m.user_set_subst_flags;
  user_set_add_libs = m.user_set_add_libs;
  user_set_subst_libs = m.user_set_subst_libs;
  user_set_compiler = m.user_set_compiler;

  static_mfs = m.static_mfs;
}

void
StaticModel::writeStaticBytecode(const string& basename) const
{
  /* Bytecode only works when there are with as many endogenous as equations.
     (e.g. the constructor of FBEGINBLOCK_ makes this assumption) */
  assert(static_cast<int>(equations.size()) == symbol_table.endo_nbr());

  // First write the .bin file
  int u_count_int {writeBytecodeBinFile(basename + "/model/bytecode/static.bin", false)};

  Bytecode::Writer code_file {basename + "/model/bytecode/static.cod"};
  auto eq_idx = views::iota(0, static_cast<int>(equations.size()));
  auto endo_idx = views::iota(0, symbol_table.endo_nbr());

  // Declare temporary terms and the (single) block
  code_file << Bytecode::FDIMST {static_cast<int>(temporary_terms_derivatives[0].size()
                                                  + temporary_terms_derivatives[1].size())}
            << Bytecode::FBEGINBLOCK {symbol_table.endo_nbr(),
                                      BlockSimulationType::solveForwardComplete,
                                      0,
                                      symbol_table.endo_nbr(),
                                      {endo_idx.begin(), endo_idx.end()},
                                      {eq_idx.begin(), eq_idx.end()},
                                      false,
                                      u_count_int,
                                      symbol_table.endo_nbr()};

  writeBytecodeHelper<false>(code_file);
}

void
StaticModel::writeStaticBlockBytecode(const string& basename) const
{
  Bytecode::Writer code_file {basename + "/model/bytecode/block/static.cod"};

  const filesystem::path bin_filename {basename + "/model/bytecode/block/static.bin"};
  ofstream bin_file {bin_filename, ios::out | ios::binary};
  if (!bin_file.is_open())
    {
      cerr << R"(Error : Can't open file ")" << bin_filename.string() << R"(" for writing)" << '\n';
      exit(EXIT_FAILURE);
    }

  // Temporary variables declaration
  code_file << Bytecode::FDIMST {static_cast<int>(blocks_temporary_terms_idxs.size())};

  temporary_terms_t temporary_terms_written;

  for (int block {0}; block < static_cast<int>(blocks.size()); block++)
    {
      const BlockSimulationType simulation_type {blocks[block].simulation_type};
      const int block_size {blocks[block].size};

      const int u_count {simulation_type == BlockSimulationType::solveBackwardComplete
                                 || simulation_type == BlockSimulationType::solveForwardComplete
                             ? writeBlockBytecodeBinFile(bin_file, block)
                             : 0};

      code_file << Bytecode::FBEGINBLOCK {blocks[block].mfs_size,
                                          simulation_type,
                                          blocks[block].first_equation,
                                          block_size,
                                          endo_idx_block2orig,
                                          eq_idx_block2orig,
                                          blocks[block].linear,
                                          u_count,
                                          block_size};

      writeBlockBytecodeHelper<false>(code_file, block, temporary_terms_written);
    }
  code_file << Bytecode::FEND {};
}

void
StaticModel::computingPass(int derivsOrder, int paramsDerivsOrder,
                           const eval_context_t& eval_context, bool no_tmp_terms, bool block,
                           bool use_dll)
{
  initializeVariablesAndEquations();

  /* In both MATLAB and Julia, tensors for higher-order derivatives are stored
     in matrices whose columns correspond to variable multi-indices. Since we
     currently are limited to 32-bit signed integers (hence 31 bits) for matrix
     indices, check that we will not overflow (see #89). Note that such a check
     is not needed for parameter derivatives, since tensors for those are not
     stored as matrices. This check is implemented at this place for symmetry
     with DynamicModel::computingPass(). */
  if (log2(symbol_table.endo_nbr()) * derivsOrder >= numeric_limits<int>::digits)
    {
      cerr << "ERROR: The derivatives matrix of the " << modelClassName()
           << " is too large. Please decrease the approximation order." << '\n';
      exit(EXIT_FAILURE);
    }

  // Compute derivatives w.r. to all endogenous
  set<int> vars;
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      int id = symbol_table.getID(SymbolType::endogenous, i);
      vars.insert(getDerivID(id, 0));
    }

  // Launch computations
  cout << "Computing " << modelClassName() << " derivatives (order " << derivsOrder << ")." << '\n'
       << flush;

  computeDerivatives(derivsOrder, vars);

  if (paramsDerivsOrder > 0)
    {
      cout << "Computing " << modelClassName() << " derivatives w.r.t. parameters (order "
           << paramsDerivsOrder << ")." << '\n'
           << flush;
      computeParamsDerivatives(paramsDerivsOrder);
    }

  computeTemporaryTerms(!use_dll, no_tmp_terms);

  if (paramsDerivsOrder > 0 && !no_tmp_terms)
    computeParamsDerivativesTemporaryTerms();

  computingPassBlock(eval_context, no_tmp_terms);
  if (!block_decomposed && block)
    {
      cerr << "ERROR: Block decomposition requested but failed. If your model does not have a "
              "steady state, you may want to try the 'no_static' option of the 'model' block."
           << '\n';
      exit(EXIT_FAILURE);
    }

  computeMCPEquationsReordering();
}

void
StaticModel::writeStaticFile(const string& basename, bool use_dll, const string& mexext,
                             const filesystem::path& matlabroot, bool julia) const
{
  filesystem::path model_dir {basename};
  model_dir /= "model";
  if (use_dll)
    {
      create_directories(model_dir / "src");
      if (block_decomposed)
        create_directories(model_dir / "src" / "block");
    }
  if (julia)
    create_directories(model_dir / "julia");
  else
    {
      auto plusfolder {packageDir(basename)};
      /* The following is not a duplicate of the same call from
         ModFile::writeMOutput(), because of planner_objective which needs its
         +objective subdirectory */
      create_directories(plusfolder);

      if (block_decomposed)
        create_directories(plusfolder / "+block");

      create_directories(plusfolder / "+debug");
    }
  create_directories(model_dir / "bytecode" / "block");

  /* PlannerObjective subclass or discretionary optimal policy models don’t
     have as many variables as equations; bytecode does not support that
     case */
  if (static_cast<int>(equations.size()) == symbol_table.endo_nbr())
    writeStaticBytecode(basename);
  if (block_decomposed)
    writeStaticBlockBytecode(basename);

  if (use_dll)
    writeModelCFiles<false>(basename, mexext, matlabroot);
  else if (julia)
    writeModelJuliaFiles<false>(basename);
  else // MATLAB/Octave
    writeModelMFiles<false>(basename);

  writeSetAuxiliaryVariablesFile<false>(basename, julia);

  if (!julia)
    {
      writeComplementarityConditionsFile<false>(basename);
      // Support for model debugging
      writeDebugModelMFiles<false>(basename);
    }
}

bool
StaticModel::exoPresentInEqs() const
{
  for (auto equation : equations)
    if (equation->hasExogenous())
      return true;
  return false;
}

void
StaticModel::writeDriverOutput(ostream& output) const
{
  output << "M_.static_tmp_nbr = [";
  for (const auto& temporary_terms_derivative : temporary_terms_derivatives)
    output << temporary_terms_derivative.size() << "; ";
  output << "];" << '\n';

  if (block_decomposed)
    writeBlockDriverOutput(output);

  writeDriverSparseIndicesHelper("static", output);

  output << "M_.static_mcp_equations_reordering = [";
  for (auto i : mcp_equations_reordering)
    output << i + 1 << "; ";
  output << "];" << '\n';
}

void
StaticModel::writeBlockDriverOutput(ostream& output) const
{
  for (int blk = 0; blk < static_cast<int>(blocks.size()); blk++)
    {
      output << "M_.block_structure_stat.block(" << blk + 1
             << ").Simulation_Type = " << static_cast<int>(blocks[blk].simulation_type) << ";"
             << '\n'
             << "M_.block_structure_stat.block(" << blk + 1 << ").endo_nbr = " << blocks[blk].size
             << ";" << '\n'
             << "M_.block_structure_stat.block(" << blk + 1 << ").mfs = " << blocks[blk].mfs_size
             << ";" << '\n'
             << "M_.block_structure_stat.block(" << blk + 1 << ").equation = [";
      for (int eq = 0; eq < blocks[blk].size; eq++)
        output << " " << getBlockEquationID(blk, eq) + 1;
      output << "];" << '\n' << "M_.block_structure_stat.block(" << blk + 1 << ").variable = [";
      for (int var = 0; var < blocks[blk].size; var++)
        output << " " << getBlockVariableID(blk, var) + 1;
      output << "];" << '\n';
    }
  output << "M_.block_structure_stat.variable_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << endo_idx_block2orig[i] + 1;
  output << "];" << '\n' << "M_.block_structure_stat.equation_reordered = [";
  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    output << " " << eq_idx_block2orig[i] + 1;
  output << "];" << '\n';

  set<pair<int, int>> row_incidence;
  for (const auto& [indices, d1] : derivatives[1])
    if (int deriv_id = indices[1]; getTypeByDerivID(deriv_id) == SymbolType::endogenous)
      {
        int eq = indices[0];
        int var {getTypeSpecificIDByDerivID(deriv_id)};
        row_incidence.emplace(eq, var);
      }
  output << "M_.block_structure_stat.incidence.sparse_IM = [" << '\n';
  for (auto [eq, var] : row_incidence)
    output << " " << eq + 1 << " " << var + 1 << ";" << '\n';
  output << "];" << '\n'
         << "M_.block_structure_stat.tmp_nbr = " << blocks_temporary_terms_idxs.size() << ";"
         << '\n';

  writeBlockDriverSparseIndicesHelper<false>(output);
}

SymbolType
StaticModel::getTypeByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < symbol_table.endo_nbr())
    return SymbolType::endogenous;
  else if (deriv_id < symbol_table.endo_nbr() + symbol_table.param_nbr())
    return SymbolType::parameter;
  else
    throw UnknownDerivIDException();
}

int
StaticModel::getLagByDerivID([[maybe_unused]] int deriv_id) const noexcept(false)
{
  return 0;
}

int
StaticModel::getSymbIDByDerivID(int deriv_id) const noexcept(false)
{
  if (deriv_id < symbol_table.endo_nbr())
    return symbol_table.getID(SymbolType::endogenous, deriv_id);
  else if (deriv_id < symbol_table.endo_nbr() + symbol_table.param_nbr())
    return symbol_table.getID(SymbolType::parameter, deriv_id - symbol_table.endo_nbr());
  else
    throw UnknownDerivIDException();
}

int
StaticModel::getTypeSpecificIDByDerivID(int deriv_id) const
{
  if (deriv_id < symbol_table.endo_nbr())
    return deriv_id;
  else if (deriv_id < symbol_table.endo_nbr() + symbol_table.param_nbr())
    return deriv_id - symbol_table.endo_nbr();
  else
    throw UnknownDerivIDException();
}

int
StaticModel::getDerivID(int symb_id, [[maybe_unused]] int lag) const noexcept(false)
{
  if (symbol_table.getType(symb_id) == SymbolType::endogenous)
    return symbol_table.getTypeSpecificID(symb_id);
  else if (symbol_table.getType(symb_id) == SymbolType::parameter)
    return symbol_table.getTypeSpecificID(symb_id) + symbol_table.endo_nbr();
  else
    /* See the special treatment in VariableNode::prepareForDerivation(),
       VariableNode::computeDerivative() and VariableNode::getChainRuleDerivative() */
    throw UnknownDerivIDException {};
}

void
StaticModel::addAllParamDerivId(set<int>& deriv_id_set)
{
  for (int i = 0; i < symbol_table.param_nbr(); i++)
    deriv_id_set.insert(i + symbol_table.endo_nbr());
}

void
StaticModel::computeChainRuleJacobian()
{
  int nb_blocks = blocks.size();

  blocks_derivatives.resize(nb_blocks);
  blocks_jacobian_sparse_column_major_order.resize(nb_blocks);
  blocks_jacobian_sparse_colptr.resize(nb_blocks);

  for (int blk = 0; blk < nb_blocks; blk++)
    {
      int nb_recursives = blocks[blk].getRecursiveSize();
      BlockSimulationType simulation_type {blocks[blk].simulation_type};

      map<int, BinaryOpNode*> recursive_vars;
      for (int i = 0; i < nb_recursives; i++)
        {
          int deriv_id = getDerivID(
              symbol_table.getID(SymbolType::endogenous, getBlockVariableID(blk, i)), 0);
          if (getBlockEquationType(blk, i) == EquationType::evaluateRenormalized)
            recursive_vars[deriv_id] = getBlockEquationRenormalizedExpr(blk, i);
          else
            recursive_vars[deriv_id] = getBlockEquationExpr(blk, i);
        }

      assert(simulation_type != BlockSimulationType::solveTwoBoundariesSimple
             && simulation_type != BlockSimulationType::solveTwoBoundariesComplete);

      int size = blocks[blk].size;
      unordered_map<expr_t, set<int>> non_null_chain_rule_derivatives;
      unordered_map<expr_t, map<int, expr_t>> chain_rule_deriv_cache;
      for (int eq = nb_recursives; eq < size; eq++)
        {
          int eq_orig = getBlockEquationID(blk, eq);
          for (int var = nb_recursives; var < size; var++)
            {
              int var_orig = getBlockVariableID(blk, var);
              expr_t d1 = equations[eq_orig]->getChainRuleDerivative(
                  getDerivID(symbol_table.getID(SymbolType::endogenous, var_orig), 0),
                  recursive_vars, non_null_chain_rule_derivatives, chain_rule_deriv_cache);
              if (d1 != Zero)
                blocks_derivatives[blk][{eq, var, 0}] = d1;
            }
        }

      // Compute the CSC representation of the Jacobian
      if (simulation_type != BlockSimulationType::evaluateForward
          && simulation_type != BlockSimulationType::evaluateBackward)
        {
          for (const auto& [indices, d1] : blocks_derivatives[blk])
            {
              auto& [eq, var, lag] {indices};
              assert(eq >= nb_recursives && var >= nb_recursives && lag == 0);
              blocks_jacobian_sparse_column_major_order[blk].try_emplace(
                  {eq - nb_recursives, var - nb_recursives}, d1);
            }
          blocks_jacobian_sparse_colptr[blk] = computeCSCColPtr(
              blocks_jacobian_sparse_column_major_order[blk], blocks[blk].mfs_size);
        }
    }
}

void
StaticModel::writeLatexFile(const string& basename, bool write_equation_tags) const
{
  writeLatexModelFile(basename, "static", ExprNodeOutputType::latexStaticModel,
                      write_equation_tags);
}

void
StaticModel::writeAuxVarInitval(ostream& output, ExprNodeOutputType output_type) const
{
  for (auto aux_equation : aux_equations)
    {
      dynamic_cast<ExprNode*>(aux_equation)->writeOutput(output, output_type);
      output << ";" << '\n';
    }
}

void
StaticModel::writeLatexAuxVarRecursiveDefinitions(ostream& output) const
{
  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;
  temporary_terms_idxs_t temporary_terms_idxs;
  for (auto aux_equation : aux_equations)
    if (aux_equation->containsExternalFunction())
      aux_equation->writeExternalFunctionOutput(output, ExprNodeOutputType::latexStaticModel,
                                                temporary_terms, temporary_terms_idxs, tef_terms);
  for (auto aux_equation : aux_equations)
    {
      output << R"(\begin{dmath})" << '\n';
      dynamic_cast<ExprNode*>(aux_equation)
          ->writeOutput(output, ExprNodeOutputType::latexStaticModel);
      output << '\n' << R"(\end{dmath})" << '\n';
    }
}

void
StaticModel::writeJsonAuxVarRecursiveDefinitions(ostream& output) const
{
  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t temporary_terms;

  for (auto aux_equation : aux_equations)
    if (aux_equation->containsExternalFunction())
      {
        vector<string> efout;
        aux_equation->writeJsonExternalFunctionOutput(efout, temporary_terms, tef_terms, false);
        for (bool printed_something {false}; const auto& it : efout)
          {
            if (exchange(printed_something, true))
              output << ", ";
            output << it;
          }
      }

  for (auto aux_equation : aux_equations)
    {
      output << R"(, {"lhs": ")";
      aux_equation->arg1->writeJsonOutput(output, temporary_terms, tef_terms, false);
      output << R"(", "rhs": ")";
      aux_equation->arg2->writeJsonOutput(output, temporary_terms, tef_terms, false);
      output << R"("})";
    }
}

void
StaticModel::writeJsonOutput(ostream& output) const
{
  output << R"("static_tmp_nbr": [)";
  for (bool printed_something {false}; const auto& tts : temporary_terms_derivatives)
    {
      if (exchange(printed_something, true))
        output << ", ";
      output << tts.size();
    }
  output << "], ";
  writeJsonSparseIndicesHelper<false>(output);
}

void
StaticModel::writeJsonComputingPassOutput(ostream& output, bool writeDetails) const
{
  auto [mlv_output, d_output] {writeJsonComputingPassOutputHelper<false>(writeDetails)};

  if (writeDetails)
    output << R"("static_model": {)";
  else
    output << R"("static_model_simple": {)";
  output << mlv_output.str();
  for (const auto& it : d_output)
    output << ", " << it.str();
  output << "}";
}

void
StaticModel::writeJsonParamsDerivatives(ostream& output, bool writeDetails) const
{
  if (!params_derivatives.size())
    return;

  auto [mlv_output, tt_output, rp_output, gp_output, rpp_output, gpp_output, hp_output,
        g3p_output] {writeJsonParamsDerivativesHelper<false>(writeDetails)};
  // g3p_output is ignored

  if (writeDetails)
    output << R"("static_model_params_derivative": {)";
  else
    output << R"("static_model_params_derivatives_simple": {)";
  output << mlv_output.str() << ", " << tt_output.str() << ", " << rp_output.str() << ", "
         << gp_output.str() << ", " << rpp_output.str() << ", " << gpp_output.str() << ", "
         << hp_output.str() << "}";
}

void
StaticModel::computeRamseyMultipliersDerivatives(int ramsey_orig_endo_nbr, bool is_matlab,
                                                 bool no_tmp_terms)
{
  // Compute derivation IDs of Lagrange multipliers
  set<int> mult_symb_ids {symbol_table.getLagrangeMultipliers()};
  vector<int> mult_deriv_ids;
  mult_deriv_ids.reserve(mult_symb_ids.size());
  for (int symb_id : mult_symb_ids)
    mult_deriv_ids.push_back(getDerivID(symb_id, 0));

  // Compute the list of aux vars for which to apply the chain rule derivation
  map<int, BinaryOpNode*> recursive_variables;
  for (auto aux_eq : aux_equations)
    {
      auto varexpr {dynamic_cast<VariableNode*>(aux_eq->arg1)};
      assert(varexpr && symbol_table.isAuxiliaryVariable(varexpr->symb_id));
      /* Determine whether the auxiliary variable has been added after the last
         Lagrange multiplier. We use the guarantee given by SymbolTable that
         symbol IDs are increasing. */
      if (varexpr->symb_id > *mult_symb_ids.crbegin())
        recursive_variables.emplace(varexpr->getDerivID(), aux_eq);
    }

  // Compute the chain rule derivatives w.r.t. multipliers
  unordered_map<expr_t, set<int>> non_null_chain_rule_derivatives;
  unordered_map<expr_t, map<int, expr_t>> cache;
  for (int eq {0}; eq < ramsey_orig_endo_nbr; eq++)
    for (int mult {0}; mult < static_cast<int>(mult_deriv_ids.size()); mult++)
      if (expr_t d {equations[eq]->getChainRuleDerivative(mult_deriv_ids[mult], recursive_variables,
                                                          non_null_chain_rule_derivatives, cache)};
          d != Zero)
        ramsey_multipliers_derivatives.try_emplace({eq, mult}, d);

  // Compute the temporary terms
  map<pair<int, int>, unordered_set<expr_t>> temp_terms_map;
  unordered_map<expr_t, pair<int, pair<int, int>>> reference_count;
  for (const auto& [row_col, d] : ramsey_multipliers_derivatives)
    d->computeTemporaryTerms({1, 0}, temp_terms_map, reference_count, is_matlab);
  /* If the user has specified the notmpterms option, clear all temporary
     terms, except those that correspond to external functions (since they are
     not optional) */
  if (no_tmp_terms)
    for (auto& it : temp_terms_map)
      erase_if(it.second, [](expr_t e) { return !dynamic_cast<AbstractExternalFunctionNode*>(e); });
  ranges::copy(temp_terms_map[{1, 0}],
               inserter(ramsey_multipliers_derivatives_temporary_terms,
                        ramsey_multipliers_derivatives_temporary_terms.begin()));
  for (int idx {0}; auto it : ramsey_multipliers_derivatives_temporary_terms)
    ramsey_multipliers_derivatives_temporary_terms_idxs[it] = idx++;

  // Compute the CSC format
  ramsey_multipliers_derivatives_sparse_colptr
      = computeCSCColPtr(ramsey_multipliers_derivatives, mult_deriv_ids.size());
}

void
StaticModel::writeDriverRamseyMultipliersDerivativesSparseIndices(ostream& output) const
{
  output << "M_.ramsey_multipliers_static_g1_sparse_rowval = int32([";
  for (auto& [row_col, d] : ramsey_multipliers_derivatives)
    output << row_col.first + 1 << ' ';
  output << "]);" << '\n' << "M_.ramsey_multipliers_static_g1_sparse_colval = int32([";
  for (auto& [row_col, d] : ramsey_multipliers_derivatives)
    output << row_col.second + 1 << ' ';
  output << "]);" << '\n' << "M_.ramsey_multipliers_static_g1_sparse_colptr = int32([";
  for (int it : ramsey_multipliers_derivatives_sparse_colptr)
    output << it + 1 << ' ';
  output << "]);" << '\n';
}

void
StaticModel::writeRamseyMultipliersDerivativesMFile(const string& basename,
                                                    int ramsey_orig_endo_nbr) const
{
  constexpr auto output_type {ExprNodeOutputType::matlabStaticModel};
  filesystem::path filename {packageDir(basename) / "ramsey_multipliers_static_g1.m"};
  ofstream output_file {filename, ios::out | ios::binary};
  if (!output_file.is_open())
    {
      cerr << "ERROR: Can't open file " << filename.string() << " for writing" << '\n';
      exit(EXIT_FAILURE);
    }

  output_file << "function g1m = ramsey_multipliers_static_g1(y, x, params, sparse_rowval, "
                 "sparse_colval, sparse_colptr)"
              << '\n'
              << "g1m_v=NaN(" << ramsey_multipliers_derivatives.size() << ",1);" << '\n';

  writeRamseyMultipliersDerivativesHelper<output_type>(output_file);

  output_file << "g1m = sparse(sparse_rowval, sparse_colval, g1m_v, " << ramsey_orig_endo_nbr
              << ", " << symbol_table.getLagrangeMultipliers().size() << ");" << '\n'
              << "end" << '\n';
  output_file.close();
}

void
StaticModel::writeRamseyMultipliersDerivativesCFile(const string& basename, const string& mexext,
                                                    const filesystem::path& matlabroot,
                                                    int ramsey_orig_endo_nbr) const
{
  constexpr auto output_type {ExprNodeOutputType::CStaticModel};
  const filesystem::path model_src_dir {filesystem::path {basename} / "model" / "src"};

  const int xlen {symbol_table.exo_nbr() + symbol_table.exo_det_nbr()};
  const int nzval {static_cast<int>(ramsey_multipliers_derivatives.size())};
  const int ncols {static_cast<int>(symbol_table.getLagrangeMultipliers().size())};

  const filesystem::path p {model_src_dir / "ramsey_multipliers_static_g1.c"};
  ofstream output {p, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << p.string() << " for writing" << '\n';
      exit(EXIT_FAILURE);
    }

  output << "#include <math.h>" << '\n'
         << '\n'
         << R"(#include "mex.h")" << '\n' // Needed for calls to external functions
         << '\n';
  writeCHelpersDefinition(output);
  writeCHelpersDeclaration(output); // Provide external definition of helpers
  output << '\n'
         << "void ramsey_multipliers_static_g1(const double *restrict y, const double *restrict x, "
            "const double *restrict params, double *restrict T, double *restrict g1m_v)"
         << '\n'
         << "{" << '\n';
  writeRamseyMultipliersDerivativesHelper<output_type>(output);
  output
      << "}" << '\n'
      << '\n'
      << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << '\n'
      << "{" << '\n'
      << "  if (nrhs != 6)" << '\n'
      << R"(    mexErrMsgTxt("Accepts exactly 6 input arguments");)" << '\n'
      << "  if (nlhs != 1)" << '\n'
      << R"(    mexErrMsgTxt("Accepts exactly 1 output argument");)" << '\n'
      << "  if (!(mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) && !mxIsSparse(prhs[0]) && "
         "mxGetNumberOfElements(prhs[0]) == "
      << symbol_table.endo_nbr() << "))" << '\n'
      << R"(    mexErrMsgTxt("y must be a real dense numeric array with )"
      << symbol_table.endo_nbr() << R"( elements");)" << '\n'
      << "  const double *restrict y = mxGetDoubles(prhs[0]);" << '\n'
      << "  if (!(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) && !mxIsSparse(prhs[1]) && "
         "mxGetNumberOfElements(prhs[1]) == "
      << xlen << "))" << '\n'
      << R"(    mexErrMsgTxt("x must be a real dense numeric array with )" << xlen
      << R"( elements");)" << '\n'
      << "  const double *restrict x = mxGetDoubles(prhs[1]);" << '\n'
      << "  if (!(mxIsDouble(prhs[2]) && !mxIsComplex(prhs[2]) && !mxIsSparse(prhs[2]) && "
         "mxGetNumberOfElements(prhs[2]) == "
      << symbol_table.param_nbr() << "))" << '\n'
      << R"(    mexErrMsgTxt("params must be a real dense numeric array with )"
      << symbol_table.param_nbr() << R"( elements");)" << '\n'
      << "  const double *restrict params = mxGetDoubles(prhs[2]);" << '\n'
      << "  if (!(mxIsInt32(prhs[3]) && mxGetNumberOfElements(prhs[3]) == " << nzval << "))" << '\n'
      << R"(    mexErrMsgTxt("sparse_rowval must be an int32 array with )" << nzval
      << R"( elements");)" << '\n'
      << "  if (!(mxIsInt32(prhs[5]) && mxGetNumberOfElements(prhs[5]) == " << ncols + 1 << "))"
      << '\n'
      << R"(    mexErrMsgTxt("sparse_colptr must be an int32 array with )" << ncols + 1
      << R"( elements");)" << '\n'
      << "  const int32_T *restrict sparse_rowval = mxGetInt32s(prhs[3]);" << '\n'
      << "  const int32_T *restrict sparse_colptr = mxGetInt32s(prhs[5]);" << '\n'
      << "  plhs[0] = mxCreateSparse(" << ramsey_orig_endo_nbr << ", " << ncols << ", " << nzval
      << ", mxREAL);" << '\n'
      << "  mwIndex *restrict ir = mxGetIr(plhs[0]), *restrict jc = mxGetJc(plhs[0]);" << '\n'
      << "  for (mwSize i = 0; i < " << nzval << "; i++)" << '\n'
      << "    *ir++ = *sparse_rowval++ - 1;" << '\n'
      << "  for (mwSize i = 0; i < " << ncols + 1 << "; i++)" << '\n'
      << "    *jc++ = *sparse_colptr++ - 1;" << '\n'
      << "  mxArray *T_mx = mxCreateDoubleMatrix("
      << ramsey_multipliers_derivatives_temporary_terms.size() << ", 1, mxREAL);" << '\n'
      << "  ramsey_multipliers_static_g1(y, x, params, mxGetDoubles(T_mx), mxGetDoubles(plhs[0]));"
      << '\n'
      << "  mxDestroyArray(T_mx);" << '\n'
      << "}" << '\n';

  output.close();

  compileMEX(packageDir(basename), "ramsey_multipliers_static_g1", mexext, {p}, matlabroot);
}

void
StaticModel::writePythonStaticFile(const string& basename, bool use_jax, bool use_numba) const
{
  filesystem::path model_dir {basename};
  model_dir /= "python";
  filesystem::create_directories(model_dir);

  string filename = model_dir.string() + "/static.py";
  ofstream output {filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << writePythonHelpersImports(use_jax, use_numba) << endl
         << writePythonHelpersImportFromInit() << endl
         << endl;

  // Write sparse indices functions
  writePythonSparseIndicesHelper("static", output, use_jax);

  // static_resid
  output << "def static_resid(y: np.ndarray, x: np.ndarray, params: np.ndarray) -> np.ndarray:"
         << endl;
  if (use_jax && use_numba)
    output << "    # Performance decorator removed for compatibility"
           << endl; // Or handle appropriately

  if (!temporary_terms_idxs.empty())
    output << "    T = " << (use_jax ? "jnp" : "np") << ".zeros(" << temporary_terms_idxs.size()
           << ")" << endl;

  temporary_terms_t temp_terms_union;
  deriv_node_temp_terms_t tef_terms;

  // Temporary terms for residuals
  writeTemporaryTerms<ExprNodeOutputType::pythonStaticModel>(
      temporary_terms_derivatives[0], temp_terms_union, temporary_terms_idxs, output, tef_terms);

  output << "    residual = " << (use_jax ? "jnp" : "np") << ".zeros(" << equations.size() << ")"
         << endl;

  // Write equations
  for (int eq = 0; eq < static_cast<int>(equations.size()); eq++)
    {
      BinaryOpNode* eq_node = equations[eq];
      expr_t lhs = eq_node->arg1;
      expr_t rhs = eq_node->arg2;
      output << "    residual[" << eq << "] = ";
      lhs->writeOutput(output, ExprNodeOutputType::pythonStaticModel, temp_terms_union,
                       temporary_terms_idxs, tef_terms);
      output << " - (";
      rhs->writeOutput(output, ExprNodeOutputType::pythonStaticModel, temp_terms_union,
                       temporary_terms_idxs, tef_terms);
      output << ")" << endl;
    }

  output << "    return residual" << endl << endl;

  // static_g1 (Jacobian)
  output << "def static_g1(y: np.ndarray, x: np.ndarray, params: np.ndarray) -> np.ndarray:"
         << endl;
  if (computed_derivs_order >= 1)
    {
      if (!temporary_terms_idxs.empty())
        output << "    T = " << (use_jax ? "jnp" : "np") << ".zeros(" << temporary_terms_idxs.size()
               << ")" << endl;

      writeTemporaryTerms<ExprNodeOutputType::pythonStaticModel>(
          temporary_terms_derivatives[1], temp_terms_union, temporary_terms_idxs, output,
          tef_terms);
    }

  output << "    g1_v = " << (use_jax ? "jnp" : "np") << ".zeros("
         << jacobian_sparse_column_major_order.size() << ")" << endl;

  int i = 0;
  for (const auto& [indices, d1] : jacobian_sparse_column_major_order)
    {
      output << "    g1_v[" << i << "] = ";
      d1->writeOutput(output, ExprNodeOutputType::pythonStaticModel, temp_terms_union,
                      temporary_terms_idxs, tef_terms);
      output << endl;
      i++;
    }
  output << "    return g1_v" << endl << endl;

  // static_g2, g3...
  for (int o = 2; o <= computed_derivs_order; o++)
    {
      output << "def static_g" << o
             << "(y: np.ndarray, x: np.ndarray, params: np.ndarray) -> np.ndarray:" << endl;
      if (!temporary_terms_idxs.empty())
        output << "    T = " << (use_jax ? "jnp" : "np") << ".zeros(" << temporary_terms_idxs.size()
               << ")" << endl;

      writeTemporaryTerms<ExprNodeOutputType::pythonStaticModel>(
          temporary_terms_derivatives[o], temp_terms_union, temporary_terms_idxs, output,
          tef_terms);

      output << "    g" << o << "_v = " << (use_jax ? "jnp" : "np") << ".zeros("
             << derivatives[o].size() << ")" << endl;

      int idx = 0;
      for (const auto& [indices, d] : derivatives[o])
        {
          output << "    g" << o << "_v[" << idx << "] = ";
          d->writeOutput(output, ExprNodeOutputType::pythonStaticModel, temp_terms_union,
                         temporary_terms_idxs, tef_terms);
          output << endl;
          idx++;
        }
      output << "    return g" << o << "_v" << endl << endl;
    }

  // Auxiliary variables
  ostringstream aux_output;
  writeAuxVarRecursiveDefinitions(aux_output, ExprNodeOutputType::pythonStaticModel);
  if (!aux_output.str().empty())
    {
      output << "def set_auxiliary_variables(y: np.ndarray, x: np.ndarray, params: np.ndarray) -> "
                "np.ndarray:"
             << endl
             << "    \"\"\"Compute auxiliary variables of the static model.\"\"\"" << endl;
      // Add indentation to each line of auxiliary equations
      istringstream aux_stream(aux_output.str());
      string line;
      while (getline(aux_stream, line))
        output << "    " << line << endl;
      output << "    return y" << endl << endl;
    }

  output.close();
}

void
StaticModel::writePythonStaticParamsDerivatives(const string& basename, bool use_jax,
                                                bool use_numba) const
{
  if (!params_derivatives.size())
    return;

  auto [tt_output, rp_output, g1p_output, rpp_output, g1pp_output, g2p_output,
        g3p_output] {writeParamsDerivativesFileHelper<ExprNodeOutputType::pythonStaticModel>()};
  // g3p_output is ignored for static model

  filesystem::path model_dir {basename};
  model_dir /= "python";
  filesystem::create_directories(model_dir);

  string filename = model_dir.string() + "/static_params_derivs.py";
  ofstream output {filename, ios::out | ios::binary};
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << writePythonHelpersImports(use_jax, use_numba) << endl
         << writePythonHelpersImportFromInit() << endl
         << endl;

  // Write main function
  output << "def static_params_derivs(y: np.ndarray, x: np.ndarray, params: np.ndarray) -> "
            "Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, "
            "np.ndarray]:"
         << endl
         << "    \"\"\"" << endl
         << "    Compute derivatives of the static model with respect to parameters." << endl
         << endl
         << "    Parameters" << endl
         << "    ----------" << endl
         << "    y : array_like" << endl
         << "        Endogenous variables" << endl
         << "    x : array_like" << endl
         << "        Exogenous variables" << endl
         << "    params : array_like" << endl
         << "        Model parameters" << endl
         << endl
         << "    Returns" << endl
         << "    -------" << endl
         << "    rp : ndarray" << endl
         << "        Jacobian of static model residuals w.r.t. parameters" << endl
         << "        Sparse format: (row, col, value) tuples" << endl
         << "    g1p : ndarray" << endl
         << "        Derivative of Jacobian w.r.t. parameters" << endl
         << "        Format: (eq, var, param, value)" << endl
         << "    rpp : ndarray" << endl
         << "        Second derivative of residuals w.r.t. parameters" << endl
         << "        Format: (eq, param1, param2, value)" << endl
         << "    g1pp : ndarray" << endl
         << "        Second derivative of Jacobian w.r.t. parameters" << endl
         << "        Format: (eq, var, param1, param2, value)" << endl
         << "    g2p : ndarray" << endl
         << "        Derivative of Hessian w.r.t. parameters" << endl
         << "        Format: (eq, var1, var2, param, value)" << endl
         << "    \"\"\"" << endl;

  // Temporary terms
  if (!params_derivs_temporary_terms_idxs.empty())
    output << "    T = " << (use_jax ? "jnp" : "np") << ".zeros("
           << params_derivs_temporary_terms_idxs.size() << ")" << endl;

  // Write temporary terms (already indented by writeTemporaryTerms)
  if (!tt_output.str().empty())
    output << tt_output.str();

  // Count non-zero elements for each derivative type
  int rp_count = params_derivatives.at({0, 1}).size();
  int g1p_count = params_derivatives.at({1, 1}).size();
  int rpp_count = params_derivatives.at({0, 2}).size();
  int g1pp_count = params_derivatives.at({1, 2}).size();
  int g2p_count = params_derivatives.at({2, 1}).size();

  // Initialize arrays
  output << endl;
  if (rp_count > 0)
    output << "    rp_i = " << (use_jax ? "jnp" : "np") << ".zeros(" << rp_count << ", dtype=int)"
           << endl
           << "    rp_j = " << (use_jax ? "jnp" : "np") << ".zeros(" << rp_count << ", dtype=int)"
           << endl
           << "    rp_v = " << (use_jax ? "jnp" : "np") << ".zeros(" << rp_count << ")" << endl;

  if (g1p_count > 0)
    output << "    g1p = " << (use_jax ? "jnp" : "np") << ".zeros((" << g1p_count << ", 4))"
           << endl;

  if (rpp_count > 0)
    output << "    rpp = " << (use_jax ? "jnp" : "np") << ".zeros((" << rpp_count << ", 4))"
           << endl;

  if (g1pp_count > 0)
    output << "    g1pp = " << (use_jax ? "jnp" : "np") << ".zeros((" << g1pp_count << ", 5))"
           << endl;

  if (g2p_count > 0)
    output << "    g2p = " << (use_jax ? "jnp" : "np") << ".zeros((" << g2p_count << ", 5))"
           << endl;

  output << endl;

  // Write derivatives with proper indentation
  auto write_python_derivs = [&output](ostringstream& deriv_output) {
    if (!deriv_output.str().empty())
      {
        istringstream deriv_stream(deriv_output.str());
        string line;
        while (getline(deriv_stream, line))
          output << "    " << line << endl;
      }
  };

  write_python_derivs(rp_output);
  write_python_derivs(g1p_output);
  write_python_derivs(rpp_output);
  write_python_derivs(g1pp_output);
  write_python_derivs(g2p_output);

  // Return statement
  output << endl << "    return rp_i, rp_j, rp_v, g1p, rpp, g1pp, g2p" << endl;

  output.close();
}
