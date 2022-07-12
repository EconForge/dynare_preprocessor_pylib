/*
 * Copyright © 2003-2022 Dynare Team
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

#ifndef _STATIC_MODEL_HH
#define _STATIC_MODEL_HH

#include <fstream>
#include <filesystem>

#include "ModelTree.hh"
#include "Bytecode.hh"

using namespace std;

class DynamicModel;

//! Stores a static model, as derived from the "model" block when leads and lags have been removed
class StaticModel : public ModelTree
{
private:
  //! Writes static model file (standard Matlab version)
  void writeStaticMFile(const string &basename) const;

  //! Writes static model file (C version)
  void writeStaticCFile(const string &basename) const;

  //! Writes static model file (Julia version)
  void writeStaticJuliaFile(const string &basename) const;

  //! Writes the main static function of block decomposed model (MATLAB version)
  void writeStaticBlockMFile(const string &basename) const;

  //! Writes the main static function of block decomposed model (C version)
  void writeStaticBlockCFile(const string &basename) const;

  //! Helper for writing a per-block static file of block decomposed model
  template<ExprNodeOutputType output_type>
  void writeStaticPerBlockHelper(int blk, ostream &output, temporary_terms_t &temporary_terms) const;

  //! Writes the per-block static files of block decomposed model (MATLAB version)
  void writeStaticPerBlockMFiles(const string &basename) const;

  //! Writes the per-block static files of block decomposed model (C version)
  void writeStaticPerBlockCFiles(const string &basename) const;

  //! Writes the code of the block-decomposed model in virtual machine bytecode
  void writeStaticBlockBytecode(const string &basename) const;

  //! Writes the code of the model in virtual machine bytecode
  void writeStaticBytecode(const string &basename) const;

  //! Adds per-block information for bytecode simulation in a separate .bin file
  void writeBlockBytecodeBinFile(const string &basename, int num,
                                 int &u_count_int, bool &file_open) const;

  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  void evaluateJacobian(const eval_context_t &eval_context, jacob_map_t *j_m, bool dynamic);

  //! Write derivative bytecode of an equation w.r. to a variable
  void writeBytecodeDerivative(BytecodeWriter &code_file, int eq, int symb_id, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const;
  //! Write chain rule derivative bytecode of an equation w.r. to a variable
  void writeBytecodeChainRuleDerivative(BytecodeWriter &code_file, int blk, int eq, int var, int lag, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs, const deriv_node_temp_terms_t &tef_terms) const;

  //! Get the type corresponding to a derivation ID
  SymbolType getTypeByDerivID(int deriv_id) const noexcept(false) override;
  //! Get the lag corresponding to a derivation ID
  int getLagByDerivID(int deriv_id) const noexcept(false) override;
  //! Get the symbol ID corresponding to a derivation ID
  int getSymbIDByDerivID(int deriv_id) const noexcept(false) override;

  int
  getJacobianCol(int deriv_id) const override
  {
    return symbol_table.getTypeSpecificID(getSymbIDByDerivID(deriv_id));
  }
  int
  getJacobianColsNbr() const override
  {
    return symbol_table.endo_nbr();
  }

  //! Computes chain rule derivatives of the Jacobian w.r. to endogenous variables
  void computeChainRuleJacobian();

  //! Helper functions for writeStaticModel
  void writeStaticModelHelper(const string &basename,
                              const string &name, const string &retvalname,
                              const string &name_tt, size_t ttlen,
                              const string &previous_tt_name,
                              const ostringstream &init_s, const ostringstream &end_s,
                              const ostringstream &s, const ostringstream &s_tt) const;
  void writeWrapperFunctions(const string &basename, const string &ending) const;

  //! Create a legacy *_static.m file for Matlab/Octave not yet using the temporary terms array interface
  void writeStaticMatlabCompatLayer(const string &name) const;

public:
  StaticModel(SymbolTable &symbol_table_arg,
              NumericalConstants &num_constants,
              ExternalFunctionsTable &external_functions_table_arg);

  StaticModel(const StaticModel &m);
  StaticModel &operator=(const StaticModel &m);

  //! Creates the static version of a dynamic model
  explicit StaticModel(const DynamicModel &m);

  //! Writes information about the static model to the driver file
  void writeDriverOutput(ostream &output, bool block) const;

  //! Execute computations (variable sorting + derivation)
  /*!
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the static files
    \param derivsOrder order of derivation with respect to endogenous
    \param paramsDerivsOrder order of derivatives w.r. to a pair (endogenous, parameter) to be computed
  */
  void computingPass(int derivsOrder, int paramsDerivsOrder, const eval_context_t &eval_context, bool no_tmp_terms, bool block);

  //! Writes static model file (+ bytecode)
  void writeStaticFile(const string &basename, bool block, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const;

  //! Write JSON Output (used by PlannerObjectiveStatement)
  void writeJsonOutput(ostream &output) const;

  //! Write JSON representation of static model
  void writeJsonComputingPassOutput(ostream &output, bool writeDetails) const;

  //! Writes file containing static parameters derivatives
  void writeJsonParamsDerivativesFile(ostream &output, bool writeDetails) const;

  //! Writes file containing static parameters derivatives
  template<bool julia>
  void writeParamsDerivativesFile(const string &basename) const;

  //! Writes LaTeX file with the equations of the static model
  void writeLatexFile(const string &basename, bool write_equation_tags) const;

  //! Writes initializations in oo_.steady_state or steady state file for the auxiliary variables
  void writeAuxVarInitval(ostream &output, ExprNodeOutputType output_type) const;

  //! Writes definition of the auxiliary variables in a .m or .jl file
  void writeSetAuxiliaryVariables(const string &basename, bool julia) const;
  void writeAuxVarRecursiveDefinitions(ostream &output, ExprNodeOutputType output_type) const;
  void writeLatexAuxVarRecursiveDefinitions(ostream &output) const;
  void writeJsonAuxVarRecursiveDefinitions(ostream &output) const;

  //! To ensure that no exogenous is present in the planner objective
  //! See #1264
  bool exoPresentInEqs() const;

  int getDerivID(int symb_id, int lag) const noexcept(false) override;
  void addAllParamDerivId(set<int> &deriv_id_set) override;
};

template<ExprNodeOutputType output_type>
void
StaticModel::writeStaticPerBlockHelper(int blk, ostream &output, temporary_terms_t &temporary_terms) const
{
  BlockSimulationType simulation_type { blocks[blk].simulation_type };
  int block_recursive_size { blocks[blk].getRecursiveSize() };

  // The equations
  deriv_node_temp_terms_t tef_terms;

  auto write_eq_tt = [&](int eq)
                     {
                       for (auto it : blocks_temporary_terms[blk][eq])
                         {
                           if (dynamic_cast<AbstractExternalFunctionNode *>(it))
                             it->writeExternalFunctionOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs, tef_terms);

                           output << "  ";
                           it->writeOutput(output, output_type, blocks_temporary_terms[blk][eq], blocks_temporary_terms_idxs, tef_terms);
                           output << '=';
                           it->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs, tef_terms);
                           temporary_terms.insert(it);
                           output << ';' << endl;
                         }
                     };

  for (int eq {0}; eq < blocks[blk].size; eq++)
    {
      write_eq_tt(eq);

      EquationType equ_type { getBlockEquationType(blk, eq) };
      BinaryOpNode *e { getBlockEquationExpr(blk, eq) };
      expr_t lhs { e->arg1 }, rhs { e->arg2 };
      switch (simulation_type)
        {
        case BlockSimulationType::evaluateBackward:
        case BlockSimulationType::evaluateForward:
          evaluation:
          if (equ_type == EquationType::evaluateRenormalized)
            {
              e = getBlockEquationRenormalizedExpr(blk, eq);
              lhs = e->arg1;
              rhs = e->arg2;
            }
          else if (equ_type != EquationType::evaluate)
            {
              cerr << "Type mismatch for equation " << getBlockEquationID(blk, eq)+1  << endl;
              exit(EXIT_FAILURE);
            }
          output << "  ";
          lhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << '=';
          rhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ';' << endl;
          break;
        case BlockSimulationType::solveBackwardSimple:
        case BlockSimulationType::solveForwardSimple:
        case BlockSimulationType::solveBackwardComplete:
        case BlockSimulationType::solveForwardComplete:
          if (eq < block_recursive_size)
            goto evaluation;
          output << "  residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq-block_recursive_size+ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=(";
          lhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ")-(";
          rhs->writeOutput(output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          output << ");" << endl;
          break;
        default:
          cerr << "Incorrect type for block " << blk+1 << endl;
          exit(EXIT_FAILURE);
        }
    }
  // The Jacobian if we have to solve the block
  if (simulation_type != BlockSimulationType::evaluateBackward
      && simulation_type != BlockSimulationType::evaluateForward)
    {
      // Write temporary terms for derivatives
      write_eq_tt(blocks[blk].size);

      ostringstream i_output, j_output, v_output;
      for (int line_counter { ARRAY_SUBSCRIPT_OFFSET(output_type) };
           const auto &[indices, d] : blocks_derivatives[blk])
        {
          const auto &[eq, var, ignore] {indices};
          i_output << "  g1_i" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << eq+1-block_recursive_size
                   << ';' << endl;
          j_output << "  g1_j" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=' << var+1-block_recursive_size
                   << ';' << endl;
          v_output << "  g1_v" << LEFT_ARRAY_SUBSCRIPT(output_type) << line_counter
                   << RIGHT_ARRAY_SUBSCRIPT(output_type) << '=';
          d->writeOutput(v_output, output_type, temporary_terms, blocks_temporary_terms_idxs);
          v_output << ';' << endl;
          line_counter++;
        }
      output << i_output.str() << j_output.str() << v_output.str();
    }
}

template<bool julia>
void
StaticModel::writeParamsDerivativesFile(const string &basename) const
{
  if (!params_derivatives.size())
    return;

  constexpr ExprNodeOutputType output_type { julia ? ExprNodeOutputType::juliaStaticModel : ExprNodeOutputType::matlabStaticModel };

  auto [tt_output, rp_output, gp_output, rpp_output, gpp_output, hp_output, g3p_output]
    { writeParamsDerivativesFileHelper<output_type>() };
  // g3p_output is ignored

  string filename { julia ? basename + "StaticParamsDerivs.jl" : packageDir(basename) + "/static_params_derivs.m" };
  ofstream paramsDerivsFile { filename, ios::out | ios::binary };
  if (!paramsDerivsFile.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  if constexpr(!julia)
    {
      paramsDerivsFile << "function [rp, gp, rpp, gpp, hp] = static_params_derivs(y, x, params)" << endl
                       << "%" << endl
                       << "% Status : Computes derivatives of the static model with respect to the parameters" << endl
                       << "%" << endl
                       << "% Inputs : " << endl
                       << "%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order" << endl
                       << "%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order" << endl
                       << "%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order" << endl
                       << "%" << endl
                       << "% Outputs:" << endl
                       << "%   rp        [M_.eq_nbr by #params] double    Jacobian matrix of static model equations with respect to parameters " << endl
                       << "%                                              Dynare may prepend or append auxiliary equations, see M_.aux_vars" << endl
                       << "%   gp        [M_.endo_nbr by M_.endo_nbr by #params] double    Derivative of the Jacobian matrix of the static model equations with respect to the parameters" << endl
                       << "%                                                           rows: variables in declaration order" << endl
                       << "%                                                           rows: equations in order of declaration" << endl
                       << "%   rpp       [#second_order_residual_terms by 4] double   Hessian matrix of second derivatives of residuals with respect to parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: number of the first parameter in derivative" << endl
                       << "%                                                              3rd column: number of the second parameter in derivative" << endl
                       << "%                                                              4th column: value of the Hessian term" << endl
                       << "%   gpp      [#second_order_Jacobian_terms by 5] double   Hessian matrix of second derivatives of the Jacobian with respect to the parameters;" << endl
                       << "%                                                              rows: respective derivative term" << endl
                       << "%                                                              1st column: equation number of the term appearing" << endl
                       << "%                                                              2nd column: column number of variable in Jacobian of the static model" << endl
                       << "%                                                              3rd column: number of the first parameter in derivative" << endl
                       << "%                                                              4th column: number of the second parameter in derivative" << endl
                       << "%                                                              5th column: value of the Hessian term" << endl
                       << "%" << endl
                       << "%" << endl
                       << "% Warning : this file is generated automatically by Dynare" << endl
                       << "%           from model file (.mod)" << endl << endl
                       << "T = NaN(" << params_derivs_temporary_terms_idxs.size() << ",1);" << endl
                       << tt_output.str()
                       << "rp = zeros(" << equations.size() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << rp_output.str()
                       << "gp = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ", "
                       << symbol_table.param_nbr() << ");" << endl
                       << gp_output.str()
                       << "if nargout >= 3" << endl
                       << "rpp = zeros(" << params_derivatives.find({ 0, 2 })->second.size() << ",4);" << endl
                       << rpp_output.str()
                       << "gpp = zeros(" << params_derivatives.find({ 1, 2 })->second.size() << ",5);" << endl
                       << gpp_output.str()
                       << "end" << endl
                       << "if nargout >= 5" << endl
                       << "hp = zeros(" << params_derivatives.find({ 2, 1 })->second.size() << ",5);" << endl
                       << hp_output.str()
                       << "end" << endl
                       << "end" << endl;
    }
  else
    paramsDerivsFile << "module " << basename << "StaticParamsDerivs" << endl
                     << "#" << endl
                     << "# NB: this file was automatically generated by Dynare" << endl
                     << "#     from " << basename << ".mod" << endl
                     << "#" << endl
                     << "export params_derivs" << endl << endl
                     << "function params_derivs(y, x, params)" << endl
		     << "@inbounds begin" << endl
                     << tt_output.str()
		     << "end" << endl
                     << "rp = zeros(" << equations.size() << ", "
                     << symbol_table.param_nbr() << ");" << endl
		     << "@inbounds begin" << endl
                     << rp_output.str()
		     << "end" << endl
                     << "gp = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ", "
                     << symbol_table.param_nbr() << ");" << endl
		     << "@inbounds begin" << endl
                     << gp_output.str()
		     << "end" << endl
                     << "rpp = zeros(" << params_derivatives.find({ 0, 2 })->second.size() << ",4);" << endl
		     << "@inbounds begin" << endl
                     << rpp_output.str()
		     << "end" << endl
                     << "gpp = zeros(" << params_derivatives.find({ 1, 2 })->second.size() << ",5);" << endl
		     << "@inbounds begin" << endl
                     << gpp_output.str()
		     << "end" << endl
                     << "hp = zeros(" << params_derivatives.find({ 2, 1 })->second.size() << ",5);" << endl
		     << "@inbounds begin" << endl
                     << hp_output.str()
		     << "end" << endl
                     << "(rp, gp, rpp, gpp, hp)" << endl
                     << "end" << endl
                     << "end" << endl;

  paramsDerivsFile.close();
}

#endif
