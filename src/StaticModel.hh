/*
 * Copyright © 2003-2023 Dynare Team
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
  /* First-order derivatives of equations w.r.t. Lagrange multipliers, using
     chain rule derivation for auxiliary variables added after the multipliers
     (so that derivatives of optimality FOCs w.r.t. multipliers with lead or
     lag ⩾ 2 are self-contained, which is required by dyn_ramsey_static.m).
     Only used if 'ramsey_model' or 'ramsey_policy' is present.
     The first index of the key is the equation number (NB: auxiliary equations
     added after the multipliers do not appear).
     The second index is the index of the Lagrange multiplier (ordered by
     increasing symbol ID) */
  SparseColumnMajorOrderMatrix ramsey_multipliers_derivatives;
  /* Column indices for the derivatives w.r.t. Lagrange multipliers in
     Compressed Sparse Column (CSC) storage (corresponds to the “jc” vector in
     MATLAB terminology) */
  vector<int> ramsey_multipliers_derivatives_sparse_colptr;
  // Temporary terms for ramsey_multipliers_derivatives
  temporary_terms_t ramsey_multipliers_derivatives_temporary_terms;
  // Stores, for each temporary term, its index in the MATLAB/Octave vector
  temporary_terms_idxs_t ramsey_multipliers_derivatives_temporary_terms_idxs;

  // Writes static model file (MATLAB/Octave version, legacy representation)
  void writeStaticMFile(const string &basename) const;

  //! Writes the code of the block-decomposed model in virtual machine bytecode
  void writeStaticBlockBytecode(const string &basename) const;

  //! Writes the code of the model in virtual machine bytecode
  void writeStaticBytecode(const string &basename) const;

  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  void evaluateJacobian(const eval_context_t &eval_context, jacob_map_t *j_m, bool dynamic);

  SymbolType getTypeByDerivID(int deriv_id) const noexcept(false) override;
  int getLagByDerivID(int deriv_id) const noexcept(false) override;
  int getSymbIDByDerivID(int deriv_id) const noexcept(false) override;
  int getTypeSpecificIDByDerivID(int deriv_id) const override;

  int
  getJacobianCol(int deriv_id, [[maybe_unused]] bool sparse) const override
  {
    return getTypeSpecificIDByDerivID(deriv_id);
  }
  int
  getJacobianColsNbr([[maybe_unused]] bool sparse) const override
  {
    return symbol_table.endo_nbr();
  }

  void computeChainRuleJacobian() override;

  /* Helper for writing MATLAB/Octave functions for residuals/derivatives and
     their temporary terms (legacy representation) */
  void writeStaticMFileHelper(const string &basename,
                              const string &name, const string &retvalname,
                              const string &name_tt, size_t ttlen,
                              const string &previous_tt_name,
                              const ostringstream &init_s, const ostringstream &end_s,
                              const ostringstream &s, const ostringstream &s_tt) const;
  /* Writes MATLAB/Octave wrapper function for computing residuals and
     derivatives at the same time (legacy representation) */
  void writeStaticMWrapperFunction(const string &basename, const string &ending) const;

  /* Create the compatibility static.m file for MATLAB/Octave not yet using the
     temporary terms array interface (legacy representation) */
  void writeStaticMCompatFile(const string &name) const;

  int
  getBlockJacobianEndoCol([[maybe_unused]] int blk, int var, [[maybe_unused]] int lag) const override
  {
    assert(var >= blocks[blk].getRecursiveSize());
    return var - blocks[blk].getRecursiveSize();
  }

  // Write the block structure of the model in the driver file
  void writeBlockDriverOutput(ostream &output) const;

  // Helper for writing ramsey_multipliers_derivatives
  template<ExprNodeOutputType output_type>
  void writeRamseyMultipliersDerivativesHelper(ostream &output) const;

protected:
  string
  modelClassName() const override
  {
    return "static model";
  }

public:
  StaticModel(SymbolTable &symbol_table_arg,
              NumericalConstants &num_constants,
              ExternalFunctionsTable &external_functions_table_arg);

  StaticModel(const StaticModel &m);
  StaticModel &operator=(const StaticModel &m);

  //! Creates the static version of a dynamic model
  explicit StaticModel(const DynamicModel &m);

  //! Writes information about the static model to the driver file
  void writeDriverOutput(ostream &output) const;

  //! Execute computations (variable sorting + derivation + block decomposition)
  /*!
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the static files
    \param derivsOrder order of derivation with respect to endogenous
    \param paramsDerivsOrder order of derivatives w.r. to a pair (endogenous, parameter) to be computed
  */
  void computingPass(int derivsOrder, int paramsDerivsOrder, const eval_context_t &eval_context, bool no_tmp_terms, bool block, bool use_dll);

  //! Writes static model file (+ bytecode)
  void writeStaticFile(const string &basename, bool use_dll, const string &mexext, const filesystem::path &matlabroot, bool julia) const;

  //! Write JSON Output (used by PlannerObjectiveStatement)
  void writeJsonOutput(ostream &output) const;

  //! Write JSON representation of static model
  void writeJsonComputingPassOutput(ostream &output, bool writeDetails) const;

  //! Write JSON params derivatives
  void writeJsonParamsDerivatives(ostream &output, bool writeDetails) const;

  //! Writes file containing static parameters derivatives
  template<bool julia>
  void writeParamsDerivativesFile(const string &basename) const;

  //! Writes LaTeX file with the equations of the static model
  void writeLatexFile(const string &basename, bool write_equation_tags) const;

  //! Writes initializations in oo_.steady_state or steady state file for the auxiliary variables
  void writeAuxVarInitval(ostream &output, ExprNodeOutputType output_type) const;

  void writeLatexAuxVarRecursiveDefinitions(ostream &output) const;
  void writeJsonAuxVarRecursiveDefinitions(ostream &output) const;

  //! To ensure that no exogenous is present in the planner objective
  //! See #1264
  bool exoPresentInEqs() const;

  int getDerivID(int symb_id, int lag) const noexcept(false) override;
  void addAllParamDerivId(set<int> &deriv_id_set) override;

  // Fills the ramsey_multipliers_derivatives structure (see the comment there)
  void computeRamseyMultipliersDerivatives(int ramsey_orig_endo_nbr, bool is_matlab, bool no_tmp_terms);

  // Writes the sparse indices of ramsey_multipliers_derivatives to the driver file
  void writeDriverRamseyMultipliersDerivativesSparseIndices(ostream &output) const;

  // Writes ramsey_multipliers_derivatives (MATLAB/Octave version)
  void writeRamseyMultipliersDerivativesMFile(const string &basename, int ramsey_orig_endo_nbr) const;

  // Writes ramsey_multipliers_derivatives (C version)
  void writeRamseyMultipliersDerivativesCFile(const string &basename, const string &mexext, const filesystem::path &matlabroot, int ramsey_orig_endo_nbr) const;
};

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

  if constexpr(!julia)
    {
      filesystem::path filename {packageDir(basename) / "static_params_derivs.m"};
      ofstream paramsDerivsFile {filename, ios::out | ios::binary};
      if (!paramsDerivsFile.is_open())
        {
          cerr << "ERROR: Can't open file " << filename.string() << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
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
                       << "rpp = zeros(" << params_derivatives.at({ 0, 2 }).size() << ",4);" << endl
                       << rpp_output.str()
                       << "gpp = zeros(" << params_derivatives.at({ 1, 2 }).size() << ",5);" << endl
                       << gpp_output.str()
                       << "end" << endl
                       << "if nargout >= 5" << endl
                       << "hp = zeros(" << params_derivatives.at({ 2, 1 }).size() << ",5);" << endl
                       << hp_output.str()
                       << "end" << endl
                       << "end" << endl;
      paramsDerivsFile.close();
    }
  else
    {
      stringstream output;
      output << "# NB: this file was automatically generated by Dynare" << endl
             << "#     from " << basename << ".mod" << endl
             << "#" << endl
             << "function static_params_derivs(y, x, params)" << endl
             << "@inbounds begin" << endl
             << tt_output.str()
             << "rp = zeros(" << equations.size() << ", "
             << symbol_table.param_nbr() << ");" << endl
             << rp_output.str()
             << "gp = zeros(" << equations.size() << ", " << symbol_table.endo_nbr() << ", "
             << symbol_table.param_nbr() << ");" << endl
             << gp_output.str()
             << "rpp = zeros(" << params_derivatives.at({ 0, 2 }).size() << ",4);" << endl
             << rpp_output.str()
             << "gpp = zeros(" << params_derivatives.at({ 1, 2 }).size() << ",5);" << endl
             << gpp_output.str()
             << "hp = zeros(" << params_derivatives.at({ 2, 1 }).size() << ",5);" << endl
             << hp_output.str()
             << "end" << endl
             << "return (rp, gp, rpp, gpp, hp)" << endl
             << "end" << endl;

      writeToFileIfModified(output, filesystem::path{basename} / "model" / "julia" / "StaticParamsDerivs.jl");
    }
}

template<ExprNodeOutputType output_type>
void
StaticModel::writeRamseyMultipliersDerivativesHelper(ostream &output) const
{
  // Write temporary terms (which includes external function stuff)
  deriv_node_temp_terms_t tef_terms;
  temporary_terms_t unused_tt_copy;
  writeTemporaryTerms<output_type>(ramsey_multipliers_derivatives_temporary_terms,
                                   unused_tt_copy,
                                   ramsey_multipliers_derivatives_temporary_terms_idxs,
                                   output, tef_terms);

  // Write chain rule derivatives
  for (int k {0};
       auto &[row_col, d] : ramsey_multipliers_derivatives)
    {
      output << "g1m_v" << LEFT_ARRAY_SUBSCRIPT(output_type)
             << k + ARRAY_SUBSCRIPT_OFFSET(output_type)
             << RIGHT_ARRAY_SUBSCRIPT(output_type) << "=";
      d->writeOutput(output, output_type, ramsey_multipliers_derivatives_temporary_terms,
                     ramsey_multipliers_derivatives_temporary_terms_idxs, tef_terms);
      output << ";" << endl;
      k++;
    }
}

#endif
