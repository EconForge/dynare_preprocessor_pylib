/*
 * Copyright © 2003-2020 Dynare Team
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _STATIC_MODEL_HH
#define _STATIC_MODEL_HH

using namespace std;

#include <fstream>
#include <filesystem>

#include "ModelTree.hh"

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

  //! Writes the static model equations and its derivatives
  void writeStaticModel(const string &basename, ostream &StaticOutput, bool use_dll, bool julia) const;

  //! Writes the main static function of block decomposed model (MATLAB version)
  void writeStaticBlockMFile(const string &basename) const;

  //! Writes the per-block static files of block decomposed model (MATLAB version)
  void writeStaticPerBlockMFiles(const string &basename) const;

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

  //! Write derivative code of an equation w.r. to a variable
  void compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const;
  //! Write chain rule derivative code of an equation w.r. to a variable
  void compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int blk, int eq, int var, int lag, const temporary_terms_t &temporary_terms, const temporary_terms_idxs_t &temporary_terms_idxs) const;

  //! Get the type corresponding to a derivation ID
  SymbolType getTypeByDerivID(int deriv_id) const noexcept(false) override;
  //! Get the lag corresponding to a derivation ID
  int getLagByDerivID(int deriv_id) const noexcept(false) override;
  //! Get the symbol ID corresponding to a derivation ID
  int getSymbIDByDerivID(int deriv_id) const noexcept(false) override;
  //! Compute the column indices of the static Jacobian
  void computeStatJacobianCols();
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

  void writeStaticModel(ostream &DynamicOutput, bool use_dll, bool julia) const;
  void writeStaticModel(const string &dynamic_basename, bool use_dll, bool julia) const;

  //! Internal helper for the copy constructor and assignment operator
  /*! Copies all the structures that contain ExprNode*, by the converting the
      pointers into their equivalent in the new tree */
  void copyHelper(const StaticModel &m);

public:
  StaticModel(SymbolTable &symbol_table_arg,
              NumericalConstants &num_constants,
              ExternalFunctionsTable &external_functions_table_arg);

  StaticModel(const StaticModel &m);
  StaticModel(StaticModel &&) = delete;
  StaticModel &operator=(const StaticModel &m);
  /* The move assignment operator is not explicitly deleted, otherwise the
     static_cast from DynamicModel does not work. However it looks like this
     operator will not be used in the end. See
     https://en.cppreference.com/w/cpp/language/copy_initialization
     With C++17, it should be possible to explicitly delete it */
  //StaticModel & operator=(StaticModel &&) = delete;

  //! Creates the static version of a dynamic model
  explicit StaticModel(const DynamicModel &m);

  //! Writes information on block decomposition when relevant
  void writeOutput(ostream &output, bool block) const;

  //! Execute computations (variable sorting + derivation)
  /*!
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the static files
    \param derivsOrder order of derivation with respect to endogenous
    \param paramsDerivsOrder order of derivatives w.r. to a pair (endogenous, parameter) to be computed
  */
  void computingPass(int derivsOrder, int paramsDerivsOrder, const eval_context_t &eval_context, bool no_tmp_terms, bool block, bool bytecode);

  //! Writes static model file
  void writeStaticFile(const string &basename, bool block, bool bytecode, bool use_dll, const string &mexext, const filesystem::path &matlabroot, const filesystem::path &dynareroot, bool julia) const;

  //! Write JSON Output (used by PlannerObjectiveStatement)
  void writeJsonOutput(ostream &output) const;

  //! Write JSON representation of static model
  void writeJsonComputingPassOutput(ostream &output, bool writeDetails) const;

  //! Writes file containing static parameters derivatives
  void writeJsonParamsDerivativesFile(ostream &output, bool writeDetails) const;

  //! Writes file containing static parameters derivatives
  void writeParamsDerivativesFile(const string &basename, bool julia) const;

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

#endif
