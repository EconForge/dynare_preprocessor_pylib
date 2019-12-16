/*
 * Copyright © 2003-2019 Dynare Team
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
  //! global temporary terms for block decomposed models
  vector<vector<temporary_terms_t>> v_temporary_terms;

  //! local temporary terms for block decomposed models
  vector<vector<temporary_terms_t>> v_temporary_terms_local;

  vector<temporary_terms_inuse_t> v_temporary_terms_inuse;

  using first_chain_rule_derivatives_t = map<tuple<int, int, int>, expr_t>;
  first_chain_rule_derivatives_t first_chain_rule_derivatives;

  //! Writes static model file (standard Matlab version)
  void writeStaticMFile(const string &basename) const;

  //! Writes static model file (C version)
  void writeStaticCFile(const string &basename) const;

  //! Writes static model file (Julia version)
  void writeStaticJuliaFile(const string &basename) const;

  //! Writes the static model equations and its derivatives
  void writeStaticModel(const string &basename, ostream &StaticOutput, bool use_dll, bool julia) const;

  //! Writes the static function calling the block to solve (Matlab version)
  void writeStaticBlockMFSFile(const string &basename) const;

  //! Writes the Block reordred structure of the model in M output
  void writeModelEquationsOrdered_M(const string &basename) const;

  //! Writes the code of the Block reordred structure of the model in virtual machine bytecode
  void writeModelEquationsCode_Block(const string &basename, map_idx_t map_idx, vector<map_idx_t> map_idx2) const;

  //! Writes the code of the model in virtual machine bytecode
  void writeModelEquationsCode(const string &basename, map_idx_t map_idx) const;

  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  void evaluateJacobian(const eval_context_t &eval_context, jacob_map_t *j_m, bool dynamic);

  map_idx_t map_idx;

  vector<map_idx_t> map_idx2;

  //! sorts the temporary terms in the blocks order
  void computeTemporaryTermsOrdered();
  //! creates a mapping from the index of temporary terms to a natural index
  void computeTemporaryTermsMapping(temporary_terms_t &temporary_terms, map_idx_t &map_idx);

  //! Write derivative code of an equation w.r. to a variable
  void compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, map_idx_t &map_idx, temporary_terms_t temporary_terms) const;
  //! Write chain rule derivative code of an equation w.r. to a variable
  void compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int var, int lag, map_idx_t &map_idx, temporary_terms_t temporary_terms) const;

  //! Get the type corresponding to a derivation ID
  SymbolType getTypeByDerivID(int deriv_id) const noexcept(false) override;
  //! Get the lag corresponding to a derivation ID
  int getLagByDerivID(int deriv_id) const noexcept(false) override;
  //! Get the symbol ID corresponding to a derivation ID
  int getSymbIDByDerivID(int deriv_id) const noexcept(false) override;
  //! Compute the column indices of the static Jacobian
  void computeStatJacobianCols();
  //! return a map on the block jacobian
  map<tuple<int, int, int, int, int>, int> get_Derivatives(int block);
  //! Computes chain rule derivatives of the Jacobian w.r. to endogenous variables
  void computeChainRuleJacobian(blocks_derivatives_t &blocks_derivatives);
  //! Collect only the first derivatives
  map<tuple<int, int, int>, expr_t> collect_first_order_derivatives_endogenous();

  //! Collecte the derivatives w.r. to endogenous of the block, to endogenous of previouys blocks and to exogenous
  void collect_block_first_order_derivatives();

  //! Indicate if the temporary terms are computed for the overall model (true) or not (false). Default value true
  bool global_temporary_terms{true};

  //! Vector describing equations: BlockSimulationType, if BlockSimulationType == EVALUATE_s then a expr_t on the new normalized equation
  equation_type_and_normalized_equation_t equation_type_and_normalized_equation;

  //! for each block contains pair< Simulation_Type, pair < Block_Size, Recursive_part_Size >>
  block_type_firstequation_size_mfs_t block_type_firstequation_size_mfs;

  //! for all blocks derivatives description
  blocks_derivatives_t blocks_derivatives;

  //! The jacobian without the elements below the cutoff
  dynamic_jacob_map_t dynamic_jacobian;

  //! Vector indicating if the block is linear in endogenous variable (true) or not (false)
  vector<bool> blocks_linear;

  //! Map the derivatives for a block tuple<lag, eq, var>
  using derivative_t = map<tuple<int, int, int>, expr_t>;
  //! Vector of derivative for each blocks
  vector<derivative_t> derivative_endo, derivative_other_endo, derivative_exo, derivative_exo_det;

  //!List for each block and for each lag-leag all the other endogenous variables and exogenous variables
  using var_t = set<int>;
  using lag_var_t = map<int, var_t>;

  //! for each block described the number of static, forward, backward and mixed variables in the block
  /*! tuple<static, forward, backward, mixed> */
  vector<tuple<int, int, int, int>> block_col_type;

  //!Maximum lead and lag for each block on endogenous of the block, endogenous of the previous blocks, exogenous and deterministic exogenous
  vector<pair<int, int>> endo_max_leadlag_block, other_endo_max_leadlag_block, exo_max_leadlag_block, exo_det_max_leadlag_block, max_leadlag_block;

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
  StaticModel & operator=(const StaticModel &m);
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

  //! Adds informations for simulation in a binary file for a block decomposed model
  void Write_Inf_To_Bin_File_Block(const string &basename, int num,
                                   int &u_count_int, bool &file_open) const;

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

  //! Return the number of blocks
  unsigned int
  getNbBlocks() const override
  {
    return (block_type_firstequation_size_mfs.size());
  };
  //! Determine the simulation type of each block
  BlockSimulationType
  getBlockSimulationType(int block_number) const override
  {
    return (get<0>(block_type_firstequation_size_mfs[block_number]));
  };
  //! Return the first equation number of a block
  unsigned int
  getBlockFirstEquation(int block_number) const override
  {
    return (get<1>(block_type_firstequation_size_mfs[block_number]));
  };
  //! Return the size of the block block_number
  unsigned int
  getBlockSize(int block_number) const override
  {
    return (get<2>(block_type_firstequation_size_mfs[block_number]));
  };
  //! Return the number of exogenous variable in the block block_number
  unsigned int
  getBlockExoSize(int block_number) const override
  {
    return 0;
  };
  //! Return the number of colums in the jacobian matrix for exogenous variable in the block block_number
  unsigned int
  getBlockExoColSize(int block_number) const override
  {
    return 0;
  }
  //! Return the number of feedback variable of the block block_number
  unsigned int
  getBlockMfs(int block_number) const override
  {
    return (get<3>(block_type_firstequation_size_mfs[block_number]));
  };
  //! Return the maximum lag in a block
  unsigned int
  getBlockMaxLag(int block_number) const override
  {
    return (block_lag_lead[block_number].first);
  };
  //! Return the maximum lead in a block
  unsigned int
  getBlockMaxLead(int block_number) const override
  {
    return (block_lag_lead[block_number].second);
  };
  //! Return the type of equation (equation_number) belonging to the block block_number
  EquationType
  getBlockEquationType(int block_number, int equation_number) const override
  {
    return (equation_type_and_normalized_equation[equation_reordered[get<1>(block_type_firstequation_size_mfs[block_number])+equation_number]].first);
  };
  //! Return true if the equation has been normalized
  bool
  isBlockEquationRenormalized(int block_number, int equation_number) const override
  {
    return (equation_type_and_normalized_equation[equation_reordered[get<1>(block_type_firstequation_size_mfs[block_number])+equation_number]].first == E_EVALUATE_S);
  };
  //! Return the expr_t of the equation equation_number belonging to the block block_number
  expr_t
  getBlockEquationExpr(int block_number, int equation_number) const override
  {
    return (equations[equation_reordered[get<1>(block_type_firstequation_size_mfs[block_number])+equation_number]]);
  };
  //! Return the expr_t of the renormalized equation equation_number belonging to the block block_number
  expr_t
  getBlockEquationRenormalizedExpr(int block_number, int equation_number) const override
  {
    return (equation_type_and_normalized_equation[equation_reordered[get<1>(block_type_firstequation_size_mfs[block_number])+equation_number]].second);
  };
  //! Return the original number of equation equation_number belonging to the block block_number
  int
  getBlockEquationID(int block_number, int equation_number) const override
  {
    return (equation_reordered[get<1>(block_type_firstequation_size_mfs[block_number])+equation_number]);
  };
  //! Return the original number of variable variable_number belonging to the block block_number
  int
  getBlockVariableID(int block_number, int variable_number) const override
  {
    return (variable_reordered[get<1>(block_type_firstequation_size_mfs[block_number])+variable_number]);
  };
  //! Return the original number of the exogenous variable varexo_number belonging to the block block_number
  int
  getBlockVariableExoID(int block_number, int variable_number) const override
  {
    return 0;
  };
  //! Return the position of equation_number in the block number belonging to the block block_number
  int
  getBlockInitialEquationID(int block_number, int equation_number) const override
  {
    return (static_cast<int>(inv_equation_reordered[equation_number]) - static_cast<int>(get<1>(block_type_firstequation_size_mfs[block_number])));
  };
  //! Return the position of variable_number in the block number belonging to the block block_number
  int
  getBlockInitialVariableID(int block_number, int variable_number) const override
  {
    return (static_cast<int>(inv_variable_reordered[variable_number]) - static_cast<int>(get<1>(block_type_firstequation_size_mfs[block_number])));
  };
  //! Return the position of variable_number in the block number belonging to the block block_number
  int
  getBlockInitialExogenousID(int block_number, int variable_number) const override
  {
    return -1;
  };
  //! Return the position of the deterministic exogenous variable_number in the block number belonging to the block block_number
  int
  getBlockInitialDetExogenousID(int block_number, int variable_number) const override
  {
    return -1;
  };
  //! Return the position of the other endogenous variable_number in the block number belonging to the block block_number
  int
  getBlockInitialOtherEndogenousID(int block_number, int variable_number) const override
  {
    return -1;
  };
};

#endif
