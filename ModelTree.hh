/*
 * Copyright (C) 2003-2009 Dynare Team
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

#ifndef _MODELTREE_HH
#define _MODELTREE_HH

using namespace std;

#include <string>
#include <vector>
#include <deque>
#include <map>
#include <ostream>

#include "DataTree.hh"

//! Vector describing equations: BlockSimulationType, if BlockSimulationType == EVALUATE_s then a NodeID on the new normalized equation
typedef vector<pair<EquationType, NodeID > > t_equation_type_and_normalized_equation;

//! Vector describing variables: max_lag in the block, max_lead in the block
typedef vector<pair< int, int> > t_lag_lead_vector;

//! for each block contains pair< pair<Simulation_Type, first_equation>, pair < Block_Size, Recursive_part_Size > >
typedef vector<pair< pair< BlockSimulationType, int>, pair<int, int> > > t_block_type_firstequation_size_mfs;

//! for a block contains derivatives pair< pair<block_equation_number, block_variable_number> , pair<lead_lag, NodeID> >
typedef vector< pair<pair<int, int>, pair< int, NodeID > > > t_block_derivatives_equation_variable_laglead_nodeid;

//! for all blocks derivatives description
typedef vector<t_block_derivatives_equation_variable_laglead_nodeid> t_blocks_derivatives;

//! Shared code for static and dynamic models
class ModelTree : public DataTree
{
  friend class DynamicModel;
  friend class StaticModel;

protected:
  //! Stores declared and generated auxiliary equations
  vector<BinaryOpNode *> equations;

  //! Only stores generated auxiliary equations, in an order meaningful for evaluation
  deque<BinaryOpNode *> aux_equations;

  //! Stores equation tags
  vector<pair<int, pair<string, string> > > equation_tags;

  //! Number of non-zero derivatives
  int NNZDerivatives[3];

  typedef map<pair<int, int>, NodeID> first_derivatives_type;
  //! First order derivatives
  /*! First index is equation number, second is variable w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Variable indices are those of the getDerivID() method.
  */
  first_derivatives_type first_derivatives;

  typedef map<pair<int, pair<int, int> >, NodeID> second_derivatives_type;
  //! Second order derivatives
  /*! First index is equation number, second and third are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only second order derivatives where var1 >= var2 (for obvious symmetry reasons).
    Variable indices are those of the getDerivID() method.
  */
  second_derivatives_type second_derivatives;

  typedef map<pair<int, pair<int, pair<int, int> > >, NodeID> third_derivatives_type;
  //! Third order derivatives
  /*! First index is equation number, second, third and fourth are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only third order derivatives where var1 >= var2 >= var3 (for obvious symmetry reasons).
    Variable indices are those of the getDerivID() method.
  */
  third_derivatives_type third_derivatives;

  //! Temporary terms (those which will be noted Txxxx)
  temporary_terms_type temporary_terms;

  //! Computes 1st derivatives
  /*! \param vars the derivation IDs w.r. to which compute the derivatives */
  void computeJacobian(const set<int> &vars);
  //! Computes 2nd derivatives
  /*! \param vars the derivation IDs w.r. to which derive the 1st derivatives */
  void computeHessian(const set<int> &vars);
  //! Computes 3rd derivatives
  /*! \param vars the derivation IDs w.r. to which derive the 2nd derivatives */
  void computeThirdDerivatives(const set<int> &vars);

  //! Write derivative of an equation w.r. to a variable
  void writeDerivative(ostream &output, int eq, int symb_id, int lag, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  //! Computes temporary terms (for all equations and derivatives)
  void computeTemporaryTerms(bool is_matlab);
  //! Writes temporary terms
  void writeTemporaryTerms(const temporary_terms_type &tt, ostream &output, ExprNodeOutputType output_type) const;
  //! Compiles temporary terms
  void compileTemporaryTerms(ostream &code_file, const temporary_terms_type &tt, map_idx_type map_idx, bool dynamic, bool steady_dynamic) const;
  //! Adds informations for simulation in a binary file
  void Write_Inf_To_Bin_File(const string &basename, int &u_count_int, bool &file_open, bool is_two_boundaries, int block_mfs) const;

  //! Writes model local variables
  /*! No temporary term is used in the output, so that local parameters declarations can be safely put before temporary terms declaration in the output files */
  void writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type) const;
  //! Compiles model equations
  void compileModelEquations(ostream &code_file, const temporary_terms_type &tt, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const;

  //! Writes LaTeX model file
  void writeLatexModelFile(const string &filename, ExprNodeOutputType output_type) const;

  //! Sparse matrix of double to store the values of the Jacobian
  /*! First index is equation number, second index is endogenous type specific ID */
  typedef map<pair<int, int>, double> jacob_map;

  //! Sparse matrix of double to store the values of the Jacobian
  /*! First index is lag, second index is equation number, third index is endogenous type specific ID */
  typedef map<pair<int, pair<int, int> >, NodeID> dynamic_jacob_map;

  //! Normalization of equations
  /*! Maps endogenous type specific IDs to equation numbers */
  vector<int> endo2eq;

  //! number of equation in the prologue and in the epilogue
  unsigned int epilogue, prologue;

  //! for each block contains pair< max_lag, max_lead>
  t_lag_lead_vector block_lag_lead;

  //! Compute the matching between endogenous and variable using the jacobian contemporaneous_jacobian
  /*!
    \param contemporaneous_jacobian Jacobian used as an incidence matrix: all elements declared in the map (even if they are zero), are used as vertices of the incidence matrix
    \return True if a complete normalization has been achieved
  */
  bool computeNormalization(const jacob_map &contemporaneous_jacobian, bool verbose);

  //! Try to compute the matching between endogenous and variable using a decreasing cutoff
  /*!
    Applied to the jacobian contemporaneous_jacobian and stop when a matching is found.
    If no matching is found using a strictly positive cutoff, then a zero cutoff is applied (i.e. use a symbolic normalization); in that case, the method adds zeros in the jacobian matrices to reflect all the edges in the symbolic incidence matrix.
    If no matching is found with a zero cutoff close to zero an error message is printout.
  */
  void computeNonSingularNormalization(jacob_map &contemporaneous_jacobian, double cutoff, jacob_map &static_jacobian, dynamic_jacob_map &dynamic_jacobian);

  //! Try to normalized each unnormalized equation (matched endogenous variable only on the LHS)
  void computeNormalizedEquations(multimap<int, int> &endo2eqs) const;
  //! Evaluate the jacobian and suppress all the elements below the cutoff
  void evaluateAndReduceJacobian(const eval_context_type &eval_context, jacob_map &contemporaneous_jacobian, jacob_map &static_jacobian, dynamic_jacob_map &dynamic_jacobian, double cutoff, bool verbose);
  //! Search the equations and variables belonging to the prologue and the epilogue of the model
  void computePrologueAndEpilogue(jacob_map &static_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered, unsigned int &prologue, unsigned int &epilogue);
  //! Determine the type of each equation of model and try to normalized the unnormalized equation using computeNormalizedEquations
  t_equation_type_and_normalized_equation equationTypeDetermination(vector<BinaryOpNode *> &equations, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, int mfs);
  //! Compute the block decomposition and for a non-recusive block find the minimum feedback set
  void computeBlockDecompositionAndFeedbackVariablesForEachBlock(jacob_map &static_jacobian, dynamic_jacob_map &dynamic_jacobian, int prologue, int epilogue, vector<int> &Index_Equ_IM, vector<int> &Index_Var_IM, vector<pair<int, int> > &blocks, t_equation_type_and_normalized_equation &Equation_Type, bool verbose_, bool select_feedback_variable, int mfs, vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered) const;
  //! Reduce the number of block merging the same type equation in the prologue and the epilogue and determine the type of each block
  t_block_type_firstequation_size_mfs reduceBlocksAndTypeDetermination(dynamic_jacob_map &dynamic_jacobian, int prologue, int epilogue, vector<pair<int, int> > &blocks, vector<BinaryOpNode *> &equations, t_equation_type_and_normalized_equation &Equation_Type, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM);
  //! Determine the maximum number of lead and lag for the endogenous variable in a bloc
  void getVariableLeadLagByBlock(dynamic_jacob_map &dynamic_jacobian, vector<int > &components_set, int nb_blck_sim, int prologue, int epilogue, t_lag_lead_vector &equation_lead_lag, t_lag_lead_vector &variable_lead_lag, vector<int> equation_reordered, vector<int> variable_reordered) const;
  //! Print an abstract of the block structure of the model
  void printBlockDecomposition(vector<pair<int, int> > blocks);
  //! Determine for each block if it is linear or not
  vector<bool> BlockLinear(t_blocks_derivatives &blocks_derivatives, vector<int> &variable_reordered);

  virtual SymbolType getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException) = 0;
  virtual int getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException) = 0;
  virtual int getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException) = 0;

  //! Determine the simulation type of each block
  virtual BlockSimulationType getBlockSimulationType(int block_number) const = 0;
  //! Return the number of blocks
  virtual unsigned int getNbBlocks() const = 0;
  //! Return the first equation number of a block
  virtual unsigned int getBlockFirstEquation(int block_number) const = 0;
  //! Return the size of the block block_number
  virtual unsigned int getBlockSize(int block_number) const = 0;
  //! Return the number of feedback variable of the block block_number
  virtual unsigned int getBlockMfs(int block_number) const = 0;
  //! Return the maximum lag in a block
  virtual unsigned int getBlockMaxLag(int block_number) const = 0;
  //! Return the maximum lead in a block
  virtual unsigned int getBlockMaxLead(int block_number) const = 0;
  //! Return the type of equation (equation_number) belonging to the block block_number
  virtual EquationType getBlockEquationType(int block_number, int equation_number) const = 0;
  //! Return true if the equation has been normalized
  virtual bool isBlockEquationRenormalized(int block_number, int equation_number) const = 0;
  //! Return the NodeID of the equation equation_number belonging to the block block_number
  virtual NodeID getBlockEquationNodeID(int block_number, int equation_number) const = 0;
  //! Return the NodeID of the renormalized equation equation_number belonging to the block block_number
  virtual NodeID getBlockEquationRenormalizedNodeID(int block_number, int equation_number) const = 0;
  //! Return the original number of equation equation_number belonging to the block block_number
  virtual int getBlockEquationID(int block_number, int equation_number) const = 0;
  //! Return the original number of variable variable_number belonging to the block block_number
  virtual int getBlockVariableID(int block_number, int variable_number) const = 0;
  //! Return the position of equation_number in the block number belonging to the block block_number
  virtual int getBlockInitialEquationID(int block_number, int equation_number) const = 0;
  //! Return the position of variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialVariableID(int block_number, int variable_number) const = 0;

public:
  ModelTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Declare a node as an equation of the model
  void addEquation(NodeID eq);
  //! Adds tags to equation number i
  void addEquationTags(int i, const string &key, const string &value);
  //! Declare a node as an auxiliary equation of the model, adding it at the end of the list of auxiliary equations
  void addAuxEquation(NodeID eq);
  //! Returns the number of equations in the model
  int equation_number() const;

  inline static std::string
  c_Equation_Type(int type)
  {
    char c_Equation_Type[4][13] =
      {
        "E_UNKNOWN   ",
        "E_EVALUATE  ",
        "E_EVALUATE_S",
        "E_SOLVE     "
      };
    return (c_Equation_Type[type]);
  };

  inline static std::string
  BlockType0(BlockType type)
  {
    switch (type)
      {
      case SIMULTANS:
        return ("SIMULTANEOUS TIME SEPARABLE  ");
        break;
      case PROLOGUE:
        return ("PROLOGUE                     ");
        break;
      case EPILOGUE:
        return ("EPILOGUE                     ");
        break;
      case SIMULTAN:
        return ("SIMULTANEOUS TIME UNSEPARABLE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };

  inline static std::string
  BlockSim(int type)
  {
    switch (type)
      {
      case EVALUATE_FORWARD:
        return ("EVALUATE FORWARD             ");
        break;
      case EVALUATE_BACKWARD:
        return ("EVALUATE BACKWARD            ");
        break;
      case SOLVE_FORWARD_SIMPLE:
        return ("SOLVE FORWARD SIMPLE         ");
        break;
      case SOLVE_BACKWARD_SIMPLE:
        return ("SOLVE BACKWARD SIMPLE        ");
        break;
      case SOLVE_TWO_BOUNDARIES_SIMPLE:
        return ("SOLVE TWO BOUNDARIES SIMPLE  ");
        break;
      case SOLVE_FORWARD_COMPLETE:
        return ("SOLVE FORWARD COMPLETE       ");
        break;
      case SOLVE_BACKWARD_COMPLETE:
        return ("SOLVE BACKWARD COMPLETE      ");
        break;
      case SOLVE_TWO_BOUNDARIES_COMPLETE:
        return ("SOLVE TWO BOUNDARIES COMPLETE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };
};

#endif
