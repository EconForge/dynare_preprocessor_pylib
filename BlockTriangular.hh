/*
 * Copyright (C) 2007-2008 Dynare Team
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

#ifndef _BLOCKTRIANGULAR_HH
#define _BLOCKTRIANGULAR_HH

#include <string>
#include "CodeInterpreter.hh"
#include "ExprNode.hh"
#include "SymbolTable.hh"
//#include "ModelNormalization.hh"
//#include "ModelBlocks.hh"
#include "IncidenceMatrix.hh"
#include "ModelTree.hh"



//! Sparse matrix of double to store the values of the Jacobian
typedef map<pair<int ,int >,double> jacob_map;

typedef vector<pair<BlockSimulationType, pair<int, int> > > t_type;

//! Vector describing equations: BlockSimulationType, if BlockSimulationType == EVALUATE_s then a NodeID on the new normalized equation
typedef vector<pair<EquationType, NodeID > > t_etype;

//! Vector describing variables: max_lag in the block, max_lead in the block
typedef vector<pair< int, int> > t_vtype;

typedef set<int> temporary_terms_inuse_type;

//! For one lead/lag of one block, stores mapping of information between original model and block-decomposed model
struct IM_compact
{
  int size, u_init, u_finish, nb_endo, nb_other_endo, size_exo, size_other_endo;
  int *u, *us, *Var, *Equ, *Var_Index, *Equ_Index, *Exogenous, *Exogenous_Index, *Equ_X, *Equ_X_Index;
  int *u_other_endo, *Var_other_endo, *Equ_other_endo, *Var_Index_other_endo, *Equ_Index_other_endo;
};

//! One block of the model
struct Block
{
  int Size, Sized, nb_exo, nb_exo_det, nb_other_endo, Nb_Recursives;
  BlockType Type;
  BlockSimulationType Simulation_Type;
  int Max_Lead, Max_Lag, Nb_Lead_Lag_Endo;
  int Max_Lag_Endo, Max_Lead_Endo;
  int Max_Lag_Other_Endo, Max_Lead_Other_Endo;
  int Max_Lag_Exo, Max_Lead_Exo;
  bool is_linear;
  int *Equation, *Own_Derivative;
  EquationType *Equation_Type;
  NodeID *Equation_Normalized;
  int *Variable, *Other_Endogenous, *Exogenous;
  temporary_terms_type **Temporary_Terms_in_Equation;
  //temporary_terms_type *Temporary_terms;
  temporary_terms_inuse_type *Temporary_InUse;
  IM_compact *IM_lead_lag;
  int Code_Start, Code_Length;
};



//! The set of all blocks of the model
struct Model_Block
{
  int Size, Periods;
  Block* Block_List;
  //int *in_Block_Equ, *in_Block_Var, *in_Equ_of_Block, *in_Var_of_Block;
};


//! Creates the incidence matrix, computes prologue & epilogue, normalizes the model and computes the block decomposition
class BlockTriangular
{
private:
  //! Find equations and endogenous variables belonging to the prologue and epilogue of the model
  void Prologue_Epilogue(bool* IM, int &prologue, int &epilogue, int n, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, bool* IM0);
  //! Allocates and fills the Model structure describing the content of each block
  void Allocate_Block(int size, int *count_Equ, int count_Block, BlockType type, BlockSimulationType SimType, Model_Block * ModelBlock, t_etype &Equation_Type, int recurs_Size);
  //! Finds a matching between equations and endogenous variables
  bool Compute_Normalization(bool *IM, int equation_number, int prologue, int epilogue, bool verbose, bool *IM0, vector<int> &Index_Var_IM) const;
  //! Decomposes into recurive blocks the non purely recursive equations and determines for each block the minimum feedback variables
  void Compute_Block_Decomposition_and_Feedback_Variables_For_Each_Block(bool *IM, int nb_var, int prologue, int epilogue, vector<int> &Index_Equ_IM, vector<int> &Index_Var_IM, vector<pair<int, int> > &blocks, t_etype &Equation_Type, bool verbose_) const;
  //! determines the type of each equation of the model (could be evaluated or need to be solved)
  t_etype Equation_Type_determination(vector<BinaryOpNode *> &equations, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM);
  //! Tries to merge the consecutive blocks in a single block and determine the type of each block: recursive, simultaneous, ...
  t_type Reduce_Blocks_and_type_determination(int prologue, int epilogue, vector<pair<int, int> > &blocks, vector<BinaryOpNode *> &equations, t_etype &Equation_Type);
  //! Compute for each variable its maximum lead and lag in its block
  t_vtype Get_Variable_LeadLag_By_Block(vector<int > &components_set, int nb_blck_sim, int prologue, int epilogue) const;
public:
  SymbolTable &symbol_table;
  /*Blocks blocks;
  Normalization normalization;*/
  IncidenceMatrix incidencematrix;
  NumericalConstants &num_const;
  DataTree *Normalized_Equation;
  BlockTriangular(SymbolTable &symbol_table_arg, NumericalConstants &num_const_arg);
  ~BlockTriangular();
  //! Frees the Model structure describing the content of each block
  void Free_Block(Model_Block* ModelBlock) const;



  void Normalize_and_BlockDecompose_Static_0_Model(jacob_map &j_m, vector<BinaryOpNode *> &equations, t_etype &V_Equation_Type, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives);
  void Normalize_and_BlockDecompose(bool* IM, Model_Block* ModelBlock, int n, int &prologue, int &epilogue, vector<int> &Index_Var_IM, vector<int> &Index_Equ_IM, bool* IM_0 , jacob_map &j_m, vector<BinaryOpNode *> &equations, t_etype &equation_simulation_type, map<pair<int, pair<int, int> >, NodeID> &first_order_endo_derivatives);
  vector<int> Index_Equ_IM;
  vector<int> Index_Var_IM;
  int prologue, epilogue;
  bool bt_verbose;
  Model_Block* ModelBlock;
  int periods;
  inline static std::string BlockType0(int type)
  {
    switch (type)
      {
      case 0:
        return ("SIMULTANEOUS TIME SEPARABLE  ");
        break;
      case 1:
        return ("PROLOGUE                     ");
        break;
      case 2:
        return ("EPILOGUE                     ");
        break;
      case 3:
        return ("SIMULTANEOUS TIME UNSEPARABLE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };
  inline static std::string BlockSim(int type)
  {
    switch (type)
      {
      case EVALUATE_FORWARD:
      //case EVALUATE_FORWARD_R:
        return ("EVALUATE FORWARD             ");
        break;
      case EVALUATE_BACKWARD:
      //case EVALUATE_BACKWARD_R:
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
  inline static std::string c_Equation_Type(int type)
  {
    char c_Equation_Type[5][13]=
    {
    "E_UNKNOWN   ",
    "E_EVALUATE  ",
    //"E_EVALUATE_R",
    "E_EVALUATE_S",
    "E_SOLVE     "
    };
    return(c_Equation_Type[type]);
  };
};
#endif
