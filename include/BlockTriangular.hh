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
#include "ModelNormalization.hh"
#include "ModelBlocks.hh"
#include "IncidenceMatrix.hh"







//! Matrix of doubles for representing jacobian
typedef map<pair<int ,int >,double> jacob_map;

typedef vector<pair<BlockSimulationType, int> > t_type;

//! Create the incidence matrix, computes prologue & epilogue, normalizes the model and computes the block decomposition
class BlockTriangular
{
  //friend class IncidenceMatrix;
public:
  const SymbolTable &symbol_table;
  BlockTriangular(const SymbolTable &symbol_table_arg);
  //BlockTriangular(const IncidenceMatrix &incidence_matrix_arg);
  //const SymbolTable &symbol_table;
  Blocks blocks;
  Normalization normalization;
  IncidenceMatrix incidencematrix;
  void Normalize_and_BlockDecompose_Static_0_Model(const jacob_map &j_m, vector<BinaryOpNode *> equations);
  bool Normalize_and_BlockDecompose(bool* IM, Model_Block* ModelBlock, int n, int* prologue, int* epilogue, simple* Index_Var_IM, simple* Index_Equ_IM, bool Do_Normalization, bool mixing, bool* IM_0 , jacob_map j_m, vector<BinaryOpNode *> equations);
  void Prologue_Epilogue(bool* IM, int* prologue, int* epilogue, int n, simple* Index_Var_IM, simple* Index_Equ_IM, bool* IM0);
  void Allocate_Block(int size, int *count_Equ, int count_Block, BlockType type, BlockSimulationType SimType, Model_Block * ModelBlock);
  void Free_Block(Model_Block* ModelBlock) const;
  t_type Reduce_Blocks_and_type_determination(int prologue, int epilogue, block_result_t* res, vector<BinaryOpNode *> equations );
  simple *Index_Equ_IM;
  simple *Index_Var_IM;
  int prologue, epilogue;
  bool bt_verbose;
  //int endo_nbr, exo_nbr;
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
      case EVALUATE_FORWARD_R:
        return ("EVALUATE FORWARD             ");
        break;
      case EVALUATE_BACKWARD:
      case EVALUATE_BACKWARD_R:
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
