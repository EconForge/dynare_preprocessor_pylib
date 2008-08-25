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

#ifndef BLOCKTRIANGULAR_H
#define BLOCKTRIANGULAR_H

#include <string>
#include "ExprNode.hh"
#include "SymbolTable.hh"
#include "ModelNormalization.hh"
#include "ModelBlocks.hh"

/*!
  \class  BlockTriangular
  \brief  Creat the incidence matrice and reorder the model's equations.
*/

#include "ExprNode.hh"

struct List_IM
{
  List_IM* pNext;
  int lead_lag;
  bool* IM;
};

typedef map<pair<int ,int >,double> jacob_map;

class BlockTriangular
{
public:
  BlockTriangular(const SymbolTable &symbol_table_arg);
  /*! The incidence matrix for each lead and lags */
  const SymbolTable &symbol_table;
  Blocks blocks;
  Normalization normalization;
  List_IM* Build_IM(int lead_lag);
  List_IM* Get_IM(int lead_lag);
  bool* bGet_IM(int lead_lag) const;
  void fill_IM(int equation, int variable_endo, int lead_lag);
  void unfill_IM(int equation, int variable_endo, int lead_lag);
  void init_incidence_matrix(int nb_endo);
  void Print_IM(int n) const;
  void Free_IM(List_IM* First_IM) const;
  void Print_SIM(bool* IM, int n) const;
  void Normalize_and_BlockDecompose_Static_0_Model(const jacob_map &j_m);
  bool Normalize_and_BlockDecompose(bool* IM, Model_Block* ModelBlock, int n, int* prologue, int* epilogue, simple* Index_Var_IM, simple* Index_Equ_IM, bool Do_Normalization, bool mixing, bool* IM_0 , jacob_map j_m);
  void Prologue_Epilogue(bool* IM, int* prologue, int* epilogue, int n, simple* Index_Var_IM, simple* Index_Equ_IM, bool* IM0);
  void swap_IM_c(bool *SIM, int pos1, int pos2, int pos3, simple* Index_Var_IM, simple* Index_Equ_IM, int n);
  void Allocate_Block(int size, int *count_Equ, int *count_Block, int type, Model_Block * ModelBlock);
  void Free_Block(Model_Block* ModelBlock) const;
  List_IM *First_IM ;
  List_IM *Last_IM ;
  simple *Index_Equ_IM;
  simple *Index_Var_IM;
  int prologue, epilogue;
  int Model_Max_Lead, Model_Max_Lag, periods;
  bool bt_verbose;
  int endo_nbr;
  Model_Block* ModelBlock;
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
      case EVALUATE_FOREWARD:
      case EVALUATE_FOREWARD_R:
        return ("EVALUATE FOREWARD            ");
        break;
      case EVALUATE_BACKWARD:
      case EVALUATE_BACKWARD_R:
        return ("EVALUATE BACKWARD            ");
        break;
      case SOLVE_FOREWARD_SIMPLE:
        return ("SOLVE FOREWARD SIMPLE        ");
        break;
      case SOLVE_BACKWARD_SIMPLE:
        return ("SOLVE BACKWARD SIMPLE        ");
        break;
      case SOLVE_TWO_BOUNDARIES_SIMPLE:
        return ("SOLVE TWO BOUNDARIES SIMPLE  ");
        break;
      case SOLVE_FOREWARD_COMPLETE:
        return ("SOLVE FOREWARD COMPLETE      ");
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
//------------------------------------------------------------------------------
#endif
