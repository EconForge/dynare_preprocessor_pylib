/*
 * Copyright © 2007-2022 Dynare Team
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

#ifndef _BYTECODE_HH
#define _BYTECODE_HH

#include <fstream>
#include <vector>
#include <utility>
#include <ios>

#include "CommonEnums.hh"

#ifdef BYTECODE_MEX
# include <dynmex.h>
# include <cstring>
#endif

using namespace std;

// The different opcodes of bytecode
enum class Tags
  {
   FLDZ, // Loads a zero onto the stack
   FLDC, // Loads a constant term onto the stack

   FDIMT, // Defines the number of temporary terms - dynamic context (the period has to be indicated)
   FDIMST, // Defines the number of temporary terms - static context (the period hasn’t to be indicated)
   FLDT, // Loads a temporary term onto the stack - dynamic context (the period has to be indicated)
   FLDST, // Loads a temporary term onto the stack - static context (the period hasn’t to be indicated)
   FSTPT, // Stores a temporary term from the stack - dynamic context (the period has to be indicated)
   FSTPST, // Stores a temporary term from the stack - static context (the period hasn’t to be indicated)

   FLDU, // Loads an element of the vector U onto the stack - dynamic context (the period has to be indicated)
   FLDSU, // Loads an element of the vector U onto the stack - static context (the period hasn’t to be indicated)
   FSTPU, // Stores an element of the vector U from the stack - dynamic context (the period has to be indicated)
   FSTPSU, // Stores an element of the vector U from the stack - static context (the period hasn’t to be indicated)

   FLDV, // Loads a variable (described in SymbolType) onto the stack - dynamic context (the period has to be indicated)
   FLDSV, // Loads a variable (described in SymbolType) onto the stack - static context (the period hasn’t to be indicated)
   FLDVS, // Loads a variable (described in SymbolType) onto the stack - dynamic context but inside the STEADY_STATE operator (the period hasn’t to be indicated)
   FSTPV, // Stores a variable (described in SymbolType) from the stack - dynamic context (the period has to be indicated)
   FSTPSV, // Stores a variable (described in SymbolType) from the stack - static context (the period hasn’t to be indicated)

   FLDR, // Loads a residual onto the stack
   FSTPR, // Stores a residual from the stack

   FSTPG, // Stores a derivative from the stack
   FSTPG2, // Stores a derivative matrix for a static model from the stack
   FSTPG3, // Stores a derivative matrix for a dynamic model from the stack

   FUNARY, // A unary operator
   FBINARY, // A binary operator
   FTRINARY, // A trinary operator

   FCUML, // Cumulates the result

   FJMPIFEVAL, // Jump if evaluate = true
   FJMP, // Jump

   FBEGINBLOCK, // Marks the beginning of a model block
   FENDBLOCK, // Marks the end of a model block
   FENDEQU, // Marks the last equation of the block; for a block that has to be solved, the derivatives appear just after this flag
   FEND, // Marks the end of the model code

   FNUMEXPR, // Stores the expression type and references

   FCALL, // Call an external function
   FLDTEF, // Loads the result of an external function onto the stack
   FSTPTEF, // Stores the result of an external function from the stack
   FLDTEFD, // Loads the result of the 1st derivative of an external function onto the stack
   FSTPTEFD, // Stores the result of the 1st derivative of an external function from the stack
   FLDTEFDD, // Loads the result of the 2nd derivative of an external function onto the stack
   FSTPTEFDD // Stores the result of the 2nd derivative of an external function from the stack
  };

enum class ExpressionType
  {
   TemporaryTerm,
   ModelEquation,
   FirstEndoDerivative,
   FirstOtherEndoDerivative,
   FirstExoDerivative,
   FirstExodetDerivative,
  };

enum class ExternalFunctionCallType
  {
   levelWithoutDerivative,
   levelWithFirstDerivative,
   levelWithFirstAndSecondDerivative,
   separatelyProvidedFirstDerivative,
   numericalFirstDerivative,
   separatelyProvidedSecondDerivative,
   numericalSecondDerivative
  };

struct Block_contain_type
{
  int Equation, Variable, Own_Derivative;
};

class BytecodeWriter;

class BytecodeInstruction
{
  template<typename B>
  friend BytecodeWriter &operator<<(BytecodeWriter &code_file, const B &instr);
protected:
  Tags op_code;
public:
  explicit BytecodeInstruction(Tags op_code_arg) : op_code{op_code_arg}
  {
  };
};

template<typename T1>
class TagWithOneArgument : public BytecodeInstruction
{
protected:
  T1 arg1;
public:
  TagWithOneArgument(Tags op_code_arg, T1 arg_arg1) : BytecodeInstruction{op_code_arg},
                                                      arg1{arg_arg1}
  {
  };
};

template<typename T1, typename T2>
class TagWithTwoArguments : public BytecodeInstruction
{
protected:
  T1 arg1;
  T2 arg2;
public:
  TagWithTwoArguments(Tags op_code_arg, T1 arg_arg1, T2 arg_arg2) :
    BytecodeInstruction{op_code_arg}, arg1{arg_arg1}, arg2{arg_arg2}
  {
  };
};

template<typename T1, typename T2, typename T3>
class TagWithThreeArguments : public BytecodeInstruction
{
protected:
  T1 arg1;
  T2 arg2;
  T3 arg3;
public:
  TagWithThreeArguments(Tags op_code_arg, T1 arg_arg1, T2 arg_arg2, T3 arg_arg3) :
    BytecodeInstruction{op_code_arg}, arg1{arg_arg1}, arg2{arg_arg2}, arg3{arg_arg3}
  {
  };
};

template<typename T1, typename T2, typename T3, typename T4>
class TagWithFourArguments : public BytecodeInstruction
{
protected:
  T1 arg1;
  T2 arg2;
  T3 arg3;
  T4 arg4;
public:
  TagWithFourArguments(Tags op_code_arg, T1 arg_arg1, T2 arg_arg2, T3 arg_arg3, T4 arg_arg4) :
    BytecodeInstruction{op_code_arg}, arg1{arg_arg1}, arg2{arg_arg2},
    arg3{move(arg_arg3)}, arg4{arg_arg4}
  {
  };
};

class FLDZ_ : public BytecodeInstruction
{
public:
  FLDZ_() : BytecodeInstruction{Tags::FLDZ}
  {
  };
};

class FEND_ : public BytecodeInstruction
{
public:
  FEND_() : BytecodeInstruction{Tags::FEND}
  {
  };
};

class FENDBLOCK_ : public BytecodeInstruction
{
public:
  FENDBLOCK_() : BytecodeInstruction{Tags::FENDBLOCK}
  {
  };
};

class FENDEQU_ : public BytecodeInstruction
{
public:
  FENDEQU_() : BytecodeInstruction{Tags::FENDEQU}
  {
  };
};

class FCUML_ : public BytecodeInstruction
{
public:
  FCUML_() : BytecodeInstruction{Tags::FCUML}
  {
  };
};

class FDIMT_ : public TagWithOneArgument<int>
{
public:
  explicit FDIMT_(int size_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FDIMT, size_arg}
  {
  };
  int
  get_size()
  {
    return arg1;
  };
};

class FDIMST_ : public TagWithOneArgument<int>
{
public:
  explicit FDIMST_(int size_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FDIMST, size_arg}
  {
  };
  int
  get_size()
  {
    return arg1;
  };
};

class FLDC_ : public TagWithOneArgument<double>
{
public:
  explicit FLDC_(double value_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FLDC, value_arg}
  {
  };
  double
  get_value()
  {
    return arg1;
  };
};

class FLDU_ : public TagWithOneArgument<int>
{
public:
  explicit FLDU_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FLDU, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDSU_ : public TagWithOneArgument<int>
{
public:
  explicit FLDSU_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FLDSU, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDR_ : public TagWithOneArgument<int>
{
public:
  explicit FLDR_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FLDR, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDT_ : public TagWithOneArgument<int>
{
public:
  explicit FLDT_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FLDT, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDST_ : public TagWithOneArgument<int>
{
public:
  explicit FLDST_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FLDST, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPT_ : public TagWithOneArgument<int>
{
public:
  explicit FSTPT_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FSTPT, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPST_ : public TagWithOneArgument<int>
{
public:
  explicit FSTPST_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FSTPST, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPR_ : public TagWithOneArgument<int>
{
public:
  explicit FSTPR_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FSTPR, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPU_ : public TagWithOneArgument<int>
{
public:
  explicit FSTPU_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FSTPU, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPSU_ : public TagWithOneArgument<int>
{
public:
  explicit FSTPSU_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FSTPSU, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPG_ : public TagWithOneArgument<int>
{
public:
  explicit FSTPG_(int pos_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FSTPG, pos_arg}
  {
  };
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPG2_ : public TagWithTwoArguments<int, int>
{
public:
  FSTPG2_(int row_arg, int col_arg) : TagWithTwoArguments::TagWithTwoArguments{Tags::FSTPG2, row_arg, col_arg}
  {
  };
  int
  get_row()
  {
    return arg1;
  };
  int
  get_col()
  {
    return arg2;
  };
};

class FSTPG3_ : public TagWithFourArguments<int, int, int, int>
{
public:
  FSTPG3_(int row_arg, int col_arg, int lag_arg, int col_pos_arg) : TagWithFourArguments::TagWithFourArguments{Tags::FSTPG3, row_arg, col_arg, lag_arg, col_pos_arg}
  {
  };
  int
  get_row()
  {
    return arg1;
  };
  int
  get_col()
  {
    return arg2;
  };
  int
  get_lag()
  {
    return arg2;
  };
  int
  get_col_pos()
  {
    return arg4;
  };
};

class FUNARY_ : public TagWithOneArgument<UnaryOpcode>
{
public:
  explicit FUNARY_(UnaryOpcode op_type_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FUNARY, op_type_arg}
  {
  };
  UnaryOpcode
  get_op_type()
  {
    return arg1;
  };
};

class FBINARY_ : public TagWithOneArgument<BinaryOpcode>
{
public:
  explicit FBINARY_(BinaryOpcode op_type_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FBINARY, op_type_arg}
  {
  };
  BinaryOpcode
  get_op_type()
  {
    return arg1;
  };
};

class FTRINARY_ : public TagWithOneArgument<TrinaryOpcode>
{
public:
  explicit FTRINARY_(TrinaryOpcode op_type_arg) : TagWithOneArgument::TagWithOneArgument{Tags::FTRINARY, op_type_arg}
  {
  };
  TrinaryOpcode
  get_op_type()
  {
    return arg1;
  };
};

class FJMPIFEVAL_ : public TagWithOneArgument<int>
{
public:
  explicit FJMPIFEVAL_(int arg_pos) : TagWithOneArgument::TagWithOneArgument{Tags::FJMPIFEVAL, arg_pos}
  {
  };
  int
  get_pos()
  {
    return arg1;
  }
};

class FJMP_ : public TagWithOneArgument<int>
{
public:
  explicit FJMP_(int arg_pos) : TagWithOneArgument::TagWithOneArgument{Tags::FJMP, arg_pos}
  {
  };
  int
  get_pos()
  {
    return arg1;
  }
};

class FLDTEF_ : public TagWithOneArgument<int>
{
public:
  explicit FLDTEF_(int number) : TagWithOneArgument::TagWithOneArgument{Tags::FLDTEF, number}
  {
  };
  int
  get_number()
  {
    return arg1;
  }
};

class FSTPTEF_ : public TagWithOneArgument<int>
{
public:
  explicit FSTPTEF_(int number) : TagWithOneArgument::TagWithOneArgument{Tags::FSTPTEF, number}
  {
  };
  int
  get_number()
  {
    return arg1;
  }
};

class FLDTEFD_ : public TagWithTwoArguments<int, int>
{
public:
  FLDTEFD_(int indx, int row) : TagWithTwoArguments::TagWithTwoArguments{Tags::FLDTEFD, indx, row}
  {
  };
  int
  get_indx()
  {
    return arg1;
  };
  int
  get_row()
  {
    return arg2;
  };
};

class FSTPTEFD_ : public TagWithTwoArguments<int, int>
{
public:
  FSTPTEFD_(int indx, int row) : TagWithTwoArguments::TagWithTwoArguments{Tags::FSTPTEFD, indx, row}
  {
  };
  int
  get_indx()
  {
    return arg1;
  };
  int
  get_row()
  {
    return arg2;
  };
};

class FLDTEFDD_ : public TagWithThreeArguments<int, int, int>
{
public:
  FLDTEFDD_(int indx, int row, int col) : TagWithThreeArguments::TagWithThreeArguments{Tags::FLDTEFDD, indx, row, col}
  {
  };
  int
  get_indx()
  {
    return arg1;
  };
  int
  get_row()
  {
    return arg2;
  };
  int
  get_col()
  {
    return arg3;
  };
};

class FSTPTEFDD_ : public TagWithThreeArguments<int, int, int>
{
public:
  FSTPTEFDD_(int indx, int row, int col) : TagWithThreeArguments::TagWithThreeArguments{Tags::FSTPTEF, indx, row, col}
  {
  };
  int
  get_indx()
  {
    return arg1;
  };
  int
  get_row()
  {
    return arg2;
  };
  int
  get_col()
  {
    return arg3;
  };
};

class FLDVS_ : public TagWithTwoArguments<SymbolType, int>
{
public:
  FLDVS_(SymbolType type_arg, int pos_arg) : TagWithTwoArguments::TagWithTwoArguments{Tags::FLDVS, type_arg, pos_arg}
  {
  };
  SymbolType
  get_type()
  {
    return arg1;
  };
  int
  get_pos()
  {
    return arg2;
  };
};

class FLDSV_ : public TagWithTwoArguments<SymbolType, int>
{
public:
  FLDSV_(SymbolType type_arg, int pos_arg) : TagWithTwoArguments::TagWithTwoArguments{Tags::FLDSV, type_arg, pos_arg}
  {
  };
  SymbolType
  get_type()
  {
    return arg1;
  };
  int
  get_pos()
  {
    return arg2;
  };
};

class FSTPSV_ : public TagWithTwoArguments<SymbolType, int>
{
public:
  FSTPSV_(SymbolType type_arg, int pos_arg) : TagWithTwoArguments::TagWithTwoArguments{Tags::FSTPSV, type_arg, pos_arg}
  {
  };
  SymbolType
  get_type()
  {
    return arg1;
  };
  int
  get_pos()
  {
    return arg2;
  };
};

class FLDV_ : public TagWithThreeArguments<SymbolType, int, int>
{
public:
  FLDV_(SymbolType type_arg, int pos_arg, int lead_lag_arg) :
    TagWithThreeArguments::TagWithThreeArguments{Tags::FLDV, type_arg, pos_arg, lead_lag_arg}
  {
  };
  SymbolType
  get_type()
  {
    return arg1;
  };
  int
  get_pos()
  {
    return arg2;
  };
  int
  get_lead_lag()
  {
    return arg3;
  };
};

class FSTPV_ : public TagWithThreeArguments<SymbolType, int, int>
{
public:
  FSTPV_(SymbolType type_arg, int pos_arg, int lead_lag_arg) :
    TagWithThreeArguments::TagWithThreeArguments{Tags::FSTPV, type_arg, pos_arg, lead_lag_arg}
  {
  };
  SymbolType
  get_type()
  {
    return arg1;
  };
  int
  get_pos()
  {
    return arg2;
  };
  int
  get_lead_lag()
  {
    return arg3;
  };
};

class FCALL_ : public BytecodeInstruction
{
  template<typename B>
  friend BytecodeWriter &operator<<(BytecodeWriter &code_file, const B &instr);
private:
  int nb_output_arguments, nb_input_arguments, indx;
  string func_name;
  string arg_func_name;
  int add_input_arguments{0}, row{0}, col{0};
  ExternalFunctionCallType call_type;
public:
  FCALL_() : BytecodeInstruction{Tags::FCALL}
  {
  };
  FCALL_(int nb_output_arguments_arg, int nb_input_arguments_arg, string func_name_arg, int indx_arg, ExternalFunctionCallType call_type_arg) :
    BytecodeInstruction{Tags::FCALL},
    nb_output_arguments{nb_output_arguments_arg},
    nb_input_arguments{nb_input_arguments_arg},
    indx{indx_arg},
    func_name{move(func_name_arg)},
    call_type{call_type_arg}
  {
  };
  string
  get_function_name()
  {
    //printf("get_function_name => func_name=%s\n",func_name.c_str());fflush(stdout);
    return func_name;
  };
  int
  get_nb_output_arguments()
  {
    return nb_output_arguments;
  };
  int
  get_nb_input_arguments()
  {
    return nb_input_arguments;
  };
  int
  get_indx()
  {
    return indx;
  };
  void
  set_arg_func_name(string arg_arg_func_name)
  {
    arg_func_name = arg_arg_func_name;
  };
  string
  get_arg_func_name()
  {
    return arg_func_name;
  };
  void
  set_nb_add_input_arguments(int arg_add_input_arguments)
  {
    add_input_arguments = arg_add_input_arguments;
  };
  int
  get_nb_add_input_arguments()
  {
    return add_input_arguments;
  };
  void
  set_row(int arg_row)
  {
    row = arg_row;
  };
  int
  get_row()
  {
    return row;
  }
  void
  set_col(int arg_col)
  {
    col = arg_col;
  };
  int
  get_col()
  {
    return col;
  };
  ExternalFunctionCallType
  get_call_type()
  {
    return call_type;
  }
#ifdef BYTECODE_MEX

  char *
  load(char *code)
  {
    op_code = Tags::FCALL; code += sizeof(op_code);
    memcpy(&nb_output_arguments, code, sizeof(nb_output_arguments)); code += sizeof(nb_output_arguments);
    memcpy(&nb_input_arguments, code, sizeof(nb_input_arguments)); code += sizeof(nb_input_arguments);
    memcpy(&indx, code, sizeof(indx)); code += sizeof(indx);
    memcpy(&add_input_arguments, code, sizeof(add_input_arguments)); code += sizeof(add_input_arguments);
    memcpy(&row, code, sizeof(row)); code += sizeof(row);
    memcpy(&col, code, sizeof(col)); code += sizeof(col);
    memcpy(&call_type, code, sizeof(call_type)); code += sizeof(call_type);
    int size;
    memcpy(&size, code, sizeof(size)); code += sizeof(size);
    char *name = static_cast<char *>(mxMalloc((size+1)*sizeof(char)));
    memcpy(name, code, size); code += size;
    name[size] = 0;
    func_name = name;
    mxFree(name);
    memcpy(&size, code, sizeof(size)); code += sizeof(size);
    name = static_cast<char *>(mxMalloc((size+1)*sizeof(char)));
    memcpy(name, code, size); code += size;
    name[size] = 0;
    arg_func_name = name;
    mxFree(name);
    return code;
  }
#endif
};

class FNUMEXPR_ : public BytecodeInstruction
{
private:
  ExpressionType expression_type;
  int equation, dvariable1, lag1;
public:
  FNUMEXPR_(const ExpressionType expression_type_arg, int equation_arg) :
    BytecodeInstruction{Tags::FNUMEXPR},
    expression_type{expression_type_arg},
    equation{equation_arg},
    dvariable1{0},
    lag1{0}
  {
  };
  FNUMEXPR_(const ExpressionType expression_type_arg, int equation_arg, int dvariable1_arg) :
    BytecodeInstruction{Tags::FNUMEXPR},
    expression_type{expression_type_arg},
    equation{equation_arg},
    dvariable1{dvariable1_arg},
    lag1{0}
  {
  };
  FNUMEXPR_(const ExpressionType expression_type_arg, int equation_arg, int dvariable1_arg, int lag1_arg) :
    BytecodeInstruction{Tags::FNUMEXPR},
    expression_type{expression_type_arg},
    equation{equation_arg},
    dvariable1{dvariable1_arg},
    lag1{lag1_arg}
  {
  };
  ExpressionType
  get_expression_type()
  {
    return expression_type;
  }
  int
  get_equation()
  {
    return equation;
  };
  int
  get_dvariable1()
  {
    return dvariable1;
  };
  int
  get_lag1()
  {
    return lag1;
  };
};

class FBEGINBLOCK_ : public BytecodeInstruction
{
  template<typename B>
  friend BytecodeWriter &operator<<(BytecodeWriter &code_file, const B &instr);
private:
  int size{0};
  BlockSimulationType type;
  vector<int> variable;
  vector<int> equation;
  vector<int> other_endogenous;
  vector<int> exogenous;
  vector<int> det_exogenous;
  bool is_linear{false};
  vector<Block_contain_type> Block_Contain_;
  int endo_nbr{0};
  int Max_Lag{0};
  int Max_Lead{0};
  int u_count_int{0};
  int nb_col_jacob{0};
  int det_exo_size, exo_size, other_endo_size;
  int nb_col_det_exo_jacob, nb_col_exo_jacob, nb_col_other_endo_jacob;
public:
  FBEGINBLOCK_() : BytecodeInstruction{Tags::FBEGINBLOCK},
                   type{BlockSimulationType::unknown}
  {
  }
  FBEGINBLOCK_(int size_arg, BlockSimulationType type_arg, int first_element, int block_size,
               const vector<int> &variable_arg, const vector<int> &equation_arg,
               bool is_linear_arg, int endo_nbr_arg, int Max_Lag_arg, int Max_Lead_arg, int &u_count_int_arg, int nb_col_jacob_arg,
               int det_exo_size_arg, int nb_col_det_exo_jacob_arg, int exo_size_arg, int nb_col_exo_jacob_arg, int other_endo_size_arg, int nb_col_other_endo_jacob_arg,
               vector<int> det_exogenous_arg, vector<int> exogenous_arg, vector<int> other_endogenous_arg) :
    BytecodeInstruction{Tags::FBEGINBLOCK},
    size{size_arg},
    type{type_arg},
    variable{variable_arg.begin()+first_element, variable_arg.begin()+(first_element+block_size)},
    equation{equation_arg.begin()+first_element, equation_arg.begin()+(first_element+block_size)},
    other_endogenous{move(other_endogenous_arg)},
    exogenous{move(exogenous_arg)},
    det_exogenous{move(det_exogenous_arg)},
    is_linear{is_linear_arg},
    endo_nbr{endo_nbr_arg},
    Max_Lag{Max_Lag_arg},
    Max_Lead{Max_Lead_arg},
    u_count_int{u_count_int_arg},
    nb_col_jacob{nb_col_jacob_arg},
    det_exo_size{det_exo_size_arg},
    exo_size{exo_size_arg},
    other_endo_size{other_endo_size_arg},
    nb_col_det_exo_jacob{nb_col_det_exo_jacob_arg},
    nb_col_exo_jacob{nb_col_exo_jacob_arg},
    nb_col_other_endo_jacob{nb_col_other_endo_jacob_arg}
  {
  }
  FBEGINBLOCK_(int size_arg, BlockSimulationType type_arg, int first_element, int block_size,
               const vector<int> &variable_arg, const vector<int> &equation_arg,
               bool is_linear_arg, int endo_nbr_arg, int Max_Lag_arg, int Max_Lead_arg, int &u_count_int_arg, int nb_col_jacob_arg) :
    BytecodeInstruction{Tags::FBEGINBLOCK},
    size{size_arg},
    type{type_arg},
    variable{variable_arg.begin()+first_element, variable_arg.begin()+(first_element+block_size)},
    equation{equation_arg.begin()+first_element, equation_arg.begin()+(first_element+block_size)},
    is_linear{is_linear_arg},
    endo_nbr{endo_nbr_arg},
    Max_Lag{Max_Lag_arg},
    Max_Lead{Max_Lead_arg},
    u_count_int{u_count_int_arg},
    nb_col_jacob{nb_col_jacob_arg},
    det_exo_size{0},
    exo_size{0},
    other_endo_size{0},
    nb_col_det_exo_jacob{0},
    nb_col_exo_jacob{0},
    nb_col_other_endo_jacob{0}
  {
  }
  int
  get_size()
  {
    return size;
  };
  BlockSimulationType
  get_type()
  {
    return type;
  };
  bool
  get_is_linear()
  {
    return is_linear;
  };
  int
  get_endo_nbr()
  {
    return endo_nbr;
  };
  int
  get_Max_Lag()
  {
    return Max_Lag;
  };
  int
  get_Max_Lead()
  {
    return Max_Lead;
  };
  int
  get_u_count_int()
  {
    return u_count_int;
  };
  vector<Block_contain_type>
  get_Block_Contain()
  {
    return Block_Contain_;
  };
  int
  get_nb_col_jacob()
  {
    return nb_col_jacob;
  };
  int
  get_exo_size()
  {
    return exo_size;
  };
  int
  get_nb_col_exo_jacob()
  {
    return nb_col_exo_jacob;
  };
  int
  get_det_exo_size()
  {
    return det_exo_size;
  };
  int
  get_nb_col_det_exo_jacob()
  {
    return nb_col_det_exo_jacob;
  };
  int
  get_other_endo_size()
  {
    return other_endo_size;
  };
  int
  get_nb_col_other_endo_jacob()
  {
    return nb_col_other_endo_jacob;
  };
  vector<int>
  get_endogenous()
  {
    return variable;
  }
  vector<int>
  get_exogenous()
  {
    return exogenous;
  }
#ifdef BYTECODE_MEX

  char *
  load(char *code)
  {
    op_code = Tags::FBEGINBLOCK; code += sizeof(op_code);
    memcpy(&size, code, sizeof(size)); code += sizeof(size);
    memcpy(&type, code, sizeof(type)); code += sizeof(type);
    for (int i = 0; i < size; i++)
      {
        Block_contain_type bc;
        memcpy(&bc.Variable, code, sizeof(bc.Variable)); code += sizeof(bc.Variable);
        memcpy(&bc.Equation, code, sizeof(bc.Equation)); code += sizeof(bc.Equation);
        Block_Contain_.push_back(bc);
      }
    if (type == BlockSimulationType::solveTwoBoundariesSimple
        || type == BlockSimulationType::solveTwoBoundariesComplete
        || type == BlockSimulationType::solveBackwardComplete
        || type == BlockSimulationType::solveForwardComplete)
      {
        memcpy(&is_linear, code, sizeof(is_linear)); code += sizeof(is_linear);
        memcpy(&endo_nbr, code, sizeof(endo_nbr)); code += sizeof(endo_nbr);
        memcpy(&Max_Lag, code, sizeof(Max_Lag)); code += sizeof(Max_Lag);
        memcpy(&Max_Lead, code, sizeof(Max_Lead)); code += sizeof(Max_Lead);
        memcpy(&u_count_int, code, sizeof(u_count_int)); code += sizeof(u_count_int);
      }
    memcpy(&nb_col_jacob, code, sizeof(nb_col_jacob)); code += sizeof(nb_col_jacob);
    memcpy(&det_exo_size, code, sizeof(det_exo_size)); code += sizeof(det_exo_size);
    memcpy(&nb_col_det_exo_jacob, code, sizeof(nb_col_det_exo_jacob)); code += sizeof(nb_col_det_exo_jacob);
    memcpy(&exo_size, code, sizeof(exo_size)); code += sizeof(exo_size);
    memcpy(&nb_col_exo_jacob, code, sizeof(nb_col_exo_jacob)); code += sizeof(nb_col_exo_jacob);
    memcpy(&other_endo_size, code, sizeof(other_endo_size)); code += sizeof(other_endo_size);
    memcpy(&nb_col_other_endo_jacob, code, sizeof(nb_col_other_endo_jacob)); code += sizeof(nb_col_other_endo_jacob);

    for (int i{0}; i < det_exo_size; i++)
      {
        int tmp_i;
        memcpy(&tmp_i, code, sizeof(tmp_i)); code += sizeof(tmp_i);
        det_exogenous.push_back(tmp_i);
      }
    for (int i{0}; i < exo_size; i++)
      {
        int tmp_i;
        memcpy(&tmp_i, code, sizeof(tmp_i)); code += sizeof(tmp_i);
        exogenous.push_back(tmp_i);
      }
    for (int i{0}; i < other_endo_size; i++)
      {
        int tmp_i;
        memcpy(&tmp_i, code, sizeof(tmp_i)); code += sizeof(tmp_i);
        other_endogenous.push_back(tmp_i);
      }
    return code;
  };
#endif
};

// Superclass of std::ofstream for writing a sequence of bytecode instructions
class BytecodeWriter : private ofstream
{
  template<typename B>
  friend BytecodeWriter &operator<<(BytecodeWriter &code_file, const B &instr);
private:
  // Stores the positions of all instructions in the byte stream
  vector<pos_type> instructions_positions;
public:
  BytecodeWriter(const string &filename);
  // Returns the number of the next instruction to be written
  int
  getInstructionCounter() const
  {
    return static_cast<int>(instructions_positions.size());
  }
  /* Overwrites an existing instruction, given its number.
     It is the responsibility of the caller to ensure that the new instruction
     occupies exactly as many bytes as the former one. */
  template<typename B>
  void
  overwriteInstruction(int instruction_number, const B &new_instruction)
  {
    seekp(instructions_positions.at(instruction_number));
    *this << new_instruction;
    instructions_positions.pop_back();
    seekp(0, ios_base::end);
  }
};

// Overloads of operator<< for writing bytecode instructions

template<typename B>
BytecodeWriter &
operator<<(BytecodeWriter &code_file, const B &instr)
{
  code_file.instructions_positions.push_back(code_file.tellp());
  code_file.write(reinterpret_cast<const char *>(&instr), sizeof(B));
  return code_file;
}

template<>
BytecodeWriter &operator<<(BytecodeWriter &code_file, const FCALL_ &instr);

template<>
BytecodeWriter &operator<<(BytecodeWriter &code_file, const FBEGINBLOCK_ &instr);


#ifdef BYTECODE_MEX
using tags_liste_t = vector<pair<Tags, void * >>;
class CodeLoad
{
private:
  char *code;
  int nb_blocks;
  vector<size_t> begin_block;
public:

  int
  get_block_number() const
  {
    return nb_blocks;
  };

  size_t
  get_begin_block(int block) const
  {
    return begin_block[block];
  }
  void *
  get_current_code() const
  {
    return code;
  };
  tags_liste_t
  get_op_code(string file_name)
  {
    tags_liste_t tags_liste;
    ifstream CompiledCode;
    streamoff Code_Size;
    CompiledCode.open(file_name + ".cod", ios::in | ios::binary| ios::ate);
    if (!CompiledCode.is_open())
      return tags_liste;
    Code_Size = CompiledCode.tellg();
    CompiledCode.seekg(ios::beg);
    code = static_cast<char *>(mxMalloc(Code_Size));
    CompiledCode.seekg(0);
    CompiledCode.read(reinterpret_cast<char *>(code), Code_Size);
    CompiledCode.close();
    nb_blocks = 0;
    bool done = false;
    int instruction = 0;
    while (!done)
      {
        switch (*reinterpret_cast<Tags *>(code))
          {
          case Tags::FLDZ:
# ifdef DEBUGL
            mexPrintf("FLDZ = %d size = %d\n", Tags::FLDZ, sizeof(FLDZ_));
# endif
            tags_liste.emplace_back(Tags::FLDZ, code);
            code += sizeof(FLDZ_);
            break;
          case Tags::FEND:
# ifdef DEBUGL
            mexPrintf("FEND\n");
# endif
            tags_liste.emplace_back(Tags::FEND, code);
            code += sizeof(FEND_);
            done = true;
            break;
          case Tags::FENDBLOCK:
# ifdef DEBUGL
            mexPrintf("FENDBLOCK\n");
# endif
            tags_liste.emplace_back(Tags::FENDBLOCK, code);
            code += sizeof(FENDBLOCK_);
            break;
          case Tags::FENDEQU:
# ifdef DEBUGL
            mexPrintf("FENDEQU\n");
# endif
            tags_liste.emplace_back(Tags::FENDEQU, code);
            code += sizeof(FENDEQU_);
            break;
          case Tags::FCUML:
# ifdef DEBUGL
            mexPrintf("FCUML\n");
# endif
            tags_liste.emplace_back(Tags::FCUML, code);
            code += sizeof(FCUML_);
            break;
          case Tags::FDIMT:
# ifdef DEBUGL
            mexPrintf("FDIMT = %d size = %d\n", Tags::FDIMT, sizeof(FDIMT_));
# endif
            tags_liste.emplace_back(Tags::FDIMT, code);
            code += sizeof(FDIMT_);
            break;
          case Tags::FDIMST:
# ifdef DEBUGL
            mexPrintf("FDIMST\n");
# endif
            tags_liste.emplace_back(Tags::FDIMST, code);
            code += sizeof(FDIMST_);
            break;
          case Tags::FNUMEXPR:
# ifdef DEBUGL
            mexPrintf("FNUMEXPR\n");
# endif
            tags_liste.emplace_back(Tags::FNUMEXPR, code);
            code += sizeof(FNUMEXPR_);
            break;
          case Tags::FLDC:
# ifdef DEBUGL
            mexPrintf("FLDC\n");
# endif
            tags_liste.emplace_back(Tags::FLDC, code);
            code += sizeof(FLDC_);
            break;
          case Tags::FLDU:
# ifdef DEBUGL
            mexPrintf("FLDU\n");
# endif
            tags_liste.emplace_back(Tags::FLDU, code);
            code += sizeof(FLDU_);
            break;
          case Tags::FLDSU:
# ifdef DEBUGL
            mexPrintf("FLDSU\n");
# endif
            tags_liste.emplace_back(Tags::FLDSU, code);
            code += sizeof(FLDSU_);
            break;
          case Tags::FLDR:
# ifdef DEBUGL
            mexPrintf("FLDR\n");
# endif
            tags_liste.emplace_back(Tags::FLDR, code);
            code += sizeof(FLDR_);
            break;
          case Tags::FLDT:
# ifdef DEBUGL
            mexPrintf("FLDT\n");
# endif
            tags_liste.emplace_back(Tags::FLDT, code);
            code += sizeof(FLDT_);
            break;
          case Tags::FLDST:
# ifdef DEBUGL
            mexPrintf("FLDST\n");
# endif
            tags_liste.emplace_back(Tags::FLDST, code);
            code += sizeof(FLDST_);
            break;
          case Tags::FSTPT:
# ifdef DEBUGL
            mexPrintf("FSTPT = %d size = %d\n", Tags::FSTPT, sizeof(FSTPT_));
# endif
            tags_liste.emplace_back(Tags::FSTPT, code);
            code += sizeof(FSTPT_);
            break;
          case Tags::FSTPST:
# ifdef DEBUGL
            mexPrintf("FSTPST\n");
# endif
            tags_liste.emplace_back(Tags::FSTPST, code);
            code += sizeof(FSTPST_);
            break;
          case Tags::FSTPR:
# ifdef DEBUGL
            mexPrintf("FSTPR\n");
# endif
            tags_liste.emplace_back(Tags::FSTPR, code);
            code += sizeof(FSTPR_);
            break;
          case Tags::FSTPU:
# ifdef DEBUGL
            mexPrintf("FSTPU\n");
# endif
            tags_liste.emplace_back(Tags::FSTPU, code);
            code += sizeof(FSTPU_);
            break;
          case Tags::FSTPSU:
# ifdef DEBUGL
            mexPrintf("FSTPSU\n");
# endif
            tags_liste.emplace_back(Tags::FSTPSU, code);
            code += sizeof(FSTPSU_);
            break;
          case Tags::FSTPG:
# ifdef DEBUGL
            mexPrintf("FSTPG\n");
# endif
            tags_liste.emplace_back(Tags::FSTPG, code);
            code += sizeof(FSTPG_);
            break;
          case Tags::FSTPG2:
# ifdef DEBUGL
            mexPrintf("FSTPG2\n");
# endif
            tags_liste.emplace_back(Tags::FSTPG2, code);
            code += sizeof(FSTPG2_);
            break;
          case Tags::FSTPG3:
# ifdef DEBUGL
            mexPrintf("FSTPG3\n");
# endif
            tags_liste.emplace_back(Tags::FSTPG3, code);
            code += sizeof(FSTPG3_);
            break;
          case Tags::FUNARY:
# ifdef DEBUGL
            mexPrintf("FUNARY\n");
# endif
            tags_liste.emplace_back(Tags::FUNARY, code);
            code += sizeof(FUNARY_);
            break;
          case Tags::FBINARY:
# ifdef DEBUGL
            mexPrintf("FBINARY\n");
# endif
            tags_liste.emplace_back(Tags::FBINARY, code);
            code += sizeof(FBINARY_);
            break;
          case Tags::FTRINARY:
# ifdef DEBUGL
            mexPrintf("FTRINARY\n");
# endif
            tags_liste.emplace_back(Tags::FTRINARY, code);
            code += sizeof(FTRINARY_);
            break;
          case Tags::FLDVS:
# ifdef DEBUGL
            mexPrintf("FLDVS\n");
# endif
            tags_liste.emplace_back(Tags::FLDVS, code);
            code += sizeof(FLDVS_);
            break;
          case Tags::FLDSV:
# ifdef DEBUGL
            mexPrintf("FLDSV\n");
# endif
            tags_liste.emplace_back(Tags::FLDSV, code);
            code += sizeof(FLDSV_);
            break;
          case Tags::FSTPSV:
# ifdef DEBUGL
            mexPrintf("FSTPSV\n");
# endif
            tags_liste.emplace_back(Tags::FSTPSV, code);
            code += sizeof(FSTPSV_);
            break;
          case Tags::FLDV:
# ifdef DEBUGL
            mexPrintf("FLDV\n");
# endif
            tags_liste.emplace_back(Tags::FLDV, code);
            code += sizeof(FLDV_);
            break;
          case Tags::FSTPV:
# ifdef DEBUGL
            mexPrintf("FSTPV\n");
# endif
            tags_liste.emplace_back(Tags::FSTPV, code);
            code += sizeof(FSTPV_);
            break;
          case Tags::FBEGINBLOCK:
# ifdef DEBUGL
            mexPrintf("FBEGINBLOCK\n");
# endif
            {
              // TODO: remove default FBEGINBLOCK_ constructor when the following is remove
              auto *fbegin_block = new FBEGINBLOCK_;

              code = fbegin_block->load(code);

              begin_block.push_back(tags_liste.size());
              tags_liste.emplace_back(Tags::FBEGINBLOCK, fbegin_block);
              nb_blocks++;
            }
            break;
          case Tags::FJMPIFEVAL:
# ifdef DEBUGL
            mexPrintf("FJMPIFEVAL\n");
# endif
            tags_liste.emplace_back(Tags::FJMPIFEVAL, code);
            code += sizeof(FJMPIFEVAL_);
            break;
          case Tags::FJMP:
# ifdef DEBUGL
            mexPrintf("FJMP\n");
# endif
            tags_liste.emplace_back(Tags::FJMP, code);
            code += sizeof(FJMP_);
            break;
          case Tags::FCALL:
            {
# ifdef DEBUGL
              mexPrintf("FCALL\n");
# endif
              // TODO: remove default FCALL_ constructor when the following is remove
              auto *fcall = new FCALL_;

              code = fcall->load(code);

              tags_liste.emplace_back(Tags::FCALL, fcall);
# ifdef DEBUGL
              mexPrintf("FCALL finish\n"); mexEvalString("drawnow;");
              mexPrintf("-- *code=%d\n", *code); mexEvalString("drawnow;");
# endif
            }
            break;
          case Tags::FLDTEF:
# ifdef DEBUGL
            mexPrintf("FLDTEF\n");
# endif
            tags_liste.emplace_back(Tags::FLDTEF, code);
            code += sizeof(FLDTEF_);
            break;
          case Tags::FSTPTEF:
# ifdef DEBUGL
            mexPrintf("FSTPTEF\n");
# endif
            tags_liste.emplace_back(Tags::FSTPTEF, code);
            code += sizeof(FSTPTEF_);
            break;
          case Tags::FLDTEFD:
# ifdef DEBUGL
            mexPrintf("FLDTEFD\n");
# endif
            tags_liste.emplace_back(Tags::FLDTEFD, code);
            code += sizeof(FLDTEFD_);
            break;
          case Tags::FSTPTEFD:
# ifdef DEBUGL
            mexPrintf("FSTPTEFD\n");
# endif
            tags_liste.emplace_back(Tags::FSTPTEFD, code);
            code += sizeof(FSTPTEFD_);
            break;
          case Tags::FLDTEFDD:
# ifdef DEBUGL
            mexPrintf("FLDTEFDD\n");
# endif
            tags_liste.emplace_back(Tags::FLDTEFDD, code);
            code += sizeof(FLDTEFDD_);
            break;
          case Tags::FSTPTEFDD:
# ifdef DEBUGL
            mexPrintf("FSTPTEFDD\n");
# endif
            tags_liste.emplace_back(Tags::FSTPTEFDD, code);
            code += sizeof(FSTPTEFDD_);
            break;
          default:
            mexPrintf("Unknown Tag value=%d code=%x\n", *code, code);
            done = true;
          }
        instruction++;
      }
    return tags_liste;
  };
};
#endif // BYTECODE_MEX

#endif // _BYTECODE_HH
