/*
 * Copyright © 2007-2023 Dynare Team
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
#include <filesystem>
#include <type_traits>

#include "CommonEnums.hh"

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

struct BytecodeInstruction
{
  const Tags op_code;
  explicit BytecodeInstruction(Tags op_code_arg) :
    op_code {op_code_arg}
  {
  }
protected:
  /* This is a base class, so the destructor should be either public+virtual or
     protected+non-virtual. We opt for the latter, because otherwise this class
     would no longer be POD; its memory representation would also include
     runtime type information, and our crude serialization technique (copying the
     whole object from memory) would thus not work. */
  ~BytecodeInstruction() = default;
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
protected:
  // See BytecodeInstruction destructor for the rationale
  ~TagWithOneArgument() = default;
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
protected:
  // See BytecodeInstruction destructor for the rationale
  ~TagWithTwoArguments() = default;
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
protected:
  // See BytecodeInstruction destructor for the rationale
  ~TagWithThreeArguments() = default;
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
protected:
  // See BytecodeInstruction destructor for the rationale
  ~TagWithFourArguments() = default;
};

class FLDZ_ final : public BytecodeInstruction
{
public:
  FLDZ_() : BytecodeInstruction{Tags::FLDZ}
  {
  };
};

class FEND_ final : public BytecodeInstruction
{
public:
  FEND_() : BytecodeInstruction{Tags::FEND}
  {
  };
};

class FENDBLOCK_ final : public BytecodeInstruction
{
public:
  FENDBLOCK_() : BytecodeInstruction{Tags::FENDBLOCK}
  {
  };
};

class FENDEQU_ final : public BytecodeInstruction
{
public:
  FENDEQU_() : BytecodeInstruction{Tags::FENDEQU}
  {
  };
};

class FDIMT_ final : public TagWithOneArgument<int>
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

class FDIMST_ final : public TagWithOneArgument<int>
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

class FLDC_ final : public TagWithOneArgument<double>
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

class FLDU_ final : public TagWithOneArgument<int>
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

class FLDSU_ final : public TagWithOneArgument<int>
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

class FLDR_ final : public TagWithOneArgument<int>
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

class FLDT_ final : public TagWithOneArgument<int>
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

class FLDST_ final : public TagWithOneArgument<int>
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

class FSTPT_ final : public TagWithOneArgument<int>
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

class FSTPST_ final : public TagWithOneArgument<int>
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

class FSTPR_ final : public TagWithOneArgument<int>
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

class FSTPU_ final : public TagWithOneArgument<int>
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

class FSTPSU_ final : public TagWithOneArgument<int>
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

class FSTPG_ final : public TagWithOneArgument<int>
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

class FSTPG2_ final : public TagWithTwoArguments<int, int>
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

class FSTPG3_ final : public TagWithFourArguments<int, int, int, int>
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

class FUNARY_ final : public TagWithOneArgument<UnaryOpcode>
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

class FBINARY_ final : public TagWithOneArgument<BinaryOpcode>
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

class FTRINARY_ final : public TagWithOneArgument<TrinaryOpcode>
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

class FJMPIFEVAL_ final : public TagWithOneArgument<int>
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

class FJMP_ final : public TagWithOneArgument<int>
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

class FLDTEF_ final : public TagWithOneArgument<int>
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

class FSTPTEF_ final : public TagWithOneArgument<int>
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

class FLDTEFD_ final : public TagWithTwoArguments<int, int>
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

class FSTPTEFD_ final : public TagWithTwoArguments<int, int>
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

class FLDTEFDD_ final : public TagWithThreeArguments<int, int, int>
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

class FSTPTEFDD_ final : public TagWithThreeArguments<int, int, int>
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

class FLDVS_ final : public TagWithTwoArguments<SymbolType, int>
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

class FLDSV_ final : public TagWithTwoArguments<SymbolType, int>
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

class FSTPSV_ final : public TagWithTwoArguments<SymbolType, int>
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

class FLDV_ final : public TagWithThreeArguments<SymbolType, int, int>
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

class FSTPV_ final : public TagWithThreeArguments<SymbolType, int, int>
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

class FCALL_ final : public BytecodeInstruction
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
  FCALL_(int nb_output_arguments_arg, int nb_input_arguments_arg, string func_name_arg, int indx_arg, ExternalFunctionCallType call_type_arg) :
    BytecodeInstruction{Tags::FCALL},
    nb_output_arguments{nb_output_arguments_arg},
    nb_input_arguments{nb_input_arguments_arg},
    indx{indx_arg},
    func_name{move(func_name_arg)},
    call_type{call_type_arg}
  {
  };
  /* Deserializing constructor.
     Updates the code pointer to point beyond the bytes read. */
  FCALL_(char *&code) :
    BytecodeInstruction{Tags::FCALL}
  {
    code += sizeof(op_code);

    auto read_member = [&code](auto &member)
    {
      member = *reinterpret_cast<add_pointer_t<decltype(member)>>(code);
      code += sizeof member;
    };

    read_member(nb_output_arguments);
    read_member(nb_input_arguments);
    read_member(indx);
    read_member(add_input_arguments);
    read_member(row);
    read_member(col);
    read_member(call_type);

    int size;
    read_member(size);
    func_name = code;
    code += size+1;

    read_member(size);
    arg_func_name = code;
    code += size+1;
  }

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
};

class FNUMEXPR_ final : public BytecodeInstruction
{
private:
  ExpressionType expression_type;
  int equation; // Equation number (non-block-specific) (or temporary term number for ExpressionType::TemporaryTerm)
  int dvariable1; // For derivatives, type-specific ID of the derivation variable
  int lag1; // For derivatives, lead/lag of the derivation variable
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

class FBEGINBLOCK_ final : public BytecodeInstruction
{
  template<typename B>
  friend BytecodeWriter &operator<<(BytecodeWriter &code_file, const B &instr);
private:
  int size{0};
  BlockSimulationType type;
  vector<int> variable;
  vector<int> equation;
  vector<int> exogenous;
  vector<int> det_exogenous;
  bool is_linear{false};
  vector<Block_contain_type> Block_Contain_;
  int u_count_int{0};
  int nb_col_jacob{0};
  int det_exo_size, exo_size;
public:
  /* Constructor when derivatives w.r.t. exogenous are present (only makes
     sense when there is no block-decomposition, since there is no provision for
     derivatives w.r.t. endogenous not belonging to the block) */
  FBEGINBLOCK_(int size_arg, BlockSimulationType type_arg, int first_element, int block_size,
               const vector<int> &variable_arg, const vector<int> &equation_arg,
               bool is_linear_arg, int u_count_int_arg, int nb_col_jacob_arg,
               int det_exo_size_arg, int exo_size_arg,
               vector<int> det_exogenous_arg, vector<int> exogenous_arg) :
    BytecodeInstruction{Tags::FBEGINBLOCK},
    size{size_arg},
    type{type_arg},
    variable{variable_arg.begin()+first_element, variable_arg.begin()+(first_element+block_size)},
    equation{equation_arg.begin()+first_element, equation_arg.begin()+(first_element+block_size)},
    exogenous{move(exogenous_arg)},
    det_exogenous{move(det_exogenous_arg)},
    is_linear{is_linear_arg},
    u_count_int{u_count_int_arg},
    nb_col_jacob{nb_col_jacob_arg},
    det_exo_size{det_exo_size_arg},
    exo_size{exo_size_arg}
  {
  }
  // Constructor when derivatives w.r.t. exogenous are absent
  FBEGINBLOCK_(int size_arg, BlockSimulationType type_arg, int first_element, int block_size,
               const vector<int> &variable_arg, const vector<int> &equation_arg,
               bool is_linear_arg, int u_count_int_arg, int nb_col_jacob_arg) :
    BytecodeInstruction{Tags::FBEGINBLOCK},
    size{size_arg},
    type{type_arg},
    variable{variable_arg.begin()+first_element, variable_arg.begin()+(first_element+block_size)},
    equation{equation_arg.begin()+first_element, equation_arg.begin()+(first_element+block_size)},
    is_linear{is_linear_arg},
    u_count_int{u_count_int_arg},
    nb_col_jacob{nb_col_jacob_arg},
    det_exo_size{0},
    exo_size{0}
  {
  }
  /* Deserializing constructor.
     Updates the code pointer to point beyond the bytes read. */
  FBEGINBLOCK_(char *&code) :
    BytecodeInstruction{Tags::FBEGINBLOCK}
  {
    code += sizeof(op_code);

    auto read_member = [&code](auto &member)
    {
      member = *reinterpret_cast<add_pointer_t<decltype(member)>>(code);
      code += sizeof member;
    };

    read_member(size);
    read_member(type);
    for (int i {0}; i < size; i++)
      {
        Block_contain_type bc;
        read_member(bc.Variable);
        read_member(bc.Equation);
        Block_Contain_.push_back(move(bc));
      }
    if (type == BlockSimulationType::solveTwoBoundariesSimple
        || type == BlockSimulationType::solveTwoBoundariesComplete
        || type == BlockSimulationType::solveBackwardComplete
        || type == BlockSimulationType::solveForwardComplete)
      {
        read_member(is_linear);
        read_member(u_count_int);
      }
    read_member(nb_col_jacob);
    read_member(det_exo_size);
    read_member(exo_size);

    for (int i {0}; i < det_exo_size; i++)
      {
        int tmp_i;
        read_member(tmp_i);
        det_exogenous.push_back(tmp_i);
      }
    for (int i {0}; i < exo_size; i++)
      {
        int tmp_i;
        read_member(tmp_i);
        exogenous.push_back(tmp_i);
      }
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
  get_det_exo_size()
  {
    return det_exo_size;
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
  BytecodeWriter(const filesystem::path &filename);
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

#endif // _BYTECODE_HH
