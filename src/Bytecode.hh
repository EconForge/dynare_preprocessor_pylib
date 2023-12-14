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

#ifndef BYTECODE_HH
#define BYTECODE_HH

#include <filesystem>
#include <fstream>
#include <ios>
#include <type_traits>
#include <utility>
#include <vector>

#include "CommonEnums.hh"

using namespace std;

namespace Bytecode
{

// The different tags encoding a bytecode instruction
enum class Tag
{
  FLDZ, // Loads a zero onto the stack
  FLDC, // Loads a constant term onto the stack

  FDIMT, // Defines the number of temporary terms - dynamic context (the period has to be indicated)
  FDIMST, // Defines the number of temporary terms - static context (the period hasn’t to be
          // indicated)
  FLDT,  // Loads a temporary term onto the stack - dynamic context (the period has to be indicated)
  FLDST, // Loads a temporary term onto the stack - static context (the period hasn’t to be
         // indicated)
  FSTPT, // Stores a temporary term from the stack - dynamic context (the period has to be
         // indicated)
  FSTPST, // Stores a temporary term from the stack - static context (the period hasn’t to be
          // indicated)

  FLDU,  // Loads an element of the vector U onto the stack - dynamic context (the period has to be
         // indicated)
  FLDSU, // Loads an element of the vector U onto the stack - static context (the period hasn’t to
         // be indicated)
  FSTPU, // Stores an element of the vector U from the stack - dynamic context (the period has to be
         // indicated)
  FSTPSU, // Stores an element of the vector U from the stack - static context (the period hasn’t to
          // be indicated)

  FLDV,  // Loads a variable (described in SymbolType) onto the stack - dynamic context (the period
         // has to be indicated)
  FLDSV, // Loads a variable (described in SymbolType) onto the stack - static context (the period
         // hasn’t to be indicated)
  FLDVS, // Loads a variable (described in SymbolType) onto the stack - dynamic context but inside
         // the STEADY_STATE operator (the period hasn’t to be indicated)
  FSTPV, // Stores a variable (described in SymbolType) from the stack - dynamic context (the period
         // has to be indicated)
  FSTPSV, // Stores a variable (described in SymbolType) from the stack - static context (the period
          // hasn’t to be indicated)

  FLDR,  // Loads a residual onto the stack
  FSTPR, // Stores a residual from the stack

  FSTPG,  // Stores a derivative from the stack
  FSTPG2, // Stores a derivative matrix for a static model from the stack
  FSTPG3, // Stores a derivative matrix for a dynamic model from the stack

  FUNARY,   // A unary operator
  FBINARY,  // A binary operator
  FTRINARY, // A trinary operator

  FJMPIFEVAL, // Jump if evaluate = true
  FJMP,       // Jump

  FBEGINBLOCK, // Marks the beginning of a model block
  FENDBLOCK,   // Marks the end of a model block
  FENDEQU,     // Marks the last equation of the block; for a block that has to be solved, the
               // derivatives appear just after this flag
  FEND,        // Marks the end of the model code

  FNUMEXPR, // Stores the expression type and references

  FCALL,    // Call an external function
  FLDTEF,   // Loads the result of an external function onto the stack
  FSTPTEF,  // Stores the result of an external function from the stack
  FLDTEFD,  // Loads the result of the 1st derivative of an external function onto the stack
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

class Writer;

struct Instruction
{
  const Tag tag;
  explicit Instruction(Tag tag_arg) : tag {tag_arg}
  {
  }

protected:
  /* This is a base class, so the destructor should be either public+virtual or
     protected+non-virtual. We opt for the latter, because otherwise this class
     would no longer be POD; its memory representation would also include
     runtime type information, and our crude serialization technique (copying the
     whole object from memory) would thus not work. */
  ~Instruction() = default;
};

template<typename T1>
class InstructionWithOneArgument : public Instruction
{
protected:
  T1 arg1;

public:
  InstructionWithOneArgument(Tag tag_arg, T1 arg_arg1) : Instruction {tag_arg}, arg1 {arg_arg1}
  {
  }

protected:
  // See Instruction destructor for the rationale
  ~InstructionWithOneArgument() = default;
};

template<typename T1, typename T2>
class InstructionWithTwoArguments : public Instruction
{
protected:
  T1 arg1;
  T2 arg2;

public:
  InstructionWithTwoArguments(Tag tag_arg, T1 arg_arg1, T2 arg_arg2) :
      Instruction {tag_arg}, arg1 {arg_arg1}, arg2 {arg_arg2}
  {
  }

protected:
  // See Instruction destructor for the rationale
  ~InstructionWithTwoArguments() = default;
};

template<typename T1, typename T2, typename T3>
class InstructionWithThreeArguments : public Instruction
{
protected:
  T1 arg1;
  T2 arg2;
  T3 arg3;

public:
  InstructionWithThreeArguments(Tag tag_arg, T1 arg_arg1, T2 arg_arg2, T3 arg_arg3) :
      Instruction {tag_arg}, arg1 {arg_arg1}, arg2 {arg_arg2}, arg3 {arg_arg3}
  {
  }

protected:
  // See Instruction destructor for the rationale
  ~InstructionWithThreeArguments() = default;
};

template<typename T1, typename T2, typename T3, typename T4>
class InstructionWithFourArguments : public Instruction
{
protected:
  T1 arg1;
  T2 arg2;
  T3 arg3;
  T4 arg4;

public:
  InstructionWithFourArguments(Tag tag_arg, T1 arg_arg1, T2 arg_arg2, T3 arg_arg3, T4 arg_arg4) :
      Instruction {tag_arg},
      arg1 {arg_arg1},
      arg2 {arg_arg2},
      arg3 {move(arg_arg3)},
      arg4 {arg_arg4}
  {
  }

protected:
  // See Instruction destructor for the rationale
  ~InstructionWithFourArguments() = default;
};

class FLDZ final : public Instruction
{
public:
  FLDZ() : Instruction {Tag::FLDZ}
  {
  }
};

class FEND final : public Instruction
{
public:
  FEND() : Instruction {Tag::FEND}
  {
  }
};

class FENDBLOCK final : public Instruction
{
public:
  FENDBLOCK() : Instruction {Tag::FENDBLOCK}
  {
  }
};

class FENDEQU final : public Instruction
{
public:
  FENDEQU() : Instruction {Tag::FENDEQU}
  {
  }
};

class FDIMT final : public InstructionWithOneArgument<int>
{
public:
  explicit FDIMT(int size_arg) : InstructionWithOneArgument {Tag::FDIMT, size_arg}
  {
  }
  int
  get_size()
  {
    return arg1;
  };
};

class FDIMST final : public InstructionWithOneArgument<int>
{
public:
  explicit FDIMST(int size_arg) : InstructionWithOneArgument {Tag::FDIMST, size_arg}
  {
  }
  int
  get_size()
  {
    return arg1;
  };
};

class FLDC final : public InstructionWithOneArgument<double>
{
public:
  explicit FLDC(double value_arg) : InstructionWithOneArgument {Tag::FLDC, value_arg}
  {
  }
  double
  get_value()
  {
    return arg1;
  };
};

class FLDU final : public InstructionWithOneArgument<int>
{
public:
  explicit FLDU(int pos_arg) : InstructionWithOneArgument {Tag::FLDU, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDSU final : public InstructionWithOneArgument<int>
{
public:
  explicit FLDSU(int pos_arg) : InstructionWithOneArgument {Tag::FLDSU, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDR final : public InstructionWithOneArgument<int>
{
public:
  explicit FLDR(int pos_arg) : InstructionWithOneArgument {Tag::FLDR, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDT final : public InstructionWithOneArgument<int>
{
public:
  explicit FLDT(int pos_arg) : InstructionWithOneArgument {Tag::FLDT, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FLDST final : public InstructionWithOneArgument<int>
{
public:
  explicit FLDST(int pos_arg) : InstructionWithOneArgument {Tag::FLDST, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPT final : public InstructionWithOneArgument<int>
{
public:
  explicit FSTPT(int pos_arg) : InstructionWithOneArgument {Tag::FSTPT, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPST final : public InstructionWithOneArgument<int>
{
public:
  explicit FSTPST(int pos_arg) : InstructionWithOneArgument {Tag::FSTPST, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPR final : public InstructionWithOneArgument<int>
{
public:
  explicit FSTPR(int pos_arg) : InstructionWithOneArgument {Tag::FSTPR, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPU final : public InstructionWithOneArgument<int>
{
public:
  explicit FSTPU(int pos_arg) : InstructionWithOneArgument {Tag::FSTPU, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPSU final : public InstructionWithOneArgument<int>
{
public:
  explicit FSTPSU(int pos_arg) : InstructionWithOneArgument {Tag::FSTPSU, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPG final : public InstructionWithOneArgument<int>
{
public:
  explicit FSTPG(int pos_arg) : InstructionWithOneArgument {Tag::FSTPG, pos_arg}
  {
  }
  int
  get_pos()
  {
    return arg1;
  };
};

class FSTPG2 final : public InstructionWithTwoArguments<int, int>
{
public:
  FSTPG2(int row_arg, int col_arg) : InstructionWithTwoArguments {Tag::FSTPG2, row_arg, col_arg}
  {
  }
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

class FSTPG3 final : public InstructionWithFourArguments<int, int, int, int>
{
public:
  FSTPG3(int row_arg, int col_arg, int lag_arg, int col_pos_arg) :
      InstructionWithFourArguments {Tag::FSTPG3, row_arg, col_arg, lag_arg, col_pos_arg}
  {
  }
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
    return arg3;
  };
  int
  get_col_pos()
  {
    return arg4;
  };
};

class FUNARY final : public InstructionWithOneArgument<UnaryOpcode>
{
public:
  explicit FUNARY(UnaryOpcode op_type_arg) : InstructionWithOneArgument {Tag::FUNARY, op_type_arg}
  {
  }
  UnaryOpcode
  get_op_type()
  {
    return arg1;
  };
};

class FBINARY final : public InstructionWithOneArgument<BinaryOpcode>
{
public:
  explicit FBINARY(BinaryOpcode op_type_arg) :
      InstructionWithOneArgument {Tag::FBINARY, op_type_arg}
  {
  }
  BinaryOpcode
  get_op_type()
  {
    return arg1;
  };
};

class FTRINARY final : public InstructionWithOneArgument<TrinaryOpcode>
{
public:
  explicit FTRINARY(TrinaryOpcode op_type_arg) :
      InstructionWithOneArgument {Tag::FTRINARY, op_type_arg}
  {
  }
  TrinaryOpcode
  get_op_type()
  {
    return arg1;
  };
};

class FJMPIFEVAL final : public InstructionWithOneArgument<int>
{
public:
  explicit FJMPIFEVAL(int arg_pos) : InstructionWithOneArgument {Tag::FJMPIFEVAL, arg_pos}
  {
  }
  int
  get_pos()
  {
    return arg1;
  }
};

class FJMP final : public InstructionWithOneArgument<int>
{
public:
  explicit FJMP(int arg_pos) : InstructionWithOneArgument {Tag::FJMP, arg_pos}
  {
  }
  int
  get_pos()
  {
    return arg1;
  }
};

class FLDTEF final : public InstructionWithOneArgument<int>
{
public:
  explicit FLDTEF(int number) : InstructionWithOneArgument {Tag::FLDTEF, number}
  {
  }
  int
  get_number()
  {
    return arg1;
  }
};

class FSTPTEF final : public InstructionWithOneArgument<int>
{
public:
  explicit FSTPTEF(int number) : InstructionWithOneArgument {Tag::FSTPTEF, number}
  {
  }
  int
  get_number()
  {
    return arg1;
  }
};

class FLDTEFD final : public InstructionWithTwoArguments<int, int>
{
public:
  FLDTEFD(int indx, int row) : InstructionWithTwoArguments {Tag::FLDTEFD, indx, row}
  {
  }
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

class FSTPTEFD final : public InstructionWithTwoArguments<int, int>
{
public:
  FSTPTEFD(int indx, int row) : InstructionWithTwoArguments {Tag::FSTPTEFD, indx, row}
  {
  }
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

class FLDTEFDD final : public InstructionWithThreeArguments<int, int, int>
{
public:
  FLDTEFDD(int indx, int row, int col) :
      InstructionWithThreeArguments {Tag::FLDTEFDD, indx, row, col}
  {
  }
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

class FSTPTEFDD final : public InstructionWithThreeArguments<int, int, int>
{
public:
  FSTPTEFDD(int indx, int row, int col) :
      InstructionWithThreeArguments {Tag::FSTPTEF, indx, row, col}
  {
  }
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

class FLDVS final : public InstructionWithTwoArguments<SymbolType, int>
{
public:
  FLDVS(SymbolType type_arg, int pos_arg) :
      InstructionWithTwoArguments {Tag::FLDVS, type_arg, pos_arg}
  {
  }
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

class FLDSV final : public InstructionWithTwoArguments<SymbolType, int>
{
public:
  FLDSV(SymbolType type_arg, int pos_arg) :
      InstructionWithTwoArguments {Tag::FLDSV, type_arg, pos_arg}
  {
  }
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

class FSTPSV final : public InstructionWithTwoArguments<SymbolType, int>
{
public:
  FSTPSV(SymbolType type_arg, int pos_arg) :
      InstructionWithTwoArguments {Tag::FSTPSV, type_arg, pos_arg}
  {
  }
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

class FLDV final : public InstructionWithThreeArguments<SymbolType, int, int>
{
public:
  FLDV(SymbolType type_arg, int pos_arg, int lead_lag_arg) :
      InstructionWithThreeArguments {Tag::FLDV, type_arg, pos_arg, lead_lag_arg}
  {
  }
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

class FSTPV final : public InstructionWithThreeArguments<SymbolType, int, int>
{
public:
  FSTPV(SymbolType type_arg, int pos_arg, int lead_lag_arg) :
      InstructionWithThreeArguments {Tag::FSTPV, type_arg, pos_arg, lead_lag_arg}
  {
  }
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

class FCALL final : public Instruction
{
  template<typename B>
  friend Writer& operator<<(Writer& code_file, const B& instr);

private:
  int nb_output_arguments, nb_input_arguments, indx;
  string func_name;
  string arg_func_name;
  int add_input_arguments {0}, row {0}, col {0};
  ExternalFunctionCallType call_type;

public:
  FCALL(int nb_output_arguments_arg, int nb_input_arguments_arg, string func_name_arg, int indx_arg,
        ExternalFunctionCallType call_type_arg) :
      Instruction {Tag::FCALL},
      nb_output_arguments {nb_output_arguments_arg},
      nb_input_arguments {nb_input_arguments_arg},
      indx {indx_arg},
      func_name {move(func_name_arg)},
      call_type {call_type_arg}
  {
  }
  /* Deserializing constructor.
     Updates the code pointer to point beyond the bytes read. */
  FCALL(char*& code) : Instruction {Tag::FCALL}
  {
    code += sizeof(tag);

    auto read_member = [&code](auto& member) {
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
    code += size + 1;

    read_member(size);
    arg_func_name = code;
    code += size + 1;
  }

  string
  get_function_name()
  {
    // printf("get_function_name => func_name=%s\n",func_name.c_str());fflush(stdout);
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
    arg_func_name = move(arg_arg_func_name);
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

class FNUMEXPR final : public Instruction
{
private:
  ExpressionType expression_type;
  int equation;   // Equation number (non-block-specific) (or temporary term number for
                  // ExpressionType::TemporaryTerm)
  int dvariable1; // For derivatives, type-specific ID of the derivation variable
  int lag1;       // For derivatives, lead/lag of the derivation variable
public:
  FNUMEXPR(const ExpressionType expression_type_arg, int equation_arg) :
      Instruction {Tag::FNUMEXPR},
      expression_type {expression_type_arg},
      equation {equation_arg},
      dvariable1 {0},
      lag1 {0}
  {
  }
  FNUMEXPR(const ExpressionType expression_type_arg, int equation_arg, int dvariable1_arg) :
      Instruction {Tag::FNUMEXPR},
      expression_type {expression_type_arg},
      equation {equation_arg},
      dvariable1 {dvariable1_arg},
      lag1 {0}
  {
  }
  FNUMEXPR(const ExpressionType expression_type_arg, int equation_arg, int dvariable1_arg,
           int lag1_arg) :
      Instruction {Tag::FNUMEXPR},
      expression_type {expression_type_arg},
      equation {equation_arg},
      dvariable1 {dvariable1_arg},
      lag1 {lag1_arg}
  {
  }
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

class FBEGINBLOCK final : public Instruction
{
  template<typename B>
  friend Writer& operator<<(Writer& code_file, const B& instr);

private:
  int size {0};
  BlockSimulationType type;
  vector<int> variable;
  vector<int> equation;
  vector<int> exogenous;
  vector<int> det_exogenous;
  bool is_linear {false};
  vector<Block_contain_type> Block_Contain_;
  int u_count_int {0};
  int nb_col_jacob {0};
  int det_exo_size, exo_size;

public:
  /* Constructor when derivatives w.r.t. exogenous are present (only makes
     sense when there is no block-decomposition, since there is no provision for
     derivatives w.r.t. endogenous not belonging to the block) */
  FBEGINBLOCK(int size_arg, BlockSimulationType type_arg, int first_element, int block_size,
              const vector<int>& variable_arg, const vector<int>& equation_arg, bool is_linear_arg,
              int u_count_int_arg, int nb_col_jacob_arg, int det_exo_size_arg, int exo_size_arg,
              vector<int> det_exogenous_arg, vector<int> exogenous_arg) :
      Instruction {Tag::FBEGINBLOCK},
      size {size_arg},
      type {type_arg},
      variable {variable_arg.begin() + first_element,
                variable_arg.begin() + (first_element + block_size)},
      equation {equation_arg.begin() + first_element,
                equation_arg.begin() + (first_element + block_size)},
      exogenous {move(exogenous_arg)},
      det_exogenous {move(det_exogenous_arg)},
      is_linear {is_linear_arg},
      u_count_int {u_count_int_arg},
      nb_col_jacob {nb_col_jacob_arg},
      det_exo_size {det_exo_size_arg},
      exo_size {exo_size_arg}
  {
  }
  // Constructor when derivatives w.r.t. exogenous are absent
  FBEGINBLOCK(int size_arg, BlockSimulationType type_arg, int first_element, int block_size,
              const vector<int>& variable_arg, const vector<int>& equation_arg, bool is_linear_arg,
              int u_count_int_arg, int nb_col_jacob_arg) :
      Instruction {Tag::FBEGINBLOCK},
      size {size_arg},
      type {type_arg},
      variable {variable_arg.begin() + first_element,
                variable_arg.begin() + (first_element + block_size)},
      equation {equation_arg.begin() + first_element,
                equation_arg.begin() + (first_element + block_size)},
      is_linear {is_linear_arg},
      u_count_int {u_count_int_arg},
      nb_col_jacob {nb_col_jacob_arg},
      det_exo_size {0},
      exo_size {0}
  {
  }
  /* Deserializing constructor.
     Updates the code pointer to point beyond the bytes read. */
  FBEGINBLOCK(char*& code) : Instruction {Tag::FBEGINBLOCK}
  {
    code += sizeof(tag);

    auto read_member = [&code](auto& member) {
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
        Block_Contain_.push_back(bc);
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
class Writer : private ofstream
{
  template<typename B>
  friend Writer& operator<<(Writer& code_file, const B& instr);

private:
  // Stores the positions of all instructions in the byte stream
  vector<pos_type> instructions_positions;

public:
  Writer(const filesystem::path& filename);
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
  overwriteInstruction(int instruction_number, const B& new_instruction)
  {
    seekp(instructions_positions.at(instruction_number));
    *this << new_instruction;
    instructions_positions.pop_back();
    seekp(0, ios_base::end);
  }
};

// Overloads of operator<< for writing bytecode instructions

template<typename B>
Writer&
operator<<(Writer& code_file, const B& instr)
{
  code_file.instructions_positions.push_back(code_file.tellp());
  code_file.write(reinterpret_cast<const char*>(&instr), sizeof(B));
  return code_file;
}

template<>
Writer& operator<<(Writer& code_file, const FCALL& instr);

template<>
Writer& operator<<(Writer& code_file, const FBEGINBLOCK& instr);

}

#endif
