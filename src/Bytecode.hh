/*
 * Copyright © 2007-2024 Dynare Team
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

#include <concepts>
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

  FSTPG,  // Stores the derivative of a simple (single-equation) block in simulate mode
  FSTPG2, // Stores the derivative matrix of a block in evaluate mode

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

template<typename T>
concept IsInstruction = derived_from<T, Instruction>;

struct FLDZ final : public Instruction
{
  FLDZ() : Instruction {Tag::FLDZ}
  {
  }
};

struct FEND final : public Instruction
{
  FEND() : Instruction {Tag::FEND}
  {
  }
};

struct FENDBLOCK final : public Instruction
{
  FENDBLOCK() : Instruction {Tag::FENDBLOCK}
  {
  }
};

struct FENDEQU final : public Instruction
{
  FENDEQU() : Instruction {Tag::FENDEQU}
  {
  }
};

struct FDIMT final : public Instruction
{
  const int size;
  explicit FDIMT(int size_arg) : Instruction {Tag::FDIMT}, size {size_arg}
  {
  }
};

struct FDIMST final : public Instruction
{
  const int size;
  explicit FDIMST(int size_arg) : Instruction {Tag::FDIMST}, size {size_arg}
  {
  }
};

struct FLDC final : public Instruction
{
  const double value;
  explicit FLDC(double value_arg) : Instruction {Tag::FLDC}, value {value_arg}
  {
  }
};

struct FLDU final : public Instruction
{
  const int pos;
  explicit FLDU(int pos_arg) : Instruction {Tag::FLDU}, pos {pos_arg}
  {
  }
};

struct FLDSU final : public Instruction
{
  const int pos;
  explicit FLDSU(int pos_arg) : Instruction {Tag::FLDSU}, pos {pos_arg}
  {
  }
};

struct FLDR final : public Instruction
{
  const int pos;
  explicit FLDR(int pos_arg) : Instruction {Tag::FLDR}, pos {pos_arg}
  {
  }
};

struct FLDT final : public Instruction
{
  const int pos;
  explicit FLDT(int pos_arg) : Instruction {Tag::FLDT}, pos {pos_arg}
  {
  }
};

struct FLDST final : public Instruction
{
  const int pos;
  explicit FLDST(int pos_arg) : Instruction {Tag::FLDST}, pos {pos_arg}
  {
  }
};

struct FSTPT final : public Instruction
{
  const int pos;
  explicit FSTPT(int pos_arg) : Instruction {Tag::FSTPT}, pos {pos_arg}
  {
  }
};

struct FSTPST final : public Instruction
{
  const int pos;
  explicit FSTPST(int pos_arg) : Instruction {Tag::FSTPST}, pos {pos_arg}
  {
  }
};

struct FSTPR final : public Instruction
{
  const int pos;
  explicit FSTPR(int pos_arg) : Instruction {Tag::FSTPR}, pos {pos_arg}
  {
  }
};

struct FSTPU final : public Instruction
{
  const int pos;
  explicit FSTPU(int pos_arg) : Instruction {Tag::FSTPU}, pos {pos_arg}
  {
  }
};

struct FSTPSU final : public Instruction
{
  const int pos;
  explicit FSTPSU(int pos_arg) : Instruction {Tag::FSTPSU}, pos {pos_arg}
  {
  }
};

struct FSTPG final : public Instruction
{
  explicit FSTPG() : Instruction {Tag::FSTPG}
  {
  }
};

struct FSTPG2 final : public Instruction
{
  const int row, col;
  FSTPG2(int row_arg, int col_arg) : Instruction {Tag::FSTPG2}, row {row_arg}, col {col_arg}
  {
  }
};

struct FUNARY final : public Instruction
{
  const UnaryOpcode op_code;
  explicit FUNARY(UnaryOpcode op_code_arg) : Instruction {Tag::FUNARY}, op_code {op_code_arg}
  {
  }
};

struct FBINARY final : public Instruction
{
  const BinaryOpcode op_code;
  explicit FBINARY(BinaryOpcode op_code_arg) : Instruction {Tag::FBINARY}, op_code {op_code_arg}
  {
  }
};

struct FTRINARY final : public Instruction
{
  const TrinaryOpcode op_code;
  explicit FTRINARY(TrinaryOpcode op_code_arg) : Instruction {Tag::FTRINARY}, op_code {op_code_arg}
  {
  }
};

struct FJMPIFEVAL final : public Instruction
{
  const int pos;
  explicit FJMPIFEVAL(int pos_arg) : Instruction {Tag::FJMPIFEVAL}, pos {pos_arg}
  {
  }
};

struct FJMP final : public Instruction
{
  const int pos;
  explicit FJMP(int pos_arg) : Instruction {Tag::FJMP}, pos {pos_arg}
  {
  }
};

struct FLDTEF final : public Instruction
{
  const int number;
  explicit FLDTEF(int number_arg) : Instruction {Tag::FLDTEF}, number {number_arg}
  {
  }
};

struct FSTPTEF final : public Instruction
{
  const int number;
  explicit FSTPTEF(int number_arg) : Instruction {Tag::FSTPTEF}, number {number_arg}
  {
  }
};

struct FLDTEFD final : public Instruction
{
  const int indx, row;
  FLDTEFD(int indx_arg, int row_arg) : Instruction {Tag::FLDTEFD}, indx {indx_arg}, row {row_arg}
  {
  }
};

struct FSTPTEFD final : public Instruction
{
  const int indx, row;
  FSTPTEFD(int indx_arg, int row_arg) : Instruction {Tag::FSTPTEFD}, indx {indx_arg}, row {row_arg}
  {
  }
};

struct FLDTEFDD final : public Instruction
{
  const int indx, row, col;
  FLDTEFDD(int indx_arg, int row_arg, int col_arg) :
      Instruction {Tag::FLDTEFDD}, indx {indx_arg}, row {row_arg}, col {col_arg}
  {
  }
};

struct FSTPTEFDD final : public Instruction
{
  const int indx, row, col;
  FSTPTEFDD(int indx_arg, int row_arg, int col_arg) :
      Instruction {Tag::FSTPTEF}, indx {indx_arg}, row {row_arg}, col {col_arg}
  {
  }
};

struct FLDVS final : public Instruction
{
  const SymbolType type;
  const int pos;
  FLDVS(SymbolType type_arg, int pos_arg) : Instruction {Tag::FLDVS}, type {type_arg}, pos {pos_arg}
  {
  }
};

struct FLDSV final : public Instruction
{
  const SymbolType type;
  const int pos;
  FLDSV(SymbolType type_arg, int pos_arg) : Instruction {Tag::FLDSV}, type {type_arg}, pos {pos_arg}
  {
  }
};

struct FSTPSV final : public Instruction
{
  const SymbolType type;
  const int pos;
  FSTPSV(SymbolType type_arg, int pos_arg) :
      Instruction {Tag::FSTPSV}, type {type_arg}, pos {pos_arg}
  {
  }
};

struct FLDV final : public Instruction
{
  const SymbolType type;
  const int pos, lead_lag;
  FLDV(SymbolType type_arg, int pos_arg, int lead_lag_arg) :
      Instruction {Tag::FLDV}, type {type_arg}, pos {pos_arg}, lead_lag {lead_lag_arg}
  {
  }
};

struct FSTPV final : public Instruction
{
  const SymbolType type;
  const int pos, lead_lag;
  FSTPV(SymbolType type_arg, int pos_arg, int lead_lag_arg) :
      Instruction {Tag::FSTPV}, type {type_arg}, pos {pos_arg}, lead_lag {lead_lag_arg}
  {
  }
};

class FCALL final : public Instruction
{
  template<IsInstruction B>
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
  }
  int
  get_nb_output_arguments()
  {
    return nb_output_arguments;
  }
  int
  get_nb_input_arguments()
  {
    return nb_input_arguments;
  }
  int
  get_indx()
  {
    return indx;
  }
  void
  set_arg_func_name(string arg_arg_func_name)
  {
    arg_func_name = move(arg_arg_func_name);
  }
  string
  get_arg_func_name()
  {
    return arg_func_name;
  }
  void
  set_nb_add_input_arguments(int arg_add_input_arguments)
  {
    add_input_arguments = arg_add_input_arguments;
  }
  int
  get_nb_add_input_arguments()
  {
    return add_input_arguments;
  }
  void
  set_row(int arg_row)
  {
    row = arg_row;
  }
  int
  get_row()
  {
    return row;
  }
  void
  set_col(int arg_col)
  {
    col = arg_col;
  }
  int
  get_col()
  {
    return col;
  }
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
  }
  int
  get_dvariable1()
  {
    return dvariable1;
  }
  int
  get_lag1()
  {
    return lag1;
  }
};

class FBEGINBLOCK final : public Instruction
{
  template<IsInstruction B>
  friend Writer& operator<<(Writer& code_file, const B& instr);

private:
  int size {0};
  BlockSimulationType type;
  vector<int> variables;
  vector<int> equations;
  vector<int> exogenous;
  vector<int> det_exogenous;
  bool is_linear {false};
  int u_count_int {0};
  int nb_col_jacob {0};
  int det_exo_size, exo_size;

public:
  /* Constructor when derivatives w.r.t. exogenous are present (only makes
     sense when there is no block-decomposition, since there is no provision for
     derivatives w.r.t. endogenous not belonging to the block) */
  FBEGINBLOCK(int size_arg, BlockSimulationType type_arg, int first_element, int block_size,
              const vector<int>& variables_arg, const vector<int>& equations_arg,
              bool is_linear_arg, int u_count_int_arg, int nb_col_jacob_arg, int det_exo_size_arg,
              int exo_size_arg, vector<int> det_exogenous_arg, vector<int> exogenous_arg) :
      Instruction {Tag::FBEGINBLOCK},
      size {size_arg},
      type {type_arg},
      variables {variables_arg.begin() + first_element,
                 variables_arg.begin() + (first_element + block_size)},
      equations {equations_arg.begin() + first_element,
                 equations_arg.begin() + (first_element + block_size)},
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
              const vector<int>& variables_arg, const vector<int>& equations_arg,
              bool is_linear_arg, int u_count_int_arg, int nb_col_jacob_arg) :
      Instruction {Tag::FBEGINBLOCK},
      size {size_arg},
      type {type_arg},
      variables {variables_arg.begin() + first_element,
                 variables_arg.begin() + (first_element + block_size)},
      equations {equations_arg.begin() + first_element,
                 equations_arg.begin() + (first_element + block_size)},
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
    variables.resize(size);
    equations.resize(size);
    for (int i {0}; i < size; i++)
      {
        read_member(variables[i]);
        read_member(equations[i]);
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
  }
  BlockSimulationType
  get_type()
  {
    return type;
  }
  bool
  get_is_linear()
  {
    return is_linear;
  }
  int
  get_u_count_int()
  {
    return u_count_int;
  }
  vector<int>
  get_variables()
  {
    return variables;
  }
  vector<int>
  get_equations()
  {
    return equations;
  }
  int
  get_nb_col_jacob()
  {
    return nb_col_jacob;
  }
  int
  get_exo_size()
  {
    return exo_size;
  }
  int
  get_det_exo_size()
  {
    return det_exo_size;
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
  template<IsInstruction B>
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
  void
  overwriteInstruction(int instruction_number, const IsInstruction auto& new_instruction)
  {
    seekp(instructions_positions.at(instruction_number));
    *this << new_instruction;
    instructions_positions.pop_back();
    seekp(0, ios_base::end);
  }
};

// Overloads of operator<< for writing bytecode instructions

template<IsInstruction B>
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
