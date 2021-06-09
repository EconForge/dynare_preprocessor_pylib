/*
 * Copyright (C) 2019-2020 Dynare Team
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

#ifndef _DIRECTIVES_HH
#define _DIRECTIVES_HH

#include "Expressions.hh"

#include <filesystem>

namespace macro
{
  class Directive : public Node
  {
    // A Parent class just for clarity
  public:
    explicit Directive(Tokenizer::location location_arg) :
      Node(move(location_arg)) { }
    // Directives can be interpreted
    virtual void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) = 0;
  };


  class TextNode : public Directive
  {
    // Class for text not interpreted by macroprocessor
    // Not a real directive node
    // Treated as such as the output is only to be interpreted
  private:
    const string text;
  public:
    TextNode(string text_arg, Tokenizer::location location_arg) :
      Directive(move(location_arg)), text{move(text_arg)} { }
    inline void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override { output << text; }
  };


  class Eval : public Directive
  {
    // Class for @{} statements
    // Not a real directive node
    // Treated as such as the output is only to be interpreted
  private:
    const ExpressionPtr expr;
  public:
    Eval(ExpressionPtr expr_arg, Tokenizer::location location_arg) :
      Directive(move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class Include : public Directive
  {
  private:
    const ExpressionPtr expr;
  public:
    Include(ExpressionPtr expr_arg, Tokenizer::location location_arg) :
      Directive(move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class IncludePath : public Directive
  {
  private:
    const ExpressionPtr expr;
  public:
    IncludePath(ExpressionPtr expr_arg, Tokenizer::location location_arg) :
      Directive(move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class Define : public Directive
  {
  private:
    const VariablePtr var;
    const FunctionPtr func;
    const ExpressionPtr value;
  public:
    Define(VariablePtr var_arg,
           ExpressionPtr value_arg,
           Tokenizer::location location_arg) :
      Directive(move(location_arg)), var{move(var_arg)}, value{move(value_arg)} { }
    Define(FunctionPtr func_arg,
           ExpressionPtr value_arg,
           Tokenizer::location location_arg) :
      Directive(move(location_arg)), func{move(func_arg)}, value{move(value_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class Echo : public Directive
  {
  private:
    const ExpressionPtr expr;
  public:
    Echo(ExpressionPtr expr_arg,
         Tokenizer::location location_arg) :
      Directive(move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class Error : public Directive
  {
  private:
    const ExpressionPtr expr;
  public:
    Error(ExpressionPtr expr_arg,
          Tokenizer::location location_arg) :
      Directive(move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class EchoMacroVars : public Directive
  {
  private:
    const bool save;
    const vector<string> vars;
  public:
    EchoMacroVars(bool save_arg,
                  Tokenizer::location location_arg) :
      Directive(move(location_arg)), save{save_arg} { }
    EchoMacroVars(bool save_arg, vector<string> vars_arg,
                  Tokenizer::location location_arg) :
      Directive(move(location_arg)), save{save_arg}, vars{move(vars_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class For : public Directive
  {
  private:
    const vector<VariablePtr> index_vec;
    const ExpressionPtr index_vals;
    const vector<DirectivePtr> statements;
  public:
    For(vector<VariablePtr> index_vec_arg,
        ExpressionPtr index_vals_arg,
        vector<DirectivePtr> statements_arg,
        Tokenizer::location location_arg) :
      Directive(move(location_arg)), index_vec{move(index_vec_arg)},
      index_vals{move(index_vals_arg)}, statements{move(statements_arg)} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  };


  class If : public Directive
  {
  protected:
    /* Every if statement and the associated body to execute are stored in a
     * pair<ExpressionPtr, vector<DirectivePtr>>, where the ExpressionPtr is the condition
     * and vector<DirectivePtr> is the series of statements to execute if the condition evaluates
     * to true.
     * The `if` statement is the first element in the vector
     * If there exist any `elseif` statements, they follow
     * If there is an `else` statement it is the last element in the vector. Its condition is true.
     */
    const vector<pair<ExpressionPtr, vector<DirectivePtr>>> expr_and_body;
    const bool ifdef, ifndef;
  public:
    If(vector<pair<ExpressionPtr, vector<DirectivePtr>>> expr_and_body_arg,
       Tokenizer::location location_arg,
       bool ifdef_arg = false, bool ifndef_arg = false) :
      Directive(move(location_arg)), expr_and_body{move(expr_and_body_arg)},
      ifdef{ifdef_arg}, ifndef{ifndef_arg} { }
    void interpret(ostream &output, Environment &env, vector<filesystem::path> &paths) override;
  protected:
    void interpretBody(const vector<DirectivePtr> &body, ostream &output,
                       Environment &env, vector<filesystem::path> &paths);
  };

  class Ifdef : public If
  {
  public:
    Ifdef(vector<pair<ExpressionPtr, vector<DirectivePtr>>> expr_and_body_arg,
          Tokenizer::location location_arg) :
      If(move(expr_and_body_arg), move(location_arg), true, false) { }
  };


  class Ifndef : public If
  {
  public:
    Ifndef(vector<pair<ExpressionPtr, vector<DirectivePtr>>> expr_and_body_arg,
           Tokenizer::location location_arg) :
      If(move(expr_and_body_arg), move(location_arg), false, true) { }
  };
}
#endif
