/*
 * Copyright (C) 2019 Dynare Team
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

#ifndef _DIRECTIVES_HH
#define _DIRECTIVES_HH

#include "Expressions.hh"

using namespace std;

namespace macro
{
  class Directive : public Node
  {
    // A Parent class just for clarity
  public:
    Directive(Environment &env_arg, const Tokenizer::location location_arg) : Node(env_arg, move(location_arg)) { }
    // Directives can be interpreted
    virtual void interpret(ostream &output, bool no_line_macro) = 0;
  };

  using DirectivePtr = shared_ptr<Directive>;


  class TextNode : public Directive
  {
    // Class for text not interpreted by macroprocessor
    // Not a real directive node
    // Treated as such as the output is only to be interpreted
  private:
    const string text;
  public:
    TextNode(const string text_arg, Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), text{move(text_arg)} { }
    inline void interpret(ostream &output, bool no_line_macro) override { output << text; }
  };


  class Eval : public Directive
  {
    // Class for @{} statements
    // Not a real directive node
    // Treated as such as the output is only to be interpreted
  private:
    const ExpressionPtr expr;
  public:
    Eval(const ExpressionPtr expr_arg, Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, bool no_line_macro) override;
  };


  class Include : public Directive
  {
  private:
    const ExpressionPtr expr;
    string name;
  public:
    Include(const ExpressionPtr expr_arg, Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, bool no_line_macro) override;
    inline string getName() const { return name; }
  };


  class IncludePath : public Directive
  {
  private:
    const ExpressionPtr expr;
    string path;
  public:
    IncludePath(const ExpressionPtr expr_arg, Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, bool no_line_macro) override;
    inline string getPath() const { return path; }
  };


  class Define : public Directive
  {
  private:
    const VariablePtr var;
    const FunctionPtr func;
    const ExpressionPtr value;
  public:
    Define(const VariablePtr var_arg,
           const ExpressionPtr value_arg,
           Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), var{move(var_arg)}, value{move(value_arg)} { }
    Define(const FunctionPtr func_arg,
           const ExpressionPtr value_arg,
           Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), func{move(func_arg)}, value{move(value_arg)} { }
    void interpret(ostream &output, bool no_line_macro) override;
  };


  class Echo : public Directive
  {
  private:
    const ExpressionPtr expr;
  public:
    Echo(const ExpressionPtr expr_arg,
         Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, bool no_line_macro) override;
  };


  class Error : public Directive
  {
  private:
    const ExpressionPtr expr;
  public:
    Error(const ExpressionPtr expr_arg,
          Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), expr{move(expr_arg)} { }
    void interpret(ostream &output, bool no_line_macro) override;
  };


  class EchoMacroVars : public Directive
  {
  private:
    bool save;
  public:
    EchoMacroVars(bool save_arg,
                  Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), save{save_arg} { }
    void interpret(ostream &output, bool no_line_macro) override;
  };


  class For : public Directive
  {
  private:
    const vector<VariablePtr> index_vec;
    const ExpressionPtr index_vals;
    vector<DirectivePtr> statements;
  public:
    For(const vector<VariablePtr> index_vec_arg,
        const ExpressionPtr index_vals_arg,
        const vector<DirectivePtr> statements_arg,
        Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), index_vec{move(index_vec_arg)},
      index_vals{move(index_vals_arg)}, statements{statements_arg} { }
    void interpret(ostream &output, bool no_line_macro) override;
  };


  class If : public Directive
  {
  protected:
    const ExpressionPtr condition;
    vector<DirectivePtr> if_statements;
    vector<DirectivePtr> else_statements;
  public:
    If(const ExpressionPtr condition_arg,
       const vector<DirectivePtr> if_statements_arg,
       Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)), condition{move(condition_arg)}, if_statements{if_statements_arg} { }
    If(const ExpressionPtr condition_arg,
       const vector<DirectivePtr> if_statements_arg,
       const vector<DirectivePtr> else_statements_arg,
       Environment &env_arg, const Tokenizer::location location_arg) :
      Directive(env_arg, move(location_arg)),
      condition{move(condition_arg)}, if_statements{if_statements_arg}, else_statements{else_statements_arg} { }
    void interpret(ostream &output, bool no_line_macro) override;
  protected:
    void loopIf(ostream &output, bool no_line_macro);
    void loopElse(ostream &output, bool no_line_macro);
  };


  class Ifdef : public If
  {
  public:
    Ifdef(const ExpressionPtr condition_arg,
          const vector<DirectivePtr> if_statements_arg,
          Environment &env_arg, const Tokenizer::location location_arg) :
      If(move(condition_arg), move(if_statements_arg), env_arg, move(location_arg)) { }
    Ifdef(const ExpressionPtr condition_arg,
          const vector<DirectivePtr> if_statements_arg,
          const vector<DirectivePtr> else_statements_arg,
          Environment &env_arg, const Tokenizer::location location_arg) :
      If(move(condition_arg), move(if_statements_arg), move(else_statements_arg), env_arg, move(location_arg)) { }
    void interpret(ostream &output, bool no_line_macro) override;
  };


  class Ifndef : public If
  {
  public:
    Ifndef(const ExpressionPtr condition_arg,
           const vector<DirectivePtr> if_statements_arg,
           Environment &env_arg, const Tokenizer::location location_arg) :
      If(move(condition_arg), move(if_statements_arg), env_arg, move(location_arg)) { }
    Ifndef(const ExpressionPtr condition_arg,
           const vector<DirectivePtr> if_statements_arg,
           const vector<DirectivePtr> else_statements_arg,
           Environment &env_arg, const Tokenizer::location location_arg) :
      If(move(condition_arg), move(if_statements_arg), move(else_statements_arg), env_arg, move(location_arg)) { }
    void interpret(ostream &output, bool no_line_macro) override;
  };
}
#endif
