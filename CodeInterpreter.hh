/*
 * Copyright (C) 2007-2009 Dynare Team
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

#ifndef _CODEINTERPRETER_HH
#define _CODEINTERPRETER_HH
//#define DEBUGL
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <vector>
#ifdef LINBCG
# include "linbcg.hh"
#endif
#ifdef BYTE_CODE
#ifndef DEBUG_EX
  #include "mex.h"
#else
  #include "mex_interface.hh"
#endif
#endif

#ifdef _MSC_VER
typedef __int8 int8_t;
typedef unsigned __int8 uint8_t;
typedef __int16 int16_t;
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
# include <stdint.h>
#endif

using namespace std;

/**
* \enum Tags
* \brief The differents flags of the bytecode
*/
enum Tags
  {
    FLDZ,         //!< Stores zero in the stack - 0
    FLDC,         //!< Stores a constant term in the stack - 1

    FDIMT,        //!< Defines the number of temporary terms - dynamic context (the period has to be indicated) - 2
    FDIMST,       //!< Defines the number of temporary terms - static context (the period hasn't to be indicated) - 3
    FLDT,         //!< Stores a temporary term in the stack - dynamic context (the period has to be indicated) - 4
    FLDST,        //!< Stores a temporary term in the stack - static context (the period hasn't to be indicated) - 5
    FSTPT,        //!< Loads a temporary term from the stack - dynamic context (the period has to be indicated) - 6
    FSTPST,       //!< Loads a temporary term from the stack - static context (the period hasn't to be indicated) - 7

    FLDU,         //!< Stores an element of the vector U in the stack - dynamic context (the period has to be indicated) - 8
    FLDSU,        //!< Stores an element of the vector U in the stack - static context (the period hasn't to be indicated) - 9
    FSTPU,        //!< Loads an element of the vector U from the stack - dynamic context (the period has to be indicated) - A
    FSTPSU,       //!< Loads an element of the vector U from the stack - static context (the period hasn't to be indicated) - B

    FLDV,         //!< Stores a variable (described in SymbolType) in the stack - dynamic context (the period has to be indicated) - C
    FLDSV,        //!< Stores a variable (described in SymbolType) in the stack - static context (the period hasn't to be indicated) - D
    FLDVS,        //!< Stores a variable (described in SymbolType) in the stack - dynamic context but inside the STEADYSTATE function (the period hasn't to be indicated) - E
    FSTPV,        //!< Loads a variable (described in SymbolType) from the stack - dynamic context (the period has to be indicated) - F
    FSTPSV,       //!< Loads a variable (described in SymbolType) from the stack - static context (the period hasn't to be indicated) - 10

    FLDR,         //!< Stores a residual in the stack - 11
    FSTPR,        //!< Loads a residual from the stack - 12

    FSTPG,        //!< Loads a derivative from the stack - 13

    FUNARY,       //!< A Unary operator - 14
    FBINARY,      //!< A binary operator - 15

    FCUML,        //!< Cumulates the result - 16

    FBEGINBLOCK,  //!< Defines the begining of a model block - 17
    FENDBLOCK,    //!< Defines the end of a model block - 18
    FENDEQU,      //!< Defines the last equation of the block. For block that has to be solved, the derivatives appear just after this flag - 19
    FEND,         //!< Defines the end of the model code - 1A

    FOK           //!< Used for debugging purpose - 1B

};


enum BlockType
  {
    SIMULTANS = 0, //!< Simultaneous time separable block
    PROLOGUE = 1,  //!< Prologue block (one equation at the beginning, later merged)
    EPILOGUE = 2,  //!< Epilogue block (one equation at the beginning, later merged)
    SIMULTAN = 3   //!< Simultaneous time unseparable block
  };

enum EquationType
  {
    E_UNKNOWN,              //!< Unknown equation type
    E_EVALUATE,             //!< Simple evaluation, normalized variable on left-hand side
    E_EVALUATE_S,           //!< Simple evaluation, normalize using the first order derivative
    E_SOLVE                 //!< No simple evaluation of the equation, it has to be solved
  };



enum BlockSimulationType
  {
    UNKNOWN,                      //!< Unknown simulation type
    EVALUATE_FORWARD,             //!< Simple evaluation, normalized variable on left-hand side, forward
    EVALUATE_BACKWARD,            //!< Simple evaluation, normalized variable on left-hand side, backward
    SOLVE_FORWARD_SIMPLE,         //!< Block of one equation, newton solver needed, forward
    SOLVE_BACKWARD_SIMPLE,        //!< Block of one equation, newton solver needed, backward
    SOLVE_TWO_BOUNDARIES_SIMPLE,  //!< Block of one equation, newton solver needed, forward & ackward
    SOLVE_FORWARD_COMPLETE,       //!< Block of several equations, newton solver needed, forward
    SOLVE_BACKWARD_COMPLETE,      //!< Block of several equations, newton solver needed, backward
    SOLVE_TWO_BOUNDARIES_COMPLETE,//!< Block of several equations, newton solver needed, forward and backwar
  };

//! Enumeration of possible symbol types
/*! Warning: do not to change existing values for 0 to 4: the values matter for homotopy_setup command */
enum SymbolType
  {
    eEndogenous = 0,               //!< Endogenous
    eExogenous = 1,                //!< Exogenous
    eExogenousDet = 2,             //!< Exogenous deterministic
    eParameter = 4,                //!< Parameter
    eModelLocalVariable = 10,      //!< Local variable whose scope is model (pound expression)
    eModFileLocalVariable = 11,    //!< Local variable whose scope is mod file (model excluded)
    eUnknownFunction = 12          //!< Function unknown to the preprocessor
  };

enum UnaryOpcode
  {
    oUminus,
    oExp,
    oLog,
    oLog10,
    oCos,
    oSin,
    oTan,
    oAcos,
    oAsin,
    oAtan,
    oCosh,
    oSinh,
    oTanh,
    oAcosh,
    oAsinh,
    oAtanh,
    oSqrt,
    oSteadyState,
    oExpectation
  };

enum BinaryOpcode
  {
    oPlus,
    oMinus,
    oTimes,
    oDivide,
    oPower,
    oEqual,
    oMax,
    oMin,
    oLess,
    oGreater,
    oLessEqual,
    oGreaterEqual,
    oEqualEqual,
    oDifferent
  };

enum TrinaryOpcode
  {
    oNormcdf
  };

  struct Block_contain_type
  {
    int Equation, Variable, Own_Derivative;
  };

#pragma pack(push, 1)
class TagWithoutArgument
{
protected:
  uint8_t op_code;
public:
  inline TagWithoutArgument(uint8_t op_code_arg) : op_code(op_code_arg) {};
  inline void write(ostream &CompileCode) {CompileCode.write(reinterpret_cast<char *>(this), sizeof(*this)); };
};

template < class T1 >
class TagWithOneArgument
{
protected:
  uint8_t op_code;
  T1 arg1;
public:
  inline TagWithOneArgument(uint8_t op_code_arg) : op_code(op_code_arg) {};
  inline TagWithOneArgument(uint8_t op_code_arg, T1 arg_arg1) : op_code(op_code_arg), arg1(arg_arg1) {};
  inline void write(ostream &CompileCode) {CompileCode.write(reinterpret_cast<char *>(this), sizeof(TagWithOneArgument)); };
};

template < class T1, class T2 >
class TagWithTwoArguments
{
protected:
  uint8_t op_code;
  T1 arg1;
  T2 arg2;
public:
  inline TagWithTwoArguments(uint8_t op_code_arg) : op_code(op_code_arg) {};
  inline TagWithTwoArguments(uint8_t op_code_arg, T1 arg_arg1, T2 arg_arg2) : op_code(op_code_arg), arg1(arg_arg1), arg2(arg_arg2) {};
  inline void write(ostream &CompileCode) {CompileCode.write(reinterpret_cast<char *>(this), sizeof(*this)); };
};

template < class T1, class T2, class T3 >
class TagWithThreeArguments
{
protected:
  uint8_t op_code;
  T1 arg1;
  T2 arg2;
  T3 arg3;
public:
  inline TagWithThreeArguments(uint8_t op_code_arg) : op_code(op_code_arg) {};
  inline TagWithThreeArguments(uint8_t op_code_arg, T1 arg_arg1, T2 arg_arg2, T3 arg_arg3) : op_code(op_code_arg), arg1(arg_arg1), arg2(arg_arg2), arg3(arg_arg3) {};
  inline void write(ostream &CompileCode) {CompileCode.write(reinterpret_cast<char *>(this), sizeof(*this)); };
};



class FLDZ_ : public TagWithoutArgument
{
public:
  inline FLDZ_() : TagWithoutArgument(FLDZ) {};
};

class FEND_ : public TagWithoutArgument
{
public:
  inline FEND_() : TagWithoutArgument(FEND) {};
};

class FENDBLOCK_ : public TagWithoutArgument
{
public:
  inline FENDBLOCK_() : TagWithoutArgument(FENDBLOCK) {};
};

class FENDEQU_ : public TagWithoutArgument
{
public:
  inline FENDEQU_() : TagWithoutArgument(FENDEQU) {};
};

class FCUML_ : public TagWithoutArgument
{
public:
  inline FCUML_() : TagWithoutArgument(FCUML) {};
};

class FDIMT_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FDIMT_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMT) {};
  inline FDIMT_(unsigned int size_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMT,size_arg) {};
  inline unsigned int get_size() {return arg1;};
};

class FDIMST_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FDIMST_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMST) {};
  inline FDIMST_(const unsigned int size_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMST,size_arg) {};
  inline unsigned int get_size() {return arg1;};
};

class FLDC_ : public TagWithOneArgument<double>
{
public:
  inline FLDC_() : TagWithOneArgument<double>::TagWithOneArgument(FLDC) {};
  inline FLDC_(const double value_arg) : TagWithOneArgument<double>::TagWithOneArgument(FLDC,value_arg) {};
  inline double get_value() {return arg1;};
};

class FLDU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FLDU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDU) {};
  inline FLDU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDU,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FLDSU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FLDSU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDSU) {};
  inline FLDSU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDSU,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FLDR_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FLDR_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDR) {};
  inline FLDR_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDR,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FLDT_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FLDT_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDT) {};
  inline FLDT_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDT,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FLDST_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FLDST_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDST) {};
  inline FLDST_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDST,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FSTPT_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FSTPT_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPT) {};
  inline FSTPT_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPT,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FSTPST_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FSTPST_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPST) {};
  inline FSTPST_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPST,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FSTPR_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FSTPR_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPR) {};
  inline FSTPR_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPR,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FSTPU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FSTPU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPU) {};
  inline FSTPU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPU,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FSTPSU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FSTPSU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPSU) {};
  inline FSTPSU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPSU,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};


class FSTPG_ : public TagWithOneArgument<unsigned int>
{
public:
  inline FSTPG_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPG,0) {};
  inline FSTPG_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPG,pos_arg) {};
  inline unsigned int get_pos() {return arg1;};
};

class FUNARY_ : public TagWithOneArgument<uint8_t>
{
public:
  inline FUNARY_() : TagWithOneArgument<uint8_t>::TagWithOneArgument(FUNARY) {};
  inline FUNARY_(uint8_t op_type_arg) : TagWithOneArgument<uint8_t>::TagWithOneArgument(FUNARY,op_type_arg) {};
  inline uint8_t get_op_type() {return arg1;};
};

class FBINARY_ : public TagWithOneArgument<uint8_t>
{
public:
  inline FBINARY_() : TagWithOneArgument<uint8_t>::TagWithOneArgument(FBINARY) {};
  inline FBINARY_(const int op_type_arg) : TagWithOneArgument<uint8_t>::TagWithOneArgument(FBINARY,op_type_arg) {};
  inline uint8_t get_op_type() {return arg1;};
};

class FOK_ : public TagWithOneArgument<int>
{
public:
  inline FOK_() : TagWithOneArgument<int>::TagWithOneArgument(FOK) {};
  inline FOK_(const int arg_arg) : TagWithOneArgument<int>::TagWithOneArgument(FOK,arg_arg) {};
  inline int get_arg() {return arg1;};
};

class FLDVS_ : public TagWithTwoArguments<uint8_t, unsigned int>
{
public:
  inline FLDVS_() : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDVS) {};
  inline FLDVS_(uint8_t type_arg, const unsigned int pos_arg) : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDVS, type_arg, pos_arg) {};
  inline uint8_t get_type() {return arg1;};
  inline unsigned int get_pos() {return arg2;};
};

class FLDSV_ : public TagWithTwoArguments<uint8_t, unsigned int>
{
public:
  inline FLDSV_() : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDSV) {};
  inline FLDSV_(const uint8_t type_arg, const unsigned int pos_arg) :
            TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDSV, type_arg, pos_arg) {};
  inline uint8_t get_type() {return arg1;};
  inline unsigned int get_pos() {return arg2;};
};

class FSTPSV_ : public TagWithTwoArguments<uint8_t, unsigned int>
{
public:
  inline FSTPSV_() : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FSTPSV) {};
  inline FSTPSV_(const uint8_t type_arg, const unsigned int pos_arg) :
            TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FSTPSV, type_arg, pos_arg) {};
  inline uint8_t get_type() {return arg1;};
  inline unsigned int get_pos() {return arg2;};
};



class FLDV_ : public TagWithThreeArguments<uint8_t, unsigned int, int>
{
public:
  inline FLDV_() : TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FLDV) {};
  inline FLDV_(const int type_arg, const unsigned int pos_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FLDV, type_arg, pos_arg, 0) {};
  inline FLDV_(const int type_arg, const unsigned int pos_arg, const int lead_lag_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FLDV, type_arg, pos_arg, lead_lag_arg) {};
  inline uint8_t get_type() {return arg1;};
  inline unsigned int get_pos() {return arg2;};
  inline int get_lead_lag() {return arg3;};
};

class FSTPV_ : public TagWithThreeArguments<uint8_t, unsigned int, int>
{
public:
  inline FSTPV_() : TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FSTPV) {};
  inline FSTPV_(const int type_arg, const unsigned int pos_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FSTPV, type_arg, pos_arg, 0) {};
  inline FSTPV_(const int type_arg, const unsigned int pos_arg, const int lead_lag_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FSTPV, type_arg, pos_arg, lead_lag_arg) {};
  inline uint8_t get_type() {return arg1;};
  inline unsigned int get_pos() {return arg2;};
  inline int get_lead_lag() {return arg3;};
};


class FBEGINBLOCK_
{
private:
  uint8_t op_code;
  int size;
  uint8_t type;
  int *variable;
  int *equation;
  int *own_derivatives;
  bool is_linear;
  vector<Block_contain_type> Block_Contain_;
  int endo_nbr;
  int Max_Lag;
  int Max_Lead;
  int u_count_int;
public:
  inline FBEGINBLOCK_(){ op_code = FBEGINBLOCK; size = 0; type = UNKNOWN; variable = NULL; equation = NULL; own_derivatives = NULL;
           is_linear = false; endo_nbr = 0; Max_Lag = 0; Max_Lead = 0; u_count_int = 0;};
  inline FBEGINBLOCK_(const int size_arg, const BlockSimulationType type_arg, int *variable_arg, int *equation_arg, int *own_derivatives_arg,
                      bool is_linear_arg, int endo_nbr_arg, int Max_Lag_arg, int Max_Lead_arg, int u_count_int_arg)
         { op_code = FBEGINBLOCK; size = size_arg; type = type_arg; variable = variable_arg; equation = equation_arg; own_derivatives = own_derivatives_arg;
           is_linear = is_linear_arg; endo_nbr = endo_nbr_arg; Max_Lag = Max_Lag_arg; Max_Lead = Max_Lead_arg; u_count_int = u_count_int_arg;/*Block_Contain.clear();*/};
  inline unsigned int get_size() { return size;};
  inline uint8_t get_type() { return type;};
  inline bool get_is_linear() { return is_linear;};
  inline int get_endo_nbr() { return endo_nbr;};
  inline int get_Max_Lag() { return Max_Lag;};
  inline int get_Max_Lead() { return Max_Lead;};
  inline int get_u_count_int() { return u_count_int;};
  inline vector<Block_contain_type> get_Block_Contain() {return Block_Contain_;};
  inline void write(ostream &CompileCode)
   {
     CompileCode.write(reinterpret_cast<char *>(&op_code), sizeof(op_code));
     CompileCode.write(reinterpret_cast<char *>(&size), sizeof(size));
     CompileCode.write(reinterpret_cast<char *>(&type), sizeof(type));
     for(int i = 0; i<size ;i++)
       {
         CompileCode.write(reinterpret_cast<char *>(&variable[i]), sizeof(variable[0]));
         CompileCode.write(reinterpret_cast<char *>(&equation[i]), sizeof(equation[0]));
         CompileCode.write(reinterpret_cast<char *>(&own_derivatives[i]), sizeof(own_derivatives[0]));
       }
     if (type==SOLVE_TWO_BOUNDARIES_SIMPLE || type==SOLVE_TWO_BOUNDARIES_COMPLETE ||
         type==SOLVE_BACKWARD_COMPLETE || type==SOLVE_FORWARD_COMPLETE)
       {
         CompileCode.write(reinterpret_cast<char *>(&is_linear), sizeof(is_linear));
         CompileCode.write(reinterpret_cast<char *>(&endo_nbr), sizeof(endo_nbr));
         CompileCode.write(reinterpret_cast<char *>(&Max_Lag), sizeof(Max_Lag));
         CompileCode.write(reinterpret_cast<char *>(&Max_Lead), sizeof(Max_Lead));
         CompileCode.write(reinterpret_cast<char *>(&u_count_int), sizeof(u_count_int));
       }
   };
#ifdef BYTE_CODE

  inline uint8_t* load(uint8_t *code)
    {
      op_code = FBEGINBLOCK; code += sizeof(op_code);
      memcpy(&size, code, sizeof(size)); code += sizeof(size);
      memcpy(&type, code, sizeof(type)); code += sizeof(type);
      for(int i = 0; i<size ;i++)
        {
          Block_contain_type bc;
          memcpy(&bc.Variable, code, sizeof(bc.Variable)); code += sizeof(bc.Variable);
          memcpy(&bc.Equation, code, sizeof(bc.Equation)); code += sizeof(bc.Equation);
          memcpy(&bc.Own_Derivative, code, sizeof(bc.Own_Derivative)); code += sizeof(bc.Own_Derivative);
          Block_Contain_.push_back(bc);
        }
       if (type==SOLVE_TWO_BOUNDARIES_SIMPLE || type==SOLVE_TWO_BOUNDARIES_COMPLETE ||
           type==SOLVE_BACKWARD_COMPLETE || type==SOLVE_FORWARD_COMPLETE)
        {
          memcpy(&is_linear, code, sizeof(is_linear)); code += sizeof(is_linear);
          memcpy(&endo_nbr, code, sizeof(endo_nbr)); code += sizeof(endo_nbr);
          memcpy(&Max_Lag, code, sizeof(Max_Lag)); code += sizeof(Max_Lag);
          memcpy(&Max_Lead, code, sizeof(Max_Lead)); code += sizeof(Max_Lead);
          memcpy(&u_count_int, code, sizeof(u_count_int)); code += sizeof(u_count_int);
        }
      return code;
    };
#endif
};



#ifdef BYTE_CODE
typedef vector<pair<Tags, void* > > tags_liste_type;
class CodeLoad
{
private:
  uint8_t *code;
public:

  inline void* get_current_code() {return code;};
  inline tags_liste_type get_op_code(string file_name)
   {
     tags_liste_type tags_liste;
     ifstream CompiledCode;
     streamoff Code_Size;
     CompiledCode.open((file_name + ".cod").c_str(),std::ios::in | std::ios::binary| std::ios::ate);
     if (!CompiledCode.is_open())
       {
         return tags_liste;
       }
     Code_Size=CompiledCode.tellg();
     CompiledCode.seekg(std::ios::beg);
     code=(uint8_t*)mxMalloc(Code_Size);
     CompiledCode.seekg(0);
     CompiledCode.read(reinterpret_cast<char *>(code), Code_Size);
     CompiledCode.close();
     bool done = false;
     while(!done)
       {
         switch (*code)
           {
             case FLDZ:
#ifdef DEBUGL
               mexPrintf("FLDZ = %d size = %d\n",FLDZ, sizeof(FLDZ_));
#endif
               tags_liste.push_back(make_pair(FLDZ, code));
               code += sizeof(FLDZ_);
               break;
             case FEND:
#ifdef DEBUGL
               mexPrintf("FEND\n");
#endif
               tags_liste.push_back(make_pair(FEND, code));
               code += sizeof(FEND_);
               done = true;
               break;
             case FENDBLOCK:
#ifdef DEBUGL
               mexPrintf("FENDBLOCK\n");
#endif
               tags_liste.push_back(make_pair(FENDBLOCK, code));
               code += sizeof(FENDBLOCK_);
               break;
             case FENDEQU:
#ifdef DEBUGL
               mexPrintf("FENDEQU\n");
#endif
               tags_liste.push_back(make_pair(FENDEQU, code));
               code += sizeof(FENDEQU_);
               break;
             case FCUML:
#ifdef DEBUGL
               mexPrintf("FCUML\n");
#endif
               tags_liste.push_back(make_pair(FCUML, code));
               code += sizeof(FCUML_);
               break;
             case FDIMT:
#ifdef DEBUGL
               mexPrintf("FDIMT = %d size = %d\n",FDIMT, sizeof(FDIMT_));
#endif
               tags_liste.push_back(make_pair(FDIMT, code));
               code += sizeof(FDIMT_);
               break;
             case FDIMST:
#ifdef DEBUGL
               mexPrintf("FDIMST\n");
#endif
               tags_liste.push_back(make_pair(FDIMST, code));
               code += sizeof(FDIMST_);
               break;
             case FLDC:
#ifdef DEBUGL
               mexPrintf("FLDC\n");
#endif
               tags_liste.push_back(make_pair(FLDC, code));
               code += sizeof(FLDC_);
               break;
             case FLDU:
#ifdef DEBUGL
               mexPrintf("FLDU\n");
#endif
               tags_liste.push_back(make_pair(FLDU, code));
               code += sizeof(FLDU_);
               break;
            case FLDSU:
#ifdef DEBUGL
               mexPrintf("FLDSU\n");
#endif
               tags_liste.push_back(make_pair(FLDSU, code));
               code += sizeof(FLDSU_);
               break;
            case FLDR:
#ifdef DEBUGL
               mexPrintf("FLDR\n");
#endif
               tags_liste.push_back(make_pair(FLDR, code));
               code += sizeof(FLDR_);
               break;
            case FLDT:
#ifdef DEBUGL
               mexPrintf("FLDT\n");
#endif
               tags_liste.push_back(make_pair(FLDT, code));
               code += sizeof(FLDT_);
               break;
            case FLDST:
#ifdef DEBUGL
               mexPrintf("FLDST\n");
#endif
               tags_liste.push_back(make_pair(FLDST, code));
               code += sizeof(FLDST_);
               break;
            case FSTPT:
#ifdef DEBUGL
               mexPrintf("FSTPT = %d size = %d\n",FSTPT, sizeof(FSTPT_));
#endif
               tags_liste.push_back(make_pair(FSTPT, code));
               code += sizeof(FSTPT_);
               break;
            case FSTPST:
#ifdef DEBUGL
               mexPrintf("FSTPST\n");
#endif
               tags_liste.push_back(make_pair(FSTPST, code));
               code += sizeof(FSTPST_);
               break;
            case FSTPR:
#ifdef DEBUGL
              mexPrintf("FSTPR\n");
#endif
               tags_liste.push_back(make_pair(FSTPR, code));
               code += sizeof(FSTPR_);
               break;
            case FSTPU:
#ifdef DEBUGL
               mexPrintf("FSTPU\n");
#endif
               tags_liste.push_back(make_pair(FSTPU, code));
               code += sizeof(FSTPU_);
               break;
            case FSTPSU:
#ifdef DEBUGL
               mexPrintf("FSTPSU\n");
#endif
               tags_liste.push_back(make_pair(FSTPSU, code));
               code += sizeof(FSTPSU_);
               break;
            case FSTPG:
#ifdef DEBUGL
               mexPrintf("FSTPG\n");
#endif
               tags_liste.push_back(make_pair(FSTPG, code));
               code += sizeof(FSTPG_);
               break;
            case FUNARY:
#ifdef DEBUGL
               mexPrintf("FUNARY\n");
#endif
               tags_liste.push_back(make_pair(FUNARY, code));
               code += sizeof(FUNARY_);
               break;
            case FBINARY:
#ifdef DEBUGL
               mexPrintf("FBINARY\n");
#endif
               tags_liste.push_back(make_pair(FBINARY, code));
               code += sizeof(FBINARY_);
               break;
            case FOK:
#ifdef DEBUGL
               mexPrintf("FOK\n");
#endif
               tags_liste.push_back(make_pair(FOK, code));
               code += sizeof(FOK_);
               break;
            case FLDVS:
#ifdef DEBUGL
               mexPrintf("FLDVS\n");
#endif
               tags_liste.push_back(make_pair(FLDVS, code));
               code += sizeof(FLDVS_);
               break;
            case FLDSV:
#ifdef DEBUGL
               mexPrintf("FLDSV\n");
#endif
               tags_liste.push_back(make_pair(FLDSV, code));
               code += sizeof(FLDSV_);
               break;
            case FSTPSV:
#ifdef DEBUGL
               mexPrintf("FSTPSV\n");
#endif
               tags_liste.push_back(make_pair(FSTPSV, code));
               code += sizeof(FSTPSV_);
               break;
            case FLDV:
#ifdef DEBUGL
               mexPrintf("FLDV\n");
#endif
               tags_liste.push_back(make_pair(FLDV, code));
               code += sizeof(FLDV_);
               break;
            case FSTPV:
#ifdef DEBUGL
               mexPrintf("FSTPV\n");
#endif
               tags_liste.push_back(make_pair(FSTPV, code));
               code += sizeof(FSTPV_);
               break;
            case FBEGINBLOCK:
#ifdef DEBUGL
               mexPrintf("FBEGINBLOCK\n");
#endif
               {
                 FBEGINBLOCK_ *fbegin_block = new FBEGINBLOCK_;

                 code = fbegin_block->load(code);

                 tags_liste.push_back(make_pair(FBEGINBLOCK, fbegin_block));
               }
               break;
            default:
              mexPrintf("Unknown Tag value=%d code=%x\n",*code, code);
              done = true;
           }
       }
     return tags_liste;
    };
};
#endif
#pragma pack(pop)
#endif

