/*
 * Copyright (C) 2003-2008 Dynare Team
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

%skeleton "lalr1.cc"
%require "2.3"
%defines

%{
using namespace std;

class ParsingDriver;

#include "ExprNode.hh"

/* Little hack: we redefine the macro which computes the locations, because
   we need to access the location from within the parsing driver for error
   and warning messages. */
#define YYLLOC_DEFAULT(Current, Rhs, N)                 \
  do {                                                  \
    if (N)                                              \
      {                                                 \
        (Current).begin = (Rhs)[1].begin;               \
        (Current).end   = (Rhs)[N].end;                 \
      }                                                 \
    else                                                \
      {                                                 \
        (Current).begin = (Current).end = (Rhs)[0].end;	\
      }                                                 \
    driver.location = (Current);                        \
  } while(false)

%}

%name-prefix="Dynare"

%parse-param { ParsingDriver &driver }
%lex-param { ParsingDriver &driver }

%locations
%initial-action
{
  // Initialize the location filenames to NULL
  @$.begin.filename = @$.end.filename = NULL;
}

%debug
%error-verbose

%union
{
  string *string_val;
  NodeID node_val;
};

%{
#include "ParsingDriver.hh"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex
%}

%token AR AUTOCORR
%token BAYESIAN_IRF BETA_PDF BICGSTAB
%token BVAR_DENSITY BVAR_FORECAST
%token BVAR_PRIOR_DECAY BVAR_PRIOR_FLAT BVAR_PRIOR_LAMBDA
%token BVAR_PRIOR_MU BVAR_PRIOR_OMEGA BVAR_PRIOR_TAU BVAR_PRIOR_TRAIN
%token BVAR_REPLIC
%token CALIB CALIB_VAR CHECK CONF_SIG CONSTANT CORR COVAR CUTOFF
%token DATAFILE DROP DSAMPLE DYNASAVE DYNATYPE
%token END ENDVAL EQUAL ESTIMATION ESTIMATED_PARAMS ESTIMATED_PARAMS_BOUNDS ESTIMATED_PARAMS_INIT
%token FILENAME FILTER_STEP_AHEAD FILTERED_VARS FIRST_OBS
%token <string_val> FLOAT_NUMBER
%token FORECAST
%token GAMMA_PDF GAUSSIAN_ELIMINATION GCC_COMPILER GMRES GRAPH
%token HISTVAL HP_FILTER HP_NGRID
%token INITVAL INITVAL_FILE
%token <string_val> INT_NUMBER
%token INV_GAMMA1_PDF INV_GAMMA2_PDF IRF
%token KALMAN_ALGO KALMAN_TOL
%token LAPLACE LCC_COMPILER LIK_ALGO LIK_INIT LINEAR LOAD_MH_FILE LOGLINEAR LU
%token MARKOWITZ MARGINAL_DENSITY MAX
%token METHOD MH_DROP MH_INIT_SCALE MH_JSCALE MH_MODE MH_NBLOCKS MH_REPLIC MH_RECOVER MIN
%token MODE_CHECK MODE_COMPUTE MODE_FILE MODEL MODEL_COMPARISON MODEL_INFO MSHOCKS
%token MODIFIEDHARMONICMEAN MOMENTS_VARENDO DIFFUSE_FILTER
%token <string_val> NAME
%token NO_COMPILER NOBS NOCONSTANT NOCORR NODIAGNOSTIC NOFUNCTIONS
%token NOGRAPH NOMOMENTS NOPRINT NORMAL_PDF
%token OBSERVATION_TRENDS OPTIM OPTIM_WEIGHTS ORDER OSR OSR_PARAMS
%token PARAMETERS PERIODS PLANNER_OBJECTIVE PREFILTER PRESAMPLE
%token PRINT PRIOR_TRUNC PRIOR_ANALYSIS POSTERIOR_ANALYSIS
%token <string_val> QUOTED_STRING
%token QZ_CRITERIUM
%token RELATIVE_IRF REPLIC RPLOT
%token SHOCKS SIGMA_E SIMUL SIMUL_ALGO SIMUL_SEED SMOOTHER SOLVE_ALGO
%token SPARSE SPARSE_DLL STDERR STEADY STOCH_SIMUL
%token TEX RAMSEY_POLICY PLANNER_DISCOUNT
%token <string_val> TEX_NAME
%token UNIFORM_PDF UNIT_ROOT_VARS USE_DLL
%token VALUES VAR VAREXO VAREXO_DET VAROBS
%token XLS_SHEET XLS_RANGE
%token NORMCDF
%token HOMOTOPY_SETUP HOMOTOPY_MODE HOMOTOPY_STEPS
%left COMMA
%left EQUAL_EQUAL EXCLAMATION_EQUAL
%left LESS GREATER LESS_EQUAL GREATER_EQUAL
%left PLUS MINUS
%left TIMES DIVIDE
%left UMINUS UPLUS
%nonassoc POWER
%token EXP LOG LN LOG10 SIN COS TAN ASIN ACOS ATAN SINH COSH TANH ASINH ACOSH ATANH SQRT
/* GSA analysis */
%token DYNARE_SENSITIVITY IDENTIFICATION MORRIS STAB REDFORM PPRIOR PRIOR_RANGE PPOST ILPTAU GLUE MORRIS_NLIV
%token MORRIS_NTRA NSAM LOAD_REDFORM LOAD_RMSE LOAD_STAB ALPHA2_STAB KSSTAT LOGTRANS_REDFORM THRESHOLD_REDFORM
%token KSSTAT_REDFORM ALPHA2_REDFORM NAMENDO NAMLAGENDO NAMEXO RMSE LIK_ONLY VAR_RMSE PFILT_RMSE ISTART_RMSE
%token ALPHA_RMSE ALPHA2_RMSE TRANS_IDENT
/* end of GSA analysis*/

%type <node_val> expression
%type <node_val> equation hand_side model_var
%type <string_val> signed_float signed_integer prior
%type <string_val> filename
%type <string_val> value value1 vec_int_elem vec_int_1 vec_int
%type <string_val> vec_value_1 vec_value
%type <string_val> calib_arg2 range number

%%

%start statement_list;

statement_list : statement
               | statement_list statement
               ;

statement : declaration
          | periods
          | cutoff
          | markowitz
          | model
          | initval
          | initval_file
          | endval
          | histval
          | init_param
          | shocks
          | mshocks
          | sigma_e
          | steady
          | check
          | simul
          | stoch_simul
          | estimation
          | prior_analysis
          | posterior_analysis
          | estimated_params
          | estimated_params_bounds
          | estimated_params_init
          | varobs
          | observation_trends
          | unit_root_vars
          | dsample
          | rplot
          | optim_weights
          | osr_params
          | osr
          | calib_var
          | calib
          | dynatype
          | dynasave
          | model_comparison
          | model_info
          | planner_objective
          | ramsey_policy
          | bvar_density
          | bvar_forecast
          | dynare_sensitivity
          | homotopy_setup
          | forecast
          ;

declaration : parameters
            | var
            | varexo
            | varexo_det
            ;


dsample : DSAMPLE INT_NUMBER ';'
          { driver.dsample($2); }
        | DSAMPLE INT_NUMBER INT_NUMBER ';'
          { driver.dsample($2, $3); }
        ;

rplot : RPLOT symbol_list ';' { driver.rplot(); };

var : VAR var_list ';';

varexo : VAREXO varexo_list ';';

varexo_det : VAREXO_DET varexo_det_list ';';

parameters : PARAMETERS parameter_list ';';

var_list : var_list NAME
           { driver.declare_endogenous($2); }
         | var_list COMMA NAME
           { driver.declare_endogenous($3); }
         | NAME
           { driver.declare_endogenous($1); }
         | var_list NAME TEX_NAME
           { driver.declare_endogenous($2, $3); }
         | var_list COMMA NAME TEX_NAME
           { driver.declare_endogenous($3, $4); }
         | NAME TEX_NAME
           { driver.declare_endogenous($1, $2); }
         ;

varexo_list : varexo_list NAME
              { driver.declare_exogenous($2); }
            | varexo_list COMMA NAME
              { driver.declare_exogenous($3); }
            | NAME
              { driver.declare_exogenous($1); }
            | varexo_list NAME TEX_NAME
              { driver.declare_exogenous($2, $3); }
            | varexo_list COMMA NAME TEX_NAME
              { driver.declare_exogenous($3, $4); }
            | NAME TEX_NAME
              { driver.declare_exogenous($1, $2); }
            ;

varexo_det_list : varexo_det_list NAME
                  { driver.declare_exogenous_det($2); }
                | varexo_det_list COMMA NAME
                  { driver.declare_exogenous_det($3); }
                | NAME
                  { driver.declare_exogenous_det($1); }
                | varexo_det_list NAME TEX_NAME
                  { driver.declare_exogenous_det($2, $3); }
                | varexo_det_list COMMA NAME TEX_NAME
                  { driver.declare_exogenous_det($3, $4); }
                | NAME TEX_NAME
                   { driver.declare_exogenous_det($1, $2); }
                ;

parameter_list : parameter_list NAME
                 { driver.declare_parameter($2); }
               | parameter_list COMMA NAME
                 { driver.declare_parameter($3); }
               | NAME
                 { driver.declare_parameter($1); }
               | parameter_list NAME TEX_NAME
                 { driver.declare_parameter($2, $3); }
               | parameter_list COMMA NAME TEX_NAME
                 { driver.declare_parameter($3, $4); }
               | NAME TEX_NAME
                 { driver.declare_parameter($1, $2); }
               ;

periods : PERIODS INT_NUMBER ';'
          { driver.periods($2); }
        | PERIODS EQUAL INT_NUMBER ';'
          { driver.periods($3); }
        ;

cutoff : CUTOFF FLOAT_NUMBER ';'
         { driver.cutoff($2); }
       | CUTOFF EQUAL FLOAT_NUMBER ';'
         { driver.cutoff($3); }
       ;

markowitz : MARKOWITZ FLOAT_NUMBER ';'
            { driver.markowitz($2); }
          | MARKOWITZ EQUAL FLOAT_NUMBER ';'
            { driver.markowitz($3); }
          ;


init_param : NAME EQUAL expression ';' { driver.init_param($1, $3); };

expression : '(' expression ')'
             { $$ = $2;}
           | NAME
             { $$ = driver.add_expression_variable($1); }
           | FLOAT_NUMBER
             { $$ = driver.add_constant($1); }
           | INT_NUMBER
             { $$ = driver.add_constant($1); }
           | expression PLUS expression
             { $$ = driver.add_plus($1, $3); }
           | expression MINUS expression
             { $$ = driver.add_minus($1, $3); }
           | expression DIVIDE expression
             { $$ = driver.add_divide($1, $3); }
           | expression TIMES expression
             { $$ = driver.add_times($1, $3); }
           | expression POWER expression
             { $$ = driver.add_power($1, $3); }
           | expression LESS expression
             { $$ = driver.add_less($1, $3); }
           | expression GREATER expression
             { $$ = driver.add_greater($1, $3); }
           | expression LESS_EQUAL expression
             { $$ = driver.add_less_equal($1, $3); }
           | expression GREATER_EQUAL expression
             { $$ = driver.add_greater_equal($1, $3); }
           | expression EQUAL_EQUAL expression
             { $$ = driver.add_equal_equal($1, $3); }
           | expression EXCLAMATION_EQUAL expression
             { $$ = driver.add_different($1, $3); }
           | MINUS expression %prec UMINUS
             { $$ = driver.add_uminus($2); }
           | PLUS expression %prec UPLUS
             { $$ = $2; }
           | EXP '(' expression ')'
             { $$ = driver.add_exp($3); }
           | LOG '(' expression ')'
             { $$ = driver.add_log($3); }
           | LN '(' expression ')'
             { $$ = driver.add_log($3); }
           | LOG10 '(' expression ')'
             { $$ = driver.add_log10($3); }
           | SIN '(' expression ')'
             { $$ = driver.add_sin($3); }
           | COS '(' expression ')'
             { $$ = driver.add_cos($3); }
           | TAN '(' expression ')'
             { $$ = driver.add_tan($3); }
           | ASIN '(' expression ')'
             { $$ = driver.add_asin($3); }
           | ACOS '(' expression ')'
             { $$ = driver.add_acos($3); }
           | ATAN '(' expression ')'
             { $$ = driver.add_atan($3); }
           | SQRT '(' expression ')'
             { $$ = driver.add_sqrt($3); }
           | MAX '(' expression COMMA expression ')'
             { $$ = driver.add_max($3, $5); }
           | MIN '(' expression COMMA expression ')'
             { $$ = driver.add_min($3, $5); }
           | NAME '(' comma_expression ')'
             { $$ = driver.add_unknown_function($1); }
           | NORMCDF '(' expression COMMA expression COMMA expression ')'
             { $$ = driver.add_normcdf($3, $5, $7); }
           | NORMCDF '(' expression ')'
             { $$ = driver.add_normcdf($3); }
           ;

comma_expression : expression
                   { driver.add_unknown_function_arg($1); }
                   | comma_expression COMMA expression
                   { driver.add_unknown_function_arg($3); }
                 ;

initval : INITVAL ';' initval_list END
          { driver.end_initval(); }

initval_file : INITVAL_FILE '(' FILENAME EQUAL NAME ')' ';'
               { driver.initval_file($5); }
             ;

endval : ENDVAL ';' initval_list END { driver.end_endval(); };

initval_list : initval_list initval_elem
             | initval_elem
             ;

initval_elem : NAME EQUAL expression ';' { driver.init_val($1, $3); };

histval : HISTVAL ';' histval_list END { driver.end_histval(); };

histval_list : histval_list histval_elem
             | histval_elem
             ;

histval_elem : NAME '(' signed_integer ')' EQUAL expression ';' { driver.hist_val($1, $3, $6); };

model_sparse_dll_options_list : model_sparse_dll_options_list COMMA model_sparse_dll_options
                          | model_sparse_dll_options
                          ;

model_sparse_options_list : model_sparse_options_list COMMA model_sparse_common_options
                          | model_sparse_common_options
                          ;

model_sparse_dll_options : model_compiler_options
                         | model_sparse_common_options
                         ;

model_compiler_options : LCC_COMPILER
                       { driver.init_compiler(0); }
                     | GCC_COMPILER
                       { driver.init_compiler(1); }
                     | NO_COMPILER
                       { driver.init_compiler(2); }
                      ;

model_sparse_common_options : o_cutoff
                     | o_markowitz
                     ;

model : MODEL ';' { driver.begin_model(); }
        equation_list END { driver.reset_data_tree(); }
      | MODEL '(' o_linear ')' ';' { driver.begin_model(); }
        equation_list END { driver.reset_data_tree(); }
      | MODEL '(' USE_DLL ')' ';' { driver.begin_model(); driver.use_dll(); }
        equation_list END { driver.reset_data_tree(); }
      | MODEL '(' SPARSE_DLL COMMA model_sparse_dll_options_list ')'
        { driver.begin_model(); driver.sparse_dll(); } ';'
        equation_list END { driver.reset_data_tree(); }
      | MODEL '(' SPARSE_DLL ')' { driver.begin_model(); driver.sparse_dll(); } ';'
        equation_list END { driver.reset_data_tree(); }
      | MODEL '(' SPARSE COMMA model_sparse_options_list ')'
        { driver.begin_model(); driver.sparse(); } ';'
        equation_list END { driver.reset_data_tree(); }
      | MODEL '(' SPARSE ')' { driver.begin_model(); driver.sparse(); } ';'
        equation_list END { driver.reset_data_tree(); }
      ;

equation_list : equation_list equation
              | equation_list pound_expression
              | equation
              | pound_expression
              ;

equation : hand_side EQUAL hand_side ';'
           { $$ = driver.add_model_equal($1, $3); }
         | hand_side ';'
           { $$ = driver.add_model_equal_with_zero_rhs($1); }
         ;

hand_side : '(' hand_side ')'
            { $$ = $2;}
          | model_var
          | FLOAT_NUMBER
            { $$ = driver.add_constant($1); }
          | INT_NUMBER
            { $$ = driver.add_constant($1); }
          | hand_side PLUS hand_side
            { $$ = driver.add_plus($1, $3); }
          | hand_side MINUS hand_side
            { $$ = driver.add_minus($1, $3); }
          | hand_side DIVIDE hand_side
            { $$ = driver.add_divide($1, $3); }
          | hand_side TIMES hand_side
            { $$ = driver.add_times($1, $3); }
          | hand_side LESS hand_side
            { $$ = driver.add_less($1, $3); }
          | hand_side GREATER hand_side
            { $$ = driver.add_greater($1, $3); }
          | hand_side LESS_EQUAL hand_side
            { $$ = driver.add_less_equal($1, $3); }
          | hand_side GREATER_EQUAL hand_side
            { $$ = driver.add_greater_equal($1, $3); }
          | hand_side EQUAL_EQUAL hand_side
            { $$ = driver.add_equal_equal($1, $3); }
          | hand_side EXCLAMATION_EQUAL hand_side
            { $$ = driver.add_different($1, $3); }
          | hand_side POWER hand_side
            { $$ = driver.add_power($1, $3); }
          | MINUS hand_side %prec UMINUS
            { $$ = driver.add_uminus($2); }
          | PLUS hand_side
            { $$ = $2; }
          | EXP '(' hand_side ')'
            { $$ = driver.add_exp($3); }
          | LOG '(' hand_side ')'
            { $$ = driver.add_log($3); }
          | LN '(' hand_side ')'
            { $$ = driver.add_log($3); }
          | LOG10 '(' hand_side ')'
            { $$ = driver.add_log10($3); }
          | SIN '(' hand_side ')'
            { $$ = driver.add_sin($3); }
          | COS '(' hand_side ')'
            { $$ = driver.add_cos($3); }
          | TAN '(' hand_side ')'
            { $$ = driver.add_tan($3); }
          | ASIN '(' hand_side ')'
            { $$ = driver.add_asin($3); }
          | ACOS '(' hand_side ')'
            { $$ = driver.add_acos($3); }
          | ATAN '(' hand_side ')'
            { $$ = driver.add_atan($3); }
          | SQRT '(' hand_side ')'
            { $$ = driver.add_sqrt($3); }
          | MAX '(' hand_side COMMA hand_side ')'
             { $$ = driver.add_max($3, $5); }
          | MIN '(' hand_side COMMA hand_side ')'
             { $$ = driver.add_min($3, $5); }
          | NORMCDF '(' hand_side COMMA hand_side COMMA hand_side ')'
             { $$ = driver.add_normcdf($3, $5, $7); }
          | NORMCDF '(' hand_side ')'
             { $$ = driver.add_normcdf($3); }
          ;

pound_expression: '#' NAME EQUAL hand_side ';'
                  { driver.declare_and_init_model_local_variable($2, $4); };

model_var : NAME
            { $$ = driver.add_model_variable($1); }
          | NAME '(' signed_integer ')'
            { $$ = driver.add_model_variable($1, $3); }
          ;

shocks : SHOCKS ';' shock_list END { driver.end_shocks(); };

mshocks : MSHOCKS ';' shock_list END { driver.end_mshocks(); };

shock_list : shock_list shock_elem
           | shock_elem
           ;

shock_elem : VAR NAME ';' PERIODS period_list ';' VALUES value_list ';'
             { driver.add_det_shock($2); }
           | VAR NAME ';' STDERR expression ';'
             { driver.add_stderr_shock($2, $5); }
           | VAR NAME EQUAL expression ';'
             { driver.add_var_shock($2, $4); }
           | VAR NAME COMMA NAME EQUAL expression ';'
             { driver.add_covar_shock($2, $4, $6); }
           | CORR NAME COMMA NAME EQUAL expression ';'
             { driver.add_correl_shock($2, $4, $6); }
           ;

period_list : period_list COMMA INT_NUMBER
              { driver.add_period($3); }
            | period_list INT_NUMBER
              { driver.add_period($2); }
            | period_list COMMA INT_NUMBER ':' INT_NUMBER
              { driver.add_period($3, $5); }
            | period_list INT_NUMBER ':' INT_NUMBER
              { driver.add_period($2, $4); }
            | INT_NUMBER ':' INT_NUMBER
              { driver.add_period($1, $3); }
            | INT_NUMBER
              { driver.add_period($1); }
            ;


sigma_e : SIGMA_E EQUAL '[' triangular_matrix ']' ';' { driver.do_sigma_e(); };

value_list
 	:  value_list COMMA expression
    {driver.add_value($3);}
  |  value_list number
    {driver.add_value($2);}
	| expression
    {driver.add_value($1);}
	;

triangular_matrix : triangular_matrix ';' triangular_row
                    { driver.end_of_row(); }
                  | triangular_row
                    { driver.end_of_row(); }
                  ;

triangular_row : triangular_row COMMA '(' expression ')'
                 { driver.add_to_row($4); }
               | triangular_row COMMA FLOAT_NUMBER
                 { driver.add_to_row_const($3); }
               | triangular_row COMMA INT_NUMBER
                 { driver.add_to_row_const($3); }
               | triangular_row '(' expression ')'
                 { driver.add_to_row($3); }
               | triangular_row FLOAT_NUMBER
                 { driver.add_to_row_const($2); }
               | triangular_row INT_NUMBER
                 { driver.add_to_row_const($2); }
               | '(' expression ')'
                 { driver.add_to_row($2); }
               | FLOAT_NUMBER
                 { driver.add_to_row_const($1); }
               | INT_NUMBER
                 { driver.add_to_row_const($1); }
               ;

steady : STEADY ';'
         { driver.steady(); }
       | STEADY '(' steady_options_list ')' ';'
         { driver.steady(); }
       ;

steady_options_list : steady_options_list COMMA steady_options
                    | steady_options
                    ;

steady_options : o_solve_algo
               | o_homotopy_mode
               | o_homotopy_steps
               ;

check : CHECK ';'
        { driver.check(); }
      | CHECK '(' check_options_list ')' ';'
        { driver.check(); }
      ;

check_options_list : check_options_list COMMA check_options
                   | check_options
                   ;

check_options : o_solve_algo;

model_info : MODEL_INFO ';'
             { driver.model_info(); }
           ;

simul : SIMUL ';'
        { driver.simulate(); }
      | SIMUL '(' simul_options_list ')' ';'
        { driver.simulate(); }
      ;

simul_options_list : simul_options_list COMMA simul_options
                   | simul_options
                   ;

simul_options : o_periods
              | o_datafile
              | o_method
              ;

stoch_simul : STOCH_SIMUL ';'
              { driver.stoch_simul(); }
            | STOCH_SIMUL '(' stoch_simul_options_list ')' ';'
              { driver.stoch_simul(); }
            | STOCH_SIMUL symbol_list ';'
              { driver.stoch_simul(); }
            | STOCH_SIMUL '(' stoch_simul_options_list ')' symbol_list ';'
              { driver.stoch_simul(); }
            ;

stoch_simul_options_list : stoch_simul_options_list COMMA stoch_simul_options
                         | stoch_simul_options
                         ;

stoch_simul_options : o_solve_algo
                    | o_simul_algo
                    | o_linear
                    | o_order
                    | o_replic
                    | o_drop
                    | o_ar
                    | o_nocorr
                    | o_nofunctions
                    | o_nomoments
                    | o_nograph
                    | o_irf
                    | o_relative_irf
                    | o_hp_filter
                    | o_hp_ngrid
                    | o_periods
                    | o_cutoff
                    | o_simul
                    | o_simul_seed
                    | o_qz_criterium
                    | o_print
                    | o_noprint
                    ;

symbol_list : symbol_list NAME
               { driver.add_in_symbol_list($2); }
             | symbol_list COMMA NAME
               { driver.add_in_symbol_list($3); }
             | NAME
               { driver.add_in_symbol_list($1); }
             ;

symbol_list_ext : symbol_list
                | ':'
                  {
                    string *colon = new string(":");
                    driver.add_in_symbol_list(colon);
                  }
                ;

signed_integer : PLUS INT_NUMBER
                 { $$ = $2; }
               | MINUS INT_NUMBER
                { $2->insert(0, "-"); $$ = $2; }
               | INT_NUMBER
                { $$ = $1; }
               ;

signed_float : PLUS FLOAT_NUMBER
               { $$ = $2; }
             | MINUS FLOAT_NUMBER
               { $2->insert(0, "-"); $$ = $2; }
             | FLOAT_NUMBER
               { $$ = $1; }
             | signed_integer
             ;

estimated_params : ESTIMATED_PARAMS ';' estimated_list END { driver.estimated_params(); };

estimated_list : estimated_list estimated_elem
                 { driver.add_estimated_params_element(); }
               | estimated_elem
                 { driver.add_estimated_params_element(); }
               ;

estimated_elem : estimated_elem1 COMMA estimated_elem2 ';';

estimated_elem1 : STDERR NAME
                  {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *$2;
                    delete $2;
                  }
                | NAME
                  {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *$1;
                    delete $1;
                  }
                | CORR NAME COMMA NAME
                  {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *$2;
                    driver.estim_params.name2 = *$4;
                    delete $2;
                    delete $4;
                  }
                ;

estimated_elem2 : prior COMMA estimated_elem3
                  {
                    driver.estim_params.prior = *$1;
                    delete $1;
                  }
                | value COMMA prior COMMA estimated_elem3
                  {
                    driver.estim_params.init_val = *$1;
                    driver.estim_params.prior = *$3;
                    delete $1;
                    delete $3;
                  }
                | value COMMA value COMMA value COMMA prior COMMA estimated_elem3
                  {
                    driver.estim_params.init_val = *$1;
                    driver.estim_params.low_bound = *$3;
                    driver.estim_params.up_bound = *$5;
                    driver.estim_params.prior = *$7;
                    delete $1;
                    delete $3;
                    delete $5;
                    delete $7;
                  }
                | value
                  {
                    driver.estim_params.init_val = *$1;
                    delete $1;
                  }
                | value COMMA value COMMA value
                  {
                    driver.estim_params.init_val = *$1;
                    driver.estim_params.low_bound = *$3;
                    driver.estim_params.up_bound = *$5;
                    delete $1;
                    delete $3;
                    delete $5;
                  }
                ;

estimated_elem3 : value COMMA value
                  {
                    driver.estim_params.mean = *$1;
                    driver.estim_params.std = *$3;
                    delete $1;
                    delete $3;
                  }
                | value COMMA value COMMA value
                  {
                    driver.estim_params.mean = *$1;
                    driver.estim_params.std = *$3;
                    driver.estim_params.p3 = *$5;
                    delete $1;
                    delete $3;
                    delete $5;
                  }
                | value COMMA value COMMA value COMMA value
                  {
                    driver.estim_params.mean = *$1;
                    driver.estim_params.std = *$3;
                    driver.estim_params.p3 = *$5;
                    driver.estim_params.p4 = *$7;
                    delete $1;
                    delete $3;
                    delete $5;
                    delete $7;
                  }
                | value COMMA value COMMA value COMMA value COMMA value
                  {
                    driver.estim_params.mean = *$1;
                    driver.estim_params.std = *$3;
                    driver.estim_params.p3 = *$5;
                    driver.estim_params.p4 = *$7;
                    driver.estim_params.jscale = *$9;
                    delete $1;
                    delete $3;
                    delete $5;
                    delete $7;
                    delete $9;
                  }
                ;

estimated_params_init : ESTIMATED_PARAMS_INIT ';' estimated_init_list END
                        { driver.estimated_params_init(); };

estimated_init_list : estimated_init_list estimated_init_elem
                      { driver.add_estimated_params_element(); }
                    | estimated_init_elem
                      { driver.add_estimated_params_element(); }
                    ;

estimated_init_elem : STDERR NAME COMMA value ';'
                      {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *$2;
                        driver.estim_params.init_val = *$4;
                        delete $2;
                        delete $4;
                      }
                    | CORR NAME COMMA NAME COMMA value ';'
                      {
                        driver.estim_params.type = 3;
                        driver.estim_params.name = *$2;
                        driver.estim_params.name2 = *$4;
                        driver.estim_params.init_val = *$6;
                        delete $2;
                        delete $4;
                        delete $6;
                      }
                    | NAME COMMA value ';'
                      {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *$1;
                        driver.estim_params.init_val = *$3;
                        delete $1;
                        delete $3;
                      }
                    ;

estimated_params_bounds : ESTIMATED_PARAMS_BOUNDS ';' estimated_bounds_list END
                          { driver.estimated_params_bounds(); };

estimated_bounds_list : estimated_bounds_list estimated_bounds_elem
                        { driver.add_estimated_params_element(); }
                      | estimated_bounds_elem
                        { driver.add_estimated_params_element(); }
                      ;

estimated_bounds_elem : STDERR NAME COMMA value COMMA value ';'
                        {
                          driver.estim_params.type = 1;
                          driver.estim_params.name = *$2;
                          driver.estim_params.low_bound = *$4;
                          driver.estim_params.up_bound = *$6;
                          delete $2;
                          delete $4;
                          delete $6;
                        }
                      | CORR NAME COMMA NAME COMMA value COMMA value ';'
                        {
                          driver.estim_params.type = 3;
                          driver.estim_params.name = *$2;
                          driver.estim_params.name2 = *$4;
                          driver.estim_params.low_bound = *$6;
                          driver.estim_params.up_bound = *$8;
                          delete $2;
                          delete $4;
                          delete $6;
                          delete $8;
                        }
                      | NAME COMMA value COMMA value ';'
                        {
                          driver.estim_params.type = 2;
                          driver.estim_params.name = *$1;
                          driver.estim_params.low_bound = *$3;
                          driver.estim_params.up_bound = *$5;
                          delete $1;
                          delete $3;
                          delete $5;
                        }
                      ;

prior : BETA_PDF
        { $$ = new string("1"); }
      | GAMMA_PDF
        { $$ = new string("2"); }
      | NORMAL_PDF
        { $$ = new string("3"); }
      | INV_GAMMA1_PDF
        { $$ = new string("4"); }
      | UNIFORM_PDF
        { $$ = new string("5"); }
      | INV_GAMMA2_PDF
        { $$ = new string("6"); }
      ;

value : { $$ = new string("NaN"); }
      | value1
      ;

value1 : INT_NUMBER
       | FLOAT_NUMBER
       | NAME
       | MINUS INT_NUMBER
         { $2->insert(0, "-"); $$ = $2; }
       | MINUS FLOAT_NUMBER
         { $2->insert(0, "-"); $$ = $2; }
       ;

estimation : ESTIMATION ';'
             { driver.run_estimation(); }
           | ESTIMATION '(' estimation_options_list ')' ';'
             { driver.run_estimation(); }
           | ESTIMATION symbol_list ';'
             { driver.run_estimation(); }
           | ESTIMATION '(' estimation_options_list ')' symbol_list ';'
             { driver.run_estimation(); }
           ;

estimation_options_list : estimation_options_list COMMA estimation_options
                        | estimation_options
                        ;

estimation_options : o_datafile
                   | o_nobs
                   | o_first_obs
                   | o_prefilter
                   | o_presample
                   | o_lik_algo
                   | o_lik_init
                   | o_nograph
                   | o_conf_sig
                   | o_mh_replic
                   | o_mh_drop
                   | o_mh_jscale
                   | o_optim
                   | o_mh_init_scale
                   | o_mode_file
                   | o_mode_compute
                   | o_mode_check
                   | o_prior_trunc
                   | o_mh_mode
                   | o_mh_nblcks
                   | o_load_mh_file
                   | o_loglinear
                   | o_nodiagnostic
                   | o_bayesian_irf
                   | o_irf
                   | o_tex
                   | o_forecast
                   | o_smoother
                   | o_moments_varendo
                   | o_filtered_vars
                   | o_kalman_algo
                   | o_kalman_tol
                   | o_xls_sheet
                   | o_xls_range
                   | o_filter_step_ahead
                   | o_solve_algo
                   | o_constant
                   | o_noconstant
                   | o_mh_recover
                   | o_diffuse_filter
                   | o_order
                   ;

prior_analysis : PRIOR_ANALYSIS '(' prior_posterior_options_list ')' ';'
                 { driver.run_prior_analysis(); }
               | PRIOR_ANALYSIS '(' prior_posterior_options_list ')' symbol_list ';'
                 { driver.run_prior_analysis(); }
               ;

prior_posterior_options_list : prior_posterior_options_list COMMA prior_posterior_options
                             | prior_posterior_options
                             ;

prior_posterior_options : o_nograph
                        | o_conf_sig
                        | o_prior_trunc
                        | o_bayesian_irf
                        | o_irf
                        | o_tex
                        | o_forecast
                        | o_smoother
                        | o_moments_varendo
                        | o_filtered_vars
                        | o_xls_sheet
                        | o_xls_range
                        | o_filter_step_ahead
                        ;

posterior_analysis : POSTERIOR_ANALYSIS '(' prior_posterior_options_list ')' ';'
                     { driver.run_posterior_analysis(); }
                   | POSTERIOR_ANALYSIS '(' prior_posterior_options_list ')' symbol_list ';'
                     { driver.run_posterior_analysis(); }
                   ;

list_optim_option : '\'' NAME '\'' COMMA '\'' NAME '\''
                    { driver.optim_options_string($2, $6); }
                  | '\'' NAME '\'' COMMA value
                    { driver.optim_options_num($2, $5); }
                  ;

optim_options : list_optim_option
              | optim_options COMMA list_optim_option;
              ;

varobs : VAROBS symbol_list ';' { driver.set_varobs(); };

observation_trends : OBSERVATION_TRENDS ';' trend_list END { driver.set_trends(); };

trend_list : trend_list trend_element
           | trend_element
           ;

trend_element :  NAME '(' expression ')' ';' { driver.set_trend_element($1, $3); };

unit_root_vars : UNIT_ROOT_VARS symbol_list ';' { driver.set_unit_root_vars(); };

optim_weights : OPTIM_WEIGHTS ';' optim_weights_list END { driver.optim_weights(); };

optim_weights_list : optim_weights_list NAME expression ';'
                     { driver.set_optim_weights($2, $3); }
                   | optim_weights_list NAME COMMA NAME expression ';'
                     { driver.set_optim_weights($2, $4, $5); }
                   | NAME expression ';'
                     { driver.set_optim_weights($1, $2); }
                   | NAME COMMA NAME expression ';'
                     { driver.set_optim_weights($1, $3, $4); }
                   ;

osr_params : OSR_PARAMS symbol_list ';' { driver.set_osr_params(); };

osr : OSR ';'
      { driver.run_osr(); }
    | OSR '(' stoch_simul_options ')' ';'
      { driver.run_osr(); }
    | OSR symbol_list ';'
      { driver.run_osr(); }
    | OSR '(' stoch_simul_options ')' symbol_list ';'
      {driver.run_osr(); }
    ;

calib_var : CALIB_VAR ';' calib_var_list END { driver.run_calib_var(); };

calib_var_list : calib_var_list calib_arg1
               | calib_arg1
               ;

calib_arg1 : NAME calib_arg2 EQUAL expression ';'
             { driver.set_calib_var($1, $2, $4); }
           | NAME COMMA NAME calib_arg2 EQUAL expression ';'
             { driver.set_calib_covar($1, $3, $4, $6); }
           | AUTOCORR NAME '(' INT_NUMBER ')' calib_arg2 EQUAL expression ';'
             { driver.set_calib_ac($2, $4, $6, $8); }
           ;

calib_arg2 : { $$ = new string("1"); }
           | '(' INT_NUMBER ')'
             { $$ = $2; }
           | '(' FLOAT_NUMBER ')'
             { $$ = $2; }
           ;

calib : CALIB ';'
        { driver.run_calib(0); }
      | CALIB '(' COVAR ')' ';'
        { driver.run_calib(1); }
      ;

dynatype : DYNATYPE '(' NAME ')'';'
           { driver.run_dynatype($3); }
         | DYNATYPE '(' NAME ')' symbol_list ';'
           { driver.run_dynatype($3); }
         | DYNATYPE NAME ';'
           { driver.run_dynatype($2); }
         | DYNATYPE '(' NAME '.' NAME ')'';'
           { driver.run_dynatype($3, $5); }
         | DYNATYPE '(' NAME '.' NAME ')' symbol_list ';'
           { driver.run_dynatype($3, $5); }
         | DYNATYPE NAME '.' NAME ';'
           { driver.run_dynatype($2,$4); }
         ;

dynasave : DYNASAVE '(' NAME ')'';'
           { driver.run_dynasave($3); }
         | DYNASAVE '(' NAME ')' symbol_list ';'
           { driver.run_dynasave($3); }
         | DYNASAVE NAME ';'
           { driver.run_dynasave($2); }
         | DYNASAVE '(' NAME '.' NAME ')'';'
           { driver.run_dynasave($3, $5); }
         | DYNASAVE '(' NAME '.' NAME ')' symbol_list ';'
           { driver.run_dynasave($3, $5); }
         | DYNASAVE NAME '.' NAME ';'
           { driver.run_dynasave($2, $4); }
         ;

model_comparison : MODEL_COMPARISON mc_filename_list ';'
                   { driver.run_model_comparison(); }
                 | MODEL_COMPARISON '(' o_marginal_density ')' mc_filename_list ';'
                   { driver.run_model_comparison(); }
                 ;

filename : NAME
           { $$ = $1; }
         | QUOTED_STRING
           { $$ = $1; }
         ;

mc_filename_list : filename
                   { driver.add_mc_filename($1); }
                 | filename '(' value ')'
                   { driver.add_mc_filename($1, $3); }
                 | mc_filename_list filename
                   { driver.add_mc_filename($2); }
                 | mc_filename_list filename '(' value ')'
                   { driver.add_mc_filename($2, $4); }
                 | mc_filename_list COMMA filename
                   { driver.add_mc_filename($3); }
                 | mc_filename_list COMMA filename '(' value ')'
                   { driver.add_mc_filename($3, $5); }
                 ;

planner_objective : PLANNER_OBJECTIVE { driver.begin_planner_objective(); }
                    hand_side { driver.end_planner_objective($3); } ';';

ramsey_policy : RAMSEY_POLICY ';'
                { driver.ramsey_policy(); }
              | RAMSEY_POLICY '(' ramsey_policy_options_list ')' ';'
                { driver.ramsey_policy(); }
              | RAMSEY_POLICY symbol_list ';'
                { driver.ramsey_policy(); }
              | RAMSEY_POLICY '(' ramsey_policy_options_list ')' symbol_list ';'
                { driver.ramsey_policy(); }
              ;

ramsey_policy_options_list : ramsey_policy_options_list COMMA ramsey_policy_options
                           | ramsey_policy_options
                           ;

ramsey_policy_options : stoch_simul_options
                      | o_planner_discount
                      ;

bvar_prior_option : o_bvar_prior_tau
                  | o_bvar_prior_decay
                  | o_bvar_prior_lambda
                  | o_bvar_prior_mu
                  | o_bvar_prior_omega
                  | o_bvar_prior_flat
                  | o_bvar_prior_train
                  ;

bvar_common_option : bvar_prior_option
                   | o_datafile
                   | o_xls_sheet
                   | o_xls_range
                   | o_first_obs
                   | o_presample
                   | o_nobs
                   | o_prefilter
                   | o_constant
                   | o_noconstant
                   ;

bvar_density_options_list : bvar_common_option COMMA bvar_density_options_list
                          | bvar_common_option
                          ;

bvar_density : BVAR_DENSITY INT_NUMBER ';'
               { driver.bvar_density($2); }
             | BVAR_DENSITY '(' bvar_density_options_list ')' INT_NUMBER ';'
               { driver.bvar_density($5); }
             ;

bvar_forecast_option : bvar_common_option
                     | o_forecast
                     | o_conf_sig
                     | o_bvar_replic
                     ;

bvar_forecast_options_list : bvar_forecast_option COMMA bvar_forecast_options_list
                           | bvar_forecast_option
                           ;

bvar_forecast : BVAR_FORECAST INT_NUMBER ';'
                { driver.bvar_forecast($2); }
              | BVAR_FORECAST '(' bvar_forecast_options_list ')' INT_NUMBER ';'
                { driver.bvar_forecast($5); }
              ;

dynare_sensitivity : DYNARE_SENSITIVITY ';'
                     { driver.dynare_sensitivity(); }
                   | DYNARE_SENSITIVITY '(' dynare_sensitivity_options_list ')' ';'
                     { driver.dynare_sensitivity(); }
                   ;

dynare_sensitivity_options_list : dynare_sensitivity_option COMMA dynare_sensitivity_options_list
                                | dynare_sensitivity_option
                                ;

dynare_sensitivity_option : o_gsa_identification
                          | o_gsa_morris
                          | o_gsa_stab
                          | o_gsa_redform
                          | o_gsa_pprior
                          | o_gsa_prior_range
                          | o_gsa_ppost
                          | o_gsa_ilptau
                          | o_gsa_glue
                          | o_gsa_morris_nliv
                          | o_gsa_morris_ntra
                          | o_gsa_nsam
                          | o_gsa_load_redform
                          | o_gsa_load_rmse
                          | o_gsa_load_stab
                          | o_gsa_alpha2_stab
                          | o_gsa_ksstat
                          | o_gsa_logtrans_redform
                          | o_gsa_ksstat_redform
                          | o_gsa_alpha2_redform
                          | o_gsa_rmse
                          | o_gsa_lik_only
                          | o_gsa_pfilt_rmse
                          | o_gsa_istart_rmse
                          | o_gsa_alpha_rmse
                          | o_gsa_alpha2_rmse
                          | o_gsa_threshold_redform
                          | o_gsa_namendo
                          | o_gsa_namexo
                          | o_gsa_namlagendo
                          | o_gsa_var_rmse
                          | o_datafile
                          | o_nobs
                          | o_first_obs
                          | o_prefilter
                          | o_presample
                          | o_nograph
                          | o_conf_sig
                          | o_loglinear
                          | o_mode_file
                          | o_gsa_trans_ident
                          ;
 

homotopy_setup: HOMOTOPY_SETUP ';' homotopy_list END 
               { driver.end_homotopy();};

homotopy_list : homotopy_item
              | homotopy_list homotopy_item
              ;

homotopy_item : NAME COMMA expression COMMA expression ';'
                { driver.homotopy_val($1, $3, $5);}
              | NAME COMMA expression ';'
                { driver.homotopy_val($1, NULL, $3);}
              ;

forecast: FORECAST ';' {driver.forecast();}
          | FORECAST '(' forecast_options ')' ';' {driver.forecast();}
          | FORECAST symbol_list ';' {driver.forecast();}
          | FORECAST '(' forecast_options ')' symbol_list ';' {driver.forecast();}
          ;

forecast_options: forecast_option
          | forecast_options COMMA forecast_option
          ;

forecast_option: o_periods
          | o_conf_sig
          ;

number : INT_NUMBER
       | FLOAT_NUMBER
       ;

o_solve_algo : SOLVE_ALGO EQUAL INT_NUMBER { driver.option_num("solve_algo", $3); };
o_simul_algo : SIMUL_ALGO EQUAL INT_NUMBER { driver.option_num("simul_algo", $3); };
o_linear : LINEAR { driver.linear(); };
o_order : ORDER EQUAL INT_NUMBER { driver.option_num("order", $3); };
o_replic : REPLIC EQUAL INT_NUMBER { driver.option_num("replic", $3); };
o_drop : DROP EQUAL INT_NUMBER { driver.option_num("drop", $3); };
o_ar : AR EQUAL INT_NUMBER { driver.option_num("ar", $3); };
o_nocorr : NOCORR { driver.option_num("nocorr", "1"); };
o_nofunctions : NOFUNCTIONS { driver.option_num("nofunctions", "1"); };
o_nomoments : NOMOMENTS { driver.option_num("nomoments", "1"); };
o_irf : IRF EQUAL INT_NUMBER { driver.option_num("irf", $3); };
o_hp_filter : HP_FILTER EQUAL INT_NUMBER { driver.option_num("hp_filter", $3); };
o_hp_ngrid : HP_NGRID EQUAL INT_NUMBER { driver.option_num("hp_ngrid", $3); };
o_periods : PERIODS EQUAL INT_NUMBER
            { driver.option_num("periods", $3); driver.option_num("simul", "1"); };
o_cutoff : CUTOFF EQUAL number { driver.option_num("cutoff", $3); }
o_method : METHOD EQUAL INT_NUMBER { driver.option_num("simulation_method",$3);}
           | METHOD EQUAL LU { driver.option_num("simulation_method", "0"); }
           | METHOD EQUAL GAUSSIAN_ELIMINATION { driver.option_num("simulation_method", "1"); }
           | METHOD EQUAL GMRES { driver.option_num("simulation_method", "2"); }
           | METHOD EQUAL BICGSTAB { driver.option_num("simulation_method", "3"); };
o_markowitz : MARKOWITZ EQUAL number { driver.option_num("markowitz", $3); };
o_simul : SIMUL { driver.option_num("simul", "1"); };
o_simul_seed : SIMUL_SEED EQUAL INT_NUMBER { driver.option_num("simul_seed", $3)} ;
o_qz_criterium : QZ_CRITERIUM EQUAL number { driver.option_num("qz_criterium", $3) };
o_datafile : DATAFILE EQUAL NAME { driver.option_str("datafile", $3); };
o_nobs : NOBS EQUAL vec_int
         { driver.option_num("nobs", $3); }
       | NOBS EQUAL INT_NUMBER
         { driver.option_num("nobs", $3); }
       ;
o_first_obs : FIRST_OBS EQUAL INT_NUMBER { driver.option_num("first_obs", $3); };
o_prefilter : PREFILTER EQUAL INT_NUMBER { driver.option_num("prefilter", $3); };
o_presample : PRESAMPLE EQUAL INT_NUMBER { driver.option_num("presample", $3); };
o_lik_algo : LIK_ALGO EQUAL INT_NUMBER { driver.option_num("lik_algo", $3); };
o_lik_init : LIK_INIT EQUAL INT_NUMBER { driver.option_num("lik_init", $3); };
o_nograph : NOGRAPH
            { driver.option_num("nograph","1"); };
          | GRAPH
            { driver.option_num("nograph", "0"); }
          ;
o_conf_sig : CONF_SIG EQUAL number { driver.option_num("conf_sig", $3); };
o_mh_replic : MH_REPLIC EQUAL INT_NUMBER { driver.option_num("mh_replic", $3); };
o_mh_drop : MH_DROP EQUAL number { driver.option_num("mh_drop", $3); };
o_mh_jscale : MH_JSCALE EQUAL number { driver.option_num("mh_jscale", $3); };
o_optim : OPTIM  EQUAL '(' optim_options ')';
o_mh_init_scale : MH_INIT_SCALE EQUAL number { driver.option_num("mh_init_scale", $3); };
o_mode_file : MODE_FILE EQUAL NAME { driver.option_str("mode_file", $3); };
o_mode_compute : MODE_COMPUTE EQUAL INT_NUMBER { driver.option_num("mode_compute", $3); };
o_mode_check : MODE_CHECK { driver.option_num("mode_check", "1"); };
o_prior_trunc : PRIOR_TRUNC EQUAL number { driver.option_num("prior_trunc", $3); };
o_mh_mode : MH_MODE EQUAL INT_NUMBER { driver.option_num("mh_mode", $3); };
o_mh_nblcks : MH_NBLOCKS EQUAL INT_NUMBER { driver.option_num("mh_nblck", $3); };
o_load_mh_file : LOAD_MH_FILE { driver.option_num("load_mh_file", "1"); };
o_loglinear : LOGLINEAR { driver.option_num("loglinear", "1"); };
o_nodiagnostic : NODIAGNOSTIC { driver.option_num("nodiagnostic", "1"); };
o_bayesian_irf : BAYESIAN_IRF { driver.option_num("bayesian_irf", "1"); };
o_tex : TEX { driver.option_num("TeX", "1"); };
o_forecast : FORECAST EQUAL INT_NUMBER { driver.option_num("forecast", $3); };
o_smoother : SMOOTHER { driver.option_num("smoother", "1"); };
o_moments_varendo : MOMENTS_VARENDO { driver.option_num("moments_varendo", "1"); };
o_filtered_vars : FILTERED_VARS { driver.option_num("filtered_vars", "1"); };
o_relative_irf : RELATIVE_IRF { driver.option_num("relative_irf", "1"); };
o_kalman_algo : KALMAN_ALGO EQUAL INT_NUMBER { driver.option_num("kalman_algo", $3); };
o_kalman_tol : KALMAN_TOL EQUAL INT_NUMBER { driver.option_num("kalman_tol", $3); };
o_marginal_density : MARGINAL_DENSITY EQUAL LAPLACE
                     { driver.option_str("mc_marginal_density", "laplace"); }
                   | MARGINAL_DENSITY EQUAL MODIFIEDHARMONICMEAN
                     { driver.option_str("mc_marginal_density", "modifiedharmonicmean"); }
                   ;
o_print : PRINT { driver.option_num("noprint", "0"); };
o_noprint : NOPRINT { driver.option_num("noprint", "1"); };
o_xls_sheet : XLS_SHEET EQUAL NAME { driver.option_str("xls_sheet", $3); };
o_xls_range : XLS_RANGE EQUAL range { driver.option_str("xls_range", $3); };
o_filter_step_ahead : FILTER_STEP_AHEAD EQUAL vec_int { driver.option_num("filter_step_ahead", $3); };
o_constant : CONSTANT { driver.option_num("noconstant", "0"); };
o_noconstant : NOCONSTANT { driver.option_num("noconstant", "1"); };
o_mh_recover : MH_RECOVER { driver.option_num("mh_recover", "1"); };
o_diffuse_filter: DIFFUSE_FILTER {driver.option_num("diffuse_filter", "1"); };
o_planner_discount : PLANNER_DISCOUNT EQUAL number { driver.option_num("planner_discount",$3); };

o_bvar_prior_tau : BVAR_PRIOR_TAU EQUAL signed_float { driver.option_num("bvar_prior_tau", $3); };
o_bvar_prior_decay : BVAR_PRIOR_DECAY EQUAL number { driver.option_num("bvar_prior_decay", $3); };
o_bvar_prior_lambda : BVAR_PRIOR_LAMBDA EQUAL signed_float { driver.option_num("bvar_prior_lambda", $3); };
o_bvar_prior_mu : BVAR_PRIOR_MU EQUAL number { driver.option_num("bvar_prior_mu", $3); };
o_bvar_prior_omega : BVAR_PRIOR_OMEGA EQUAL INT_NUMBER { driver.option_num("bvar_prior_omega", $3); };
o_bvar_prior_flat : BVAR_PRIOR_FLAT { driver.option_num("bvar_prior_flat", "1"); };
o_bvar_prior_train : BVAR_PRIOR_TRAIN EQUAL INT_NUMBER { driver.option_num("bvar_prior_train", $3); };
o_bvar_replic : BVAR_REPLIC EQUAL INT_NUMBER { driver.option_num("bvar_replic", $3); };

o_gsa_identification : IDENTIFICATION EQUAL INT_NUMBER { driver.option_num("identification", $3); }; /*not in doc */
o_gsa_morris : MORRIS EQUAL INT_NUMBER { driver.option_num("morris", $3); };
o_gsa_stab : STAB  EQUAL INT_NUMBER { driver.option_num("stab", $3); };
o_gsa_redform : REDFORM  EQUAL INT_NUMBER { driver.option_num("redform", $3); };
o_gsa_pprior : PPRIOR  EQUAL INT_NUMBER { driver.option_num("pprior", $3); };
o_gsa_prior_range : PRIOR_RANGE  EQUAL INT_NUMBER { driver.option_num("prior_range", $3); };
o_gsa_ppost : PPOST  EQUAL INT_NUMBER { driver.option_num("ppost", $3); };
o_gsa_ilptau : ILPTAU EQUAL INT_NUMBER { driver.option_num("ilptau", $3); };
o_gsa_glue : GLUE EQUAL INT_NUMBER { driver.option_num("glue", $3); };
o_gsa_morris_nliv : MORRIS_NLIV EQUAL INT_NUMBER { driver.option_num("morris_nliv", $3); };
o_gsa_morris_ntra : MORRIS_NTRA EQUAL INT_NUMBER { driver.option_num("morris_ntra", $3); };
o_gsa_nsam : NSAM EQUAL INT_NUMBER { driver.option_num("Nsam", $3); }; /* not in doc ??*/
o_gsa_load_redform : LOAD_REDFORM EQUAL INT_NUMBER { driver.option_num("load_redform", $3); };
o_gsa_load_rmse : LOAD_RMSE EQUAL INT_NUMBER { driver.option_num("load_rmse", $3); };
o_gsa_load_stab : LOAD_STAB EQUAL INT_NUMBER { driver.option_num("load_stab", $3); };
o_gsa_alpha2_stab : ALPHA2_STAB EQUAL number { driver.option_num("alpha2_stab", $3); };
o_gsa_ksstat : KSSTAT EQUAL number { driver.option_num("ksstat", $3); };
o_gsa_logtrans_redform : LOGTRANS_REDFORM EQUAL INT_NUMBER { driver.option_num("logtrans_redform", $3); };
o_gsa_threshold_redform : THRESHOLD_REDFORM EQUAL vec_value { driver.option_num("threshold_redform",$3); };
o_gsa_ksstat_redform : KSSTAT_REDFORM EQUAL number { driver.option_num("ksstat_redfrom", $3); };
o_gsa_alpha2_redform : ALPHA2_REDFORM EQUAL number { driver.option_num("alpha2_redform", $3); };
o_gsa_namendo : NAMENDO EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("namendo"); };
o_gsa_namlagendo : NAMLAGENDO EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("namlagendo"); };
o_gsa_namexo : NAMEXO EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("namexo"); };
o_gsa_rmse : RMSE EQUAL INT_NUMBER { driver.option_num("rmse", $3); };
o_gsa_lik_only : LIK_ONLY EQUAL INT_NUMBER { driver.option_num("lik_only", $3); };
o_gsa_var_rmse : VAR_RMSE EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("var_rmse"); };
o_gsa_pfilt_rmse : PFILT_RMSE EQUAL number { driver.option_num("pfilt_rmse", $3); };
o_gsa_istart_rmse : ISTART_RMSE EQUAL INT_NUMBER { driver.option_num("istart_rmse", $3); };
o_gsa_alpha_rmse : ALPHA_RMSE EQUAL number { driver.option_num("alpha_rmse", $3); };
o_gsa_alpha2_rmse : ALPHA2_RMSE EQUAL number { driver.option_num("alpha2_rmse", $3); };
o_gsa_trans_ident : TRANS_IDENT EQUAL INT_NUMBER { driver.option_num("trans_ident", $3); };

o_homotopy_mode : HOMOTOPY_MODE EQUAL INT_NUMBER {driver.option_num("homotopy_mode",$3); };
o_homotopy_steps : HOMOTOPY_STEPS EQUAL INT_NUMBER {driver.option_num("homotopy_steps",$3); };

range : NAME ':' NAME
        {
          $1->append(":");
          $1->append(*$3);
          delete $3;
          $$ = $1;
        };

vec_int_elem : INT_NUMBER
             | INT_NUMBER ':' INT_NUMBER
               {
                 $1->append(":");
                 $1->append(*$3);
                 delete $3;
                 $$ = $1;
               }
             ;

vec_int_1 : '[' vec_int_elem
            { $2->insert(0, "["); $$ = $2;}
          | vec_int_1 vec_int_elem
            {
              $1->append(" ");
              $1->append(*$2);
              delete $2;
              $$ = $1;
            }
          ;

vec_int : vec_int_1 ']' { $1->append("]"); $$ = $1; };


vec_value_1 : '[' value1
            { $2->insert(0, "["); $$ = $2;}
          | vec_value_1 value1
            {
              $1->append(" ");
              $1->append(*$2);
              delete $2;
              $$ = $1;
            }
          ;

vec_value : vec_value_1 ']' { $1->append("]"); $$ = $1; };

%%

void
Dynare::parser::error(const Dynare::parser::location_type &l,
                      const string &m)
{
  driver.error(l, m);
}

/*
  Local variables:
  mode: C++
  End:
*/
