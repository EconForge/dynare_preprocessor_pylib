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

%{
using namespace std;

#include "ParsingDriver.hh"
#include "DynareBison.hh"

// Announce to Flex the prototype we want for lexing function
#define YY_DECL                                                \
  Dynare::parser::token_type                                   \
    DynareFlex::lex(Dynare::parser::semantic_type *yylval,     \
                    Dynare::parser::location_type *yylloc,     \
                    ParsingDriver &driver)
 
// Shortcut to access tokens defined by Bison
typedef Dynare::parser::token token;

/* By default yylex returns int, we use token_type.
   Unfortunately yyterminate by default returns 0, which is
   not of token_type.  */
#define yyterminate() return Dynare::parser::token_type (0);

int comment_caller, line_caller;
/* Particular value : when sigma_e command is found
 this flag is set to 1, when command finished it is set to 0
 */
int sigma_e = 0;
%}

%option c++

%option prefix="Dynare"

%option case-insensitive noyywrap nounput batch debug never-interactive

%x COMMENT
%x DYNARE_STATEMENT
%x DYNARE_BLOCK
%x NATIVE
%x LINE1
%x LINE2
%x LINE3

%{
// Increments location counter for every token read
#define YY_USER_ACTION yylloc->columns(yyleng);
%}
%%
 /* Code put at the beginning of yylex() */
%{
  // Reset location before reading token
  yylloc->step();
%}

 /* Rules for matching $line directives */
<*>^@#line\ \"  { line_caller = YYSTATE; BEGIN(LINE1); }
<LINE1>[^\"]*   {
                  if (yylloc->begin.filename)
                    delete yylloc->begin.filename;
                  yylloc->begin.filename = yylloc->end.filename = new string(yytext);
                  BEGIN(LINE2);
                }
<LINE2>\"       BEGIN(LINE3);
<LINE3>[0-9]+   {
                  yylloc->begin.line = yylloc->end.line = atoi(yytext) - 1;
                  BEGIN(line_caller);
                }

 /* spaces, tabs and carriage returns are ignored */
<*>[ \t\r\f]+  { yylloc->step(); }
<*>[\n]+       { yylloc->lines(yyleng); yylloc->step(); }

 /* Comments */
<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK>["%"].*
<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK>["/"]["/"].*
<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK>"/*"   {comment_caller = YY_START; BEGIN COMMENT;}

<COMMENT>"*/"        {BEGIN comment_caller;}
<COMMENT>.

 /* Begin of a Dynare statement */
<INITIAL>var {BEGIN DYNARE_STATEMENT; return token::VAR;}
<INITIAL>varexo {BEGIN DYNARE_STATEMENT; return token::VAREXO;}
<INITIAL>varexo_det {BEGIN DYNARE_STATEMENT; return token::VAREXO_DET;}
<INITIAL>parameters {BEGIN DYNARE_STATEMENT; return token::PARAMETERS;}
<INITIAL>periods 	{BEGIN DYNARE_STATEMENT; return token::PERIODS;}
<INITIAL>cutoff 	{BEGIN DYNARE_STATEMENT; return token::CUTOFF;}
<INITIAL>markowitz 	{BEGIN DYNARE_STATEMENT; return token::MARKOWITZ;}
<INITIAL>estimation {BEGIN DYNARE_STATEMENT; return token::ESTIMATION;}
<INITIAL>prior_analysis {BEGIN DYNARE_STATEMENT; return token::PRIOR_ANALYSIS;}
<INITIAL>posterior_analysis {BEGIN DYNARE_STATEMENT; return token::POSTERIOR_ANALYSIS;}
<INITIAL>varobs 	{BEGIN DYNARE_STATEMENT; return token::VAROBS;}
<INITIAL>unit_root_vars	{BEGIN DYNARE_STATEMENT; return token::UNIT_ROOT_VARS;}
<INITIAL>rplot	 	{BEGIN DYNARE_STATEMENT; return token::RPLOT;}
<INITIAL>osr_params 	{BEGIN DYNARE_STATEMENT; return token::OSR_PARAMS;}
<INITIAL>osr	 	{BEGIN DYNARE_STATEMENT; return token::OSR;}
<INITIAL>dynatype	{BEGIN DYNARE_STATEMENT; return token::DYNATYPE;}
<INITIAL>dynasave 	{BEGIN DYNARE_STATEMENT; return token::DYNASAVE;}
<INITIAL>model_comparison 	{BEGIN DYNARE_STATEMENT; return token::MODEL_COMPARISON;}

<INITIAL>steady {BEGIN DYNARE_STATEMENT; return token::STEADY;}
<INITIAL>check {BEGIN DYNARE_STATEMENT; return token::CHECK;}
<INITIAL>simul {BEGIN DYNARE_STATEMENT; return token::SIMUL;}
<INITIAL>stoch_simul {BEGIN DYNARE_STATEMENT; return token::STOCH_SIMUL;}
<INITIAL>dsample {BEGIN DYNARE_STATEMENT; return token::DSAMPLE;}
<INITIAL>Sigma_e {BEGIN DYNARE_STATEMENT; sigma_e = 1; return token::SIGMA_E;}
<INITIAL>calib {BEGIN DYNARE_STATEMENT; return token::CALIB;}
<INITIAL>planner_objective {BEGIN DYNARE_STATEMENT; return token::PLANNER_OBJECTIVE;}
<INITIAL>ramsey_policy {BEGIN DYNARE_STATEMENT; return token::RAMSEY_POLICY;}

<INITIAL>bvar_density {BEGIN DYNARE_STATEMENT; return token::BVAR_DENSITY; }
<INITIAL>bvar_forecast {BEGIN DYNARE_STATEMENT; return token::BVAR_FORECAST; }
<INITIAL>dynare_sensitivity {BEGIN DYNARE_STATEMENT; return token::DYNARE_SENSITIVITY;}
<INITIAL>initval_file {BEGIN DYNARE_STATEMENT; return token::INITVAL_FILE;}
 /* End of a Dynare statement */

<DYNARE_STATEMENT>; {
  if (!sigma_e)
    BEGIN INITIAL;
  return Dynare::parser::token_type (yytext[0]);
}


 /* Begin of a Dynare block */
<INITIAL>model {BEGIN DYNARE_BLOCK; return token::MODEL;}
<INITIAL>initval {BEGIN DYNARE_BLOCK; return token::INITVAL;}
<INITIAL>endval {BEGIN DYNARE_BLOCK; return token::ENDVAL;}
<INITIAL>histval {BEGIN DYNARE_BLOCK; return token::HISTVAL;}
<INITIAL>shocks {BEGIN DYNARE_BLOCK; return token::SHOCKS;}
<INITIAL>estimated_params {BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS;}
 /* priors is an alias for estimated_params */
<INITIAL>priors {BEGIN DYNARE_BLOCK;return token::ESTIMATED_PARAMS;}
<INITIAL>estimated_params_init 		{BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS_INIT;}
<INITIAL>estimated_params_bounds 	{BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS_BOUNDS;}
<INITIAL>observation_trends {BEGIN DYNARE_BLOCK; return token::OBSERVATION_TRENDS;}
<INITIAL>optim_weights {BEGIN DYNARE_BLOCK; return token::OPTIM_WEIGHTS;}
<INITIAL>calib_var 	{BEGIN DYNARE_BLOCK; return token::CALIB_VAR;}
<INITIAL>homotopy_setup {BEGIN DYNARE_BLOCK; return token::HOMOTOPY_SETUP;}

 /* End of a Dynare block */
<DYNARE_BLOCK>end[ \t\n]*; 	{BEGIN INITIAL; return token::END;}

 /* Inside  of a Dynare statement */
<DYNARE_STATEMENT>datafile 		{return token::DATAFILE;}
<DYNARE_STATEMENT>method    {return token::METHOD;}
<DYNARE_STATEMENT>nobs 			{return token::NOBS;}
<DYNARE_STATEMENT>first_obs 		{return token::FIRST_OBS;}
<DYNARE_STATEMENT>prefilter 		{return token::PREFILTER;}
<DYNARE_STATEMENT>presample 		{return token::PRESAMPLE;}
<DYNARE_STATEMENT>lik_algo  		{return token::LIK_ALGO;}
<DYNARE_STATEMENT>lik_init  		{return token::LIK_INIT;}
<DYNARE_STATEMENT>graph   		{return token::GRAPH;}
<DYNARE_STATEMENT>nograph   		{return token::NOGRAPH;}
<DYNARE_STATEMENT>print   		{return token::PRINT;}
<DYNARE_STATEMENT>noprint   		{return token::NOPRINT;}
<DYNARE_STATEMENT>conf_sig  		{return token::CONF_SIG;}
<DYNARE_STATEMENT>mh_replic 		{return token::MH_REPLIC;}
<DYNARE_STATEMENT>mh_drop   		{return token::MH_DROP;}
<DYNARE_STATEMENT>mh_jscale   		{return token::MH_JSCALE;}
<DYNARE_STATEMENT>mh_init_scale 	{return token::MH_INIT_SCALE;}
<DYNARE_STATEMENT>mode_file 		{return token::MODE_FILE;}
<DYNARE_STATEMENT>mode_compute 	{return token::MODE_COMPUTE;}
<DYNARE_STATEMENT>mode_check 		{return token::MODE_CHECK;}
<DYNARE_STATEMENT>prior_trunc 	{return token::PRIOR_TRUNC;}
<DYNARE_STATEMENT>mh_mode 		{return token::MH_MODE;}
<DYNARE_STATEMENT>mh_nblocks 		{return token::MH_NBLOCKS;}
<DYNARE_STATEMENT>load_mh_file 	{return token::LOAD_MH_FILE;}
<DYNARE_STATEMENT>loglinear 		{return token::LOGLINEAR;}
<DYNARE_STATEMENT>nodiagnostic 	{return token::NODIAGNOSTIC;}
<DYNARE_STATEMENT>kalman_algo 	{return token::KALMAN_ALGO;}
<DYNARE_STATEMENT>kalman_tol 	{return token::KALMAN_TOL;}
<DYNARE_STATEMENT>forecast 	{return token::FORECAST;}
<DYNARE_STATEMENT>smoother 	{return token::SMOOTHER;}
<DYNARE_STATEMENT>bayesian_irf 	{return token::BAYESIAN_IRF;}
<DYNARE_STATEMENT>moments_varendo {return token::MOMENTS_VARENDO;}
<DYNARE_STATEMENT>filtered_vars	{return token::FILTERED_VARS;}
<DYNARE_STATEMENT>filter_step_ahead	{return token::FILTER_STEP_AHEAD;}
<DYNARE_STATEMENT>relative_irf 	{return token::RELATIVE_IRF;}
<DYNARE_STATEMENT>tex		{return token::TEX;}
<DYNARE_STATEMENT>nomoments	{return token::NOMOMENTS;}
<DYNARE_STATEMENT>corr		{return token::CORR;}
<DYNARE_STATEMENT>nocorr	{return token::NOCORR;}
<DYNARE_STATEMENT>optim		{return token::OPTIM;}
<DYNARE_STATEMENT>periods	{return token::PERIODS;}
<DYNARE_STATEMENT>cutoff	{return token::CUTOFF;}
<DYNARE_STATEMENT>markowitz	{return token::MARKOWITZ;}
<DYNARE_STATEMENT>model_comparison_approximation {return token::MODEL_COMPARISON;}
<DYNARE_STATEMENT>laplace {return token::LAPLACE;}
<DYNARE_STATEMENT>modifiedharmonicmean {return token::MODIFIEDHARMONICMEAN;}
<DYNARE_STATEMENT>constant	{return token::CONSTANT;}
<DYNARE_STATEMENT>noconstant	{return token::NOCONSTANT;}
<DYNARE_STATEMENT>covar {return token::COVAR;}
<DYNARE_STATEMENT>filename {return token::FILENAME;}

<DYNARE_STATEMENT>bvar_prior_tau { return token::BVAR_PRIOR_TAU; }
<DYNARE_STATEMENT>bvar_prior_decay { return token::BVAR_PRIOR_DECAY; }
<DYNARE_STATEMENT>bvar_prior_lambda { return token::BVAR_PRIOR_LAMBDA; }
<DYNARE_STATEMENT>bvar_prior_mu { return token::BVAR_PRIOR_MU; }
<DYNARE_STATEMENT>bvar_prior_omega { return token::BVAR_PRIOR_OMEGA; }
<DYNARE_STATEMENT>bvar_prior_flat { return token::BVAR_PRIOR_FLAT; }
<DYNARE_STATEMENT>bvar_prior_train { return token::BVAR_PRIOR_TRAIN; }
<DYNARE_STATEMENT>bvar_replic { return token::BVAR_REPLIC; }

<DYNARE_STATEMENT>homotopy_mode {return token::HOMOTOPY_MODE; }
<DYNARE_STATEMENT>homotopy_steps {return token::HOMOTOPY_STEPS; }

<DYNARE_STATEMENT>[\$][^$]*[\$] {
  strtok(yytext+1, "$");
  yylval->string_val = new string(yytext + 1);
  return token::TEX_NAME;
}

 /* Inside a Dynare block */
<DYNARE_BLOCK>var {return token::VAR;}
<DYNARE_BLOCK>stderr {return token::STDERR;}
<DYNARE_BLOCK>values {return token::VALUES;}
<DYNARE_BLOCK>corr {return token::CORR;}
<DYNARE_BLOCK>periods {return token::PERIODS;}
<DYNARE_BLOCK>cutoff {return token::CUTOFF;}
<DYNARE_BLOCK>markowitz {return token::MARKOWITZ;}
<DYNARE_BLOCK>gamma_pdf {return token::GAMMA_PDF;}
<DYNARE_BLOCK>beta_pdf {return token::BETA_PDF;}
<DYNARE_BLOCK>normal_pdf {return token::NORMAL_PDF;}
<DYNARE_BLOCK>inv_gamma_pdf {return token::INV_GAMMA_PDF;}
<DYNARE_BLOCK>inv_gamma1_pdf {return token::INV_GAMMA_PDF;}
<DYNARE_BLOCK>inv_gamma2_pdf {return token::INV_GAMMA_PDF;}
<DYNARE_BLOCK>uniform_pdf {return token::UNIFORM_PDF;}

<DYNARE_BLOCK>; {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_BLOCK># {return Dynare::parser::token_type (yytext[0]);}

<DYNARE_BLOCK>autocorr {return token::AUTOCORR;}

 /* Inside Dynare statement */
<DYNARE_STATEMENT>solve_algo {return token::SOLVE_ALGO;}
<DYNARE_STATEMENT>dr_algo {return token::DR_ALGO;}
<DYNARE_STATEMENT>simul_algo {return token::SIMUL_ALGO;}
<DYNARE_STATEMENT>drop {return token::DROP;}
<DYNARE_STATEMENT>order {return token::ORDER;}
<DYNARE_STATEMENT>replic {return token::REPLIC;}
<DYNARE_STATEMENT>ar {return token::AR;}
<DYNARE_STATEMENT>nofunctions {return token::NOFUNCTIONS;}
<DYNARE_STATEMENT>irf {return token::IRF;}
<DYNARE_STATEMENT>hp_filter {return token::HP_FILTER;}
<DYNARE_STATEMENT>hp_ngrid {return token::HP_NGRID;}
<DYNARE_STATEMENT>simul_seed {return token::SIMUL_SEED;}
<DYNARE_STATEMENT>qz_criterium {return token::QZ_CRITERIUM;}
<DYNARE_STATEMENT>simul {return token::SIMUL;}
<DYNARE_STATEMENT>xls_sheet {return token::XLS_SHEET;}
<DYNARE_STATEMENT>xls_range {return token::XLS_RANGE;}
<DYNARE_STATEMENT>mh_recover {return token::MH_RECOVER;}
<DYNARE_STATEMENT>planner_discount {return token::PLANNER_DISCOUNT;}


<DYNARE_STATEMENT>[\.] {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT>[\\] {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT>[\'] {return Dynare::parser::token_type (yytext[0]);}

<DYNARE_STATEMENT,DYNARE_BLOCK>use_dll {return token::USE_DLL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>lu {return token::LU;}
<DYNARE_STATEMENT,DYNARE_BLOCK>gaussian_elimination {return token::GAUSSIAN_ELIMINATION;}
<DYNARE_STATEMENT,DYNARE_BLOCK>gmres {return token::GMRES;}
<DYNARE_STATEMENT,DYNARE_BLOCK>bicgstab {return token::BICGSTAB;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sparse {return token::SPARSE;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sparse_dll {return token::SPARSE_DLL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>gcc_compiler {return token::GCC_COMPILER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>lcc_compiler {return token::LCC_COMPILER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>no_compiler {return token::NO_COMPILER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>linear {return token::LINEAR;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[,] {return token::COMMA;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[:] {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\(\)] {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\[] {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\]] {
  if (sigma_e)
    sigma_e = 0;
  return Dynare::parser::token_type (yytext[0]);
}
<DYNARE_STATEMENT,DYNARE_BLOCK>[+] {return token::PLUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[-] {return token::MINUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[*] {return token::TIMES;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[/] {return token::DIVIDE;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[=] {return token::EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[<] {return token::LESS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[>] {return token::GREATER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>">=" {return token::GREATER_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>"<=" {return token::LESS_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>"==" {return token::EQUAL_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>"!=" {return token::EXCLAMATION_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\^] {return token::POWER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>exp {return token::EXP;}
<DYNARE_STATEMENT,DYNARE_BLOCK>log {return token::LOG;}
<DYNARE_STATEMENT,DYNARE_BLOCK>log10 {return token::LOG10;}
<DYNARE_STATEMENT,DYNARE_BLOCK>ln {return token::LN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sin {return token::SIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>cos {return token::COS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>tan {return token::TAN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>asin {return token::ASIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>acos {return token::ACOS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>atan {return token::ATAN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sinh {return token::SINH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>cosh {return token::COSH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>tanh {return token::TANH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>asinh {return token::ASINH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>acosh {return token::ACOSH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>atanh {return token::ATANH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sqrt {return token::SQRT;}
<DYNARE_STATEMENT,DYNARE_BLOCK>max {return token::MAX;}
<DYNARE_STATEMENT,DYNARE_BLOCK>min {return token::MIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>normcdf {return token::NORMCDF;}

 /* options for GSA module by Marco Ratto */
<DYNARE_STATEMENT>identification {return token::IDENTIFICATION;}
<DYNARE_STATEMENT>morris {return token::MORRIS;}
<DYNARE_STATEMENT>stab {return token::STAB;}
<DYNARE_STATEMENT>redform {return token::REDFORM;}
<DYNARE_STATEMENT>pprior {return token::PPRIOR;}
<DYNARE_STATEMENT>prior_range {return token::PRIOR_RANGE;}
<DYNARE_STATEMENT>ppost {return token::PPOST;}
<DYNARE_STATEMENT>ilptau {return token::ILPTAU;}
<DYNARE_STATEMENT>glue {return token::GLUE;}
<DYNARE_STATEMENT>morris_nliv {return token::MORRIS_NLIV;}
<DYNARE_STATEMENT>morris_ntra {return token::MORRIS_NTRA;}
<DYNARE_STATEMENT>Nsam {return token::NSAM;}
<DYNARE_STATEMENT>load_redform {return token::LOAD_REDFORM;}
<DYNARE_STATEMENT>load_rmse {return token::LOAD_RMSE;}
<DYNARE_STATEMENT>load_stab {return token::LOAD_STAB;}
<DYNARE_STATEMENT>alpha2_stab {return token::ALPHA2_STAB;}
<DYNARE_STATEMENT>ksstat {return token::KSSTAT;}
<DYNARE_STATEMENT>logtrans_redform {return token::LOGTRANS_REDFORM;}
<DYNARE_STATEMENT>threshold_redform {return token::THRESHOLD_REDFORM;}
<DYNARE_STATEMENT>ksstat_redform {return token::KSSTAT_REDFORM;}
<DYNARE_STATEMENT>alpha2_redform {return token::ALPHA2_REDFORM;}
<DYNARE_STATEMENT>namendo {return token::NAMENDO;}
<DYNARE_STATEMENT>namlagendo {return token::NAMLAGENDO;}
<DYNARE_STATEMENT>namexo {return token::NAMEXO;}
<DYNARE_STATEMENT>rmse {return token::RMSE;}
<DYNARE_STATEMENT>lik_only {return token::LIK_ONLY;}
<DYNARE_STATEMENT>var_rmse {return token::VAR_RMSE;}
<DYNARE_STATEMENT>pfilt_rmse {return token::PFILT_RMSE;}
<DYNARE_STATEMENT>istart_rmse {return token::ISTART_RMSE;}
<DYNARE_STATEMENT>alpha_rmse {return token::ALPHA_RMSE;}
<DYNARE_STATEMENT>alpha2_rmse {return token::ALPHA2_RMSE;}
<DYNARE_STATEMENT>trans_ident {return token::TRANS_IDENT;}
 /* end of GSA options */

<DYNARE_STATEMENT,DYNARE_BLOCK>[A-Za-z_][A-Za-z0-9_]* {
  yylval->string_val = new string(yytext);
  return token::NAME;
}

<DYNARE_STATEMENT,DYNARE_BLOCK>((([0-9]*\.[0-9]+)|([0-9]+\.))([edED][-+]?[0-9]+)?)|([0-9]+[edED][-+]?[0-9]+) {
  yylval->string_val = new string(yytext);
  return token::FLOAT_NUMBER;
}

<DYNARE_STATEMENT,DYNARE_BLOCK>[0-9]+ {
  yylval->string_val = new string(yytext);
  return token::INT_NUMBER;
}

 /* an instruction starting with a recognized symbol (which is not a modfile local variable)
    is passed as NAME,
    otherwise it is a native statement until the end of the line
 */
<INITIAL>[A-Za-z_][A-Za-z0-9_]* {
  if (driver.symbol_exists_and_is_not_modfile_local_variable(yytext))
    {
      BEGIN DYNARE_STATEMENT;
      yylval->string_val = new string(yytext);
      return token::NAME;
    }
  else
    {
      /* Enter a native block */
      BEGIN NATIVE;
      yyless(0);
    }
}

 /* Enter a native block */
<INITIAL>. { BEGIN NATIVE; yyless(0); }

 /* Add the native statement */
<NATIVE>.* { driver.add_native(yytext); BEGIN INITIAL; }

<*><<EOF>> {
             if (yylloc->begin.filename)
               delete yylloc->begin.filename;
             yyterminate();
           }

<*>.      { driver.error(*yylloc, "character unrecognized by lexer"); }
%%

DynareFlex::DynareFlex(istream* in, ostream* out)
  : DynareFlexLexer(in, out)
{
}

/* This implementation of DynareFlexLexer::yylex() is required to fill the
 * vtable of the class DynareFlexLexer. We define the scanner's main yylex
 * function via YY_DECL to reside in the DynareFlex class instead. */

#ifdef yylex
# undef yylex
#endif

int
DynareFlexLexer::yylex()
{
  cerr << "DynareFlexLexer::yylex() has been called, that should never happen!" << endl;
  exit(-1);
}

/*
  Local variables:
  mode: C++
  End:
*/
