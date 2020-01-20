/* -*- C++ -*- */
/*
 * Copyright © 2003-2019 Dynare Team
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

#include <cstring>
#include "ParsingDriver.hh"

// Announce to Flex the prototype we want for lexing function
#define YY_DECL                                                \
  Dynare::parser::token_type                                   \
    DynareFlex::lex(Dynare::parser::semantic_type *yylval,     \
                    Dynare::parser::location_type *yylloc,     \
                    ParsingDriver &driver)

// Shortcut to access tokens defined by Bison
using token = Dynare::parser::token;

/* By default yylex returns int, we use token_type.
   Unfortunately yyterminate by default returns 0, which is
   not of token_type.  */
#define yyterminate() return Dynare::parser::token_type (0);

int comment_caller, line_caller;
/* Particular value : when sigma_e command is found
 this flag is set to 1, when command finished it is set to 0
 */
int sigma_e = 0;
string eofbuff;
%}

%option c++

%option prefix="Dynare"

%option case-insensitive noyywrap batch debug never-interactive

 /* NB: if new start conditions are defined, add them in the line for <<EOF>> */
%x COMMENT
%x DYNARE_STATEMENT
%x DYNARE_BLOCK
%x VERBATIM_BLOCK
%x NATIVE
%x NATIVE_COMMENT
%x DATES_STATEMENT
%x LINE1
%x LINE2
%x LINE3

%{
// Increments location counter for every token read
#define YY_USER_ACTION location_increment(yylloc, yytext);
%}

DATE -?[0-9]+([ya]|m([1-9]|1[0-2])|q[1-4]|w([1-9]{1}|[1-4][0-9]|5[0-2]))

%%
 /* Code put at the beginning of yylex() */
%{
  // Reset location before reading token
  yylloc->step();
%}

 /* Rules for matching $line directives */
<*>^@#line\ \"  { line_caller = YYSTATE; BEGIN(LINE1); }
<LINE1>[^\"]*   {
                  filename = string(yytext);
                  BEGIN(LINE2);
                }
<LINE2>\"       BEGIN(LINE3);
<LINE3>[0-9]+   {
                  yylloc->begin.line = yylloc->end.line = atoi(yytext) - 1;
                  BEGIN(line_caller);
                }

 /* spaces, tabs and carriage returns are ignored */
<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK,COMMENT,DATES_STATEMENT,LINE1,LINE2,LINE3>[[:space:]]+  { yylloc->step(); }

 /* Comments */
<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK,DATES_STATEMENT>%.*
<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK,DATES_STATEMENT>"//".*
<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK,DATES_STATEMENT>"/*"   {comment_caller = YY_START; BEGIN COMMENT;}

<COMMENT>"*/"        {BEGIN comment_caller;}
<COMMENT>.

 /* Begin of a Dynare statement */
<INITIAL>var {BEGIN DYNARE_STATEMENT; return token::VAR;}
<INITIAL>varexo {BEGIN DYNARE_STATEMENT; return token::VAREXO;}
<INITIAL>varexo_det {BEGIN DYNARE_STATEMENT; return token::VAREXO_DET;}
<INITIAL>trend_var {BEGIN DYNARE_STATEMENT; return token::TREND_VAR;}
<INITIAL>log_trend_var {BEGIN DYNARE_STATEMENT; return token::LOG_TREND_VAR;}
<INITIAL>predetermined_variables {BEGIN DYNARE_STATEMENT; return token::PREDETERMINED_VARIABLES;}
<INITIAL>parameters {BEGIN DYNARE_STATEMENT; return token::PARAMETERS;}
<INITIAL>model_local_variable {BEGIN DYNARE_STATEMENT; return token::MODEL_LOCAL_VARIABLE;}
<INITIAL>periods 	{BEGIN DYNARE_STATEMENT; return token::PERIODS;}
<INITIAL>model_info {BEGIN DYNARE_STATEMENT; return token::MODEL_INFO;}
<INITIAL>estimation {BEGIN DYNARE_STATEMENT; return token::ESTIMATION;}
<INITIAL>var_estimation {BEGIN DYNARE_STATEMENT; return token::VAR_ESTIMATION;}
<INITIAL>set_time {BEGIN DYNARE_STATEMENT; return token::SET_TIME;}
<INITIAL>data {BEGIN DYNARE_STATEMENT; return token::DATA;}
<INITIAL>varobs 	{BEGIN DYNARE_STATEMENT; return token::VAROBS;}
<INITIAL>varexobs 	{BEGIN DYNARE_STATEMENT; return token::VAREXOBS;}
<INITIAL>unit_root_vars	{BEGIN DYNARE_STATEMENT; return token::UNIT_ROOT_VARS;}
<INITIAL>rplot	 	{BEGIN DYNARE_STATEMENT; return token::RPLOT;}
<INITIAL>osr_params 	{BEGIN DYNARE_STATEMENT; return token::OSR_PARAMS;}
<INITIAL>osr	 	{BEGIN DYNARE_STATEMENT; return token::OSR;}
<INITIAL>dynatype	{BEGIN DYNARE_STATEMENT; return token::DYNATYPE;}
<INITIAL>dynasave 	{BEGIN DYNARE_STATEMENT; return token::DYNASAVE;}
<INITIAL>model_comparison 	{BEGIN DYNARE_STATEMENT; return token::MODEL_COMPARISON;}
<INITIAL>change_type  {BEGIN DYNARE_STATEMENT; return token::CHANGE_TYPE;}
<INITIAL>load_params_and_steady_state  {BEGIN DYNARE_STATEMENT; return token::LOAD_PARAMS_AND_STEADY_STATE;}
<INITIAL>save_params_and_steady_state  {BEGIN DYNARE_STATEMENT; return token::SAVE_PARAMS_AND_STEADY_STATE;}
<INITIAL>write_latex_dynamic_model  {BEGIN DYNARE_STATEMENT; return token::WRITE_LATEX_DYNAMIC_MODEL;}
<INITIAL>write_latex_static_model  {BEGIN DYNARE_STATEMENT; return token::WRITE_LATEX_STATIC_MODEL;}
<INITIAL>write_latex_original_model  {BEGIN DYNARE_STATEMENT; return token::WRITE_LATEX_ORIGINAL_MODEL;}
<INITIAL>write_latex_steady_state_model  {BEGIN DYNARE_STATEMENT; return token::WRITE_LATEX_STEADY_STATE_MODEL;}

<INITIAL>steady {BEGIN DYNARE_STATEMENT; return token::STEADY;}
<INITIAL>check {BEGIN DYNARE_STATEMENT; return token::CHECK;}
<INITIAL>simul {BEGIN DYNARE_STATEMENT; return token::SIMUL;}
<INITIAL>stoch_simul {BEGIN DYNARE_STATEMENT; return token::STOCH_SIMUL;}
<INITIAL>var_model {BEGIN DYNARE_STATEMENT; return token::VAR_MODEL;}
<INITIAL>trend_component_model {BEGIN DYNARE_STATEMENT; return token::TREND_COMPONENT_MODEL;}
<INITIAL>var_expectation_model {BEGIN DYNARE_STATEMENT; return token::VAR_EXPECTATION_MODEL;}
<INITIAL>pac_model {BEGIN DYNARE_STATEMENT; return token::PAC_MODEL;}
<INITIAL>dsample {BEGIN DYNARE_STATEMENT; return token::DSAMPLE;}
<INITIAL>Sigma_e {BEGIN DYNARE_STATEMENT; sigma_e = 1; return token::SIGMA_E;}
<INITIAL>planner_objective {BEGIN DYNARE_STATEMENT; return token::PLANNER_OBJECTIVE;}
<INITIAL>ramsey_model {BEGIN DYNARE_STATEMENT; return token::RAMSEY_MODEL;}
<INITIAL>ramsey_policy {BEGIN DYNARE_STATEMENT; return token::RAMSEY_POLICY;}
<INITIAL>evaluate_planner_objective {BEGIN DYNARE_STATEMENT; return token::EVALUATE_PLANNER_OBJECTIVE;}
<INITIAL>discretionary_policy {BEGIN DYNARE_STATEMENT; return token::DISCRETIONARY_POLICY;}
<INITIAL>identification {BEGIN DYNARE_STATEMENT; return token::IDENTIFICATION;}

<INITIAL>bvar_density {BEGIN DYNARE_STATEMENT; return token::BVAR_DENSITY; }
<INITIAL>bvar_forecast {BEGIN DYNARE_STATEMENT; return token::BVAR_FORECAST; }
<INITIAL>dynare_sensitivity {BEGIN DYNARE_STATEMENT; return token::DYNARE_SENSITIVITY;}
<INITIAL>initval_file {BEGIN DYNARE_STATEMENT; return token::INITVAL_FILE;}
<INITIAL>histval_file {BEGIN DYNARE_STATEMENT; return token::HISTVAL_FILE;}
<INITIAL>forecast {BEGIN DYNARE_STATEMENT; return token::FORECAST;}
<INITIAL>shock_decomposition {BEGIN DYNARE_STATEMENT; return token::SHOCK_DECOMPOSITION;}
<INITIAL>realtime_shock_decomposition {BEGIN DYNARE_STATEMENT; return token::REALTIME_SHOCK_DECOMPOSITION;}
<INITIAL>plot_shock_decomposition {BEGIN DYNARE_STATEMENT; return token::PLOT_SHOCK_DECOMPOSITION;}
<INITIAL>initial_condition_decomposition {BEGIN DYNARE_STATEMENT; return token::INITIAL_CONDITION_DECOMPOSITION;}
<INITIAL>squeeze_shock_decomposition {BEGIN DYNARE_STATEMENT; return token::SQUEEZE_SHOCK_DECOMPOSITION;}
<INITIAL>sbvar {BEGIN DYNARE_STATEMENT; return token::SBVAR;}
<INITIAL>ms_estimation {BEGIN DYNARE_STATEMENT; return token::MS_ESTIMATION;}
<INITIAL>ms_simulation {BEGIN DYNARE_STATEMENT; return token::MS_SIMULATION;}
<INITIAL>ms_compute_mdd {BEGIN DYNARE_STATEMENT; return token::MS_COMPUTE_MDD;}
<INITIAL>ms_compute_probabilities {BEGIN DYNARE_STATEMENT; return token::MS_COMPUTE_PROBABILITIES;}
<INITIAL>ms_forecast {BEGIN DYNARE_STATEMENT; return token::MS_FORECAST;}
<INITIAL>ms_irf {BEGIN DYNARE_STATEMENT; return token::MS_IRF;}
<INITIAL>ms_variance_decomposition {BEGIN DYNARE_STATEMENT; return token::MS_VARIANCE_DECOMPOSITION;}
<INITIAL>conditional_forecast {BEGIN DYNARE_STATEMENT; return token::CONDITIONAL_FORECAST;}
<INITIAL>plot_conditional_forecast {BEGIN DYNARE_STATEMENT; return token::PLOT_CONDITIONAL_FORECAST;}
<INITIAL>gmm_estimation {BEGIN DYNARE_STATEMENT; return token::GMM_ESTIMATION;}
<INITIAL>smm_estimation {BEGIN DYNARE_STATEMENT; return token::SMM_ESTIMATION;}

<INITIAL>markov_switching {BEGIN DYNARE_STATEMENT; return token::MARKOV_SWITCHING;}
<INITIAL>svar {BEGIN DYNARE_STATEMENT; return token::SVAR;}
<INITIAL>svar_global_identification_check {BEGIN DYNARE_STATEMENT; return token::SVAR_GLOBAL_IDENTIFICATION_CHECK;}
<INITIAL>external_function {BEGIN DYNARE_STATEMENT; return token::EXTERNAL_FUNCTION;}
 /* End of a Dynare statement */
<INITIAL>calib_smoother { BEGIN DYNARE_STATEMENT; return token::CALIB_SMOOTHER; } 
<INITIAL>model_diagnostics {BEGIN DYNARE_STATEMENT; return token::MODEL_DIAGNOSTICS;}
<INITIAL>extended_path {BEGIN DYNARE_STATEMENT; return token::EXTENDED_PATH;}
<INITIAL>smoother2histval {BEGIN DYNARE_STATEMENT; return token::SMOOTHER2HISTVAL;}
<INITIAL>perfect_foresight_setup {BEGIN DYNARE_STATEMENT; return token::PERFECT_FORESIGHT_SETUP;}
<INITIAL>perfect_foresight_solver {BEGIN DYNARE_STATEMENT; return token::PERFECT_FORESIGHT_SOLVER;}
<INITIAL>det_cond_forecast {BEGIN DYNARE_STATEMENT; return token::DET_COND_FORECAST;}
<INITIAL>compilation_setup {BEGIN DYNARE_STATEMENT; return token::COMPILATION_SETUP;}

<DYNARE_STATEMENT>; {
  if (!sigma_e)
    BEGIN INITIAL;
  return Dynare::parser::token_type (yytext[0]);
}


 /* Begin of a Dynare block */
<INITIAL>model {BEGIN DYNARE_BLOCK; return token::MODEL;}
<INITIAL>steady_state_model {BEGIN DYNARE_BLOCK; return token::STEADY_STATE_MODEL;}
<INITIAL>initval {BEGIN DYNARE_BLOCK; return token::INITVAL;}
<INITIAL>endval {BEGIN DYNARE_BLOCK; return token::ENDVAL;}
<INITIAL>histval {BEGIN DYNARE_BLOCK; return token::HISTVAL;}
<INITIAL>shocks {BEGIN DYNARE_BLOCK; return token::SHOCKS;}
<INITIAL>shock_groups {BEGIN DYNARE_BLOCK; return token::SHOCK_GROUPS;}
<INITIAL>init2shocks {BEGIN DYNARE_BLOCK; return token::INIT2SHOCKS;}
<INITIAL>mshocks {BEGIN DYNARE_BLOCK; return token::MSHOCKS;}
<INITIAL>estimated_params {BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS;}
<INITIAL>epilogue {BEGIN DYNARE_BLOCK; return token::EPILOGUE;}
 /* priors is an alias for estimated_params */
<INITIAL>priors {BEGIN DYNARE_BLOCK;return token::ESTIMATED_PARAMS;}
<INITIAL>estimated_params_init 		{BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS_INIT;}
<INITIAL>estimated_params_bounds 	{BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS_BOUNDS;}
<INITIAL>osr_params_bounds              {BEGIN DYNARE_BLOCK; return token::OSR_PARAMS_BOUNDS;}
<INITIAL>observation_trends {BEGIN DYNARE_BLOCK; return token::OBSERVATION_TRENDS;}
<INITIAL>optim_weights {BEGIN DYNARE_BLOCK; return token::OPTIM_WEIGHTS;}
<INITIAL>homotopy_setup {BEGIN DYNARE_BLOCK; return token::HOMOTOPY_SETUP;}
<INITIAL>conditional_forecast_paths {BEGIN DYNARE_BLOCK; return token::CONDITIONAL_FORECAST_PATHS;}
<INITIAL>svar_identification {BEGIN DYNARE_BLOCK; return token::SVAR_IDENTIFICATION;}
<INITIAL>moment_calibration {BEGIN DYNARE_BLOCK; return token::MOMENT_CALIBRATION;}
<INITIAL>irf_calibration {BEGIN DYNARE_BLOCK; return token::IRF_CALIBRATION;}
<INITIAL>ramsey_constraints {BEGIN DYNARE_BLOCK; return token::RAMSEY_CONSTRAINTS;}
<INITIAL>restrictions {BEGIN DYNARE_BLOCK; return token::RESTRICTIONS;}
<INITIAL>generate_irfs {BEGIN DYNARE_BLOCK; return token::GENERATE_IRFS;}

 /* For the semicolon after an "end" keyword */
<INITIAL>; {return Dynare::parser::token_type (yytext[0]);}

 /* End of a Dynare block */
<DYNARE_BLOCK>end 	{BEGIN INITIAL; return token::END;}

<DYNARE_STATEMENT>subsamples {return token::SUBSAMPLES;}
<DYNARE_STATEMENT>options {return token::OPTIONS;}
<DYNARE_STATEMENT>prior {
  yylval->build<string>(yytext);
  return token::PRIOR;
}
<INITIAL>std {BEGIN DYNARE_STATEMENT; return token::STD;}
<INITIAL>corr {BEGIN DYNARE_STATEMENT; return token::CORR;}
<DYNARE_STATEMENT>function {return token::FUNCTION;}
<DYNARE_STATEMENT>sampling_draws {return token::SAMPLING_DRAWS;}
<INITIAL>prior_function {BEGIN DYNARE_STATEMENT; return token::PRIOR_FUNCTION;}
<INITIAL>posterior_function {BEGIN DYNARE_STATEMENT; return token::POSTERIOR_FUNCTION;}

 /* Inside  of a Dynare statement */
<DYNARE_STATEMENT>{DATE} {
                           char *yycopy = strdup(yytext);
                           char *uput = yycopy + yyleng;
                           unput(')');
                           unput('\'');
                           while (uput > yycopy)
                             unput(*--uput);
                           unput('\'');
                           unput('(');
                           unput('s');
                           unput('e');
                           unput('t');
                           unput('a');
                           unput('d');
                           free( yycopy );
                         }
<DYNARE_STATEMENT>${DATE} { yylloc->step();
#if (YY_FLEX_MAJOR_VERSION > 2) || (YY_FLEX_MAJOR_VERSION == 2 && YY_FLEX_MINOR_VERSION >= 6)
                            yyout << yytext + 1;
#else
                            *yyout << yytext + 1;
#endif
                          }
<DYNARE_STATEMENT>dates  {dates_parens_nb=0; BEGIN DATES_STATEMENT; yylval->build<string>("dates");}
<DYNARE_STATEMENT>file                  {return token::FILE;}
<DYNARE_STATEMENT>datafile 		{return token::DATAFILE;}
<DYNARE_STATEMENT>dirname       {return token::DIRNAME;}
<DYNARE_STATEMENT>nobs 			{return token::NOBS;}
<DYNARE_STATEMENT>last_obs 		{return token::LAST_OBS;}
<DYNARE_STATEMENT>first_obs 		{return token::FIRST_OBS;}
<DYNARE_STATEMENT>mean                  {return token::MEAN;}
<DYNARE_STATEMENT>stdev                 {return token::STDEV;}
<DYNARE_STATEMENT>truncate              {return token::TRUNCATE;}
<DYNARE_STATEMENT>domain                {return token::DOMAINN;}
<DYNARE_STATEMENT>variance              {return token::VARIANCE;}
<DYNARE_STATEMENT>mode                  {return token::MODE;}
<DYNARE_STATEMENT>interval              {return token::INTERVAL;}
<DYNARE_STATEMENT>shape                 {return token::SHAPE;}
<DYNARE_STATEMENT>shift                 {return token::SHIFT;}
<DYNARE_STATEMENT>bounds                {return token::BOUNDS;}
<DYNARE_STATEMENT>init                  {return token::INIT;}
<DYNARE_STATEMENT>jscale                {return token::JSCALE;}
<DYNARE_STATEMENT>prefilter 		{return token::PREFILTER;}
<DYNARE_STATEMENT>presample 		{return token::PRESAMPLE;}
<DYNARE_STATEMENT>lik_algo  		{return token::LIK_ALGO;}
<DYNARE_STATEMENT>lik_init  		{return token::LIK_INIT;}
<DYNARE_STATEMENT>taper_steps       {return token::TAPER_STEPS;}
<DYNARE_STATEMENT>geweke_interval   {return token::GEWEKE_INTERVAL;}
<DYNARE_STATEMENT>raftery_lewis_qrs {return token::RAFTERY_LEWIS_QRS;}
<DYNARE_STATEMENT>raftery_lewis_diagnostics {return token::RAFTERY_LEWIS_DIAGNOSTICS;}
<DYNARE_STATEMENT>graph   		{return token::GRAPH;}
<DYNARE_STATEMENT>nograph   		{return token::NOGRAPH;}
<DYNARE_STATEMENT>posterior_graph   		{return token::POSTERIOR_GRAPH;}
<DYNARE_STATEMENT>posterior_nograph   		{return token::POSTERIOR_NOGRAPH;}
<DYNARE_STATEMENT>nodisplay     {return token::NODISPLAY;}
<DYNARE_STATEMENT>graph_format  {return token::GRAPH_FORMAT;}
<DYNARE_STATEMENT>eps  {yylval->build<string>(yytext); return token::EPS;}
<DYNARE_STATEMENT>pdf  {yylval->build<string>(yytext); return token::PDF;}
<DYNARE_STATEMENT>fig  {yylval->build<string>(yytext); return token::FIG;}
<DYNARE_STATEMENT>none  {yylval->build<string>(yytext); return token::NONE;}
<DYNARE_STATEMENT>print   		{return token::PRINT;}
<DYNARE_STATEMENT>noprint   		{return token::NOPRINT;}
<DYNARE_STATEMENT>conf_sig  		{return token::CONF_SIG;}
<DYNARE_STATEMENT>mh_conf_sig  		{return token::MH_CONF_SIG;}
<DYNARE_STATEMENT>mh_replic 		{return token::MH_REPLIC;}
<DYNARE_STATEMENT>mh_drop   		{return token::MH_DROP;}
<DYNARE_STATEMENT>mh_jscale   		{return token::MH_JSCALE;}
<DYNARE_STATEMENT>mh_init_scale 	{return token::MH_INIT_SCALE;}
<DYNARE_STATEMENT>mh_tune_jscale   	{return token::MH_TUNE_JSCALE;}
<DYNARE_STATEMENT>mode_file 		{return token::MODE_FILE;}
<DYNARE_STATEMENT>mode_compute 	{return token::MODE_COMPUTE;}
<DYNARE_STATEMENT>mode_check 		{return token::MODE_CHECK;}
<DYNARE_STATEMENT>mode_check_neighbourhood_size 		{return token::MODE_CHECK_NEIGHBOURHOOD_SIZE;}
<DYNARE_STATEMENT>mode_check_symmetric_plots 		{return token::MODE_CHECK_SYMMETRIC_PLOTS;}
<DYNARE_STATEMENT>mode_check_number_of_points 		{return token::MODE_CHECK_NUMBER_OF_POINTS;}
<DYNARE_STATEMENT>prior_trunc 	{return token::PRIOR_TRUNC;}
<DYNARE_STATEMENT>mh_mode 		{return token::MH_MODE;}
<DYNARE_STATEMENT>mh_nblocks 		{return token::MH_NBLOCKS;}
<DYNARE_STATEMENT>load_mh_file 	{return token::LOAD_MH_FILE;}
<DYNARE_STATEMENT>load_results_after_load_mh 	{return token::LOAD_RESULTS_AFTER_LOAD_MH;}
<DYNARE_STATEMENT>loglinear 		{return token::LOGLINEAR;}
<DYNARE_STATEMENT>linear_approximation 		{return token::LINEAR_APPROXIMATION;}
<DYNARE_STATEMENT>logdata 	{return token::LOGDATA;}
<DYNARE_STATEMENT>nodiagnostic 	{return token::NODIAGNOSTIC;}
<DYNARE_STATEMENT>kalman_algo 	{return token::KALMAN_ALGO;}
<DYNARE_STATEMENT>fast_kalman_filter {return token::FAST_KALMAN_FILTER;}
<DYNARE_STATEMENT>kalman_tol 	{return token::KALMAN_TOL;}
<DYNARE_STATEMENT>diffuse_kalman_tol 	{return token::DIFFUSE_KALMAN_TOL;}
<DYNARE_STATEMENT>forecast 	{return token::FORECAST;}
<DYNARE_STATEMENT>smoother 	{return token::SMOOTHER;}
<DYNARE_STATEMENT>bayesian_irf 	{return token::BAYESIAN_IRF;}
<DYNARE_STATEMENT>dsge_var 	{return token::DSGE_VAR;}
<DYNARE_STATEMENT>dsge_varlag 	{return token::DSGE_VARLAG;}
<DYNARE_STATEMENT>moments_varendo {return token::MOMENTS_VARENDO;}
<DYNARE_STATEMENT>contemporaneous_correlation {return token::CONTEMPORANEOUS_CORRELATION;}
<DYNARE_STATEMENT>posterior_max_subsample_draws	{return token::POSTERIOR_MAX_SUBSAMPLE_DRAWS;}
<DYNARE_STATEMENT>filtered_vars	{return token::FILTERED_VARS;}
<DYNARE_STATEMENT>filter_step_ahead	{return token::FILTER_STEP_AHEAD;}
<DYNARE_STATEMENT,DYNARE_BLOCK>relative_irf {return token::RELATIVE_IRF;}
<DYNARE_STATEMENT>tex		{return token::TEX;}
<DYNARE_STATEMENT>nomoments	{return token::NOMOMENTS;}
<DYNARE_STATEMENT>std		{return token::STD;}
<DYNARE_STATEMENT>corr		{return token::CORR;}
<DYNARE_STATEMENT>nocorr	{return token::NOCORR;}
<DYNARE_STATEMENT>optim		{return token::OPTIM;}
<DYNARE_STATEMENT>periods	{return token::PERIODS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>model_name	{return token::MODEL_NAME;}
<DYNARE_STATEMENT>auxiliary_model_name    {return token::AUXILIARY_MODEL_NAME;}
<DYNARE_STATEMENT>endogenous_terminal_period 	{return token::ENDOGENOUS_TERMINAL_PERIOD;}
<DYNARE_STATEMENT>sub_draws	{return token::SUB_DRAWS;}
<DYNARE_STATEMENT>minimal_solving_periods {return token::MINIMAL_SOLVING_PERIODS;}
<DYNARE_STATEMENT>markowitz	{return token::MARKOWITZ;}
<DYNARE_STATEMENT>marginal_density {return token::MARGINAL_DENSITY;}
<DYNARE_STATEMENT>laplace       {return token::LAPLACE;}
<DYNARE_STATEMENT>modifiedharmonicmean {return token::MODIFIEDHARMONICMEAN;}
<DYNARE_STATEMENT>constant	{return token::CONSTANT;}
<DYNARE_STATEMENT>noconstant	{return token::NOCONSTANT;}
<DYNARE_STATEMENT>filename      {return token::FILENAME;}
<DYNARE_STATEMENT>diffuse_filter {return token::DIFFUSE_FILTER;}
<DYNARE_STATEMENT>plot_priors   {return token::PLOT_PRIORS;}
<DYNARE_STATEMENT>aim_solver {return token::AIM_SOLVER;}
<DYNARE_STATEMENT>partial_information {return token::PARTIAL_INFORMATION;}
<DYNARE_STATEMENT>conditional_variance_decomposition {return token::CONDITIONAL_VARIANCE_DECOMPOSITION;}
<DYNARE_STATEMENT>name {return token::EXT_FUNC_NAME;}
<DYNARE_STATEMENT>nargs {return token::EXT_FUNC_NARGS;}
<DYNARE_STATEMENT>first_deriv_provided {return token::FIRST_DERIV_PROVIDED;}
<DYNARE_STATEMENT>second_deriv_provided {return token::SECOND_DERIV_PROVIDED;}
<DYNARE_STATEMENT>freq {return token::FREQ;}
<DYNARE_STATEMENT>monthly {return token::MONTHLY; }
<DYNARE_STATEMENT>quarterly {return token::QUARTERLY; }
<DYNARE_STATEMENT>initial_year {return token::INITIAL_YEAR;}
<DYNARE_STATEMENT>initial_subperiod {return token::INITIAL_SUBPERIOD;}
<DYNARE_STATEMENT>final_year {return token::FINAL_YEAR;}
<DYNARE_STATEMENT>final_subperiod {return token::FINAL_SUBPERIOD;}
<DYNARE_STATEMENT>vlist {return token::VLIST;}
<DYNARE_STATEMENT>vlistlog {return token::VLISTLOG;}
<DYNARE_STATEMENT>vlistper {return token::VLISTPER;}
<DYNARE_STATEMENT>keep_kalman_algo_if_singularity_is_detected {return token::KEEP_KALMAN_ALGO_IF_SINGULARITY_IS_DETECTED;}
<DYNARE_STATEMENT>restriction_fname {return token::RESTRICTION_FNAME;}
<DYNARE_STATEMENT>nlags {return token::NLAGS;}
<DYNARE_STATEMENT>restrictions {return token::RESTRICTIONS;}
<DYNARE_BLOCK>crossequations {return token::CROSSEQUATIONS;}
<DYNARE_BLOCK>covariance {return token::COVARIANCE;}
<DYNARE_BLOCK>adl {return token::ADL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>diff {return token::DIFF;}
<DYNARE_STATEMENT>cross_restrictions {return token::CROSS_RESTRICTIONS;}
<DYNARE_STATEMENT>contemp_reduced_form {return token::CONTEMP_REDUCED_FORM;}
<DYNARE_STATEMENT>real_pseudo_forecast {return token::REAL_PSEUDO_FORECAST;}
<DYNARE_STATEMENT>no_bayesian_prior {return token::NO_BAYESIAN_PRIOR;}
<DYNARE_STATEMENT>dummy_obs {return token::DUMMY_OBS;}
<DYNARE_STATEMENT>spectral_density {return token::SPECTRAL_DENSITY;}
<DYNARE_STATEMENT>nstates {return token::NSTATES;}
<DYNARE_STATEMENT>indxscalesstates {return token::INDXSCALESSTATES;}
<DYNARE_STATEMENT>fixed_point {return token::FIXED_POINT;}
<DYNARE_STATEMENT>doubling {return token::DOUBLING;}
<DYNARE_STATEMENT>plot_init_date {return token::PLOT_INIT_DATE;}
<DYNARE_STATEMENT>plot_end_date {return token::PLOT_END_DATE;}
<DYNARE_STATEMENT>square_root_solver {return token::SQUARE_ROOT_SOLVER;}
<DYNARE_STATEMENT>cycle_reduction {return token::CYCLE_REDUCTION;}
<DYNARE_STATEMENT>logarithmic_reduction {return token::LOGARITHMIC_REDUCTION;}
<DYNARE_STATEMENT>use_univariate_filters_if_singularity_is_detected {return token::USE_UNIVARIATE_FILTERS_IF_SINGULARITY_IS_DETECTED;}
<DYNARE_STATEMENT>hybrid {return token::HYBRID;}
<DYNARE_STATEMENT>default {return token::DEFAULT;}
<DYNARE_STATEMENT>init2shocks {return token::INIT2SHOCKS;}

<DYNARE_STATEMENT>number_of_particles {return token::NUMBER_OF_PARTICLES;}
<DYNARE_STATEMENT>resampling {return token::RESAMPLING;}
<DYNARE_STATEMENT>systematic {return token::SYSTEMATIC;}
<DYNARE_STATEMENT>generic {return token::GENERIC;}
<DYNARE_STATEMENT>resampling_threshold {return token::RESAMPLING_THRESHOLD;}
<DYNARE_STATEMENT>resampling_method {return token::RESAMPLING_METHOD;}
<DYNARE_STATEMENT>kitagawa {return token::KITAGAWA;}
<DYNARE_STATEMENT>smooth {return token::SMOOTH;}
<DYNARE_STATEMENT>stratified {return token::STRATIFIED;}
<DYNARE_STATEMENT>cpf_weights {return token::CPF_WEIGHTS;}
<DYNARE_STATEMENT>amisanotristani {return token::AMISANOTRISTANI;}
<DYNARE_STATEMENT>murrayjonesparslow {return token::MURRAYJONESPARSLOW;}
<DYNARE_STATEMENT>filter_algorithm {return token::FILTER_ALGORITHM;}
<DYNARE_STATEMENT>nonlinear_filter_initialization {return token::NONLINEAR_FILTER_INITIALIZATION;}
<DYNARE_STATEMENT>proposal_approximation {return token::PROPOSAL_APPROXIMATION;}
<DYNARE_STATEMENT>cubature {return token::CUBATURE;}
<DYNARE_STATEMENT>unscented {return token::UNSCENTED;}
<DYNARE_STATEMENT>montecarlo {return token::MONTECARLO;}
<DYNARE_STATEMENT>distribution_approximation {return token::DISTRIBUTION_APPROXIMATION;}
<DYNARE_STATEMENT>proposal_distribution {return token::PROPOSAL_DISTRIBUTION;}
<DYNARE_STATEMENT>no_posterior_kernel_density {return token::NO_POSTERIOR_KERNEL_DENSITY;}
<DYNARE_STATEMENT>rescale_prediction_error_covariance {return token::RESCALE_PREDICTION_ERROR_COVARIANCE;}
<DYNARE_STATEMENT>use_penalized_objective_for_hessian {return token::USE_PENALIZED_OBJECTIVE_FOR_HESSIAN;}
<DYNARE_STATEMENT>expression    {return token::EXPRESSION;}

<DYNARE_STATEMENT>alpha {
  yylval->build<string>(yytext);
  return token::ALPHA;
}
<DYNARE_STATEMENT>beta {
  yylval->build<string>(yytext);
  return token::BETA;
}
<DYNARE_STATEMENT>gamma {
  yylval->build<string>(yytext);
  return token::GAMMA;
}
<DYNARE_STATEMENT>inv_gamma {
  yylval->build<string>(yytext);
  return token::INV_GAMMA;
}
<DYNARE_STATEMENT>inv_gamma1 {
  yylval->build<string>(yytext);
  return token::INV_GAMMA1;
}
<DYNARE_STATEMENT>inv_gamma2 {
  yylval->build<string>(yytext);
  return token::INV_GAMMA2;
}
<DYNARE_STATEMENT>dirichlet {
  yylval->build<string>(yytext);
  return token::DIRICHLET;
}
<DYNARE_STATEMENT>weibull {
  yylval->build<string>(yytext);
  return token::WEIBULL;
}
<DYNARE_STATEMENT>normal {
  yylval->build<string>(yytext);
  return token::NORMAL;
}
<DYNARE_STATEMENT>uniform {
  yylval->build<string>(yytext);
  return token::UNIFORM;
}
<DYNARE_STATEMENT>gsig2_lmdm {return token::GSIG2_LMDM;}
<DYNARE_STATEMENT>specification {return token::SPECIFICATION;}
<DYNARE_STATEMENT>sims_zha {return token::SIMS_ZHA;}
<DYNARE_STATEMENT>q_diag {return token::Q_DIAG;}
<DYNARE_STATEMENT>flat_prior {return token::FLAT_PRIOR;}
<DYNARE_STATEMENT>ncsk {return token::NCSK;}
<DYNARE_STATEMENT>nstd {return token::NSTD;}
<DYNARE_STATEMENT>ninv {
  yylval->build<string>(yytext);
  return token::NINV;
}
<DYNARE_STATEMENT>indxparr {return token::INDXPARR;}
<DYNARE_STATEMENT>indxovr {return token::INDXOVR;}
<DYNARE_STATEMENT>aband {
  yylval->build<string>(yytext);
  return token::ABAND;
}
<DYNARE_STATEMENT>write_equation_tags {return token::WRITE_EQUATION_TAGS;}
<DYNARE_STATEMENT>eqtags {return token::EQTAGS;}
<DYNARE_STATEMENT>targets {return token::TARGETS;}
<DYNARE_STATEMENT>indxap {return token::INDXAP;}
<DYNARE_STATEMENT>apband {return token::APBAND;}
<DYNARE_STATEMENT>indximf {return token::INDXIMF;}
<DYNARE_STATEMENT>indxfore {return token::INDXFORE;}
<DYNARE_STATEMENT>foreband {return token::FOREBAND;}
<DYNARE_STATEMENT>indxgforehat {return token::INDXGFOREHAT;}
<DYNARE_STATEMENT>indxgimfhat {return token::INDXGIMFHAT;}
<DYNARE_STATEMENT>indxestima {return token::INDXESTIMA;}
<DYNARE_STATEMENT>indxgdls {return token::INDXGDLS;}
<DYNARE_STATEMENT>eq_ms {return token::EQ_MS;}
<DYNARE_STATEMENT>cms {
  yylval->build<string>(yytext);
  return token::CMS;
}
<DYNARE_STATEMENT>ncms {
  yylval->build<string>(yytext);
  return token::NCMS;
}
<DYNARE_STATEMENT>eq_cms {return token::EQ_CMS;}
<DYNARE_STATEMENT>tlindx {return token::TLINDX;}
<DYNARE_STATEMENT>tlnumber {return token::TLNUMBER;}
<DYNARE_STATEMENT>cnum {
  yylval->build<string>(yytext);
  return token::CNUM;
}
<DYNARE_STATEMENT>nodecomposition {return token::NODECOMPOSITION;};
<DYNARE_BLOCK>use_calibration {return token::USE_CALIBRATION;}
<DYNARE_STATEMENT>output_file_tag {return token::OUTPUT_FILE_TAG;}
<DYNARE_STATEMENT>file_tag {return token::FILE_TAG;};
<DYNARE_STATEMENT>no_create_init {return token::NO_CREATE_INIT;};
<DYNARE_STATEMENT>simulation_file_tag {return token::SIMULATION_FILE_TAG;};
<DYNARE_STATEMENT>horizon {return token::HORIZON;}
<DYNARE_STATEMENT>parameter_uncertainty {return token::PARAMETER_UNCERTAINTY;}
<DYNARE_STATEMENT>no_error_bands {return token::NO_ERROR_BANDS;}
<DYNARE_STATEMENT>error_band_percentiles {return token::ERROR_BAND_PERCENTILES;}
<DYNARE_STATEMENT>shock_draws {return token::SHOCK_DRAWS;}
<DYNARE_STATEMENT>shocks_per_parameter {return token::SHOCKS_PER_PARAMETER;}
<DYNARE_STATEMENT>thinning_factor {return token::THINNING_FACTOR;}
<DYNARE_STATEMENT>free_parameters {return token::FREE_PARAMETERS;}
<DYNARE_STATEMENT>median {return token::MEDIAN;}
<DYNARE_STATEMENT>regime {return token::REGIME;}
<DYNARE_STATEMENT>regimes {return token::REGIMES;}
<DYNARE_STATEMENT>data_obs_nbr {return token::DATA_OBS_NBR;}
<DYNARE_STATEMENT>filtered_probabilities {return token::FILTERED_PROBABILITIES;}
<DYNARE_STATEMENT>real_time_smoothed {return token::REAL_TIME_SMOOTHED;}
<DYNARE_STATEMENT>proposal_type {return token::PROPOSAL_TYPE;}
<DYNARE_STATEMENT>proposal_lower_bound {return token::PROPOSAL_LOWER_BOUND;}
<DYNARE_STATEMENT>proposal_upper_bound {return token::PROPOSAL_UPPER_BOUND;}
<DYNARE_STATEMENT>proposal_draws {return token::PROPOSAL_DRAWS;}
<DYNARE_STATEMENT>use_mean_center {return token::USE_MEAN_CENTER;}
<DYNARE_STATEMENT>adaptive_mh_draws {return token::ADAPTIVE_MH_DRAWS;}
<DYNARE_STATEMENT>coefficients_prior_hyperparameters {return token::COEFFICIENTS_PRIOR_HYPERPARAMETERS;}
<DYNARE_STATEMENT>convergence_starting_value {return token::CONVERGENCE_STARTING_VALUE;}
<DYNARE_STATEMENT>convergence_ending_value {return token::CONVERGENCE_ENDING_VALUE;}
<DYNARE_STATEMENT>convergence_increment_value {return token::CONVERGENCE_INCREMENT_VALUE;}
<DYNARE_STATEMENT>max_iterations_starting_value {return token::MAX_ITERATIONS_STARTING_VALUE;}
<DYNARE_STATEMENT>max_iterations_increment_value {return token::MAX_ITERATIONS_INCREMENT_VALUE;}
<DYNARE_STATEMENT>max_block_iterations {return token::MAX_BLOCK_ITERATIONS;}
<DYNARE_STATEMENT>max_repeated_optimization_runs {return token::MAX_REPEATED_OPTIMIZATION_RUNS;}
<DYNARE_STATEMENT>maxit {return token::MAXIT;}
<DYNARE_STATEMENT>function_convergence_criterion {return token::FUNCTION_CONVERGENCE_CRITERION;}
<DYNARE_STATEMENT>parameter_convergence_criterion {return token::PARAMETER_CONVERGENCE_CRITERION;}
<DYNARE_STATEMENT>number_of_large_perturbations {return token::NUMBER_OF_LARGE_PERTURBATIONS;}
<DYNARE_STATEMENT>number_of_small_perturbations {return token::NUMBER_OF_SMALL_PERTURBATIONS;}
<DYNARE_STATEMENT>number_of_posterior_draws_after_perturbation {return token::NUMBER_OF_POSTERIOR_DRAWS_AFTER_PERTURBATION;}
<DYNARE_STATEMENT>max_number_of_stages {return token::MAX_NUMBER_OF_STAGES;}
<DYNARE_STATEMENT>random_function_convergence_criterion {return token::RANDOM_FUNCTION_CONVERGENCE_CRITERION;}
<DYNARE_STATEMENT>random_parameter_convergence_criterion {return token::RANDOM_PARAMETER_CONVERGENCE_CRITERION;}
<DYNARE_STATEMENT>tolf {return token::TOLF;}
<DYNARE_STATEMENT>tolx {return token::TOLX;}
<DYNARE_STATEMENT>opt_algo {return token::OPT_ALGO;}
<DYNARE_STATEMENT>add_flags {return token::ADD_FLAGS;}
<DYNARE_STATEMENT>substitute_flags {return token::SUBSTITUTE_FLAGS;}
<DYNARE_STATEMENT>add_libs {return token::ADD_LIBS;}
<DYNARE_STATEMENT>substitute_libs {return token::SUBSTITUTE_LIBS;}
<DYNARE_STATEMENT>compiler {return token::COMPILER;}
<DYNARE_STATEMENT>instruments {return token::INSTRUMENTS;}
<DYNARE_STATEMENT>hessian  {
  yylval->build<string>(yytext);
  return token::HESSIAN;
}
<DYNARE_STATEMENT>prior_variance  {
  yylval->build<string>(yytext);
  return token::PRIOR_VARIANCE;
}
<DYNARE_STATEMENT>identity_matrix  {
  yylval->build<string>(yytext);
  return token::IDENTITY_MATRIX;
}
<DYNARE_STATEMENT>mcmc_jumping_covariance {return token::MCMC_JUMPING_COVARIANCE;}

 /* These four (var, varexo, varexo_det, parameters) are for change_type */
<DYNARE_STATEMENT>var { return token::VAR; }
<DYNARE_STATEMENT>varexo { return token::VAREXO; }
<DYNARE_STATEMENT>varexo_det { return token::VAREXO_DET; }
<DYNARE_STATEMENT>parameters { return token::PARAMETERS; }
<DYNARE_STATEMENT>predetermined_variables { return token::PREDETERMINED_VARIABLES; }

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
<DYNARE_STATEMENT>homotopy_force_continue {return token::HOMOTOPY_FORCE_CONTINUE;}
<DYNARE_STATEMENT>nocheck {return token::NOCHECK; }

<DYNARE_STATEMENT>controlled_varexo {return token::CONTROLLED_VAREXO; }
<DYNARE_STATEMENT>parameter_set {return token::PARAMETER_SET; }
<DYNARE_STATEMENT>init_state {return token::INIT_STATE; }
<DYNARE_STATEMENT>fast_realtime {return token::FAST_REALTIME; }
<DYNARE_STATEMENT>save_realtime {return token::SAVE_REALTIME;}
<DYNARE_STATEMENT>detail_plot {return token::DETAIL_PLOT;}
<DYNARE_STATEMENT>flip {return token::FLIP;}
<DYNARE_STATEMENT>interactive {return token::INTERACTIVE;}
<DYNARE_STATEMENT>screen_shocks {return token::SCREEN_SHOCKS;}
<DYNARE_STATEMENT>steadystate {return token::STEADYSTATE;}
<DYNARE_STATEMENT>type {return token::TYPE;}
<DYNARE_STATEMENT>qoq {return token::QOQ; }
<DYNARE_STATEMENT>yoy {return token::YOY; }
<DYNARE_STATEMENT>aoa {return token::AOA; }
<DYNARE_STATEMENT>unconditional {return token::UNCONDITIONAL; }
<DYNARE_STATEMENT>conditional {return token::CONDITIONAL; }
<DYNARE_STATEMENT>fig_name {return token::FIG_NAME;}
<DYNARE_STATEMENT>write_xls {return token::WRITE_XLS;}
<DYNARE_STATEMENT>realtime {return token::REALTIME;}
<DYNARE_STATEMENT>vintage {return token::VINTAGE;}
<DYNARE_STATEMENT>prior_mode {return token::PRIOR_MODE; }
<DYNARE_STATEMENT>prior_mean {return token::PRIOR_MEAN; }
<DYNARE_STATEMENT>posterior_mode {return token::POSTERIOR_MODE; }
<DYNARE_STATEMENT>posterior_mean {return token::POSTERIOR_MEAN; }
<DYNARE_STATEMENT>posterior_median {return token::POSTERIOR_MEDIAN; }
<DYNARE_STATEMENT>mle_mode {return token::MLE_MODE; }
<DYNARE_STATEMENT>k_order_solver {return token::K_ORDER_SOLVER; }
<DYNARE_STATEMENT>filter_covariance {return token::FILTER_COVARIANCE; }
<DYNARE_STATEMENT>filter_decomposition {return token::FILTER_DECOMPOSITION; }
<DYNARE_STATEMENT>smoothed_state_uncertainty {return token::SMOOTHED_STATE_UNCERTAINTY; }
<DYNARE_STATEMENT>selected_variables_only {return token::SELECTED_VARIABLES_ONLY; }
<DYNARE_STATEMENT>pruning {return token::PRUNING; }
<DYNARE_STATEMENT>save_draws {return token::SAVE_DRAWS; }
<DYNARE_STATEMENT>deflator {return token::DEFLATOR;}
<DYNARE_STATEMENT>log_deflator {return token::LOG_DEFLATOR;}
<DYNARE_STATEMENT>epilogue {return token::EPILOGUE;}
<DYNARE_STATEMENT>growth_factor {return token::GROWTH_FACTOR;}
<DYNARE_STATEMENT>log_growth_factor {return token::LOG_GROWTH_FACTOR;}
<DYNARE_STATEMENT>growth {return token::GROWTH;}
<DYNARE_STATEMENT>cova_compute {return token::COVA_COMPUTE;}
<DYNARE_STATEMENT>discretionary_tol {return token::DISCRETIONARY_TOL;}
<DYNARE_STATEMENT>analytic_derivation {return token::ANALYTIC_DERIVATION;}
<DYNARE_STATEMENT>analytic_derivation_mode {return token::ANALYTIC_DERIVATION_MODE;}
<DYNARE_STATEMENT>solver_periods {return token::SOLVER_PERIODS;}
<DYNARE_STATEMENT>endogenous_prior {return token::ENDOGENOUS_PRIOR;}
<DYNARE_STATEMENT>consider_all_endogenous {return token::CONSIDER_ALL_ENDOGENOUS;}
<DYNARE_STATEMENT>consider_only_observed {return token::CONSIDER_ONLY_OBSERVED;}
<DYNARE_STATEMENT>infile {return token::INFILE;}
<DYNARE_STATEMENT>invars {return token::INVARS;}
<DYNARE_STATEMENT>period {return token::PERIOD;}
<DYNARE_STATEMENT>outfile {return token::OUTFILE;}
<DYNARE_STATEMENT>outvars {return token::OUTVARS;}
<DYNARE_STATEMENT>huge_number {return token::HUGE_NUMBER;}
<DYNARE_STATEMENT>dr_display_tol {return token::DR_DISPLAY_TOL;}
<DYNARE_STATEMENT>posterior_sampling_method {return token::POSTERIOR_SAMPLING_METHOD;}
<DYNARE_STATEMENT>posterior_sampler_options {return token::POSTERIOR_SAMPLER_OPTIONS;}
<DYNARE_STATEMENT>silent_optimizer {return token::SILENT_OPTIMIZER;}
<DYNARE_STATEMENT>lmmcp {return token::LMMCP;}
<DYNARE_STATEMENT>occbin {return token::OCCBIN;}
<DYNARE_STATEMENT>centered_moments {return token::CENTERED_MOMENTS; }
<DYNARE_STATEMENT>autolag {return token::AUTOLAG; }
<DYNARE_STATEMENT>recursive_order_estimation {return token::RECURSIVE_ORDER_ESTIMATION; }
<DYNARE_STATEMENT>bartlett_kernel_lag {return token::BARTLETT_KERNEL_LAG; }
<DYNARE_STATEMENT>optimal {
  yylval->build<string>(yytext);
  return token::OPTIMAL;
}
<DYNARE_STATEMENT>diagonal  {
  yylval->build<string>(yytext);
  return token::DIAGONAL;
}
<DYNARE_STATEMENT>weighting_matrix {return token::WEIGHTING_MATRIX; }
<DYNARE_STATEMENT>penalized_estimator {return token::PENALIZED_ESTIMATOR; }
<DYNARE_STATEMENT>verbose {return token::VERBOSE; }
<DYNARE_STATEMENT>simulation_multiple {return token::SIMULATION_MULTIPLE; }
<DYNARE_STATEMENT>seed {return token::SEED; }
<DYNARE_STATEMENT>bounded_shock_support {return token::BOUNDED_SHOCK_SUPPORT; }
<DYNARE_STATEMENT>analytical_girf {return token::ANALYTICAL_GIRF; }
<DYNARE_STATEMENT>irf_in_percent {return token::IRF_IN_PERCENT; }
<DYNARE_STATEMENT>emas_girf {return token::EMAS_GIRF; }
<DYNARE_STATEMENT>emas_drop {return token::EMAS_DROP; }
<DYNARE_STATEMENT>emas_tolf {return token::EMAS_TOLF; }
<DYNARE_STATEMENT>emas_max_iter {return token::EMAS_MAX_ITER; }
<DYNARE_STATEMENT>variable {return token::VARIABLE;}
<DYNARE_STATEMENT>no_identification_strength {return token::NO_IDENTIFICATION_STRENGTH;}
<DYNARE_STATEMENT>no_identification_reducedform {return token::NO_IDENTIFICATION_REDUCEDFORM;}
<DYNARE_STATEMENT>no_identification_moments {return token::NO_IDENTIFICATION_MOMENTS;}
<DYNARE_STATEMENT>no_identification_minimal {return token::NO_IDENTIFICATION_MINIMAL;}
<DYNARE_STATEMENT>no_identification_spectrum {return token::NO_IDENTIFICATION_SPECTRUM;}
<DYNARE_STATEMENT>normalize_jacobians {return token::NORMALIZE_JACOBIANS;}
<DYNARE_STATEMENT>grid_nbr {return token::GRID_NBR;}
<DYNARE_STATEMENT>tol_rank {return token::TOL_RANK;}
<DYNARE_STATEMENT>tol_deriv {return token::TOL_DERIV;}
<DYNARE_STATEMENT>tol_sv {return token::TOL_SV;}
<DYNARE_STATEMENT>checks_via_subsets {return token::CHECKS_VIA_SUBSETS;}
<DYNARE_STATEMENT>max_dim_subsets_groups {return token::MAX_DIM_SUBSETS_GROUPS;}
<DYNARE_STATEMENT>max_nrows {return token::MAX_NROWS;}
<DYNARE_STATEMENT>with_epilogue {return token::WITH_EPILOGUE;}

<DYNARE_STATEMENT>\$[^$]*\$ {
  strtok(yytext+1, "$");
  yylval->build<string>(yytext + 1);
  return token::TEX_NAME;
}

 /* Inside a Dynare block */
<DYNARE_BLOCK>var {return token::VAR;}
<DYNARE_BLOCK>stderr {return token::STDERR;}
<DYNARE_BLOCK>values {return token::VALUES;}
<DYNARE_BLOCK>corr {return token::CORR;}
<DYNARE_BLOCK>periods {return token::PERIODS;}
<DYNARE_BLOCK>cutoff {return token::CUTOFF;}
<DYNARE_BLOCK>mfs	{return token::MFS;}
<DYNARE_BLOCK>balanced_growth_test_tol {return token::BALANCED_GROWTH_TEST_TOL;}
<DYNARE_BLOCK>gamma_pdf {return token::GAMMA_PDF;}
<DYNARE_BLOCK>beta_pdf {return token::BETA_PDF;}
<DYNARE_BLOCK>normal_pdf {return token::NORMAL_PDF;}
<DYNARE_BLOCK>inv_gamma_pdf {return token::INV_GAMMA_PDF;}
<DYNARE_BLOCK>inv_gamma1_pdf {return token::INV_GAMMA1_PDF;}
<DYNARE_BLOCK>inv_gamma2_pdf {return token::INV_GAMMA2_PDF;}
<DYNARE_BLOCK>uniform_pdf {return token::UNIFORM_PDF;}
<DYNARE_BLOCK>weibull_pdf {return token::WEIBULL_PDF;}
<DYNARE_BLOCK>dsge_prior_weight {return token::DSGE_PRIOR_WEIGHT;}

<DYNARE_BLOCK>; {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_BLOCK># {return Dynare::parser::token_type (yytext[0]);}

<DYNARE_BLOCK>restriction {return token::RESTRICTION;}

 /* Inside Dynare statement */
<DYNARE_STATEMENT>solve_algo {return token::SOLVE_ALGO;}
<DYNARE_STATEMENT>dr_algo {return token::DR_ALGO;}
<DYNARE_STATEMENT>simul_algo {return token::SIMUL_ALGO;}
<DYNARE_STATEMENT>stack_solve_algo {return token::STACK_SOLVE_ALGO;}
<DYNARE_STATEMENT>robust_lin_solve {return token::ROBUST_LIN_SOLVE;}
<DYNARE_STATEMENT>drop {return token::DROP;}
<DYNARE_STATEMENT>order {return token::ORDER;}
<DYNARE_STATEMENT>sylvester {return token::SYLVESTER;}
<DYNARE_STATEMENT>lyapunov {return token::LYAPUNOV;}
<DYNARE_STATEMENT>dr {
  yylval->build<string>(yytext);
  return token::DR;
 }
<DYNARE_STATEMENT>sylvester_fixed_point_tol {return token::SYLVESTER_FIXED_POINT_TOL;}
<DYNARE_STATEMENT>lyapunov_fixed_point_tol {return token::LYAPUNOV_FIXED_POINT_TOL;}
<DYNARE_STATEMENT>lyapunov_doubling_tol {return token::LYAPUNOV_DOUBLING_TOL;}
<DYNARE_STATEMENT>dr_cycle_reduction_tol {return token::DR_CYCLE_REDUCTION_TOL;}
<DYNARE_STATEMENT>dr_logarithmic_reduction_tol {return token::DR_LOGARITHMIC_REDUCTION_TOL;}
<DYNARE_STATEMENT>dr_logarithmic_reduction_maxiter {return token::DR_LOGARITHMIC_REDUCTION_MAXITER;}
<DYNARE_STATEMENT>replic {return token::REPLIC;}
<DYNARE_STATEMENT>ar {return token::AR;}
<DYNARE_STATEMENT>nofunctions {return token::NOFUNCTIONS;}
<DYNARE_STATEMENT>irf {return token::IRF;}
<DYNARE_STATEMENT>irf_shocks {return token::IRF_SHOCKS;}
<DYNARE_STATEMENT>hp_filter {return token::HP_FILTER;}
<DYNARE_STATEMENT>one_sided_hp_filter {return token::ONE_SIDED_HP_FILTER;}
<DYNARE_STATEMENT>bandpass_filter {return token::BANDPASS_FILTER;}
<DYNARE_STATEMENT>hp_ngrid {return token::HP_NGRID;}
<DYNARE_STATEMENT>filtered_theoretical_moments_grid {return token::FILTERED_THEORETICAL_MOMENTS_GRID;}
<DYNARE_STATEMENT>simul_seed {return token::SIMUL_SEED;}
<DYNARE_STATEMENT>qz_criterium {return token::QZ_CRITERIUM;}
<DYNARE_STATEMENT>qz_zero_threshold {return token::QZ_ZERO_THRESHOLD;}
<DYNARE_STATEMENT>simul {return token::SIMUL;}
<DYNARE_STATEMENT>simul_replic {return token::SIMUL_REPLIC;}
<DYNARE_STATEMENT>xls_sheet {return token::XLS_SHEET;}
<DYNARE_STATEMENT>xls_range {return token::XLS_RANGE;}
<DYNARE_STATEMENT>series {return token::SERIES;}
<DYNARE_STATEMENT>mh_recover {return token::MH_RECOVER;}
<DYNARE_STATEMENT>planner_discount {return token::PLANNER_DISCOUNT;}
<DYNARE_STATEMENT>planner_discount_latex_name {return token::PLANNER_DISCOUNT_LATEX_NAME;}
<DYNARE_STATEMENT>calibration {return token::CALIBRATION;}
<DYNARE_STATEMENT>irf_plot_threshold {return token::IRF_PLOT_THRESHOLD;}
<DYNARE_STATEMENT>no_homotopy {return token::NO_HOMOTOPY;}

<DYNARE_BLOCK>stderr_multiples {return token::STDERR_MULTIPLES;}
<DYNARE_BLOCK>diagonal_only {return token::DIAGONAL_ONLY;}
<DYNARE_BLOCK>equation {return token::EQUATION;}
<DYNARE_BLOCK>exclusion {return token::EXCLUSION;}
<DYNARE_BLOCK>lag {return token::LAG;}
<DYNARE_BLOCK>coeff {return token::COEFF;}
<DYNARE_BLOCK>overwrite {return token::OVERWRITE;}
<DYNARE_STATEMENT,DYNARE_BLOCK>upper_cholesky {return token::UPPER_CHOLESKY;}
<DYNARE_STATEMENT,DYNARE_BLOCK>lower_cholesky {return token::LOWER_CHOLESKY;}
<DYNARE_STATEMENT>chain {return token::CHAIN;}
<DYNARE_STATEMENT>number_of_lags {return token::NUMBER_OF_LAGS;}
<DYNARE_STATEMENT>number_of_regimes {return token::NUMBER_OF_REGIMES;}
<DYNARE_STATEMENT>duration {return token::DURATION;}
<DYNARE_STATEMENT>coefficients {return token::COEFFICIENTS;}
<DYNARE_STATEMENT>variances {return token::VARIANCES;}
<DYNARE_STATEMENT>equations {return token::EQUATIONS;}

<DYNARE_STATEMENT,DYNARE_BLOCK>\. {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT>\\ {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT>\' {return Dynare::parser::token_type (yytext[0]);}

<DYNARE_BLOCK>use_dll {return token::USE_DLL;}
<DYNARE_BLOCK>block {return token::BLOCK;}
<DYNARE_BLOCK>bytecode {return token::BYTECODE;}
<DYNARE_BLOCK>linear_decomposition {return token::LINEAR_DECOMPOSITION;}
<DYNARE_BLOCK>all_values_required {return token::ALL_VALUES_REQUIRED;}
<DYNARE_BLOCK>no_static {return token::NO_STATIC;}
<DYNARE_BLOCK>differentiate_forward_vars {return token::DIFFERENTIATE_FORWARD_VARS;}
<DYNARE_BLOCK>parallel_local_files {return token::PARALLEL_LOCAL_FILES;}

<DYNARE_STATEMENT,DYNARE_BLOCK>linear {return token::LINEAR;}

<DYNARE_STATEMENT,DYNARE_BLOCK>, {return token::COMMA;}
<DYNARE_STATEMENT,DYNARE_BLOCK>: {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\(\)] {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT,DYNARE_BLOCK>\[ {return Dynare::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT,DYNARE_BLOCK>\] {
  if (sigma_e)
    sigma_e = 0;
  return Dynare::parser::token_type (yytext[0]);
}
<DYNARE_STATEMENT,DYNARE_BLOCK>\+ {return token::PLUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>-  {return token::MINUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>\* {return token::TIMES;}
<DYNARE_STATEMENT,DYNARE_BLOCK>\/ {return token::DIVIDE;}
<DYNARE_STATEMENT,DYNARE_BLOCK>= {return token::EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>< {return token::LESS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>> {return token::GREATER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>>= {return token::GREATER_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK><= {return token::LESS_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>== {return token::EQUAL_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>!= {return token::EXCLAMATION_EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>\^ {return token::POWER;}
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
<DYNARE_STATEMENT,DYNARE_BLOCK>sqrt {return token::SQRT;}
<DYNARE_STATEMENT,DYNARE_BLOCK>cbrt {return token::CBRT;}
<DYNARE_STATEMENT,DYNARE_BLOCK>max {return token::MAX;}
<DYNARE_STATEMENT,DYNARE_BLOCK>min {return token::MIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>abs {return token::ABS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sign {return token::SIGN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>normcdf {return token::NORMCDF;}
<DYNARE_STATEMENT,DYNARE_BLOCK>normpdf {return token::NORMPDF;}
<DYNARE_STATEMENT,DYNARE_BLOCK>erf {return token::ERF;}
<DYNARE_STATEMENT,DYNARE_BLOCK>steady_state {return token::STEADY_STATE;}
<DYNARE_STATEMENT,DYNARE_BLOCK>expectation {return token::EXPECTATION;}
<DYNARE_BLOCK>var_expectation {return token::VAR_EXPECTATION;}
<DYNARE_BLOCK>pac_expectation {return token::PAC_EXPECTATION;}
<DYNARE_STATEMENT>discount {return token::DISCOUNT;}
<DYNARE_STATEMENT>steady_state_growth {return token::STEADY_STATE_GROWTH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>varobs {return token::VAROBS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>varexobs {return token::VAREXOBS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>nan {return token::NAN_CONSTANT;}
<DYNARE_STATEMENT,DYNARE_BLOCK>inf {return token::INF_CONSTANT;}
<DYNARE_STATEMENT,DYNARE_BLOCK>constants {return token::CONSTANTS;}

 /* options for GSA module by Marco Ratto */
<DYNARE_STATEMENT>identification {return token::IDENTIFICATION;}
<DYNARE_STATEMENT>morris {return token::MORRIS;}
<DYNARE_STATEMENT>stab {return token::STAB;}
<DYNARE_STATEMENT>redform {return token::REDFORM;}
<DYNARE_STATEMENT>pprior {return token::PPRIOR;}
<DYNARE_STATEMENT>prior_range {return token::PRIOR_RANGE;}
<DYNARE_STATEMENT>ppost {return token::PPOST;}
<DYNARE_STATEMENT>ilptau {return token::ILPTAU;}
<DYNARE_STATEMENT>morris_nliv {return token::MORRIS_NLIV;}
<DYNARE_STATEMENT>morris_ntra {return token::MORRIS_NTRA;}
<DYNARE_STATEMENT>Nsam {return token::NSAM;}
<DYNARE_STATEMENT>load_redform {return token::LOAD_REDFORM;}
<DYNARE_STATEMENT>load_rmse {return token::LOAD_RMSE;}
<DYNARE_STATEMENT>load_stab {return token::LOAD_STAB;}
<DYNARE_STATEMENT>alpha2_stab {return token::ALPHA2_STAB;}
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
<DYNARE_STATEMENT>load_ident_files {return token::LOAD_IDENT_FILES;}
<DYNARE_STATEMENT>useautocorr {return token::USEAUTOCORR;}
<DYNARE_STATEMENT>neighborhood_width {return token::NEIGHBORHOOD_WIDTH;}
<DYNARE_STATEMENT>pvalue_ks {return token::PVALUE_KS;}
<DYNARE_STATEMENT>pvalue_corr {return token::PVALUE_CORR;}
 /* end of GSA options */

 /* For identification() statement */
<DYNARE_STATEMENT>prior_mc {return token::PRIOR_MC;}
<DYNARE_STATEMENT>advanced {return token::ADVANCED;}
<DYNARE_STATEMENT>max_dim_cova_group {return token::MAX_DIM_COVA_GROUP;}
<DYNARE_STATEMENT>gsa_sample_file {return token::GSA_SAMPLE_FILE;}

<DYNARE_STATEMENT>use_shock_groups {return token::USE_SHOCK_GROUPS;}
<DYNARE_STATEMENT>colormap {return token::COLORMAP;}

<DYNARE_STATEMENT,DYNARE_BLOCK>[a-z_][a-z0-9_]* {
  yylval->build<string>(yytext);
  return token::NAME;
}

<DYNARE_STATEMENT,DYNARE_BLOCK>((([0-9]*\.[0-9]+)|([0-9]+\.))([ed][-+]?[0-9]+)?)|([0-9]+[ed][-+]?[0-9]+) {
  yylval->build<string>(yytext);
  return token::FLOAT_NUMBER;
}

<DYNARE_STATEMENT,DYNARE_BLOCK>[0-9]+ {
  yylval->build<string>(yytext);
  return token::INT_NUMBER;
}

<DATES_STATEMENT>\( { yylval->as<string>().append(yytext); dates_parens_nb++; }
<DATES_STATEMENT>\) {
                      yylval->as<string>().append(yytext);
                      if (--dates_parens_nb == 0)
                      {
                        BEGIN DYNARE_STATEMENT;
                        return token::DATES;
                      }
                    }
<DATES_STATEMENT>.  { yylval->as<string>().append(yytext); }

<DYNARE_BLOCK>\|e { return token::PIPE_E; }
<DYNARE_BLOCK>\|x { return token::PIPE_X; }
<DYNARE_BLOCK>\|p { return token::PIPE_P; }

<DYNARE_STATEMENT,DYNARE_BLOCK>\'[^\']+\' {
  yylval->build<string>(yytext + 1).pop_back();
  return token::QUOTED_STRING;
}


 /* Verbatim Block */
<INITIAL>verbatim[[:space:]]*;   {
                                   BEGIN VERBATIM_BLOCK;
                                 }
<VERBATIM_BLOCK>end[[:space:]]*; {
                                   BEGIN INITIAL;
                                 }
<VERBATIM_BLOCK>\n      {
                          if (strlen(yytext) > 1)
                             driver.add_verbatim_remove_charset(yytext, "\n");
                        }
<VERBATIM_BLOCK>.       { yymore(); }
<VERBATIM_BLOCK><<EOF>> {
                          driver.add_verbatim(eofbuff);
                          yyterminate();
                        }


 /* An instruction starting with a recognized symbol (which is not a modfile local
    or an external function) is passed as NAME, otherwise it is a native statement
    until the end of the line.
    We exclude modfile local vars because the user may want to modify their value
    using a Matlab assignment statement.
    We also exclude external functions because the user may have used a Matlab matrix
    element in initval (in which case Dynare recognizes the matrix name as an external
    function symbol), and may want to modify the matrix later with Matlab statements.
 */
<INITIAL>[a-z_][a-z0-9_]* {
  if (driver.symbol_exists_and_is_not_modfile_local_or_external_function(yytext))
    {
      BEGIN DYNARE_STATEMENT;
      yylval->build<string>(yytext);
      return token::NAME;
    }
  else
    {
      /* Enter a native block */
      BEGIN NATIVE;
      yyless(0);
    }
}

 /*
    For joint prior statement, match [symbol, symbol, ...]
    If no match, begin native and push everything back on stack

    We produce SYMBOL_VEC in Flex (instead of matching `'[' symbol_list ']'`
    in Bison because the pattern also matches potential native statements
    (e.g. function returns from a MATLAB/Octave function). Hence, we need to
    be able to back out of the statement if we realize it's a native statement
    and move to the NATIVE context
 */
<INITIAL>\[([[:space:]]*[a-z_][a-z0-9_]*[[:space:]]*,{1}[[:space:]]*)*([[:space:]]*[a-z_][a-z0-9_]*[[:space:]]*){1}\] {
  string yytextcpy = string(yytext);
  yytextcpy.erase(remove(yytextcpy.begin(), yytextcpy.end(), '['), yytextcpy.end());
  yytextcpy.erase(remove(yytextcpy.begin(), yytextcpy.end(), ']'), yytextcpy.end());
  yytextcpy.erase(remove(yytextcpy.begin(), yytextcpy.end(), ' '), yytextcpy.end());
  istringstream ss(yytextcpy);
  string token;
  vector<string> val;

  bool dynare_statement = true;

  while(getline(ss, token, ','))
    if (driver.symbol_exists_and_is_not_modfile_local_or_external_function(token))
      val.push_back(token);
    else
      {
        BEGIN NATIVE;
        yyless(0);
        dynare_statement = false;
        break;
      }
  if (dynare_statement)
    {
      BEGIN DYNARE_STATEMENT;
      yylval->build<vector<string>>(val);
      return token::SYMBOL_VEC;
    }
}

 /* Enter a native block */
<INITIAL>. { BEGIN NATIVE; yyless(0); }

 /* Add the native statement */
<NATIVE>{
  [^/%*\n\.\'\"]*             |
  \'                          |
  \'[^\'\n]*\'                |
  \"[^\"\n]*\"                |
  \.{1,2}                     |
  \*                          |
  \/                          { yymore(); eofbuff = string(yytext); }
  \.{3,}[[:space:]]*\n        { driver.add_native_remove_charset(yytext, "\n"); }
  \n                          {
                                if (strlen(yytext) > 1)
                                  driver.add_native_remove_charset(yytext, "\n");
                                BEGIN INITIAL;
                              }
  <<EOF>>                     {
                                driver.add_native(eofbuff);
                                yyterminate();
                              }
  \.{3,}[[:space:]]*%.*\n     |
  %[^\n]*                     { driver.add_native_remove_charset(yytext, "%"); }
  \.{3,}[[:space:]]*"//".*\n  |
  "//"[^\n]*                  { driver.add_native_remove_charset(yytext, "//"); }
  \.{3,}[[:space:]]*"/*"      {
                                driver.add_native_remove_charset(yytext, "/*");
                                BEGIN NATIVE_COMMENT;
                              }
  "/*"                        {
                                driver.add_native_remove_charset(yytext, "/*");
                                comment_caller = NATIVE;
                                BEGIN COMMENT;
                              }
}

<NATIVE_COMMENT>"*/"[[:space:]]*\n   { BEGIN NATIVE; }
<NATIVE_COMMENT>.

<INITIAL,DYNARE_STATEMENT,DYNARE_BLOCK,COMMENT,DATES_STATEMENT,LINE1,LINE2,LINE3,NATIVE_COMMENT><<EOF>> { yyterminate(); }

<*>.      { driver.error(*yylloc, "character unrecognized by lexer"); }
%%

DynareFlex::DynareFlex(istream* in, ostream* out)
  : DynareFlexLexer{in, out}
{
}

void
DynareFlex::location_increment(Dynare::parser::location_type *yylloc, const char *yytext)
{
  while (*yytext != 0)
    if (*yytext++ == '\n')
      yylloc->lines(1);
    else
      yylloc->columns(1);
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
  exit(EXIT_FAILURE);
}
