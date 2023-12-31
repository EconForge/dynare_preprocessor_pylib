// -*- C++ -*-
/*
 * Copyright © 2003-2023 Dynare Team
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

%language "c++"
%require "3.2"
%defines
%define api.value.type variant
%define api.namespace {Dynare}
%define api.location.file "DynareBisonLocation.hh"
%define parse.assert
%define parse.error verbose
%define parse.trace

%code top {
class ParsingDriver;
}

%code requires {
#include "CommonEnums.hh"
#include "ExprNode.hh"
}

%param { ParsingDriver &driver }

%locations
%initial-action
{
  // Initialize the locations' filenames to the filename maintained by the lexer
  @$.begin.filename = @$.end.filename = &(driver.lexer->filename);
}

%code {
/* Little hack: we redefine the macro which computes the locations, because
   we need to access the location from within the parsing driver for error
   and warning messages. */
#define YYLLOC_DEFAULT(Current, Rhs, N)                         \
  do {                                                          \
    if (N)                                                      \
      {                                                         \
        (Current).begin = YYRHSLOC(Rhs, 1).begin;               \
        (Current).end   = YYRHSLOC(Rhs, N).end;                 \
      }                                                         \
    else                                                        \
      {                                                         \
        (Current).begin = (Current).end = YYRHSLOC(Rhs, 0).end;	\
      }                                                         \
    driver.location = (Current);                                \
  } while(false)

#include "ParsingDriver.hh"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex
}

%token AIM_SOLVER ANALYTIC_DERIVATION ANALYTIC_DERIVATION_MODE AR POSTERIOR_SAMPLING_METHOD
%token BALANCED_GROWTH_TEST_TOL BAYESIAN_IRF BETA_PDF BLOCK USE_CALIBRATION SILENT_OPTIMIZER
%token BVAR_DENSITY BVAR_FORECAST BVAR_IRF NODECOMPOSITION DR_DISPLAY_TOL HUGE_NUMBER FIG_NAME WRITE_XLS
%token BVAR_PRIOR_DECAY BVAR_PRIOR_FLAT BVAR_PRIOR_LAMBDA INTERACTIVE SCREEN_SHOCKS STEADYSTATE
%token BVAR_PRIOR_MU BVAR_PRIOR_OMEGA BVAR_PRIOR_TAU BVAR_PRIOR_TRAIN DETAIL_PLOT TYPE
%token BVAR_REPLIC BYTECODE ALL_VALUES_REQUIRED PROPOSAL_DISTRIBUTION REALTIME VINTAGE
%token CALIB_SMOOTHER CHANGE_TYPE CHECK CONDITIONAL_FORECAST CONDITIONAL_FORECAST_PATHS CONF_SIG CONSTANT CONTROLLED_VAREXO CORR CUTOFF CYCLE_REDUCTION LOGARITHMIC_REDUCTION
%token COMMA CONSIDER_ALL_ENDOGENOUS CONSIDER_ALL_ENDOGENOUS_AND_AUXILIARY CONSIDER_ONLY_OBSERVED INITIAL_CONDITION_DECOMPOSITION
%token DATAFILE FILE SERIES DOUBLING DR_CYCLE_REDUCTION_TOL DR_LOGARITHMIC_REDUCTION_TOL DR_LOGARITHMIC_REDUCTION_MAXITER DR_ALGO DROP DSAMPLE DYNASAVE DYNATYPE CALIBRATION DIFFERENTIATE_FORWARD_VARS
%token END ENDVAL EQUAL ESTIMATION ESTIMATED_PARAMS ESTIMATED_PARAMS_BOUNDS ESTIMATED_PARAMS_INIT EXTENDED_PATH ENDOGENOUS_PRIOR EXPRESSION
%token FILENAME DIRNAME FILTER_STEP_AHEAD FILTERED_VARS FIRST_OBS FIRST_SIMULATION_PERIOD LAST_SIMULATION_PERIOD LAST_OBS
%token SET_TIME OSR_PARAMS_BOUNDS KEEP_KALMAN_ALGO_IF_SINGULARITY_IS_DETECTED
%token <string> FALSE FLOAT_NUMBER DATES
%token DEFAULT FIXED_POINT FLIP OPT_ALGO COMPILATION_SETUP COMPILER ADD_FLAGS SUBSTITUTE_FLAGS ADD_LIBS SUBSTITUTE_LIBS
%token FORECAST K_ORDER_SOLVER INSTRUMENTS SHIFT MEAN STDEV VARIANCE MODE INTERVAL SHAPE DOMAINN
%token GAMMA_PDF GRAPH GRAPH_FORMAT CONDITIONAL_VARIANCE_DECOMPOSITION NOCHECK STD
%token HISTVAL HISTVAL_FILE HOMOTOPY_SETUP HOMOTOPY_MODE HOMOTOPY_STEPS HOMOTOPY_FORCE_CONTINUE HP_FILTER HP_NGRID FILTERED_THEORETICAL_MOMENTS_GRID HYBRID ONE_SIDED_HP_FILTER
%token IDENTIFICATION INF_CONSTANT INITVAL INITVAL_FILE BOUNDS JSCALE INIT INFILE INVARS
%token <string> INT_NUMBER
%token CONDITIONAL_LIKELIHOOD
%token INV_GAMMA_PDF INV_GAMMA1_PDF INV_GAMMA2_PDF IRF IRF_SHOCKS IRF_PLOT_THRESHOLD IRF_CALIBRATION
%token FAST_KALMAN_FILTER KALMAN_ALGO KALMAN_TOL DIFFUSE_KALMAN_TOL SCHUR_VEC_TOL SUBSAMPLES OPTIONS TOLF TOLX PLOT_INIT_DATE PLOT_END_DATE
%token LAPLACE LIK_INIT LINEAR LINEAR_DECOMPOSITION LOAD_IDENT_FILES LOAD_MH_FILE LOAD_RESULTS_AFTER_LOAD_MH LOAD_PARAMS_AND_STEADY_STATE LOGLINEAR LOGDATA LYAPUNOV LINEAR_APPROXIMATION
%token LYAPUNOV_COMPLEX_THRESHOLD LYAPUNOV_FIXED_POINT_TOL LYAPUNOV_DOUBLING_TOL LOG_DEFLATOR LOG_TREND_VAR LOG_GROWTH_FACTOR
%token MATCHED_MOMENTS MARKOWITZ MARGINAL_DENSITY MAX MAXIT
%token MFS MH_CONF_SIG MH_DROP MH_INIT_SCALE MH_INIT_SCALE_FACTOR MH_JSCALE MH_TUNE_JSCALE MH_TUNE_GUESS MH_POSTERIOR_MODE_ESTIMATION MH_NBLOCKS MH_REPLIC MH_RECOVER MH_INITIALIZE_FROM_PREVIOUS_MCMC MH_INITIALIZE_FROM_PREVIOUS_MCMC_DIRECTORY MH_INITIALIZE_FROM_PREVIOUS_MCMC_RECORD MH_INITIALIZE_FROM_PREVIOUS_MCMC_PRIOR
%token POSTERIOR_MAX_SUBSAMPLE_DRAWS MIN MINIMAL_SOLVING_PERIODS
%token MODE_CHECK MODE_CHECK_NEIGHBOURHOOD_SIZE MODE_CHECK_SYMMETRIC_PLOTS MODE_CHECK_NUMBER_OF_POINTS MODE_COMPUTE MODE_FILE MODEL MODEL_COMPARISON MODEL_INFO MSHOCKS ABS SIGN
%token MODEL_DIAGNOSTICS MODIFIEDHARMONICMEAN MOMENTS_VARENDO CONTEMPORANEOUS_CORRELATION DIFFUSE_FILTER SUB_DRAWS TAPER_STEPS GEWEKE_INTERVAL RAFTERY_LEWIS_QRS RAFTERY_LEWIS_DIAGNOSTICS BROOKS_GELMAN_PLOTROWS MCMC_JUMPING_COVARIANCE MOMENT_CALIBRATION
%token NUMBER_OF_PARTICLES RESAMPLING SYSTEMATIC GENERIC RESAMPLING_THRESHOLD RESAMPLING_METHOD KITAGAWA STRATIFIED SMOOTH
%token CPF_WEIGHTS AMISANOTRISTANI MURRAYJONESPARSLOW WRITE_EQUATION_TAGS FILTER_INITIAL_STATE
%token NONLINEAR_FILTER_INITIALIZATION FILTER_ALGORITHM PROPOSAL_APPROXIMATION CUBATURE UNSCENTED MONTECARLO DISTRIBUTION_APPROXIMATION
%token <string> NAME
%token USE_PENALIZED_OBJECTIVE_FOR_HESSIAN INIT_STATE FAST_REALTIME RESCALE_PREDICTION_ERROR_COVARIANCE GENERATE_IRFS
%token NAN_CONSTANT NO_STATIC NOBS NOCONSTANT NODISPLAY NOCORR NODIAGNOSTIC NOFUNCTIONS NO_HOMOTOPY
%token NOGRAPH POSTERIOR_NOGRAPH POSTERIOR_GRAPH NOMOMENTS NOPRINT NORMAL_PDF SAVE_DRAWS MODEL_NAME STDERR_MULTIPLES DIAGONAL_ONLY
%token DETERMINISTIC_TRENDS OBSERVATION_TRENDS OPTIM OPTIM_WEIGHTS ORDER OSR OSR_PARAMS MAX_DIM_COVA_GROUP ADVANCED OUTFILE OUTVARS OVERWRITE DISCOUNT
%token PARALLEL_LOCAL_FILES PARAMETERS PARAMETER_SET PARTIAL_INFORMATION PERIODS PERIOD PLANNER_OBJECTIVE PLOT_CONDITIONAL_FORECAST PLOT_PRIORS PREFILTER PRESAMPLE
%token PERFECT_FORESIGHT_SETUP PERFECT_FORESIGHT_SOLVER NO_POSTERIOR_KERNEL_DENSITY FUNCTION
%token PERFECT_FORESIGHT_WITH_EXPECTATION_ERRORS_SETUP PERFECT_FORESIGHT_WITH_EXPECTATION_ERRORS_SOLVER
%token PRINT PRIOR_MC PRIOR_TRUNC PRIOR_MODE PRIOR_MEAN POSTERIOR_MODE POSTERIOR_MEAN POSTERIOR_MEDIAN MLE_MODE PRUNING PARTICLE_FILTER_OPTIONS
%token <string> QUOTED_STRING
%token QZ_CRITERIUM QZ_ZERO_THRESHOLD DSGE_VAR DSGE_VARLAG DSGE_PRIOR_WEIGHT TRUNCATE PIPE_E PIPE_X PIPE_P
%token RELATIVE_IRF REPLIC SIMUL_REPLIC RPLOT SAVE_PARAMS_AND_STEADY_STATE PARAMETER_UNCERTAINTY TARGETS
%token SHOCKS HETEROSKEDASTIC_SHOCKS SHOCK_DECOMPOSITION SHOCK_GROUPS USE_SHOCK_GROUPS SIMUL SIMUL_ALGO SIMUL_SEED ENDOGENOUS_TERMINAL_PERIOD
%token SMOOTHER SMOOTHER2HISTVAL SQUARE_ROOT_SOLVER STACK_SOLVE_ALGO STEADY_STATE_MODEL SOLVE_ALGO SOLVER_PERIODS ROBUST_LIN_SOLVE
%token STDERR STEADY STOCH_SIMUL SYLVESTER SYLVESTER_FIXED_POINT_TOL REGIMES REGIME REALTIME_SHOCK_DECOMPOSITION CONDITIONAL UNCONDITIONAL
%token TEX RAMSEY_MODEL RAMSEY_POLICY RAMSEY_CONSTRAINTS PLANNER_DISCOUNT PLANNER_DISCOUNT_LATEX_NAME
%token DISCRETIONARY_POLICY DISCRETIONARY_TOL EVALUATE_PLANNER_OBJECTIVE
%token OCCBIN_SETUP OCCBIN_SOLVER OCCBIN_WRITE_REGIMES OCCBIN_GRAPH SIMUL_MAXIT LIKELIHOOD_MAXIT SMOOTHER_MAXIT SIMUL_PERIODS LIKELIHOOD_PERIODS SMOOTHER_PERIODS
%token SIMUL_CURB_RETRENCH LIKELIHOOD_CURB_RETRENCH SMOOTHER_CURB_RETRENCH SIMUL_CHECK_AHEAD_PERIODS SIMUL_MAX_CHECK_AHEAD_PERIODS SIMUL_RESET_CHECK_AHEAD_PERIODS
%token LIKELIHOOD_CHECK_AHEAD_PERIODS LIKELIHOOD_MAX_CHECK_AHEAD_PERIODS SMOOTHER_CHECK_AHEAD_PERIODS SMOOTHER_MAX_CHECK_AHEAD_PERIODS
%token SIMUL_DEBUG SMOOTHER_DEBUG SIMUL_PERIODIC_SOLUTION LIKELIHOOD_PERIODIC_SOLUTION SMOOTHER_PERIODIC_SOLUTION
%token LIKELIHOOD_INVERSION_FILTER SMOOTHER_INVERSION_FILTER FILTER_USE_RELEXATION
%token LIKELIHOOD_PIECEWISE_KALMAN_FILTER SMOOTHER_PIECEWISE_KALMAN_FILTER LIKELIHOOD_MAX_KALMAN_ITERATIONS
%token <string> TEX_NAME TRUE BIND RELAX ERROR_BIND ERROR_RELAX
%token UNIFORM_PDF UNIT_ROOT_VARS USE_DLL USEAUTOCORR GSA_SAMPLE_FILE USE_UNIVARIATE_FILTERS_IF_SINGULARITY_IS_DETECTED
%token VALUES SCALES VAR VAREXO VAREXO_DET VARIABLE VAROBS VAREXOBS PREDETERMINED_VARIABLES VAR_EXPECTATION VAR_EXPECTATION_MODEL PLOT_SHOCK_DECOMPOSITION MODEL_LOCAL_VARIABLE
%token WRITE_LATEX_DYNAMIC_MODEL WRITE_LATEX_STATIC_MODEL WRITE_LATEX_ORIGINAL_MODEL WRITE_LATEX_STEADY_STATE_MODEL
%token XLS_SHEET XLS_RANGE LMMCP BANDPASS_FILTER COLORMAP VAR_MODEL PAC_MODEL QOQ YOY AOA PAC_EXPECTATION TREND_COMPONENT_MODEL
%left EQUAL_EQUAL EXCLAMATION_EQUAL
%left LESS GREATER LESS_EQUAL GREATER_EQUAL
%left PLUS MINUS
%left TIMES DIVIDE
%precedence UNARY
%nonassoc POWER
%token EXP LOG LN LOG10 SIN COS TAN ASIN ACOS ATAN SINH COSH TANH ASINH ACOSH ATANH ERF ERFC DIFF ADL AUXILIARY_MODEL_NAME
%token SQRT CBRT NORMCDF NORMPDF STEADY_STATE EXPECTATION
/* GSA analysis */
%token DYNARE_SENSITIVITY MORRIS STAB REDFORM PPRIOR PRIOR_RANGE PPOST ILPTAU MORRIS_NLIV
%token MORRIS_NTRA NSAM LOAD_REDFORM LOAD_RMSE LOAD_STAB ALPHA2_STAB LOGTRANS_REDFORM THRESHOLD_REDFORM
%token KSSTAT_REDFORM ALPHA2_REDFORM NAMENDO NAMLAGENDO NAMEXO RMSE LIK_ONLY VAR_RMSE PFILT_RMSE ISTART_RMSE
%token ALPHA_RMSE ALPHA2_RMSE
/* end of GSA analysis*/
%token FREQ INITIAL_YEAR INITIAL_SUBPERIOD FINAL_YEAR FINAL_SUBPERIOD DATA VLIST
%token VLISTLOG VLISTPER SPECTRAL_DENSITY INIT2SHOCKS
%token RESTRICTION RESTRICTION_FNAME CROSS_RESTRICTIONS NLAGS CONTEMP_REDUCED_FORM REAL_PSEUDO_FORECAST
%token DUMMY_OBS NSTATES INDXSCALESSTATES NO_BAYESIAN_PRIOR SPECIFICATION SIMS_ZHA
%token <string> ALPHA BETA ABAND NINV CMS NCMS CNUM GAMMA INV_GAMMA INV_GAMMA1 INV_GAMMA2 NORMAL UNIFORM EPS PDF FIG DR NONE PRIOR PRIOR_VARIANCE HESSIAN IDENTITY_MATRIX DIRICHLET DIAGONAL OPTIMAL
%token GSIG2_LMDM Q_DIAG FLAT_PRIOR NCSK NSTD WEIBULL WEIBULL_PDF
%token INDXPARR INDXOVR INDXAP APBAND INDXIMF INDXFORE FOREBAND INDXGFOREHAT INDXGIMFHAT
%token INDXESTIMA INDXGDLS EQ_MS FILTER_COVARIANCE UPDATED_COVARIANCE FILTER_DECOMPOSITION SMOOTHED_STATE_UNCERTAINTY SMOOTHER_REDUX
%token EQ_CMS TLINDX TLNUMBER RESTRICTIONS POSTERIOR_SAMPLER_OPTIONS
%token OUTPUT_FILE_TAG HORIZON
%token SBVAR TREND_VAR DEFLATOR GROWTH_FACTOR MS_IRF MS_VARIANCE_DECOMPOSITION GROWTH
%token MS_ESTIMATION MS_SIMULATION MS_COMPUTE_MDD MS_COMPUTE_PROBABILITIES MS_FORECAST
%token SVAR_IDENTIFICATION EQUATION EXCLUSION LAG UPPER_CHOLESKY LOWER_CHOLESKY MONTHLY QUARTERLY
%token MARKOV_SWITCHING CHAIN DURATION NUMBER_OF_REGIMES NUMBER_OF_LAGS EPILOGUE
%token SVAR SVAR_GLOBAL_IDENTIFICATION_CHECK COEFF COEFFICIENTS VARIANCES CONSTANTS EQUATIONS
%token EXTERNAL_FUNCTION EXT_FUNC_NAME EXT_FUNC_NARGS FIRST_DERIV_PROVIDED SECOND_DERIV_PROVIDED
%token SELECTED_VARIABLES_ONLY COVA_COMPUTE SIMULATION_FILE_TAG FILE_TAG
%token NO_ERROR_BANDS ERROR_BAND_PERCENTILES SHOCKS_PER_PARAMETER NO_CREATE_INIT
%token SHOCK_DRAWS FREE_PARAMETERS MEDIAN DATA_OBS_NBR NEIGHBORHOOD_WIDTH PVALUE_KS PVALUE_CORR
%token FILTERED_PROBABILITIES REAL_TIME_SMOOTHED PRIOR_FUNCTION POSTERIOR_FUNCTION SAMPLING_DRAWS
%token PROPOSAL_TYPE PROPOSAL_UPPER_BOUND PROPOSAL_LOWER_BOUND PROPOSAL_DRAWS USE_MEAN_CENTER
%token ADAPTIVE_MH_DRAWS THINNING_FACTOR COEFFICIENTS_PRIOR_HYPERPARAMETERS
%token CONVERGENCE_STARTING_VALUE CONVERGENCE_ENDING_VALUE CONVERGENCE_INCREMENT_VALUE
%token MAX_ITERATIONS_STARTING_VALUE MAX_ITERATIONS_INCREMENT_VALUE MAX_BLOCK_ITERATIONS
%token MAX_REPEATED_OPTIMIZATION_RUNS FUNCTION_CONVERGENCE_CRITERION SAVE_REALTIME
%token PARAMETER_CONVERGENCE_CRITERION NUMBER_OF_LARGE_PERTURBATIONS NUMBER_OF_SMALL_PERTURBATIONS
%token NUMBER_OF_POSTERIOR_DRAWS_AFTER_PERTURBATION MAX_NUMBER_OF_STAGES
%token RANDOM_FUNCTION_CONVERGENCE_CRITERION RANDOM_PARAMETER_CONVERGENCE_CRITERION NO_INIT_ESTIMATION_CHECK_FIRST_OBS
%token HETEROSKEDASTIC_FILTER TIME_SHIFT STRUCTURAL CONSTANT_SIMULATION_LENGTH
%token SURPRISE OCCBIN_CONSTRAINTS
%token PAC_TARGET_INFO COMPONENT TARGET AUXNAME AUXNAME_TARGET_NONSTATIONARY PAC_TARGET_NONSTATIONARY
%token <string> KIND LL DL DD ADD MULTIPLY
/* Method of Moments */
%token GMM SMM IRF_MATCHING
%token METHOD_OF_MOMENTS MOM_METHOD SIMULATION_METHOD
%token BARTLETT_KERNEL_LAG WEIGHTING_MATRIX WEIGHTING_MATRIX_SCALING_FACTOR ANALYTIC_STANDARD_ERRORS ANALYTIC_JACOBIAN PENALIZED_ESTIMATOR VERBOSE
%token SIMULATION_MULTIPLE MOM_SEED SEED BOUNDED_SHOCK_SUPPORT ADDITIONAL_OPTIMIZER_STEPS MOM_SE_TOLX SE_TOLX MOM_BURNIN BURNIN
%token IRF_MATCHING_FILE ADD_TINY_NUMBER_TO_CHOLESKY
%token EQTAGS
%token ANALYTICAL_GIRF IRF_IN_PERCENT EMAS_GIRF EMAS_DROP EMAS_TOLF EMAS_MAX_ITER
%token NO_IDENTIFICATION_STRENGTH NO_IDENTIFICATION_REDUCEDFORM NO_IDENTIFICATION_MOMENTS
%token NO_IDENTIFICATION_MINIMAL NO_IDENTIFICATION_SPECTRUM NORMALIZE_JACOBIANS GRID_NBR
%token TOL_RANK TOL_DERIV TOL_SV CHECKS_VIA_SUBSETS MAX_DIM_SUBSETS_GROUPS
%token MAX_NROWS SQUEEZE_SHOCK_DECOMPOSITION WITH_EPILOGUE MODEL_REMOVE MODEL_REPLACE MODEL_OPTIONS
%token VAR_REMOVE ESTIMATED_PARAMS_REMOVE BLOCK_STATIC BLOCK_DYNAMIC INCIDENCE RESID NON_ZERO LEARNT_IN PLUS_EQUAL TIMES_EQUAL
%token FSOLVE_OPTIONS
%token ENDVAL_STEADY STEADY_SOLVE_ALGO STEADY_MAXIT STEADY_TOLF STEADY_TOLX STEADY_MARKOWITZ
%token HOMOTOPY_MAX_COMPLETION_SHARE HOMOTOPY_MIN_STEP_SIZE HOMOTOPY_INITIAL_STEP_SIZE HOMOTOPY_STEP_SIZE_INCREASE_SUCCESS_COUNT
%token HOMOTOPY_LINEARIZATION_FALLBACK HOMOTOPY_MARGINAL_LINEARIZATION_FALLBACK

%token <vector<string>> SYMBOL_VEC

%type <expr_t> expression expression_or_empty
%type <expr_t> equation hand_side
%type <string> non_negative_number signed_number signed_integer date_str
%type <string> filename symbol namespace_qualified_filename namespace_qualified_symbol
%type <string> date_expr signed_inf signed_number_w_inf range
%type <string> integer_range signed_integer_range boolean
%type <string> name_value_pair name_value_pair_list
%type <string> name_value_pair_with_boolean name_value_pair_with_boolean_list
%type <string> name_value_pair_with_suboptions name_value_pair_with_suboptions_list
%type <SymbolType> change_type_arg
%type <vector<string>> vec_str vec_str_1
%type <vector<string>> vec_value vec_value_1 vec_value_w_inf vec_value_w_inf_1
%type <vector<vector<string>>> vec_of_vec_value vec_of_vec_value_1
%type <vector<string>> symbol_list symbol_list_or_wildcard
%type <vector<int>> vec_int_elem vec_int_1 vec_int vec_int_number
%type <PriorDistributions> prior_pdf prior_distribution
%type <pair<expr_t,expr_t>> calibration_range
%type <pair<string,string>> partition_elem subsamples_eq_opt integer_range_w_inf
%type <vector<pair<string,string>>> partition partition_1 tag_pair_list_for_selection symbol_list_with_tex
%type <tuple<string,string,string,string>> prior_eq_opt options_eq_opt
%type <vector<pair<int, int>>> period_list
%type <vector<expr_t>> matched_moments_list value_list
%type <tuple<string, BinaryOpNode *, BinaryOpNode *, expr_t, expr_t>> occbin_constraints_regime
%type <vector<tuple<string, BinaryOpNode *, BinaryOpNode *, expr_t, expr_t>>> occbin_constraints_regimes_list
%type <map<string, expr_t>> occbin_constraints_regime_options_list
%type <pair<string, expr_t>> occbin_constraints_regime_option
%type <PacTargetKind> pac_target_kind
%type <vector<tuple<string, string, vector<pair<string, string>>>>> symbol_list_with_tex_and_partition
%%

%start statement_list;

statement_list : statement
               | statement_list statement
               ;

statement : parameters
          | var
          | varexo
          | varexo_det
          | predetermined_variables
          | model_local_variable
          | change_type
          | model
          | initval
          | initval_file
          | endval
          | histval
          | init_param
          | shocks
          | mshocks
          | heteroskedastic_shocks
          | steady
          | check
          | simul
          | stoch_simul
          | estimation
          | estimated_params
          | estimated_params_bounds
          | estimated_params_init
          | estimated_params_remove
          | set_time
          | data
          | epilogue
          | var_model
          | pac_model
          | trend_component_model
          | prior
          | prior_eq
          | subsamples
          | subsamples_eq
          | options
          | options_eq
          | varobs
          | deterministic_trends
          | observation_trends
          | filter_initial_state
          | varexobs
          | unit_root_vars
          | dsample
          | rplot
          | optim_weights
          | osr_params
          | osr_params_bounds
          | osr
          | dynatype
          | dynasave
          | model_comparison
          | model_info
          | planner_objective
          | ramsey_model
          | ramsey_policy
          | ramsey_constraints
          | evaluate_planner_objective
          | occbin_setup
          | occbin_solver
          | occbin_write_regimes
          | occbin_graph
          | discretionary_policy
          | bvar_density
          | bvar_forecast
          | bvar_irf
          | sbvar
          | dynare_sensitivity
          | homotopy_setup
          | forecast
          | load_params_and_steady_state
          | save_params_and_steady_state
          | identification
          | write_latex_dynamic_model
          | write_latex_static_model
          | write_latex_original_model
          | write_latex_steady_state_model
          | shock_decomposition
          | realtime_shock_decomposition
          | plot_shock_decomposition
          | initial_condition_decomposition
          | squeeze_shock_decomposition
          | conditional_forecast
          | conditional_forecast_paths
          | plot_conditional_forecast
          | svar_identification
          | svar_global_identification_check
          | markov_switching
          | svar
          | external_function
          | steady_state_model
          | trend_var
          | generate_irfs
          | log_trend_var
          | ms_estimation
          | ms_simulation
          | ms_compute_mdd
          | ms_compute_probabilities
          | ms_forecast
          | ms_irf
          | ms_variance_decomposition
          | calib_smoother
          | extended_path
          | model_diagnostics
          | moment_calibration
          | irf_calibration
          | smoother2histval
          | histval_file
          | perfect_foresight_setup
          | perfect_foresight_solver
          | perfect_foresight_with_expectation_errors_setup
          | perfect_foresight_with_expectation_errors_solver
          | prior_function
          | posterior_function
          | method_of_moments
          | shock_groups
          | init2shocks
          | var_expectation_model
          | compilation_setup
          | matched_moments
          | occbin_constraints
          | model_remove
          | model_replace
          | model_options
          | var_remove
          | pac_target_info
          | resid
          ;

dsample : DSAMPLE INT_NUMBER ';'
          { driver.dsample($2); }
        | DSAMPLE INT_NUMBER INT_NUMBER ';'
          { driver.dsample($2, $3); }
        ;

symbol_list : symbol_list symbol
              {
                $$ = $1;
                $$.push_back($2);
              }
            | symbol_list COMMA symbol
              {
                $$ = $1;
                $$.push_back($3);
              }
            | symbol
              { $$ = { $1 }; }
            ;

symbol_list_or_wildcard : symbol_list
                        | ':'
                          { $$ = { ":" }; }
                        ;

symbol_list_with_tex : symbol_list_with_tex symbol
                       {
                         $$ = $1;
                         $$.emplace_back($2, "");
                       }
                     | symbol_list_with_tex COMMA symbol
                       {
                         $$ = $1;
                         $$.emplace_back($3, "");
                       }
                     | symbol
                       { $$ = { { $1, "" }}; }
                     | symbol_list_with_tex symbol TEX_NAME
                       {
                         $$ = $1;
                         $$.emplace_back($2, $3);
                       }
                     | symbol_list_with_tex COMMA symbol TEX_NAME
                       {
                         $$ = $1;
                         $$.emplace_back($3, $4);
                       }
                     | symbol TEX_NAME
                       { $$ = { { $1, $2 } }; }
                     ;

partition_elem : symbol EQUAL QUOTED_STRING
                 { $$ = { $1, $3 }; }

partition_1 : '(' partition_elem
              { $$ = { $2 }; }
            | '(' COMMA partition_elem
              { $$ = { $3 }; }
            | partition_1 partition_elem
              {
                $$ = $1;
                $$.push_back($2);
              }
            | partition_1 COMMA partition_elem
              {
                $$ = $1;
                $$.push_back($3);
              }
            ;

partition : partition_1 ')'
          | partition_1 COMMA ')'
          ;

symbol_list_with_tex_and_partition : symbol_list_with_tex_and_partition symbol
                                     {
                                       $$ = $1;
                                       $$.emplace_back($2, "", vector<pair<string,string>>{});
                                     }
                                   | symbol_list_with_tex_and_partition COMMA symbol
                                     {
                                       $$ = $1;
                                       $$.emplace_back($3, "", vector<pair<string,string>>{});
                                     }
                                   | symbol
                                     { $$ = { { $1, "", {} }}; }
                                   | symbol_list_with_tex_and_partition symbol partition
                                     {
                                       $$ = $1;
                                       $$.emplace_back($2, "", $3);
                                     }
                                   | symbol_list_with_tex_and_partition COMMA symbol partition
                                     {
                                       $$ = $1;
                                       $$.emplace_back($3, "", $4);
                                     }
                                   | symbol partition
                                     { $$ = { { $1, "", $2 }}; }
                                   | symbol_list_with_tex_and_partition symbol TEX_NAME
                                     {
                                       $$ = $1;
                                       $$.emplace_back($2, $3, vector<pair<string,string>>{});
                                     }
                                   | symbol_list_with_tex_and_partition COMMA symbol TEX_NAME
                                     {
                                       $$ = $1;
                                       $$.emplace_back($3, $4, vector<pair<string,string>>{});
                                     }
                                   | symbol TEX_NAME
                                     { $$ = { { $1, $2, {} }}; }
                                   | symbol_list_with_tex_and_partition symbol TEX_NAME partition
                                     {
                                       $$ = $1;
                                       $$.emplace_back($2, $3, $4);
                                     }
                                   | symbol_list_with_tex_and_partition COMMA symbol TEX_NAME partition
                                     {
                                       $$ = $1;
                                       $$.emplace_back($3, $4, $5);
                                     }
                                   | symbol TEX_NAME partition
                                     { $$ = { { $1, $2, $3 }}; }
                                   ;

rplot : RPLOT symbol_list ';' { driver.rplot($2); };

trend_var : TREND_VAR '(' GROWTH_FACTOR EQUAL { driver.begin_model(); } hand_side ')' symbol_list_with_tex ';'
            { driver.end_trend_var(false, $6, $8); }
          ;

log_trend_var : LOG_TREND_VAR '(' LOG_GROWTH_FACTOR EQUAL { driver.begin_model(); } hand_side ')' symbol_list_with_tex ';'
                { driver.end_trend_var(true, $6, $8); }
              ;

var : VAR symbol_list_with_tex_and_partition ';'
      { driver.var($2, false); }
    | VAR '(' LOG ')' symbol_list_with_tex_and_partition ';'
      { driver.var($5, true); }
    | VAR '(' DEFLATOR EQUAL { driver.begin_model(); } hand_side ')' symbol_list_with_tex_and_partition ';'
      { driver.end_nonstationary_var(false, $6, $8, false); }
    | VAR '(' LOG COMMA DEFLATOR EQUAL { driver.begin_model(); } hand_side ')' symbol_list_with_tex_and_partition ';'
      { driver.end_nonstationary_var(false, $8, $10, true); }
    | VAR '(' LOG_DEFLATOR EQUAL { driver.begin_model(); } hand_side ')' symbol_list_with_tex_and_partition ';'
      { driver.end_nonstationary_var(true, $6, $8, false); }
    /* The case LOG + LOG_DEFLATOR is omitted, because it does not make much sense
       from an economic point of view (amounts to taking the log two times) */
    ;

var_remove : VAR_REMOVE symbol_list ';' { driver.var_remove($2); };

var_model : VAR_MODEL '(' var_model_options_list ')' ';' { driver.var_model(); }
          ;

var_model_options_list : var_model_options_list COMMA var_model_options
                       | var_model_options
                       ;

var_model_options : o_var_name
                  | o_var_eq_tags
                  | o_var_structural
                  ;

trend_component_model : TREND_COMPONENT_MODEL '('  trend_component_model_options_list ')' ';' { driver.trend_component_model(); }
                      ;

trend_component_model_options_list : trend_component_model_options_list COMMA trend_component_model_options
                                   | trend_component_model_options
                                   ;

trend_component_model_options : o_trend_component_model_name
                              | o_trend_component_model_targets
                              | o_trend_component_model_eq_tags
                              ;

pac_model : PAC_MODEL '(' { driver.begin_pac_model(); } pac_model_options_list ')' ';' { driver.pac_model(); };

pac_model_options_list : pac_model_options_list COMMA pac_model_options
                       | pac_model_options
                       ;

pac_model_options : o_pac_name
                  | o_pac_aux_model_name
                  | o_pac_discount
                  | o_pac_growth
                  | o_pac_auxname
                  | o_pac_kind
                  ;

var_expectation_model : VAR_EXPECTATION_MODEL '(' var_expectation_model_options_list ')' ';'
                        { driver.var_expectation_model(); }
                      ;

var_expectation_model_options_list : var_expectation_model_option
                                   | var_expectation_model_options_list COMMA var_expectation_model_option
                                   ;


var_expectation_model_option : VARIABLE EQUAL symbol
                               { driver.option_str("variable", $3); }
                             | EXPRESSION EQUAL { driver.begin_model(); } hand_side
                               {
                                 driver.var_expectation_model_expression = $4;
                                 driver.reset_data_tree();
                               }
                             | AUXILIARY_MODEL_NAME EQUAL symbol
                               { driver.option_str("auxiliary_model_name", $3); }
                             | HORIZON EQUAL INT_NUMBER
                               { driver.option_num("horizon", $3); }
                             | HORIZON EQUAL integer_range_w_inf
                               { driver.option_num("horizon", "[ " + $3.first + ' ' + $3.second + " ]"); }
                             | MODEL_NAME EQUAL symbol
                               { driver.option_str("model_name", $3); }
                             | DISCOUNT EQUAL expression
                               { driver.var_expectation_model_discount = $3; }
                             | TIME_SHIFT EQUAL signed_integer
                               { driver.option_num("time_shift", $3); }
                             ;

varexo : VAREXO symbol_list_with_tex_and_partition ';'
         { driver.varexo($2); }
       ;

varexo_det : VAREXO_DET symbol_list_with_tex_and_partition ';'
             { driver.varexo_det($2); }
           ;

predetermined_variables : PREDETERMINED_VARIABLES symbol_list ';'
                          { driver.predetermined_variables($2); }
                        ;

parameters : PARAMETERS symbol_list_with_tex_and_partition ';'
             { driver.parameters($2); }
           ;

model_local_variable : MODEL_LOCAL_VARIABLE symbol_list_with_tex ';'
                       { driver.model_local_variable($2); }
                     ;

change_type : CHANGE_TYPE '(' change_type_arg ')' symbol_list ';'
              { driver.change_type($3, $5); }
            ;

change_type_arg : PARAMETERS
                  { $$ = SymbolType::parameter; }
                | VAR
                  { $$ = SymbolType::endogenous; }
                | VAREXO
                  { $$ = SymbolType::exogenous; }
                | VAREXO_DET
                  { $$ = SymbolType::exogenousDet; }
                ;

init_param : symbol EQUAL expression ';' { driver.init_param($1, $3); };

expression : '(' expression ')'
             { $$ = $2;}
           | namespace_qualified_symbol
             { $$ = driver.add_expression_variable($1); }
           | non_negative_number
             { $$ = driver.add_non_negative_constant($1); }
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
           | MINUS expression %prec UNARY
             { $$ = driver.add_uminus($2); }
           | PLUS expression %prec UNARY
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
           | SINH '(' expression ')'
             { $$ = driver.add_sinh($3); }
           | COSH '(' expression ')'
             { $$ = driver.add_cosh($3); }
           | TANH '(' expression ')'
             { $$ = driver.add_tanh($3); }
           | ASINH '(' expression ')'
             { $$ = driver.add_asinh($3); }
           | ACOSH '(' expression ')'
             { $$ = driver.add_acosh($3); }
           | ATANH '(' expression ')'
             { $$ = driver.add_atanh($3); }
           | SQRT '(' expression ')'
             { $$ = driver.add_sqrt($3); }
           | CBRT '(' expression ')'
             { $$ = driver.add_cbrt($3); }
           | ABS '(' expression ')'
             { $$ = driver.add_abs($3); }
           | SIGN '(' expression ')'
             { $$ = driver.add_sign($3); }
           | MAX '(' expression COMMA expression ')'
             { $$ = driver.add_max($3, $5); }
           | MIN '(' expression COMMA expression ')'
             { $$ = driver.add_min($3, $5); }
           | namespace_qualified_symbol { driver.push_external_function_arg_vector_onto_stack(); } '(' comma_expression ')'
             { $$ = driver.add_model_var_or_external_function($1, false); }
           | NORMCDF '(' expression COMMA expression COMMA expression ')'
             { $$ = driver.add_normcdf($3, $5, $7); }
           | NORMCDF '(' expression ')'
             { $$ = driver.add_normcdf($3); }
           | NORMPDF '(' expression COMMA expression COMMA expression ')'
             { $$ = driver.add_normpdf($3, $5, $7); }
           | NORMPDF '(' expression ')'
             { $$ = driver.add_normpdf($3); }
           | ERF '(' expression ')'
             { $$ = driver.add_erf($3); }
           | ERFC '(' expression ')'
             { $$ = driver.add_erfc($3); }
           | NAN_CONSTANT
             { $$ = driver.add_nan_constant(); }
           | INF_CONSTANT
             { $$ = driver.add_inf_constant(); }
           ;

comma_expression : expression
                   { driver.add_external_function_arg($1); }
                 | comma_expression COMMA expression
                   { driver.add_external_function_arg($3); }
                 ;

expression_or_empty : %empty
                      { $$ = driver.add_nan_constant(); }
                    | expression
	            ;

initval : INITVAL ';' initval_list END ';'
          { driver.end_initval(false); }
        | INITVAL '(' ALL_VALUES_REQUIRED ')' ';' initval_list END ';'
          { driver.end_initval(true); }
        ;

initval_list : initval_list initval_elem
             | initval_elem
             ;

initval_elem : symbol EQUAL expression ';' { driver.init_val($1, $3); };

histval_file : HISTVAL_FILE '(' h_options_list ')' ';'
              { driver.histval_file();};

initval_file : INITVAL_FILE '(' h_options_list ')' ';'
              { driver.initval_file();};

h_options_list: h_options_list COMMA h_options
               | h_options
               ;

h_options: o_filename
          | o_datafile
          | o_first_obs
          | o_data_first_obs
          | o_first_simulation_period
          | o_date_first_simulation_period
          | o_last_simulation_period
          | o_date_last_simulation_period
          | o_last_obs
          | o_data_last_obs
          | o_nobs
          | o_series
          ;

endval : ENDVAL ';' endval_list END ';'
         { driver.end_endval(false); }
       | ENDVAL '(' ALL_VALUES_REQUIRED ')' ';' endval_list END ';'
         { driver.end_endval(true); }
       | ENDVAL '(' LEARNT_IN EQUAL INT_NUMBER ')' ';' endval_list END ';'
         { driver.end_endval_learnt_in($5); }
       ;

endval_list : endval_list endval_elem
            | endval_elem
            ;

endval_elem : symbol EQUAL expression ';'
              { driver.end_val(EndValLearntInStatement::LearntEndValType::level, $1, $3); };
            | symbol PLUS_EQUAL expression ';'
              { driver.end_val(EndValLearntInStatement::LearntEndValType::add, $1, $3); };
            | symbol TIMES_EQUAL expression ';'
              { driver.end_val(EndValLearntInStatement::LearntEndValType::multiply, $1, $3); };
            ;

histval : HISTVAL ';' histval_list END ';'
          { driver.end_histval(false); };
        | HISTVAL '(' ALL_VALUES_REQUIRED ')' ';' histval_list END ';'
          { driver.end_histval(true); }
        ;

histval_list : histval_list histval_elem
             | histval_elem
             ;

histval_elem : symbol '(' signed_integer ')' EQUAL expression ';' { driver.hist_val($1, $3, $6); };

epilogue : EPILOGUE ';' { driver.begin_epilogue(); }
           epilogue_equation_list END ';' { driver.end_epilogue(); }
         ;

epilogue_equation_list : epilogue_equation_list epilogue_equation
                       | epilogue_equation
                       ;

epilogue_equation : NAME { driver.add_epilogue_variable($1); } EQUAL expression ';'
                    { driver.add_epilogue_equal($1, $4); }
                  ;

compilation_setup : COMPILATION_SETUP '(' compilation_setup_options_list ')' ';' { };

compilation_setup_options_list : compilation_setup_options_list COMMA compilation_setup_option
                               | compilation_setup_option
                               ;

compilation_setup_option : SUBSTITUTE_FLAGS EQUAL QUOTED_STRING
                           { driver.compilation_setup_substitute_flags($3); }
                         | ADD_FLAGS EQUAL QUOTED_STRING
                           { driver.compilation_setup_add_flags($3); }
                         | SUBSTITUTE_LIBS EQUAL QUOTED_STRING
                           { driver.compilation_setup_substitute_libs($3); }
                         | ADD_LIBS EQUAL QUOTED_STRING
                           { driver.compilation_setup_add_libs($3); }
                         | COMPILER EQUAL QUOTED_STRING
                           { driver.compilation_setup_compiler($3); }
                         ;

matched_moments : MATCHED_MOMENTS ';' { driver.begin_matched_moments(); }
                  matched_moments_list END ';' { driver.end_matched_moments($4); }
                ;

matched_moments_list : hand_side ';'
                       { $$ = { $1 }; }
                     | matched_moments_list hand_side ';'
                       {
                         $$ = $1;
                         $$.push_back($2);
                       }
                     ;

occbin_constraints : OCCBIN_CONSTRAINTS ';' { driver.begin_occbin_constraints(); }
                     occbin_constraints_regimes_list END ';' { driver.end_occbin_constraints($4); }
                   ;


occbin_constraints_regimes_list : occbin_constraints_regime
                                  { $$ = { $1 }; }
                                | occbin_constraints_regimes_list occbin_constraints_regime
                                  {
                                    $$ = $1;
                                    $$.push_back($2);
                                  }
                                ;

occbin_constraints_regime : NAME QUOTED_STRING ';' occbin_constraints_regime_options_list
                            {
                              BinaryOpNode *bind = dynamic_cast<BinaryOpNode *>($4["bind"]);
                              if (bind && !(bind->op_code == BinaryOpcode::less
                                            || bind->op_code == BinaryOpcode::greater
                                            || bind->op_code == BinaryOpcode::lessEqual
                                            || bind->op_code == BinaryOpcode::greaterEqual))
                                driver.error("The 'bind' expression must be an inequality constraint");
                              BinaryOpNode *relax = dynamic_cast<BinaryOpNode *>($4["relax"]);
                              if (relax && !(relax->op_code == BinaryOpcode::less
                                             || relax->op_code == BinaryOpcode::greater
                                             || relax->op_code == BinaryOpcode::lessEqual
                                             || relax->op_code == BinaryOpcode::greaterEqual))
                                driver.error("The 'relax' expression must be an inequality constraint");
                              $$ = { $2, bind, relax, $4["error_bind"], $4["error_relax"] };
                            }
                          ;

occbin_constraints_regime_options_list : occbin_constraints_regime_option
                                         { $$ = { $1 }; }
                                       | occbin_constraints_regime_options_list occbin_constraints_regime_option
                                         {
                                           $$ = $1;
                                           auto [it, success] = $$.insert($2);
                                           if (!success)
                                             driver.error("The '" + $2.first + "' clause is declared multiple times");
                                         }
                                       ;

occbin_constraints_regime_option : BIND hand_side ';'
                                   { $$ = { "bind", $2 }; }
                                 | RELAX hand_side ';'
                                   { $$ = { "relax", $2 }; }
                                 | ERROR_BIND hand_side ';'
                                   { $$ = { "error_bind", $2 }; }
                                 | ERROR_RELAX hand_side ';'
                                   { $$ = { "error_relax", $2 }; }
                                 ;

pac_target_info : PAC_TARGET_INFO '(' symbol ')' ';'
                  { driver.begin_pac_target_info($3); }
                  pac_target_info_statement_list
                  END ';'
                  { driver.end_pac_target_info(); }
                  ;

pac_target_info_statement_list : pac_target_info_statement
                               | pac_target_info_statement_list pac_target_info_statement
                               ;

pac_target_info_statement : TARGET hand_side ';'
                            { driver.set_pac_target_info_target($2); }
                          | AUXNAME_TARGET_NONSTATIONARY symbol ';'
                            { driver.set_pac_target_info_auxname_target_nonstationary($2); }
                          | pac_target_info_component
                          ;

pac_target_info_component : COMPONENT hand_side ';'
                            pac_target_info_component_list
                            { driver.add_pac_target_info_component($2); }
                          ;

pac_target_info_component_list : pac_target_info_component_elem
                               | pac_target_info_component_list pac_target_info_component_elem
                               ;

pac_target_info_component_elem : GROWTH hand_side ';'
                                 { driver.set_pac_target_info_component_growth($2); }
                               | AUXNAME symbol ';'
                                 { driver.set_pac_target_info_component_auxname($2); }
                               | KIND pac_target_kind ';'
                                 { driver.set_pac_target_info_component_kind($2); }
                               ;

pac_target_kind : LL
                  { $$ = PacTargetKind::ll; }
                | DL
                  { $$ = PacTargetKind::dl; }
                | DD
                  { $$ = PacTargetKind::dd; }
                ;

resid : RESID ';'
        { driver.resid(); }
      | RESID '(' o_non_zero ')' ';'
        { driver.resid(); }
      ;


/* The tokens below must be accepted in both DYNARE_STATEMENT and DYNARE_BLOCK
   states in the lexer, because of model block and model_options statement */
model_option : BLOCK { driver.block(); }
             | o_cutoff
             | o_mfs
             | BYTECODE { driver.bytecode(); }
             | USE_DLL { driver.use_dll(); }
             | NO_STATIC { driver.no_static();}
             | DIFFERENTIATE_FORWARD_VARS { driver.differentiate_forward_vars_all(); }
             | DIFFERENTIATE_FORWARD_VARS EQUAL '(' symbol_list ')' { driver.differentiate_forward_vars_some($4); }
             | o_linear
             | PARALLEL_LOCAL_FILES EQUAL '(' parallel_local_filename_list ')'
             | BALANCED_GROWTH_TEST_TOL EQUAL non_negative_number { driver.balanced_growth_test_tol($3); }
             ;

model_options_list : model_options_list COMMA model_option
                   | model_option
                   ;

model : MODEL ';' { driver.begin_model(); }
        equation_list END ';' { driver.end_model(); }
      | MODEL '(' model_options_list ')' ';' { driver.begin_model(); }
        equation_list END ';' { driver.end_model(); }
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
         | '[' tags_list ']' hand_side EQUAL hand_side ';'
           { $$ = driver.add_model_equal($4, $6); }
         | '[' tags_list ']' hand_side ';'
           { $$ = driver.add_model_equal_with_zero_rhs($4); }
         ;

tags_list : tags_list COMMA tag_pair
          | tag_pair
          ;

tag_pair : symbol EQUAL QUOTED_STRING
           { driver.add_equation_tags($1, $3); }
         | symbol
           { driver.add_equation_tags($1, ""); }
         ;

hand_side : '(' hand_side ')'
            { $$ = $2;}
          | namespace_qualified_symbol
            { $$ = driver.add_model_variable($1); }
          | symbol PIPE_E
            { $$ = driver.declare_or_change_type(SymbolType::endogenous, $1); }
          | symbol PIPE_X
            { $$ = driver.declare_or_change_type(SymbolType::exogenous, $1); }
          | symbol PIPE_P
            { $$ = driver.declare_or_change_type(SymbolType::parameter, $1); }
          | non_negative_number
            { $$ = driver.add_non_negative_constant($1); }
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
          | EXPECTATION '(' signed_integer ')''(' hand_side ')'
	    { $$ = driver.add_expectation($3, $6); }
          | VAR_EXPECTATION '(' symbol ')'
            { $$ = driver.add_var_expectation($3); }
          | PAC_EXPECTATION '(' symbol ')'
            { $$ = driver.add_pac_expectation($3); }
          | PAC_TARGET_NONSTATIONARY '(' symbol ')'
            { $$ = driver.add_pac_target_nonstationary($3); }
          | MINUS hand_side %prec UNARY
            { $$ = driver.add_uminus($2); }
          | PLUS hand_side %prec UNARY
            { $$ = $2; }
          | EXP '(' hand_side ')'
            { $$ = driver.add_exp($3); }
          | DIFF '(' hand_side ')'
            { $$ = driver.add_diff($3); }
          | ADL '(' hand_side COMMA QUOTED_STRING ')'
            { $$ = driver.add_adl($3, $5, "1"); }
          | ADL '(' hand_side COMMA QUOTED_STRING COMMA INT_NUMBER ')'
            { $$ = driver.add_adl($3, $5, $7); }
          | ADL '(' hand_side COMMA QUOTED_STRING COMMA vec_int ')'
            { $$ = driver.add_adl($3, $5, $7); }
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
          | SINH '(' hand_side ')'
            { $$ = driver.add_sinh($3); }
          | COSH '(' hand_side ')'
            { $$ = driver.add_cosh($3); }
          | TANH '(' hand_side ')'
            { $$ = driver.add_tanh($3); }
          | ASINH '(' hand_side ')'
            { $$ = driver.add_asinh($3); }
          | ACOSH '(' hand_side ')'
            { $$ = driver.add_acosh($3); }
          | ATANH '(' hand_side ')'
            { $$ = driver.add_atanh($3); }
          | SQRT '(' hand_side ')'
            { $$ = driver.add_sqrt($3); }
          | CBRT '(' hand_side ')'
             { $$ = driver.add_cbrt($3); }
          | ABS '(' hand_side ')'
            { $$ = driver.add_abs($3); }
          | SIGN '(' hand_side ')'
            { $$ = driver.add_sign($3); }
          | MAX '(' hand_side COMMA hand_side ')'
            { $$ = driver.add_max($3, $5); }
          | MIN '(' hand_side COMMA hand_side ')'
            { $$ = driver.add_min($3, $5); }
          | namespace_qualified_symbol { driver.push_external_function_arg_vector_onto_stack(); } '(' comma_hand_side ')'
            { $$ = driver.add_model_var_or_external_function($1, true); }
          | NORMCDF '(' hand_side COMMA hand_side COMMA hand_side ')'
            { $$ = driver.add_normcdf($3, $5, $7); }
          | NORMCDF '(' hand_side ')'
            { $$ = driver.add_normcdf($3); }
          | NORMPDF '(' hand_side COMMA hand_side COMMA hand_side ')'
            { $$ = driver.add_normpdf($3, $5, $7); }
          | NORMPDF '(' hand_side ')'
            { $$ = driver.add_normpdf($3); }
          | ERF '(' hand_side ')'
            { $$ = driver.add_erf($3); }
          | ERFC '(' hand_side ')'
            { $$ = driver.add_erfc($3); }
          | STEADY_STATE '(' hand_side ')'
            { $$ = driver.add_steady_state($3); }
          ;

comma_hand_side : hand_side
                  { driver.add_external_function_arg($1); }
                | comma_hand_side COMMA hand_side
                  { driver.add_external_function_arg($3); }
                ;

pound_expression: '#' symbol EQUAL hand_side ';'
                  { driver.declare_and_init_model_local_variable($2, $4); };

model_remove : MODEL_REMOVE '(' tag_pair_list_for_selection ')' ';'
               { driver.model_remove($3); };

model_replace : MODEL_REPLACE '(' tag_pair_list_for_selection ')' ';'
               { driver.begin_model_replace($3); }
               equation_list END ';'
               { driver.end_model(); };

model_options : MODEL_OPTIONS '(' model_options_list ')' ';'

tag_pair_list_for_selection : QUOTED_STRING
                              { $$ = { { "name", $1 } }; }
                            | symbol EQUAL QUOTED_STRING
                              { $$ = { { $1, $3 } }; }
                            | tag_pair_list_for_selection COMMA QUOTED_STRING
                              {
                                $$ = $1;
                                $$.emplace_back("name", $3);
                              }
                            | tag_pair_list_for_selection COMMA symbol EQUAL QUOTED_STRING
                              {
                                $$ = $1;
                                $$.emplace_back($3, $5);
                              }
                            ;

shocks : SHOCKS ';' shock_list END ';' { driver.end_shocks(false); }
       | SHOCKS '(' OVERWRITE ')' ';' shock_list END ';' { driver.end_shocks(true); }
       | SHOCKS '(' OVERWRITE ')' ';'  END ';' { driver.end_shocks(true); }
       | SHOCKS '(' SURPRISE ')' ';' det_shock_list END ';' { driver.end_shocks_surprise(false); }
       | SHOCKS '(' SURPRISE COMMA OVERWRITE ')' ';' det_shock_list END ';' { driver.end_shocks_surprise(true); }
       | SHOCKS '(' OVERWRITE COMMA SURPRISE ')' ';' det_shock_list END ';' { driver.end_shocks_surprise(true); }
       | SHOCKS '(' LEARNT_IN EQUAL INT_NUMBER ')' ';' det_shock_list END ';' { driver.end_shocks_learnt_in($5, false); }
       | SHOCKS '(' LEARNT_IN EQUAL INT_NUMBER COMMA OVERWRITE ')' ';' det_shock_list END ';' { driver.end_shocks_learnt_in($5, true); }
       | SHOCKS '(' OVERWRITE COMMA LEARNT_IN EQUAL INT_NUMBER ')' ';' det_shock_list END ';' { driver.end_shocks_learnt_in($7, true); }
       ;

shock_list : shock_list shock_elem
           | shock_elem
           ;

shock_elem : det_shock_elem
           | VAR symbol ';' STDERR expression ';'
             { driver.add_stderr_shock($2, $5); }
           | VAR symbol EQUAL expression ';'
             { driver.add_var_shock($2, $4); }
           | VAR symbol COMMA symbol EQUAL expression ';'
             { driver.add_covar_shock($2, $4, $6); }
           | CORR symbol COMMA symbol EQUAL expression ';'
             { driver.add_correl_shock($2, $4, $6); }
           ;

det_shock_elem : VAR symbol ';' PERIODS period_list ';' VALUES value_list ';'
                 { driver.add_det_shock($2, $5, $8, ParsingDriver::DetShockType::standard); }
               | VAR symbol ';' PERIODS period_list ';' ADD value_list ';'
                 { driver.add_det_shock($2, $5, $8, ParsingDriver::DetShockType::add); }
               | VAR symbol ';' PERIODS period_list ';' MULTIPLY value_list ';'
                 { driver.add_det_shock($2, $5, $8, ParsingDriver::DetShockType::multiply); }
               ;

det_shock_list : det_shock_list det_shock_elem
               | det_shock_elem
               ;

heteroskedastic_shocks : HETEROSKEDASTIC_SHOCKS ';' heteroskedastic_shock_list END ';'
                         { driver.end_heteroskedastic_shocks(false); }
                       | HETEROSKEDASTIC_SHOCKS '(' OVERWRITE ')' ';' heteroskedastic_shock_list END ';'
                         { driver.end_heteroskedastic_shocks(true); }
                       | HETEROSKEDASTIC_SHOCKS '(' OVERWRITE ')' ';'  END ';'
                         { driver.end_heteroskedastic_shocks(true); }
                       ;

heteroskedastic_shock_list : heteroskedastic_shock_list heteroskedastic_shock_elem
                           | heteroskedastic_shock_elem
                           ;

heteroskedastic_shock_elem : VAR symbol ';' PERIODS period_list ';' VALUES value_list ';'
                             { driver.add_heteroskedastic_shock($2, $5, $8, false); }
                           | VAR symbol ';' PERIODS period_list ';' SCALES value_list ';'
                             { driver.add_heteroskedastic_shock($2, $5, $8, true); }
                           ;

svar_identification : SVAR_IDENTIFICATION {driver.begin_svar_identification();} ';' svar_identification_list END ';'
                      { driver.end_svar_identification(); }
                    ;

svar_identification_list : svar_identification_list svar_identification_elem
                         | svar_identification_elem
                         ;

svar_identification_elem : EXCLUSION LAG INT_NUMBER ';' svar_equation_list
                           { driver.combine_lag_and_restriction($3); }
                         | EXCLUSION CONSTANTS ';'
                           { driver.add_constants_exclusion(); }
                         | RESTRICTION EQUATION INT_NUMBER COMMA
			 { driver.add_restriction_equation_nbr($3);}
                           restriction_expression EQUAL
                                {driver.add_restriction_equal();}
                           restriction_expression ';'
                         | UPPER_CHOLESKY ';'
                           { driver.add_upper_cholesky(); }
                         | LOWER_CHOLESKY ';'
                           { driver.add_lower_cholesky(); }
                         ;

svar_equation_list : svar_equation_list EQUATION INT_NUMBER COMMA symbol_list ';'
                     { driver.add_restriction_in_equation($3, $5); }
                   | EQUATION INT_NUMBER COMMA symbol_list ';'
                     { driver.add_restriction_in_equation($2, $4); }
                   ;

restriction_expression : expression {driver.check_restriction_expression_constant($1);}
                       | restriction_expression_1
                       ;

restriction_expression_1 : restriction_elem_expression
                         | restriction_expression_1 restriction_elem_expression
                         ;

restriction_elem_expression : COEFF '(' symbol COMMA INT_NUMBER ')'
                                 { driver.add_positive_restriction_element($3,$5);}
                            | PLUS COEFF '(' symbol COMMA INT_NUMBER ')'
		                 { driver.add_positive_restriction_element($4,$6);}
                            | MINUS COEFF '(' symbol COMMA INT_NUMBER ')'
		                 { driver.add_negative_restriction_element($4,$6);}
                            | expression TIMES COEFF '(' symbol COMMA INT_NUMBER ')'
                                 { driver.add_positive_restriction_element($1,$5,$7);}
                            ;

svar_global_identification_check: SVAR_GLOBAL_IDENTIFICATION_CHECK ';'
                                  { driver.add_svar_global_identification_check(); }
                                ;

markov_switching : MARKOV_SWITCHING '(' ms_options_list ')' ';'
                   { driver.markov_switching(); }
                 ;

ms_options_list : ms_options_list COMMA ms_options
                | ms_options
                ;

ms_options : o_chain
           | o_duration
           | o_restrictions
           | o_number_of_regimes
           | o_number_of_lags
           | o_parameters
           ;

svar : SVAR '(' svar_options_list ')' ';'
       { driver.svar(); }
     ;

svar_options_list : svar_options_list COMMA svar_options
                  | svar_options
                  ;

svar_options : o_coefficients
             | o_variances
             | o_equations
             | o_chain
             ;

mshocks : MSHOCKS ';' mshock_list END ';' { driver.end_mshocks(false); }
        | MSHOCKS '(' OVERWRITE ')' ';' mshock_list END ';' { driver.end_mshocks(true); }
        ;

mshock_list : mshock_list det_shock_elem
            | det_shock_elem
            ;

period_list : period_list COMMA INT_NUMBER
              {
                $$ = $1;
                int p = stoi($3);
                $$.emplace_back(p, p);
              }
            | period_list INT_NUMBER
              {
                $$ = $1;
                int p = stoi($2);
                $$.emplace_back(p, p);
              }
            | period_list COMMA INT_NUMBER ':' INT_NUMBER
              {
                $$ = $1;
                int p1 = stoi($3), p2 = stoi($5);
                if (p1 > p2)
                  driver.error("Can't have first period index greater than second index in range specification");
                $$.emplace_back(p1, p2);
              }
            | period_list INT_NUMBER ':' INT_NUMBER
              {
                $$ = $1;
                int p1 = stoi($2), p2 = stoi($4);
                if (p1 > p2)
                  driver.error("Can't have first period index greater than second index in range specification");
                $$.emplace_back(p1, p2);
              }
            | INT_NUMBER ':' INT_NUMBER
              {
                int p1 = stoi($1), p2 = stoi($3);
                if (p1 > p2)
                  driver.error("Can't have first period index greater than second index in range specification");
                $$ = { { p1, p2 } };
              }
            | INT_NUMBER
              {
                int p = stoi($1);
                $$ = { { p, p } };
              }
            ;

value_list : value_list COMMA '(' expression ')'
             {
               $$ = $1;
               $$.push_back($4);
             }
           | value_list '(' expression ')'
             {
               $$ = $1;
               $$.push_back($3);
             }
           | '(' expression ')'
             { $$ = { $2 }; }
           | value_list COMMA signed_number
             {
               $$ = $1;
               $$.push_back($3.at(0) == '-' ?
                            driver.add_uminus(driver.add_non_negative_constant($3.substr(1))) :
                            driver.add_non_negative_constant($3));
             }
           | value_list signed_number
             {
               $$ = $1;
               $$.push_back($2.at(0) == '-' ?
                            driver.add_uminus(driver.add_non_negative_constant($2.substr(1))) :
                            driver.add_non_negative_constant($2));
             }
           | signed_number
             {
               $$ = { $1.at(0) == '-' ?
                      driver.add_uminus(driver.add_non_negative_constant($1.substr(1))) :
                      driver.add_non_negative_constant($1) };
             }
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
               | o_homotopy_force_continue
               | o_markowitz
               | o_steady_maxit
               | o_nocheck
               | o_steady_tolf
               | o_steady_tolx
               | o_fsolve_options
               ;

check : CHECK ';'
        { driver.check(); }
      | CHECK '(' check_options_list ')' ';'
        { driver.check(); }
      ;

check_options_list : check_options_list COMMA check_options
                   | check_options
                   ;

check_options : steady_options
	      | o_qz_zero_threshold
	      ;

model_info : MODEL_INFO ';'
             { driver.model_info(); }
            | MODEL_INFO '(' model_info_options_list ')' ';'
              { driver.model_info(); }
           ;

model_info_options_list : model_info_options_list COMMA model_info_options
                   | model_info_options
                   ;

model_info_options : o_block_static
                   | o_block_dynamic
                   | o_incidence
                   ;

perfect_foresight_setup : PERFECT_FORESIGHT_SETUP ';'
                          { driver.perfect_foresight_setup(); }
                        | PERFECT_FORESIGHT_SETUP '(' perfect_foresight_setup_options_list ')' ';'
                          { driver.perfect_foresight_setup(); }
                        ;

perfect_foresight_setup_options_list : perfect_foresight_setup_options_list COMMA perfect_foresight_setup_options
                                     | perfect_foresight_setup_options
                                     ;

perfect_foresight_setup_options : o_periods
                                | o_datafile
                                | o_endval_steady
                                ;

perfect_foresight_solver : PERFECT_FORESIGHT_SOLVER ';'
                          { driver.perfect_foresight_solver(); }
                         | PERFECT_FORESIGHT_SOLVER '(' perfect_foresight_solver_options_list ')' ';'
                          { driver.perfect_foresight_solver(); }
                         ;

perfect_foresight_solver_options_list : perfect_foresight_solver_options_list COMMA perfect_foresight_solver_options
                                     | perfect_foresight_solver_options
                                     ;

perfect_foresight_solver_options : o_stack_solve_algo
                                 | o_markowitz
                                 | o_minimal_solving_periods
                                 | o_simul_maxit
                                 | o_endogenous_terminal_period
                                 | o_linear_approximation
                                 | o_no_homotopy
                                 | o_solve_algo
                                 | o_robust_lin_solve
                                 | o_lmmcp
                                 | o_pf_tolf
                                 | o_pf_tolx
                                 | o_noprint
                                 | o_print
                                 | o_pf_steady_solve_algo
                                 | o_pf_steady_maxit
                                 | o_pf_steady_tolf
                                 | o_pf_steady_tolx
                                 | o_pf_steady_markowitz
                                 | o_homotopy_max_completion_share
                                 | o_homotopy_min_step_size
                                 | o_homotopy_initial_step_size
                                 | o_homotopy_step_size_increase_success_count
                                 | o_homotopy_linearization_fallback
                                 | o_homotopy_marginal_linearization_fallback
                                 ;

perfect_foresight_with_expectation_errors_setup : PERFECT_FORESIGHT_WITH_EXPECTATION_ERRORS_SETUP ';'
                                                  { driver.perfect_foresight_with_expectation_errors_setup(); }
                                                | PERFECT_FORESIGHT_WITH_EXPECTATION_ERRORS_SETUP '(' perfect_foresight_with_expectation_errors_setup_options_list ')' ';'
                                                  { driver.perfect_foresight_with_expectation_errors_setup(); }
                                                ;

perfect_foresight_with_expectation_errors_setup_options_list : perfect_foresight_with_expectation_errors_setup_options_list COMMA perfect_foresight_with_expectation_errors_setup_options
                                                             | perfect_foresight_with_expectation_errors_setup_options
                                                             ;

perfect_foresight_with_expectation_errors_setup_options : o_periods
                                                        | o_datafile
                                                        ;

perfect_foresight_with_expectation_errors_solver : PERFECT_FORESIGHT_WITH_EXPECTATION_ERRORS_SOLVER ';'
                                                  { driver.perfect_foresight_with_expectation_errors_solver(); }
                                                 | PERFECT_FORESIGHT_WITH_EXPECTATION_ERRORS_SOLVER '(' perfect_foresight_with_expectation_errors_solver_options_list ')' ';'
                                                  { driver.perfect_foresight_with_expectation_errors_solver(); }
                                                 ;

perfect_foresight_with_expectation_errors_solver_options_list : perfect_foresight_with_expectation_errors_solver_options_list COMMA perfect_foresight_with_expectation_errors_solver_options
                                                              | perfect_foresight_with_expectation_errors_solver_options
                                                              ;

perfect_foresight_with_expectation_errors_solver_options : o_pfwee_constant_simulation_length
                                                         | perfect_foresight_solver_options
                                                         ;

method_of_moments : METHOD_OF_MOMENTS ';'
                    { driver.method_of_moments(); }
                  | METHOD_OF_MOMENTS '(' method_of_moments_options_list ')' ';'
                    { driver.method_of_moments(); }
                  ;

method_of_moments_options_list : method_of_moments_option COMMA method_of_moments_options_list
                               | method_of_moments_option
                               ;

method_of_moments_option : o_add_tiny_number_to_cholesky
                         | o_additional_optimizer_steps
                         | o_aim_solver
                         | o_analytic_jacobian
                         | o_analytic_standard_errors
                         | o_bartlett_kernel_lag
                         | o_bounded_shock_support
                         | o_brooks_gelman_plotrows
                         | o_cova_compute
                         | o_datafile
                         | o_dirname
                         | o_dr
                         | o_dr_cycle_reduction_tol
                         | o_dr_logarithmic_reduction_maxiter
                         | o_dr_logarithmic_reduction_tol
                         | o_drop
                         | o_first_obs
                         | o_geweke_interval
                         | o_graph_format
                         | o_huge_number
                         | o_irf_matching_file
                         | o_k_order_solver
                         | o_load_mh_file
                         | o_load_results_after_load_mh
                         | o_logdata
                         | o_lyapunov
                         | o_lyapunov_complex_threshold
                         | o_lyapunov_doubling_tol
                         | o_lyapunov_fixed_point_tol
                         | o_mcmc_jumping_covariance
                         | o_mh_conf_sig
                         | o_mh_drop
                         | o_mh_init_scale_factor
                         | o_mh_initialize_from_previous_mcmc
                         | o_mh_initialize_from_previous_mcmc_directory
                         | o_mh_initialize_from_previous_mcmc_prior
                         | o_mh_initialize_from_previous_mcmc_record
                         | o_mh_jscale
                         | o_mh_nblocks
                         | o_mh_posterior_mode_estimation
                         | o_mh_recover
                         | o_mh_replic
                         | o_mh_tune_guess
                         | o_mh_tune_jscale
                         | o_mom_burnin
                         | o_mom_method
                         | o_mom_seed
                         | o_mom_se_tolx
                         | o_mode_check
                         | o_mode_check_neighbourhood_size
                         | o_mode_check_number_of_points
                         | o_mode_check_symmetric_plots
                         | o_mode_compute
                         | o_mode_file
                         | o_nobs
                         | o_no_posterior_kernel_density
                         | o_nodiagnostic
                         | o_nodisplay
                         | o_nograph
                         | o_noprint
                         | o_optim
                         | o_order
                         | o_penalized_estimator
                         | o_plot_priors
                         | o_posterior_max_subsample_draws
                         | o_posterior_sampler_options
                         | o_posterior_sampling_method
                         | o_prefilter
                         | o_prior_trunc
                         | o_pruning
                         | o_qz_criterium
                         | o_qz_zero_threshold
                         | o_raftery_lewis_diagnostics
                         | o_raftery_lewis_qrs
                         | o_relative_irf
                         | o_replic
                         | o_schur_vec_tol
                         | o_silent_optimizer
                         | o_simulation_method
                         | o_simulation_multiple
                         | o_sub_draws
                         | o_sylvester
                         | o_sylvester_fixed_point_tol
                         | o_taper_steps
                         | o_tex
                         | o_use_penalized_objective_for_hessian
                         | o_verbose
                         | o_weighting_matrix
                         | o_weighting_matrix_scaling_factor
                         | o_xls_range
                         | o_xls_sheet
                         ;

prior_function : PRIOR_FUNCTION '(' prior_posterior_function_options_list ')' ';'
                { driver.prior_posterior_function(true); }
               ;

posterior_function : POSTERIOR_FUNCTION '(' prior_posterior_function_options_list ')' ';'
                    { driver.prior_posterior_function(false); }
                   ;

prior_posterior_function_options_list : prior_posterior_function_options_list COMMA prior_posterior_function_options
                                      | prior_posterior_function_options
                                      ;

prior_posterior_function_options : o_function
                                 | o_sampling_draws
                                 ;

simul : SIMUL ';'
        { driver.simul(); }
      | SIMUL '(' simul_options_list ')' ';'
        { driver.simul(); }
      ;

simul_options_list : simul_options_list COMMA simul_options
                   | simul_options
                   ;

simul_options : perfect_foresight_setup_options
              | perfect_foresight_solver_options
              ;

external_function : EXTERNAL_FUNCTION '(' external_function_options_list ')' ';'
                    { driver.external_function(); }
                  ;

external_function_options_list : external_function_options_list COMMA external_function_options
                               | external_function_options
                               ;

external_function_options : o_ext_func_name
                          | o_ext_func_nargs
                          | o_first_deriv_provided
                          | o_second_deriv_provided
                          ;

stoch_simul : STOCH_SIMUL ';'
              { driver.stoch_simul({}); }
            | STOCH_SIMUL '(' stoch_simul_options_list ')' ';'
              { driver.stoch_simul({}); }
            | STOCH_SIMUL symbol_list ';'
              { driver.stoch_simul($2); }
            | STOCH_SIMUL '(' stoch_simul_options_list ')' symbol_list ';'
              { driver.stoch_simul($5); }
            ;

stoch_simul_options_list : stoch_simul_options_list COMMA stoch_simul_options
                         | stoch_simul_options
                         ;

stoch_simul_primary_options : o_dr_algo
                            | o_solve_algo
                            | o_simul_algo
                            | o_order
                            | o_replic
                            | o_drop
                            | o_ar
                            | o_nocorr
                            | o_contemporaneous_correlation
                            | o_nofunctions
                            | o_nomoments
                            | o_nograph
                            | o_nodisplay
                            | o_graph_format
                            | o_irf
                            | o_irf_shocks
                            | o_relative_irf
                            | o_analytical_girf
                            | o_irf_in_percent
                            | o_emas_girf
                            | o_emas_drop
                            | o_emas_tolf
                            | o_emas_max_iter
                            | o_stderr_multiples
                            | o_diagonal_only
                            | o_hp_filter
                            | o_hp_ngrid
                            | o_filtered_theoretical_moments_grid
                            | o_periods
                            | o_simul
                            | o_simul_seed
                            | o_simul_replic
                            | o_qz_criterium
                            | o_qz_zero_threshold
                            | o_print
                            | o_noprint
                            | o_aim_solver
                            | o_partial_information
                            | o_conditional_variance_decomposition
                            | o_k_order_solver
                            | o_pruning
                            | o_sylvester
                            | o_sylvester_fixed_point_tol
                            | o_dr
                            | o_dr_cycle_reduction_tol
                            | o_dr_logarithmic_reduction_tol
                            | o_dr_logarithmic_reduction_maxiter
                            | o_irf_plot_threshold
                            | o_dr_display_tol
                            | o_tex
                            ;

stoch_simul_options : stoch_simul_primary_options
                    | o_loglinear
                    | o_nodecomposition
                    | o_spectral_density
                    | o_bandpass_filter
                    | o_one_sided_hp_filter
                    ;

signed_integer : PLUS INT_NUMBER
                 { $$ = $2; }
               | MINUS INT_NUMBER
                 {
                   $$ = $2;
                   $$.insert(0, "-");
                 }
               | INT_NUMBER
               ;

non_negative_number : INT_NUMBER
                    | FLOAT_NUMBER
                    ;

signed_number : PLUS non_negative_number
               { $$ = $2; }
              | MINUS non_negative_number
               {
                 $$ = $2;
                 $$.insert(0, "-");
               }
              | non_negative_number
              ;

signed_inf : PLUS INF_CONSTANT
             { $$ = "Inf"; }
           | MINUS INF_CONSTANT
             { $$ = "-Inf"; }
           | INF_CONSTANT
             { $$ = "Inf"; }
           ;

signed_number_w_inf : signed_inf
                    | signed_number
                    ;

boolean : TRUE
        | FALSE
        ;

estimated_params : ESTIMATED_PARAMS ';' estimated_list END ';'
                   { driver.estimated_params(false); }
                 | ESTIMATED_PARAMS '(' OVERWRITE ')' ';' estimated_list END ';'
                   { driver.estimated_params(true); }
                 ;

estimated_list : estimated_list estimated_elem
                 { driver.add_estimated_params_element(); }
               | estimated_elem
                 { driver.add_estimated_params_element(); }
               ;

estimated_elem : estimated_elem1 COMMA estimated_elem2 ';'
               | estimated_elem1 ';'

estimated_elem1 : STDERR symbol
                  {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = $2;
                  }
                | symbol
                  {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = $1;
                  }
                | CORR symbol COMMA symbol
                  {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = $2;
                    driver.estim_params.name2 = $4;
                  }
                | DSGE_PRIOR_WEIGHT
                  {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = "dsge_prior_weight";
                  }
                ;

estimated_elem2 : prior_pdf COMMA estimated_elem3
                  {
                    driver.estim_params.prior = $1;
                  }
                | expression_or_empty COMMA prior_pdf COMMA estimated_elem3
                  {
                    driver.estim_params.init_val = $1;
                    driver.estim_params.prior = $3;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA prior_pdf COMMA estimated_elem3
                  {
                    driver.estim_params.init_val = $1;
                    driver.estim_params.low_bound = $3;
                    driver.estim_params.up_bound = $5;
                    driver.estim_params.prior = $7;
                  }
                | expression
                  {
                    driver.estim_params.init_val = $1;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.init_val = $1;
                    driver.estim_params.low_bound = $3;
                    driver.estim_params.up_bound = $5;
                  }
                ;

estimated_elem3 : expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                    driver.estim_params.p3 = $5;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                    driver.estim_params.p3 = $5;
                    driver.estim_params.p4 = $7;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA expression
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                    driver.estim_params.p3 = $5;
                    driver.estim_params.p4 = $7;
                    driver.estim_params.jscale = $9;
                  }
                ;

estimated_params_init : ESTIMATED_PARAMS_INIT ';' estimated_init_list END ';'
                        { driver.estimated_params_init(); }
                      | ESTIMATED_PARAMS_INIT '(' USE_CALIBRATION ')' ';' END ';'
                        { driver.estimated_params_init(true); }
                      | ESTIMATED_PARAMS_INIT '(' USE_CALIBRATION ')' ';' estimated_init_list END ';'
                        { driver.estimated_params_init(true); }
                      ;

estimated_init_list : estimated_init_list estimated_init_elem
                      { driver.add_estimated_params_element(); }
                    | estimated_init_elem
                      { driver.add_estimated_params_element(); }
                    ;

estimated_init_elem : STDERR symbol COMMA expression ';'
                      {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = $2;
                        driver.estim_params.init_val = $4;
                      }
                    | CORR symbol COMMA symbol COMMA expression ';'
                      {
                        driver.estim_params.type = 3;
                        driver.estim_params.name = $2;
                        driver.estim_params.name2 = $4;
                        driver.estim_params.init_val = $6;
                      }
                    | symbol COMMA expression ';'
                      {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = $1;
                        driver.estim_params.init_val = $3;
                      }
                    ;

estimated_params_bounds : ESTIMATED_PARAMS_BOUNDS ';' estimated_bounds_list END ';'
                          { driver.estimated_params_bounds(); };

estimated_bounds_list : estimated_bounds_list estimated_bounds_elem
                        { driver.add_estimated_params_element(); }
                      | estimated_bounds_elem
                        { driver.add_estimated_params_element(); }
                      ;

estimated_bounds_elem : STDERR symbol COMMA expression COMMA expression ';'
                        {
                          driver.estim_params.type = 1;
                          driver.estim_params.name = $2;
                          driver.estim_params.low_bound = $4;
                          driver.estim_params.up_bound = $6;
                        }
                      | CORR symbol COMMA symbol COMMA expression COMMA expression ';'
                        {
                          driver.estim_params.type = 3;
                          driver.estim_params.name = $2;
                          driver.estim_params.name2 = $4;
                          driver.estim_params.low_bound = $6;
                          driver.estim_params.up_bound = $8;
                        }
                      | symbol COMMA expression COMMA expression ';'
                        {
                          driver.estim_params.type = 2;
                          driver.estim_params.name = $1;
                          driver.estim_params.low_bound = $3;
                          driver.estim_params.up_bound = $5;
                        }
                      ;

estimated_params_remove : ESTIMATED_PARAMS_REMOVE ';' estimated_remove_list END ';'
                          { driver.estimated_params_remove(); };

estimated_remove_list : estimated_remove_elem
                      | estimated_remove_list estimated_remove_elem
                      ;

estimated_remove_elem : estimated_elem1 ';'
                        { driver.add_estimated_params_element(); }

osr_params_bounds : OSR_PARAMS_BOUNDS ';' osr_bounds_list END ';'
                    { driver.osr_params_bounds(); };

osr_bounds_list : osr_bounds_list osr_bounds_elem
                  { driver.add_osr_params_element(); }
                | osr_bounds_elem
                  { driver.add_osr_params_element(); }
                ;

osr_bounds_elem : symbol COMMA expression COMMA expression ';'
                  {
                    driver.osr_params.name = $1;
                    driver.osr_params.low_bound = $3;
                    driver.osr_params.up_bound = $5;
                  }
                ;

prior_distribution : BETA
                     { $$ = PriorDistributions::beta; }
                   | GAMMA
                     { $$ = PriorDistributions::gamma; }
                   | NORMAL
                     { $$ = PriorDistributions::normal; }
                   | INV_GAMMA
                     { $$ = PriorDistributions::invGamma; }
                   | INV_GAMMA1
                     { $$ = PriorDistributions::invGamma1; }
                   | UNIFORM
                     { $$ = PriorDistributions::uniform; }
                   | INV_GAMMA2
                     { $$ = PriorDistributions::invGamma2; }
                   | DIRICHLET
                     { $$ = PriorDistributions::dirichlet; }
                   | WEIBULL
                     { $$ = PriorDistributions::weibull; }
                   ;

prior_pdf : BETA_PDF
            { $$ = PriorDistributions::beta; }
          | GAMMA_PDF
            { $$ = PriorDistributions::gamma; }
          | NORMAL_PDF
            { $$ = PriorDistributions::normal; }
          | INV_GAMMA_PDF
            { $$ = PriorDistributions::invGamma; }
          | INV_GAMMA1_PDF
            { $$ = PriorDistributions::invGamma1; }
          | UNIFORM_PDF
            { $$ = PriorDistributions::uniform; }
          | INV_GAMMA2_PDF
            { $$ = PriorDistributions::invGamma2; }
          | WEIBULL_PDF
            { $$ = PriorDistributions::weibull; }
          ;

date_str : DATES

date_expr : date_str
          | date_expr PLUS INT_NUMBER
            { $$ = $1 + '+' + $3; }
          ;

set_time : SET_TIME '(' date_expr ')' ';'
           { driver.set_time($3); }
         ;

data : DATA '(' data_options_list ')'';'
       { driver.estimation_data(); }
     ;

data_options_list : data_options_list COMMA data_options
                  | data_options
                  ;

data_options : o_file
             | o_series
             | o_data_first_obs
             | o_data_last_obs
             | o_data_nobs
             | o_xls_sheet
             | o_xls_range
             ;

subsamples : subsamples_eq_opt '(' subsamples_name_list ')' ';'
             { driver.set_subsamples($1.first, $1.second); }
           ;

subsamples_eq : subsamples_eq_opt EQUAL subsamples_eq_opt ';'
                { driver.copy_subsamples($1.first, $1.second, $3.first, $3.second); }
              ;

subsamples_eq_opt : symbol '.' SUBSAMPLES
                    { $$ = { $1, "" }; }
                  | STD '(' symbol ')' '.' SUBSAMPLES
                    { $$ = { $3, "" }; }
                  | CORR '(' symbol COMMA symbol ')' '.' SUBSAMPLES
                    { $$ = { $3, $5 }; }
                  ;

subsamples_name_list : subsamples_name_list COMMA o_subsample_name
                     | o_subsample_name
                     ;

prior : symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = PriorDistributions::noShape; } '(' prior_options_list ')' ';'
        { driver.set_prior($1, ""); }
      | symbol '.' symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = PriorDistributions::noShape; } '(' prior_options_list ')' ';'
        { driver.set_prior($1, $3); }
      | SYMBOL_VEC '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = PriorDistributions::noShape; }  '(' joint_prior_options_list ')' ';'
        { driver.set_joint_prior($1); }
      | STD '(' symbol ')' '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = PriorDistributions::noShape; } '(' prior_options_list ')' ';'
        { driver.set_std_prior($3, ""); }
      | STD '(' symbol ')' '.' symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = PriorDistributions::noShape; } '(' prior_options_list ')' ';'
        { driver.set_std_prior($3, $6); }
      | CORR '(' symbol COMMA symbol ')' '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = PriorDistributions::noShape; } '(' prior_options_list ')' ';'
        { driver.set_corr_prior($3, $5, ""); }
      | CORR '(' symbol COMMA symbol ')' '.' symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = PriorDistributions::noShape; } '(' prior_options_list ')' ';'
        { driver.set_corr_prior($3, $5, $8); }
      ;

prior_options_list : prior_options_list COMMA prior_options
                   | prior_options
                   ;

prior_options : o_shift
              | o_mean
              | o_median
              | o_stdev
              | o_truncate
              | o_variance
              | o_mode
              | o_interval
              | o_shape
              | o_domain
              ;

joint_prior_options_list : joint_prior_options_list COMMA joint_prior_options
                         | joint_prior_options
                         ;

joint_prior_options : o_shift
                    | o_mean_vec
                    | o_median
                    | o_stdev
                    | o_truncate
                    | o_variance_mat
                    | o_mode
                    | o_interval
                    | o_shape
                    | o_domain
                    ;

prior_eq : prior_eq_opt EQUAL prior_eq_opt ';'
           { driver.copy_prior(get<0>($1), get<1>($1), get<2>($1), get<3>($1),
                               get<0>($3), get<1>($3), get<2>($3), get<3>($3)); }
         ;

prior_eq_opt : symbol '.' PRIOR
               { $$ = { "par", $1, "", "" }; }
             | symbol '.' symbol '.' PRIOR
               { $$ = { "par", $1, "", $3 }; }
             | STD '(' symbol ')' '.'  PRIOR
               { $$ = { "std", $3, "", "" }; }
             | STD '(' symbol ')' '.' symbol '.' PRIOR
               { $$ = { "std", $3, "", $6 }; }
             | CORR '(' symbol COMMA symbol ')' '.'  PRIOR
               { $$ = { "corr", $3, $5, "" }; }
             | CORR '(' symbol COMMA symbol ')' '.' symbol '.' PRIOR
               { $$ = { "corr", $3, $5, $8 }; }
             ;

options : symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_options($1, ""); }
        | symbol '.' symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_options($1, $3); }
        | STD '(' symbol ')' '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_std_options($3, ""); }
        | STD '(' symbol ')' '.' symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_std_options($3, $6); }
        | CORR '(' symbol COMMA symbol ')' '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_corr_options($3, $5, ""); }
        | CORR '(' symbol COMMA symbol ')' '.' symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_corr_options($3, $5, $8); }
        ;

options_options_list : options_options_list COMMA options_options
                     | options_options
                     ;

options_options : o_jscale
                | o_init
                | o_bounds
                ;

options_eq : options_eq_opt EQUAL options_eq_opt ';'
             { driver.copy_options(get<0>($1), get<1>($1), get<2>($1), get<3>($1),
                                   get<0>($3), get<1>($3), get<2>($3), get<3>($3)); }
           ;

options_eq_opt : symbol '.' OPTIONS
                 { $$ = { "par", $1, "", "" }; }
               | symbol '.' symbol '.' OPTIONS
                 { $$ = { "par", $1, "", $3 }; }
               | STD '(' symbol ')' '.'  OPTIONS
                 { $$ = { "std", $3, "", "" }; }
               | STD '(' symbol ')' '.' symbol '.' OPTIONS
                 { $$ = { "std", $3, "", $6 }; }
               | CORR '(' symbol COMMA symbol ')' '.'  OPTIONS
                 { $$ = { "corr", $3, $5, "" }; }
               | CORR '(' symbol COMMA symbol ')' '.' symbol '.' OPTIONS
                 { $$ = { "corr", $3, $5, $8 }; }
               ;

estimation : ESTIMATION ';'
             { driver.run_estimation({}); }
           | ESTIMATION '(' estimation_options_list ')' ';'
             { driver.run_estimation({}); }
           | ESTIMATION symbol_list ';'
             { driver.run_estimation($2); }
           | ESTIMATION '(' estimation_options_list ')' symbol_list ';'
             { driver.run_estimation($5); }
           ;

estimation_options_list : estimation_options_list COMMA estimation_options
                        | estimation_options
                        ;

estimation_options : o_datafile
                   | o_nobs
                   | o_est_first_obs
                   | o_prefilter
                   | o_presample
                   | o_lik_init
                   | o_nograph
                   | o_posterior_nograph
                   | o_nodisplay
                   | o_graph_format
                   | o_forecasts_conf_sig
                   | o_mh_conf_sig
                   | o_mh_replic
                   | o_mh_drop
                   | o_mh_jscale
                   | o_mh_tune_jscale
                   | o_mh_tune_guess
                   | o_optim
                   | o_mh_init_scale
                   | o_mh_init_scale_factor
                   | o_mode_file
                   | o_mode_compute
                   | o_additional_optimizer_steps
                   | o_mode_check
                   | o_mode_check_neighbourhood_size
                   | o_mode_check_symmetric_plots
                   | o_mode_check_number_of_points
                   | o_prior_trunc
                   | o_mh_posterior_mode_estimation
                   | o_mh_nblocks
                   | o_load_mh_file
                   | o_load_results_after_load_mh
                   | o_loglinear
                   | o_logdata
                   | o_nodecomposition
                   | o_nodiagnostic
                   | o_bayesian_irf
                   | o_relative_irf
                   | o_dsge_var
                   | o_dsge_varlag
                   | o_irf
                   | o_tex
                   | o_forecast
                   | o_smoother
                   | o_moments_varendo
                   | o_contemporaneous_correlation
                   | o_filtered_vars
                   | o_conditional_likelihood
                   | o_fast_kalman_filter
                   | o_kalman_algo
                   | o_kalman_tol
                   | o_diffuse_kalman_tol
                   | o_xls_sheet
                   | o_xls_range
                   | o_filter_step_ahead
                   | o_solve_algo
                   | o_constant
                   | o_noconstant
                   | o_mh_recover
                   | o_mh_initialize_from_previous_mcmc
                   | o_mh_initialize_from_previous_mcmc_directory
                   | o_mh_initialize_from_previous_mcmc_record
                   | o_mh_initialize_from_previous_mcmc_prior
                   | o_diffuse_filter
                   | o_plot_priors
                   | o_order
                   | o_aim_solver
                   | o_partial_information
                   | o_filter_covariance
                   | o_updated_covariance
                   | o_filter_decomposition
                   | o_smoothed_state_uncertainty
                   | o_smoother_redux
                   | o_selected_variables_only
                   | o_conditional_variance_decomposition
                   | o_cova_compute
                   | o_irf_shocks
                   | o_sub_draws
                   | o_sylvester
                   | o_sylvester_fixed_point_tol
                   | o_lyapunov
                   | o_lyapunov_fixed_point_tol
                   | o_lyapunov_doubling_tol
                   | o_dr
                   | o_dr_cycle_reduction_tol
                   | o_dr_logarithmic_reduction_tol
                   | o_dr_logarithmic_reduction_maxiter
                   | o_analytic_derivation
                   | o_ar
                   | o_endogenous_prior
                   | o_use_univariate_filters_if_singularity_is_detected
                   | o_qz_zero_threshold
                   | o_taper_steps
                   | o_geweke_interval
                   | o_raftery_lewis_qrs
                   | o_raftery_lewis_diagnostics
                   | o_brooks_gelman_plotrows
                   | o_mcmc_jumping_covariance
                   | o_irf_plot_threshold
                   | o_posterior_max_subsample_draws
                   | o_consider_all_endogenous
                   | o_consider_all_endogenous_and_auxiliary
                   | o_consider_only_observed
                   | o_number_of_particles
                   | o_particle_filter_options
                   | o_resampling
                   | o_resampling_threshold
                   | o_resampling_method
                   | o_filter_algorithm
                   | o_nonlinear_filter_initialization
                   | o_cpf_weights
                   | o_proposal_approximation
                   | o_distribution_approximation
                   | o_dirname
                   | o_huge_number
                   | o_silent_optimizer
                   | o_proposal_distribution
                   | o_no_posterior_kernel_density
                   | o_posterior_sampling_method
                   | o_posterior_sampler_options
                   | o_keep_kalman_algo_if_singularity_is_detected
                   | o_use_penalized_objective_for_hessian
                   | o_rescale_prediction_error_covariance
                   | o_analytical_girf
                   | o_irf_in_percent
                   | o_emas_girf
                   | o_emas_drop
                   | o_emas_tolf
                   | o_emas_max_iter
                   | o_stderr_multiples
                   | o_diagonal_only
                   | o_no_init_estimation_check_first_obs
                   | o_heteroskedastic_filter
                   ;

name_value_pair : QUOTED_STRING COMMA QUOTED_STRING
                  { $$ = "''" + $1 + "'',''" + $3 + "''"; }
                | QUOTED_STRING COMMA signed_number
                  { $$ = "''" + $1 + "''," + $3; }
                ;

name_value_pair_list : name_value_pair
                     | name_value_pair_list COMMA name_value_pair
                       { $$ = $1 + ',' + $3; }
                     ;

name_value_pair_with_boolean : name_value_pair
                             | QUOTED_STRING COMMA boolean
                               { $$ = "''" + $1 + "''," + $3; }
                             ;

name_value_pair_with_boolean_list : name_value_pair_with_boolean
                                  | name_value_pair_with_boolean_list COMMA name_value_pair_with_boolean
                                    { $$ = $1 + ',' + $3; }
                                  ;


name_value_pair_with_suboptions : name_value_pair
                                | QUOTED_STRING COMMA vec_str
                                  {
                                    $$ = "''" + $1 + "'',{";
                                    for (auto &it : $3)
                                      {
                                        if (&it != &($3).front())
                                          $$ += ",";
                                        $$ += "''" + it + "''";
                                      }
                                    $$ += '}';
                                  }
                                | QUOTED_STRING COMMA '(' name_value_pair_list ')'
                                  { $$ = "''" + $1 + "'',''(" + $4 + ")''"; }
                                ;

name_value_pair_with_suboptions_list : name_value_pair_with_suboptions
                                     | name_value_pair_with_suboptions_list COMMA name_value_pair_with_suboptions
                                       { $$ = $1 + ',' + $3; }
                                     ;

varobs : VAROBS { driver.check_varobs(); } varobs_list ';';

varobs_list : varobs_list symbol
              { driver.add_varobs($2); }
            | varobs_list COMMA symbol
              { driver.add_varobs($3); }
            | symbol
              { driver.add_varobs($1); }
            ;

varexobs : VAREXOBS { driver.check_varexobs(); } varexobs_list ';';

varexobs_list : varexobs_list symbol
              { driver.add_varexobs($2); }
            | varexobs_list COMMA symbol
              { driver.add_varexobs($3); }
            | symbol
              { driver.add_varexobs($1); }
            ;

deterministic_trends : DETERMINISTIC_TRENDS ';' trend_list END ';' { driver.set_deterministic_trends(); };
observation_trends : OBSERVATION_TRENDS ';' trend_list END ';' { driver.set_trends(); };

trend_list : trend_list trend_element
           | trend_element
           ;

trend_element :  symbol '(' expression ')' ';' { driver.set_trend_element($1, $3); };

filter_initial_state : FILTER_INITIAL_STATE ';' filter_initial_state_list END ';' { driver.set_filter_initial_state(); };

filter_initial_state_list : filter_initial_state_list filter_initial_state_element
                          | filter_initial_state_element
                          ;

filter_initial_state_element : symbol '(' signed_integer ')' EQUAL expression ';' { driver.set_filter_initial_state_element($1, $3, $6); };

unit_root_vars : UNIT_ROOT_VARS symbol_list ';' { driver.set_unit_root_vars(); };

optim_weights : OPTIM_WEIGHTS ';' optim_weights_list END ';' { driver.optim_weights(); };

optim_weights_list : optim_weights_list symbol expression ';'
                     { driver.set_optim_weights($2, $3); }
                   | optim_weights_list symbol COMMA symbol expression ';'
                     { driver.set_optim_weights($2, $4, $5); }
                   | symbol expression ';'
                     { driver.set_optim_weights($1, $2); }
                   | symbol COMMA symbol expression ';'
                     { driver.set_optim_weights($1, $3, $4); }
                   ;

osr_params : OSR_PARAMS symbol_list ';' { driver.set_osr_params($2); };


osr_options_list : osr_options_list COMMA osr_options
                 | osr_options
                 ;

osr_options : stoch_simul_primary_options
            | o_osr_maxit
            | o_osr_tolf
            | o_opt_algo
            | o_optim
            | o_huge_number
            | o_silent_optimizer
            | o_analytic_derivation
            | o_analytic_derivation_mode
            ;

osr : OSR ';'
      { driver.run_osr({}); }
    | OSR '(' osr_options_list ')' ';'
      { driver.run_osr({}); }
    | OSR symbol_list ';'
      { driver.run_osr($2); }
    | OSR '(' osr_options_list ')' symbol_list ';'
      {driver.run_osr($5); }
    ;

dynatype : DYNATYPE '(' filename ')' ';'
           { driver.run_dynatype($3, {}); }
         | DYNATYPE '(' filename ')' symbol_list ';'
           { driver.run_dynatype($3, $5); }
         ;

dynasave : DYNASAVE '(' filename ')' ';'
           { driver.run_dynasave($3, {}); }
         | DYNASAVE '(' filename ')' symbol_list ';'
           { driver.run_dynasave($3, $5); }
         ;

load_params_and_steady_state : LOAD_PARAMS_AND_STEADY_STATE '(' filename ')' ';'
                               { driver.run_load_params_and_steady_state($3); }
                             ;

save_params_and_steady_state : SAVE_PARAMS_AND_STEADY_STATE '(' filename ')' ';'
                               { driver.run_save_params_and_steady_state($3); }
                             ;

identification : IDENTIFICATION ';'
                 { driver.run_identification(); }
               | IDENTIFICATION '(' identification_options_list ')' ';'
                 { driver.run_identification(); }
               ;

identification_options_list : identification_option COMMA identification_options_list
                            | identification_option
                            ;

identification_option : o_ar
                      | o_useautocorr
                      | o_load_ident_files
                      | o_prior_mc
                      | o_advanced
                      | o_max_dim_cova_group
                      | o_gsa_prior_range
                      | o_periods
                      | o_replic
                      | o_gsa_sample_file
                      | o_parameter_set
                      | o_lik_init
                      | o_kalman_algo
                      | o_nograph
                      | o_nodisplay
                      | o_graph_format
                      | o_diffuse_filter
                      | o_prior_trunc
                      | o_analytic_derivation
                      | o_analytic_derivation_mode
                      | o_tex
                      | o_no_identification_strength
                      | o_no_identification_reducedform
                      | o_no_identification_moments
                      | o_no_identification_minimal
                      | o_no_identification_spectrum
                      | o_normalize_jacobians
                      | o_grid_nbr
                      | o_tol_rank
                      | o_tol_deriv
                      | o_tol_sv
                      | o_checks_via_subsets
                      | o_max_dim_subsets_groups
                      | o_order
                      | o_schur_vec_tol
                      ;

model_comparison : MODEL_COMPARISON mc_filename_list ';'
                   { driver.run_model_comparison(); }
                 | MODEL_COMPARISON '(' o_marginal_density ')' mc_filename_list ';'
                   { driver.run_model_comparison(); }
                 ;

filename : symbol
         | QUOTED_STRING
         ;

namespace_qualified_symbol : symbol
                           | namespace_qualified_symbol '.' symbol
                             { $$ = $1 + "." + $3; }
                           ;

namespace_qualified_filename : namespace_qualified_symbol
                             | QUOTED_STRING
                             ;

parallel_local_filename_list : filename
                               { driver.add_parallel_local_file($1); }
                             | parallel_local_filename_list COMMA filename
                               { driver.add_parallel_local_file($3); }
                             ;

mc_filename_list : filename
                   { driver.add_mc_filename($1); }
                 | filename '(' non_negative_number ')'
                   { driver.add_mc_filename($1, $3); }
                 | mc_filename_list filename
                   { driver.add_mc_filename($2); }
                 | mc_filename_list filename '(' non_negative_number ')'
                   { driver.add_mc_filename($2, $4); }
                 | mc_filename_list COMMA filename
                   { driver.add_mc_filename($3); }
                 | mc_filename_list COMMA filename '(' non_negative_number ')'
                   { driver.add_mc_filename($3, $5); }
                 ;

planner_objective : PLANNER_OBJECTIVE { driver.begin_planner_objective(); }
                    hand_side { driver.end_planner_objective($3); } ';';

ramsey_model : RAMSEY_MODEL ';'
                { driver.ramsey_model(); }
              | RAMSEY_MODEL '(' ramsey_model_options_list ')' ';'
                { driver.ramsey_model(); }
              ;

ramsey_policy : RAMSEY_POLICY ';'
                { driver.ramsey_policy({}); }
              | RAMSEY_POLICY '(' ramsey_policy_options_list ')' ';'
                { driver.ramsey_policy({}); }
              | RAMSEY_POLICY symbol_list ';'
                { driver.ramsey_policy($2); }
              | RAMSEY_POLICY '(' ramsey_policy_options_list ')' symbol_list ';'
                { driver.ramsey_policy($5); }
              ;

ramsey_constraints : RAMSEY_CONSTRAINTS ';' ramsey_constraints_list END ';'
                     { driver.add_ramsey_constraints_statement(); }
		   ;

ramsey_constraints_list : ramsey_constraints_list ramsey_constraint
                 | ramsey_constraint
		 ;

ramsey_constraint : NAME  LESS expression ';'
                    { driver.ramsey_constraint_add_less($1,$3); }
		  | NAME  GREATER  expression ';'
                    { driver.ramsey_constraint_add_greater($1,$3); }
		  | NAME  LESS_EQUAL expression ';'
                    { driver.ramsey_constraint_add_less_equal($1,$3); }
		  | NAME  GREATER_EQUAL  expression ';'
                    { driver.ramsey_constraint_add_greater_equal($1,$3); }
		  ;

evaluate_planner_objective : EVALUATE_PLANNER_OBJECTIVE ';'
                             { driver.evaluate_planner_objective(); }
                           | EVALUATE_PLANNER_OBJECTIVE '(' evaluate_planner_objective_options_list ')' ';'
                             { driver.evaluate_planner_objective(); }
                           ;

evaluate_planner_objective_options_list : evaluate_planner_objective_option COMMA evaluate_planner_objective_options_list
                                        | evaluate_planner_objective_option
                                        ;

evaluate_planner_objective_option : o_evaluate_planner_objective_periods
                                  | o_evaluate_planner_objective_drop
                                  ;

occbin_setup : OCCBIN_SETUP ';'
               { driver.occbin_setup(); }
             | OCCBIN_SETUP '(' occbin_setup_options_list ')' ';'
               { driver.occbin_setup(); }

occbin_setup_options_list : occbin_setup_option COMMA occbin_setup_options_list
                          | occbin_setup_option
                          ;

occbin_setup_option : o_occbin_simul_periods
                    | o_occbin_simul_maxit
                    | o_occbin_simul_curb_retrench
                    | o_occbin_simul_check_ahead_periods
                    | o_occbin_simul_max_check_ahead_periods
                    | o_occbin_simul_periodic_solution
                    | o_occbin_simul_debug
                    | o_occbin_simul_reset_check_ahead_periods
                    | o_occbin_likelihood_periods
                    | o_occbin_likelihood_maxit
                    | o_occbin_likelihood_curb_retrench
                    | o_occbin_likelihood_check_ahead_periods
                    | o_occbin_likelihood_max_check_ahead_periods
                    | o_occbin_likelihood_periodic_solution
                    | o_occbin_likelihood_max_kalman_iterations
                    | o_occbin_likelihood_inversion_filter
                    | o_occbin_likelihood_piecewise_kalman_filter
                    | o_occbin_smoother_periods
                    | o_occbin_smoother_maxit
                    | o_occbin_smoother_curb_retrench
                    | o_occbin_smoother_check_ahead_periods
                    | o_occbin_smoother_max_check_ahead_periods
                    | o_occbin_smoother_periodic_solution
                    | o_occbin_smoother_inversion_filter
                    | o_occbin_smoother_piecewise_kalman_filter
                    | o_occbin_smoother_debug
                    | o_occbin_filter_use_relaxation
                    ;

occbin_solver : OCCBIN_SOLVER ';'
                { driver.occbin_solver(); }
              | OCCBIN_SOLVER '(' occbin_solver_options_list ')' ';'
                { driver.occbin_solver(); }
              ;

occbin_solver_options_list : occbin_solver_option COMMA occbin_solver_options_list
                           | occbin_solver_option
                           ;

occbin_solver_option : o_occbin_simul_periods
                     | o_occbin_simul_maxit
                     | o_occbin_simul_curb_retrench
                     | o_occbin_simul_check_ahead_periods
                     | o_occbin_simul_max_check_ahead_periods
                     | o_occbin_simul_reset_check_ahead_periods
                     | o_occbin_simul_debug
                     | o_occbin_simul_periodic_solution
                     ;

occbin_write_regimes : OCCBIN_WRITE_REGIMES ';'
                       { driver.occbin_write_regimes(); }
                     | OCCBIN_WRITE_REGIMES '(' occbin_write_regimes_options_list ')' ';'
                       { driver.occbin_write_regimes(); }

occbin_write_regimes_options_list : occbin_write_regimes_option COMMA occbin_write_regimes_options_list
                                  | occbin_write_regimes_option
                                  ;

occbin_write_regimes_option : o_occbin_write_regimes_periods
                            | o_occbin_write_regimes_filename
                            | o_occbin_write_regimes_smoother
                            | o_occbin_write_regimes_simul
                            ;

occbin_graph : OCCBIN_GRAPH ';'
               { driver.occbin_graph({}); }
             | OCCBIN_GRAPH '(' occbin_graph_options_list ')' ';'
               { driver.occbin_graph({}); }
             | OCCBIN_GRAPH symbol_list ';'
               { driver.occbin_graph($2); }
             | OCCBIN_GRAPH '(' occbin_graph_options_list ')' symbol_list ';'
               { driver.occbin_graph($5); }
             ;

occbin_graph_options_list : occbin_graph_option COMMA occbin_graph_options_list
                          | occbin_graph_option
                          ;

occbin_graph_option : o_occbin_graph_noconstant ;

discretionary_policy : DISCRETIONARY_POLICY ';'
                       { driver.discretionary_policy({}); }
                     | DISCRETIONARY_POLICY '(' discretionary_policy_options_list ')' ';'
                       { driver.discretionary_policy({}); }
                     | DISCRETIONARY_POLICY symbol_list ';'
                       { driver.discretionary_policy($2); }
                     | DISCRETIONARY_POLICY '(' discretionary_policy_options_list ')' symbol_list ';'
                       { driver.discretionary_policy($5); }
                     ;

discretionary_policy_options_list : discretionary_policy_options_list COMMA discretionary_policy_options
                           | discretionary_policy_options
                           ;

discretionary_policy_options : ramsey_policy_options
                             | o_discretionary_tol;
                             | o_dp_maxit;
                             ;

ramsey_model_options_list : ramsey_model_options_list COMMA ramsey_model_options
                           | ramsey_model_options
                           ;

ramsey_model_options :  o_planner_discount
                      | o_planner_discount_latex_name
                      | o_instruments
                      ;

ramsey_policy_options_list : ramsey_policy_options_list COMMA ramsey_policy_options
                           | ramsey_policy_options
                           ;

ramsey_policy_options : stoch_simul_primary_options
                      | o_planner_discount
                      | o_instruments
                      ;

write_latex_dynamic_model : WRITE_LATEX_DYNAMIC_MODEL ';'
                            { driver.write_latex_dynamic_model(false); }
                          | WRITE_LATEX_DYNAMIC_MODEL '(' WRITE_EQUATION_TAGS ')' ';'
                            { driver.write_latex_dynamic_model(true); }
                          ;

write_latex_static_model : WRITE_LATEX_STATIC_MODEL ';'
                           { driver.write_latex_static_model(false); }
                          | WRITE_LATEX_STATIC_MODEL '(' WRITE_EQUATION_TAGS ')' ';'
                            { driver.write_latex_static_model(true); }
                         ;

write_latex_original_model : WRITE_LATEX_ORIGINAL_MODEL ';'
                           { driver.write_latex_original_model(false); }
                          | WRITE_LATEX_ORIGINAL_MODEL '(' WRITE_EQUATION_TAGS ')' ';'
                            { driver.write_latex_original_model(true); }
                         ;

write_latex_steady_state_model : WRITE_LATEX_STEADY_STATE_MODEL ';'
                                 { driver.write_latex_steady_state_model(); }
                               ;

shock_decomposition : SHOCK_DECOMPOSITION ';'
                      { driver.shock_decomposition({}); }
                    | SHOCK_DECOMPOSITION '(' shock_decomposition_options_list ')' ';'
                      { driver.shock_decomposition({}); }
                    | SHOCK_DECOMPOSITION symbol_list ';'
                      { driver.shock_decomposition($2); }
                    | SHOCK_DECOMPOSITION '(' shock_decomposition_options_list ')' symbol_list ';'
                      { driver.shock_decomposition($5); }
                    ;

realtime_shock_decomposition : REALTIME_SHOCK_DECOMPOSITION ';'
                               { driver.realtime_shock_decomposition({}); }
                             | REALTIME_SHOCK_DECOMPOSITION '(' realtime_shock_decomposition_options_list ')' ';'
                               { driver.realtime_shock_decomposition({}); }
                             | REALTIME_SHOCK_DECOMPOSITION symbol_list ';'
                               { driver.realtime_shock_decomposition($2); }
                             | REALTIME_SHOCK_DECOMPOSITION '(' realtime_shock_decomposition_options_list ')' symbol_list ';'
                               { driver.realtime_shock_decomposition($5); }
                             ;

plot_shock_decomposition : PLOT_SHOCK_DECOMPOSITION ';'
                           { driver.plot_shock_decomposition({}); }
                         | PLOT_SHOCK_DECOMPOSITION '(' plot_shock_decomposition_options_list ')' ';'
                           { driver.plot_shock_decomposition({}); }
                         | PLOT_SHOCK_DECOMPOSITION symbol_list ';'
                           { driver.plot_shock_decomposition($2); }
                         | PLOT_SHOCK_DECOMPOSITION '(' plot_shock_decomposition_options_list ')' symbol_list ';'
                           { driver.plot_shock_decomposition($5); }
                         ;

initial_condition_decomposition : INITIAL_CONDITION_DECOMPOSITION ';'
                                  { driver.initial_condition_decomposition({}); }
                                | INITIAL_CONDITION_DECOMPOSITION '(' initial_condition_decomposition_options_list ')' ';'
                                  { driver.initial_condition_decomposition({}); }
                                | INITIAL_CONDITION_DECOMPOSITION symbol_list ';'
                                  { driver.initial_condition_decomposition($2); }
                                | INITIAL_CONDITION_DECOMPOSITION '(' initial_condition_decomposition_options_list ')' symbol_list ';'
                                  { driver.initial_condition_decomposition($5); }
                                ;

squeeze_shock_decomposition : SQUEEZE_SHOCK_DECOMPOSITION ';'
                              { driver.squeeze_shock_decomposition({}); }
                            | SQUEEZE_SHOCK_DECOMPOSITION symbol_list ';'
                              { driver.squeeze_shock_decomposition($2); }
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
                     | o_bvar_conf_sig
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

bvar_irf : BVAR_IRF '(' INT_NUMBER COMMA QUOTED_STRING ')' ';'
                { driver.bvar_irf($3, $5); }

sbvar_option : o_datafile
             | o_freq
             | o_initial_year
             | o_initial_subperiod
             | o_final_year
             | o_final_subperiod
             | o_data
             | o_vlist
             | o_vlistlog
             | o_vlistper
             | o_restriction_fname
             | o_nlags
             | o_cross_restrictions
             | o_contemp_reduced_form
             | o_real_pseudo_forecast
             | o_no_bayesian_prior
             | o_dummy_obs
             | o_nstates
             | o_indxscalesstates
             | o_alpha
             | o_beta
             | o_gsig2_lmdm
             | o_q_diag
             | o_flat_prior
             | o_ncsk
             | o_nstd
             | o_ninv
             | o_indxparr
             | o_indxovr
             | o_aband
             | o_indxap
             | o_apband
             | o_indximf
             | o_indxfore
             | o_foreband
             | o_indxgforhat
             | o_indxgimfhat
             | o_indxestima
             | o_indxgdls
             | o_eq_ms
             | o_cms
             | o_ncms
             | o_eq_cms
             | o_tlindx
             | o_tlnumber
             | o_cnum
             | o_forecast
             | o_coefficients_prior_hyperparameters;
             ;

sbvar_options_list : sbvar_option COMMA sbvar_options_list
                   | sbvar_option
                   ;

sbvar : SBVAR ';'
        { driver.sbvar(); }
      | SBVAR '(' sbvar_options_list ')' ';'
        { driver.sbvar(); }
      ;

ms_variance_decomposition_option : o_output_file_tag
                                 | o_file_tag
                                 | o_simulation_file_tag
                                 | o_filtered_probabilities
                                 | o_no_error_bands
                                 | o_error_band_percentiles
                                 | o_shock_draws
                                 | o_shocks_per_parameter
                                 | o_thinning_factor
                                 | o_free_parameters
                                 | o_regime
                                 | o_regimes
                                 | o_parameter_uncertainty
                                 | o_horizon
                                 ;

ms_variance_decomposition_options_list : ms_variance_decomposition_option COMMA ms_variance_decomposition_options_list
                                       | ms_variance_decomposition_option
                                       ;

ms_variance_decomposition : MS_VARIANCE_DECOMPOSITION ';'
                            { driver.ms_variance_decomposition(); }
                          | MS_VARIANCE_DECOMPOSITION '(' ms_variance_decomposition_options_list ')' ';'
                            { driver.ms_variance_decomposition(); }
                          ;

ms_forecast_option : o_output_file_tag
                   | o_file_tag
                   | o_simulation_file_tag
                   | o_data_obs_nbr
                   | o_error_band_percentiles
                   | o_shock_draws
                   | o_shocks_per_parameter
                   | o_thinning_factor
                   | o_free_parameters
                   | o_median
                   | o_regime
                   | o_regimes
                   | o_parameter_uncertainty
                   | o_horizon
                   ;

ms_forecast_options_list : ms_forecast_option COMMA ms_forecast_options_list
                         | ms_forecast_option
                         ;

ms_forecast : MS_FORECAST ';'
              { driver.ms_forecast(); }
            | MS_FORECAST '(' ms_forecast_options_list ')' ';'
              { driver.ms_forecast(); }
            ;

ms_irf_option : o_output_file_tag
              | o_file_tag
              | o_simulation_file_tag
              | o_parameter_uncertainty
              | o_horizon
              | o_filtered_probabilities
              | o_error_band_percentiles
              | o_shock_draws
              | o_shocks_per_parameter
              | o_thinning_factor
              | o_free_parameters
              | o_median
              | o_regime
              | o_regimes
              ;

ms_irf_options_list : ms_irf_option COMMA ms_irf_options_list
                    | ms_irf_option
                    ;

ms_irf : MS_IRF ';'
         { driver.ms_irf({}); }
       | MS_IRF '(' ms_irf_options_list ')' ';'
         { driver.ms_irf({}); }
       | MS_IRF symbol_list ';'
         { driver.ms_irf($2); }
       | MS_IRF '(' ms_irf_options_list ')' symbol_list ';'
         { driver.ms_irf($5); }
       ;

ms_compute_probabilities_option : o_output_file_tag
                                | o_file_tag
                                | o_filtered_probabilities
                                | o_real_time_smoothed
                                ;

ms_compute_probabilities_options_list : ms_compute_probabilities_option COMMA ms_compute_probabilities_options_list
                                      | ms_compute_probabilities_option
                                      ;

ms_compute_probabilities : MS_COMPUTE_PROBABILITIES ';'
                           { driver.ms_compute_probabilities(); }
                         | MS_COMPUTE_PROBABILITIES '(' ms_compute_probabilities_options_list ')' ';'
                           { driver.ms_compute_probabilities(); }
                         ;

ms_compute_mdd_option : o_output_file_tag
                      | o_file_tag
                      | o_simulation_file_tag
                      | o_proposal_type
                      | o_proposal_lower_bound
                      | o_proposal_upper_bound
                      | o_proposal_draws
                      | o_use_mean_center
                      ;

ms_compute_mdd_options_list : ms_compute_mdd_option COMMA ms_compute_mdd_options_list
                            | ms_compute_mdd_option
                            ;

ms_compute_mdd : MS_COMPUTE_MDD ';'
                 { driver.ms_compute_mdd(); }
               | MS_COMPUTE_MDD '(' ms_compute_mdd_options_list ')' ';'
                 { driver.ms_compute_mdd(); }
               ;

ms_simulation_option : o_output_file_tag
                     | o_file_tag
                     | o_ms_mh_replic
                     | o_ms_drop
                     | o_thinning_factor
                     | o_adaptive_mh_draws
                     | o_save_draws
                     ;

ms_simulation_options_list : ms_simulation_option COMMA ms_simulation_options_list
                           | ms_simulation_option
                           ;

ms_simulation : MS_SIMULATION ';'
                { driver.ms_simulation(); }
              | MS_SIMULATION '(' ms_simulation_options_list ')' ';'
                { driver.ms_simulation(); }
              ;

ms_estimation_option : o_coefficients_prior_hyperparameters
                     | o_freq
                     | o_initial_year
                     | o_initial_subperiod
                     | o_final_year
                     | o_final_subperiod
                     | o_datafile
                     | o_xls_sheet
                     | o_xls_range
                     | o_nlags
                     | o_cross_restrictions
                     | o_contemp_reduced_form
                     | o_no_bayesian_prior
                     | o_alpha
                     | o_beta
                     | o_gsig2_lmdm
                     | o_specification
                     | o_output_file_tag
                     | o_file_tag
                     | o_no_create_init
                     | o_convergence_starting_value
                     | o_convergence_ending_value
                     | o_convergence_increment_value
                     | o_max_iterations_starting_value
                     | o_max_iterations_increment_value
                     | o_max_block_iterations
                     | o_max_repeated_optimization_runs
                     | o_function_convergence_criterion
                     | o_parameter_convergence_criterion
                     | o_number_of_large_perturbations
                     | o_number_of_small_perturbations
                     | o_number_of_posterior_draws_after_perturbation
                     | o_max_number_of_stages
                     | o_random_function_convergence_criterion
                     | o_random_parameter_convergence_criterion
                     ;

ms_estimation_options_list : ms_estimation_option COMMA ms_estimation_options_list
                           | ms_estimation_option
                           ;

ms_estimation : MS_ESTIMATION ';'
                { driver.ms_estimation(); }
              | MS_ESTIMATION '(' ms_estimation_options_list ')' ';'
                { driver.ms_estimation(); }
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
                          | o_gsa_morris_nliv
                          | o_gsa_morris_ntra
                          | o_gsa_nsam
                          | o_gsa_load_redform
                          | o_gsa_load_rmse
                          | o_gsa_load_stab
                          | o_gsa_alpha2_stab
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
                          | o_gsa_neighborhood_width
                          | o_gsa_pvalue_ks
                          | o_gsa_pvalue_corr
                          | o_datafile
                          | o_nobs
                          | o_first_obs
                          | o_prefilter
                          | o_presample
                          | o_nograph
                          | o_nodisplay
                          | o_graph_format
                          | o_forecasts_conf_sig
                          | o_mh_conf_sig
                          | o_loglinear
                          | o_mode_file
                          | o_load_ident_files
                          | o_useautocorr
                          | o_ar
                          | o_kalman_algo
                          | o_lik_init
                          | o_diffuse_filter
                          | o_analytic_derivation
                          | o_analytic_derivation_mode
                          ;

shock_decomposition_options_list : shock_decomposition_option COMMA shock_decomposition_options_list
                                 | shock_decomposition_option
                                 ;

shock_decomposition_option : o_parameter_set
                           | o_datafile
                           | o_use_shock_groups
                           | o_colormap
                           | o_shock_decomposition_nograph
                           | o_first_obs
                           | o_nobs
                           | o_init_state
                           | o_forecast_type
                           | o_shock_decomposition_with_epilogue
                           | o_prefilter
                           | o_loglinear
                           | o_diffuse_kalman_tol
                           | o_diffuse_filter
                           | o_xls_sheet
                           | o_xls_range
                           ;

realtime_shock_decomposition_options_list : realtime_shock_decomposition_option COMMA realtime_shock_decomposition_options_list
                                         | realtime_shock_decomposition_option
                                         ;

realtime_shock_decomposition_option : o_parameter_set
                                    | o_datafile
                                    | o_first_obs
                                    | o_nobs
                                    | o_use_shock_groups
                                    | o_colormap
                                    | o_shock_decomposition_nograph
                                    | o_shock_decomposition_presample
                                    | o_shock_decomposition_forecast
                                    | o_save_realtime
                                    | o_fast_realtime
                                    | o_shock_decomposition_with_epilogue
                                    ;

plot_shock_decomposition_options_list : plot_shock_decomposition_option COMMA plot_shock_decomposition_options_list
                                      | plot_shock_decomposition_option
                                      ;

plot_shock_decomposition_option : o_use_shock_groups
                                | o_colormap
                                | o_psd_nodisplay
                                | o_psd_graph_format
                                | o_psd_detail_plot
                                | o_psd_interactive
                                | o_psd_screen_shocks
                                | o_psd_steadystate
                                | o_psd_type
                                | o_psd_fig_name
                                | o_psd_write_xls
                                | o_psd_realtime
                                | o_psd_vintage
                                | o_psd_plot_init_date
                                | o_psd_plot_end_date
                                | o_psd_diff
                                | o_psd_flip
                                | o_psd_nograph
                                | o_psd_init2shocks
                                | o_psd_max_nrows
                                ;

initial_condition_decomposition_options_list : initial_condition_decomposition_option COMMA initial_condition_decomposition_options_list
                                             | initial_condition_decomposition_option
                                             ;

initial_condition_decomposition_option : o_icd_type
                                       | o_icd_detail_plot
                                       | o_icd_steadystate
                                       | o_icd_write_xls
                                       | o_icd_plot_init_date
                                       | o_icd_plot_end_date
                                       | o_icd_nodisplay
                                       | o_icd_graph_format
                                       | o_icd_fig_name
                                       | o_icd_diff
                                       | o_icd_flip
                                       | o_icd_colormap
                                       | o_icd_max_nrows
                                       | o_icd_with_epilogue
                                       ;

homotopy_setup: HOMOTOPY_SETUP ';' homotopy_list END ';'
               { driver.end_homotopy();};

homotopy_list : homotopy_item
              | homotopy_list homotopy_item
              ;

homotopy_item : symbol COMMA expression COMMA expression ';'
                { driver.homotopy_val($1, $3, $5);}
              | symbol COMMA expression ';'
                { driver.homotopy_val($1, nullptr, $3);}
              ;

forecast: FORECAST ';'
          { driver.forecast({}); }
        | FORECAST '(' forecast_options ')' ';'
          { driver.forecast({}); }
        | FORECAST symbol_list ';'
          { driver.forecast($2); }
        | FORECAST '(' forecast_options ')' symbol_list ';'
          { driver.forecast($5); }
        ;

forecast_options: forecast_option
          | forecast_options COMMA forecast_option
          ;

forecast_option: o_periods
          | o_forecasts_conf_sig
          | o_nograph
          | o_nodisplay
          | o_graph_format
          ;

conditional_forecast : CONDITIONAL_FORECAST '(' conditional_forecast_options ')' ';'
                       { driver.conditional_forecast(); }
                     ;

conditional_forecast_options : conditional_forecast_option
                             | conditional_forecast_options COMMA conditional_forecast_option
                             ;

conditional_forecast_option : o_periods
                            | o_replic
                            | o_conditional_forecast_conf_sig
                            | o_controlled_varexo
                            | o_parameter_set
                            ;

plot_conditional_forecast : PLOT_CONDITIONAL_FORECAST symbol_list ';'
                            { driver.plot_conditional_forecast(nullopt, $2); }
                          | PLOT_CONDITIONAL_FORECAST '(' PERIODS EQUAL INT_NUMBER ')' symbol_list ';'
                            { driver.plot_conditional_forecast($5, $7); }
                          ;

conditional_forecast_paths : CONDITIONAL_FORECAST_PATHS ';' conditional_forecast_paths_shock_list END ';'
                             { driver.conditional_forecast_paths(); }
                           ;

conditional_forecast_paths_shock_list : conditional_forecast_paths_shock_elem
                                      | conditional_forecast_paths_shock_list conditional_forecast_paths_shock_elem
                                      ;

conditional_forecast_paths_shock_elem : VAR symbol ';' PERIODS period_list ';' VALUES value_list ';'
                                        { driver.add_det_shock($2, $5, $8, ParsingDriver::DetShockType::conditional_forecast); }
                                      ;

steady_state_model : STEADY_STATE_MODEL ';' { driver.begin_steady_state_model(); }
                     steady_state_equation_list END ';' { driver.reset_data_tree(); }
                   ;

steady_state_equation_list : steady_state_equation_list steady_state_equation
                           | steady_state_equation
                           ;

steady_state_equation : symbol EQUAL expression ';'
                        { driver.add_steady_state_model_equal($1, $3); }
                      | '[' symbol_list ']' EQUAL expression ';'
                        { driver.add_steady_state_model_equal_multiple($2, $5); }
                      ;

calib_smoother : CALIB_SMOOTHER ';'
                 { driver.calib_smoother({}); }
               | CALIB_SMOOTHER '(' calib_smoother_options_list ')' ';'
                 { driver.calib_smoother({}); }
               | CALIB_SMOOTHER symbol_list ';'
                 { driver.calib_smoother($2); }
               | CALIB_SMOOTHER '(' calib_smoother_options_list ')' symbol_list ';'
                 { driver.calib_smoother($5); }
               ;

calib_smoother_options_list : calib_smoother_option COMMA calib_smoother_options_list
                            | calib_smoother_option
                            ;

calib_smoother_option : o_filtered_vars
                      | o_filter_step_ahead
                      | o_datafile
                      | o_prefilter
                      | o_kalman_algo
                      | o_loglinear
                      | o_first_obs
                      | o_filter_covariance
                      | o_updated_covariance
                      | o_filter_decomposition
                      | o_diffuse_kalman_tol
                      | o_diffuse_filter
                      | o_smoothed_state_uncertainty
                      | o_smoother_redux
                      | o_parameter_set
                      | o_xls_sheet
                      | o_xls_range
                      | o_heteroskedastic_filter
                      | o_nobs
                      ;

generate_irfs : GENERATE_IRFS ';' END ';'
                { driver.end_generate_irfs(); }
              | GENERATE_IRFS ';' generate_irfs_element_list END ';'
                { driver.end_generate_irfs(); }
              | GENERATE_IRFS '(' generate_irfs_options_list ')' ';' END ';'
                { driver.end_generate_irfs(); }
              | GENERATE_IRFS '(' generate_irfs_options_list ')' ';' generate_irfs_element_list END ';'
                { driver.end_generate_irfs(); }
              ;

generate_irfs_options_list : generate_irfs_option COMMA generate_irfs_options_list
                           | generate_irfs_option
                           ;

generate_irfs_option : o_stderr_multiples
                     | o_diagonal_only
                     ;

generate_irfs_element_list : generate_irfs_element_list generate_irfs_element
                           | generate_irfs_element
                           ;

generate_irfs_element : NAME COMMA generate_irfs_exog_element_list ';'
                        { driver.add_generate_irfs_element($1); }
                      ;

generate_irfs_exog_element_list : generate_irfs_exog_element_list COMMA symbol EQUAL signed_number
                                  { driver.add_generate_irfs_exog_element($3, $5); }
                                | symbol EQUAL signed_number
                                  { driver.add_generate_irfs_exog_element($1, $3); }
                                ;

extended_path : EXTENDED_PATH ';'
                { driver.extended_path(); }
              | EXTENDED_PATH '(' extended_path_options_list ')' ';'
                { driver.extended_path(); }
              ;

extended_path_options_list : extended_path_option COMMA extended_path_options_list
                           | extended_path_option
                           ;

extended_path_option : o_periods
                     | o_solver_periods
                     | o_extended_path_order
                     | o_hybrid
		     | o_lmmcp
                     ;

model_diagnostics : MODEL_DIAGNOSTICS ';'
                    { driver.model_diagnostics(); }
                  ;

calibration_range : '[' expression COMMA expression ']'
                    { $$ = { $2, $4 }; }
                  | PLUS
                    { $$ = { driver.add_non_negative_constant("0"), driver.add_inf_constant() }; }
                  | MINUS
                    { $$ = { driver.add_uminus(driver.add_inf_constant()), driver.add_non_negative_constant("0") }; }
                  ;

moment_calibration : MOMENT_CALIBRATION ';' moment_calibration_list END ';'
                     { driver.end_moment_calibration(); }
                   ;

moment_calibration_list : moment_calibration_item
                        | moment_calibration_list moment_calibration_item
                        ;

moment_calibration_item : symbol COMMA symbol COMMA calibration_range ';'
                          { driver.add_moment_calibration_item($1, $3, "0", $5); }
                        | symbol COMMA symbol '(' signed_integer ')' COMMA calibration_range ';'
                          { driver.add_moment_calibration_item($1, $3, $5, $8); }
                        | symbol COMMA symbol '(' signed_integer_range ')' COMMA calibration_range ';'
                          { driver.add_moment_calibration_item($1, $3, $5, $8); }
                        ;

irf_calibration : IRF_CALIBRATION ';' irf_calibration_list END ';'
                  { driver.end_irf_calibration(); }
                | IRF_CALIBRATION '(' o_relative_irf ')' ';' irf_calibration_list END ';'
                  { driver.end_irf_calibration(); }
                ;

irf_calibration_list : irf_calibration_item
                     | irf_calibration_list irf_calibration_item
                     ;

irf_calibration_item : symbol COMMA symbol COMMA calibration_range ';'
                       { driver.add_irf_calibration_item($1, "1", $3, $5); }
                     | symbol '(' INT_NUMBER ')' COMMA symbol COMMA calibration_range ';'
                       { driver.add_irf_calibration_item($1, $3, $6, $8); }
                     | symbol '(' integer_range ')' COMMA symbol COMMA calibration_range ';'
                       { driver.add_irf_calibration_item($1, $3, $6, $8); }
                     ;

smoother2histval : SMOOTHER2HISTVAL ';'
                   { driver.smoother2histval(); }
                 | SMOOTHER2HISTVAL '(' smoother2histval_options_list ')' ';'
                   { driver.smoother2histval(); }
                 ;

smoother2histval_options_list : smoother2histval_option COMMA smoother2histval_options_list
                              | smoother2histval_option
                              ;

smoother2histval_option : o_infile
                        | o_invars
                        | o_period
                        | o_outfile
                        | o_outvars
                        ;

shock_groups : SHOCK_GROUPS ';' shock_group_list END ';'
               { driver.end_shock_groups("default"); }
             | SHOCK_GROUPS '(' NAME EQUAL symbol ')' ';' shock_group_list END ';'
               {driver.end_shock_groups($5);}
             ;

shock_group_list : shock_group_list shock_group_element
                 | shock_group_element
                 ;

shock_group_element : symbol EQUAL shock_name_list ';' { driver.add_shock_group($1); }
                    | QUOTED_STRING EQUAL shock_name_list ';' { driver.add_shock_group($1); }
                    ;

shock_name_list : shock_name_list COMMA symbol {driver.add_shock_group_element($3);}
                | shock_name_list symbol {driver.add_shock_group_element($2);}
                | symbol {driver.add_shock_group_element($1);}
                ;

init2shocks : INIT2SHOCKS ';' init2shocks_list END ';'
              { driver.end_init2shocks("default"); }
            | INIT2SHOCKS '(' NAME EQUAL symbol ')' ';' init2shocks_list END ';'
              {driver.end_init2shocks($5);}
            ;

init2shocks_list : init2shocks_list init2shocks_element
                 | init2shocks_element
                 ;

init2shocks_element : symbol symbol ';' { driver.add_init2shocks($1, $2); }
                    | symbol COMMA symbol ';' { driver.add_init2shocks($1, $3); }
                    ;

o_dr_algo : DR_ALGO EQUAL INT_NUMBER {
                                       if ($3 == "0")
                                         driver.warning("dr_algo option is now deprecated, and may be removed in a future version of Dynare");
                                       else
                                         driver.error("dr_algo=1 option is no longer supported");
                                     };
o_solve_algo : SOLVE_ALGO EQUAL INT_NUMBER { driver.option_num("solve_algo", $3); };
o_simul_algo : SIMUL_ALGO EQUAL INT_NUMBER {
                                             if ($3 == "0")
                                               driver.warning("simul_algo option is now deprecated, and may be removed in a future version of Dynare");
                                             else
                                               driver.error("simul_algo=1 option is no longer supported");
                                           };
o_stack_solve_algo : STACK_SOLVE_ALGO EQUAL INT_NUMBER { driver.option_num("stack_solve_algo", $3); };
o_robust_lin_solve : ROBUST_LIN_SOLVE { driver.option_num("simul.robust_lin_solve", "true"); };
o_endogenous_terminal_period : ENDOGENOUS_TERMINAL_PERIOD { driver.option_num("endogenous_terminal_period", "true"); };
o_linear : LINEAR { driver.linear(); };
o_order : ORDER EQUAL INT_NUMBER { driver.option_num("order", $3); };
o_replic : REPLIC EQUAL INT_NUMBER { driver.option_num("replic", $3); };
o_drop : DROP EQUAL INT_NUMBER { driver.option_num("drop", $3); };
o_ar : AR EQUAL INT_NUMBER { driver.option_num("ar", $3); };
o_nocorr : NOCORR { driver.option_num("nocorr", "true"); };
o_nofunctions : NOFUNCTIONS { driver.option_num("nofunctions", "true"); };
o_nomoments : NOMOMENTS { driver.option_num("nomoments", "true"); };
o_irf : IRF EQUAL INT_NUMBER { driver.option_num("irf", $3); };
o_irf_shocks : IRF_SHOCKS EQUAL '(' symbol_list ')' { driver.option_symbol_list("irf_shocks", $4); };
o_hp_filter : HP_FILTER EQUAL non_negative_number { driver.option_num("hp_filter", $3); };
o_hp_ngrid : HP_NGRID EQUAL INT_NUMBER {
                                         driver.warning("The 'hp_ngrid' option is deprecated. It has been superseded by the 'filtered_theoretical_moments_grid' option.");
                                         driver.option_num("filtered_theoretical_moments_grid", $3);
                                       };
o_filtered_theoretical_moments_grid : FILTERED_THEORETICAL_MOMENTS_GRID EQUAL INT_NUMBER { driver.option_num("filtered_theoretical_moments_grid", $3); };
o_one_sided_hp_filter : ONE_SIDED_HP_FILTER EQUAL non_negative_number { driver.option_num("one_sided_hp_filter", $3); };
o_periods : PERIODS EQUAL INT_NUMBER { driver.option_num("periods", $3); };
o_solver_periods : SOLVER_PERIODS EQUAL INT_NUMBER { driver.option_num("ep.periods", $3); };
o_extended_path_order : ORDER EQUAL INT_NUMBER { driver.option_num("ep.stochastic.order", $3); };
o_hybrid : HYBRID { driver.option_num("ep.stochastic.hybrid_order", "2"); };
o_steady_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("steady.maxit", $3); };
o_simul_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("simul.maxit", $3); };
o_bandpass_filter : BANDPASS_FILTER { driver.option_num("bandpass.indicator", "true"); }
                  | BANDPASS_FILTER EQUAL vec_int
                    {
                      driver.option_num("bandpass.indicator", "true");
                      driver.option_vec_int("bandpass.passband", $3);
                    }
                   ;
o_dp_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("dp.maxit", $3); };
o_osr_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("osr.maxit", $3); };
o_osr_tolf : TOLF EQUAL non_negative_number { driver.option_num("osr.tolf", $3); };
o_pf_tolf : TOLF EQUAL non_negative_number { driver.option_num("dynatol.f", $3); };
o_pf_tolx : TOLX EQUAL non_negative_number { driver.option_num("dynatol.x", $3); };
o_steady_tolf : TOLF EQUAL non_negative_number { driver.option_num("solve_tolf", $3); };
o_steady_tolx : TOLX EQUAL non_negative_number { driver.option_num("solve_tolx", $3); };
o_opt_algo : OPT_ALGO EQUAL INT_NUMBER { driver.option_num("osr.opt_algo", $3); }
           | OPT_ALGO EQUAL filename { driver.option_str("osr.opt_algo", $3); }
           ;
o_cutoff : CUTOFF EQUAL non_negative_number { driver.cutoff($3); };
o_markowitz : MARKOWITZ EQUAL non_negative_number { driver.option_num("markowitz", $3); };

// Options of perfect_foresight_* commands that control steady state computation at terminal condition
o_pf_steady_solve_algo : STEADY_SOLVE_ALGO EQUAL INT_NUMBER { driver.option_num("simul.steady_solve_algo", $3); };
o_pf_steady_maxit : STEADY_MAXIT EQUAL INT_NUMBER { driver.option_num("simul.steady_maxit", $3); };
o_pf_steady_tolf : STEADY_TOLF EQUAL non_negative_number { driver.option_num("simul.steady_tolf", $3); };
o_pf_steady_tolx : STEADY_TOLX EQUAL non_negative_number { driver.option_num("simul.steady_tolx", $3); };
o_pf_steady_markowitz : STEADY_MARKOWITZ EQUAL non_negative_number { driver.option_num("simul.steady_markowitz", $3); };

o_minimal_solving_periods : MINIMAL_SOLVING_PERIODS EQUAL non_negative_number { driver.option_num("minimal_solving_periods", $3); };
o_mfs : MFS EQUAL INT_NUMBER { driver.mfs($3); };
o_simul : SIMUL; // Do nothing, only here for backward compatibility
o_simul_replic : SIMUL_REPLIC EQUAL INT_NUMBER { driver.option_num("simul_replic", $3); };
o_simul_seed : SIMUL_SEED EQUAL INT_NUMBER { driver.error("'simul_seed' option is no longer supported; use 'set_dynare_seed' command instead"); } ;
o_qz_criterium : QZ_CRITERIUM EQUAL non_negative_number { driver.option_num("qz_criterium", $3); };
o_qz_zero_threshold : QZ_ZERO_THRESHOLD EQUAL non_negative_number { driver.option_num("qz_zero_threshold", $3); };
o_file : FILE EQUAL filename { driver.option_str("file", $3); };
o_pac_name : MODEL_NAME EQUAL symbol { driver.option_str("pac.model_name", $3); };
o_pac_aux_model_name : AUXILIARY_MODEL_NAME EQUAL symbol { driver.option_str("pac.aux_model_name", $3); };
o_pac_discount : DISCOUNT EQUAL symbol { driver.option_str("pac.discount", $3); };
o_pac_growth : GROWTH { driver.begin_model(); } EQUAL hand_side { driver.set_pac_growth($4); };
o_pac_auxname : AUXNAME EQUAL symbol { driver.set_pac_auxname($3); };
o_pac_kind : KIND EQUAL pac_target_kind { driver.set_pac_kind($3); };
o_var_name : MODEL_NAME EQUAL symbol { driver.option_str("var.model_name", $3); };
o_series : SERIES EQUAL symbol { driver.option_str("series", $3); };
o_datafile : DATAFILE EQUAL filename { driver.option_str("datafile", $3); };
o_filename : FILENAME EQUAL filename { driver.option_str("filename", $3); };
o_var_eq_tags : EQTAGS EQUAL vec_str { driver.option_vec_str("var.eqtags", $3); }
o_var_structural : STRUCTURAL { driver.option_num("var.structural", "true"); }
o_dirname : DIRNAME EQUAL filename { driver.option_str("dirname", $3); };
o_huge_number : HUGE_NUMBER EQUAL non_negative_number { driver.option_num("huge_number", $3); };
o_nobs : NOBS EQUAL vec_int
         { driver.option_vec_int("nobs", $3); }
       | NOBS EQUAL vec_int_number
         { driver.option_vec_int("nobs", $3); }
       ;
o_trend_component_model_name : MODEL_NAME EQUAL symbol { driver.option_str("trend_component.name", $3); };
o_trend_component_model_targets : TARGETS EQUAL vec_str { driver.option_vec_str("trend_component.targets", $3); }
o_trend_component_model_eq_tags : EQTAGS EQUAL vec_str { driver.option_vec_str("trend_component.eqtags", $3); }
o_conditional_variance_decomposition : CONDITIONAL_VARIANCE_DECOMPOSITION EQUAL vec_int
                                       { driver.option_vec_int("conditional_variance_decomposition", $3); }
                                     | CONDITIONAL_VARIANCE_DECOMPOSITION EQUAL vec_int_number
                                       { driver.option_vec_int("conditional_variance_decomposition", $3); }
                                     ;
o_est_first_obs : FIRST_OBS EQUAL vec_int
                  { driver.option_vec_int("first_obs", $3); }
                | FIRST_OBS EQUAL vec_int_number
                  { driver.option_vec_int("first_obs", $3); }
                ;
o_posterior_sampling_method : POSTERIOR_SAMPLING_METHOD EQUAL QUOTED_STRING
                              { driver.option_str("posterior_sampler_options.posterior_sampling_method", $3); } ;
o_first_obs : FIRST_OBS EQUAL INT_NUMBER { driver.option_num("first_obs", $3); };
o_data_first_obs : FIRST_OBS EQUAL date_expr { driver.option_date("firstobs", $3); } ;
o_first_simulation_period : FIRST_SIMULATION_PERIOD EQUAL INT_NUMBER { driver.option_num("first_simulation_period", $3); };
o_date_first_simulation_period : FIRST_SIMULATION_PERIOD EQUAL date_expr { driver.option_date("firstsimulationperiod", $3); } ;
o_last_simulation_period : LAST_SIMULATION_PERIOD EQUAL INT_NUMBER { driver.option_num("last_simulation_period", $3); };
o_date_last_simulation_period : LAST_SIMULATION_PERIOD EQUAL date_expr { driver.option_date("lastsimulationperiod", $3); } ;
o_last_obs : LAST_OBS EQUAL INT_NUMBER { driver.option_num("last_obs", $3); };
o_data_last_obs : LAST_OBS EQUAL date_expr { driver.option_date("lastobs", $3); } ;
o_keep_kalman_algo_if_singularity_is_detected : KEEP_KALMAN_ALGO_IF_SINGULARITY_IS_DETECTED { driver.option_num("kalman.keep_kalman_algo_if_singularity_is_detected", "true"); } ;
o_data_nobs : NOBS EQUAL INT_NUMBER { driver.option_num("nobs", $3); };
o_shift : SHIFT EQUAL signed_number { driver.option_num("shift", $3); };
o_shape : SHAPE EQUAL prior_distribution { driver.prior_shape = $3; };
o_mode : MODE EQUAL signed_number { driver.option_num("mode", $3); };
o_mean : MEAN EQUAL signed_number { driver.option_num("mean", $3); };
o_mean_vec : MEAN EQUAL vec_value { driver.option_vec_value("mean", $3); };
o_truncate : TRUNCATE EQUAL vec_value { driver.option_vec_value("truncate", $3); };
o_stdev : STDEV EQUAL non_negative_number { driver.option_num("stdev", $3); };
o_jscale : JSCALE EQUAL non_negative_number { driver.option_num("jscale", $3); };
o_init : INIT EQUAL signed_number { driver.option_num("init", $3); };
o_bounds : BOUNDS EQUAL vec_value_w_inf { driver.option_vec_value("bounds", $3); };
o_domain : DOMAINN EQUAL vec_value { driver.option_vec_value("domain", $3); };
o_interval : INTERVAL EQUAL vec_value { driver.option_vec_value("interval", $3); };
o_variance : VARIANCE EQUAL expression { driver.set_prior_variance($3); }
o_variance_mat : VARIANCE EQUAL vec_of_vec_value { driver.option_vec_of_vec_value("variance",$3); }
o_prefilter : PREFILTER EQUAL INT_NUMBER { driver.option_num("prefilter", $3); };
o_presample : PRESAMPLE EQUAL INT_NUMBER { driver.option_num("presample", $3); };
o_lik_init : LIK_INIT EQUAL INT_NUMBER { driver.option_num("lik_init", $3); };
o_nograph : NOGRAPH
            { driver.option_num("nograph", "true"); }
          | GRAPH
            { driver.option_num("nograph", "false"); }
          ;
o_posterior_nograph : POSTERIOR_NOGRAPH
            { driver.option_num("no_graph.posterior", "true"); }
          | POSTERIOR_GRAPH
            { driver.option_num("no_graph.posterior", "false"); }
          ;
o_psd_nograph : NOGRAPH { driver.option_num("no_graph.plot_shock_decomposition", "true"); }
o_shock_decomposition_nograph : NOGRAPH { driver.option_num("no_graph.shock_decomposition", "true"); }
o_init_state : INIT_STATE EQUAL INT_NUMBER { driver.option_num("shock_decomp.init_state", $3); };
o_forecast_type : FORECAST EQUAL UNCONDITIONAL
                  { driver.option_str("shock_decomp.forecast_type", "unconditional"); }
                | FORECAST EQUAL CONDITIONAL
                  { driver.option_str("shock_decomp.forecast_type", "conditional"); }
o_shock_decomposition_presample : PRESAMPLE EQUAL INT_NUMBER { driver.option_num("shock_decomp.presample", $3); };
o_shock_decomposition_forecast : FORECAST EQUAL INT_NUMBER { driver.option_num("shock_decomp.forecast", $3); };
o_save_realtime : SAVE_REALTIME EQUAL vec_int { driver.option_vec_int("shock_decomp.save_realtime", $3); };
o_fast_realtime : FAST_REALTIME EQUAL vec_int
                  { driver.option_vec_int("shock_decomp.fast_realtime", $3); }
                | FAST_REALTIME EQUAL vec_int_number
                  { driver.option_vec_int("shock_decomp.fast_realtime", $3); }
                ;
o_nodisplay : NODISPLAY { driver.option_num("nodisplay", "true"); };
o_icd_nodisplay : NODISPLAY { driver.option_num("initial_condition_decomp.nodisplay", "true"); };
o_psd_nodisplay : NODISPLAY { driver.option_num("plot_shock_decomp.nodisplay", "true"); };
o_psd_init2shocks : INIT2SHOCKS { driver.option_str("plot_shock_decomp.init2shocks", "default"); }
                  | INIT2SHOCKS EQUAL symbol { driver.option_str("plot_shock_decomp.init2shocks", $3); }
                  ;
o_icd_max_nrows : MAX_NROWS EQUAL INT_NUMBER { driver.option_num("initial_condition_decomp.max_nrows", $3); };
o_psd_max_nrows : MAX_NROWS EQUAL INT_NUMBER { driver.option_num("plot_shock_decomp.max_nrows", $3); };
o_graph_format : GRAPH_FORMAT EQUAL allowed_graph_formats
                 { driver.process_graph_format_option(); }
               | GRAPH_FORMAT EQUAL '(' list_allowed_graph_formats ')'
                 { driver.process_graph_format_option(); }
               ;
o_icd_graph_format : GRAPH_FORMAT EQUAL allowed_graph_formats
                     { driver.initial_condition_decomp_process_graph_format_option(); }
                   | GRAPH_FORMAT EQUAL '(' list_allowed_graph_formats ')'
                     { driver.initial_condition_decomp_process_graph_format_option(); }
                   ;
o_psd_graph_format : GRAPH_FORMAT EQUAL allowed_graph_formats
                     { driver.plot_shock_decomp_process_graph_format_option(); }
                   | GRAPH_FORMAT EQUAL '(' list_allowed_graph_formats ')'
                     { driver.plot_shock_decomp_process_graph_format_option(); }
                   ;
o_shock_decomposition_with_epilogue : WITH_EPILOGUE { driver.option_num("shock_decomp.with_epilogue", "true"); };
o_icd_with_epilogue : WITH_EPILOGUE { driver.option_num("initial_condition_decomp.with_epilogue", "true"); };
allowed_graph_formats : EPS
                        { driver.add_graph_format("eps"); }
                      | FIG
                        { driver.add_graph_format("fig"); }
                      | PDF
                        { driver.add_graph_format("pdf"); }
                      | NONE
                        { driver.add_graph_format("none"); }
                      ;
list_allowed_graph_formats : allowed_graph_formats
                           | list_allowed_graph_formats COMMA allowed_graph_formats
                           ;

o_subsample_name : symbol EQUAL date_expr ':' date_expr
                   { driver.set_subsample_name_equal_to_date_range($1, $3, $5); }
                 ;
o_bvar_conf_sig : CONF_SIG EQUAL non_negative_number { driver.option_num("bvar.conf_sig", $3); };
o_forecasts_conf_sig : CONF_SIG EQUAL non_negative_number { driver.option_num("forecasts.conf_sig", $3); };
o_conditional_forecast_conf_sig : CONF_SIG EQUAL non_negative_number { driver.option_num("conditional_forecast.conf_sig", $3); };
o_mh_conf_sig : MH_CONF_SIG EQUAL non_negative_number { driver.option_num("mh_conf_sig", $3); };
o_mh_replic : MH_REPLIC EQUAL INT_NUMBER { driver.option_num("mh_replic", $3); };
o_posterior_max_subsample_draws : POSTERIOR_MAX_SUBSAMPLE_DRAWS EQUAL INT_NUMBER { driver.option_num("posterior_max_subsample_draws", $3); };
o_mh_drop : MH_DROP EQUAL non_negative_number { driver.option_num("mh_drop", $3); };
o_mh_jscale : MH_JSCALE EQUAL non_negative_number { driver.option_num("mh_jscale", $3); };
o_mh_tune_jscale : MH_TUNE_JSCALE EQUAL non_negative_number
                 { driver.option_num("mh_tune_jscale.target", $3); driver.option_num("mh_tune_jscale.status", "true");}
                 | MH_TUNE_JSCALE {driver.option_num("mh_tune_jscale.status", "true");};
o_mh_tune_guess : MH_TUNE_GUESS EQUAL non_negative_number { driver.option_num("mh_tune_jscale.guess", $3); };
o_optim : OPTIM  EQUAL '(' name_value_pair_with_boolean_list ')' { driver.option_str("optim_opt", $4); };
o_posterior_sampler_options : POSTERIOR_SAMPLER_OPTIONS EQUAL '(' name_value_pair_with_suboptions_list ')' { driver.option_str("posterior_sampler_options.sampling_opt", $4); };
o_proposal_distribution : PROPOSAL_DISTRIBUTION EQUAL symbol { driver.option_str("posterior_sampler_options.posterior_sampling_method.proposal_distribution", $3); };
o_no_posterior_kernel_density : NO_POSTERIOR_KERNEL_DENSITY
                             { driver.option_num("estimation.moments_posterior_density.indicator", "false"); }
                           ;
o_mh_init_scale : MH_INIT_SCALE EQUAL non_negative_number { driver.option_num("mh_init_scale", $3); };
o_mh_init_scale_factor : MH_INIT_SCALE_FACTOR EQUAL non_negative_number { driver.option_num("mh_init_scale_factor", $3); };
o_mode_file : MODE_FILE EQUAL filename { driver.option_str("mode_file", $3); };
o_mode_compute : MODE_COMPUTE EQUAL INT_NUMBER { driver.option_num("mode_compute", $3); };
               | MODE_COMPUTE EQUAL symbol { driver.option_str("mode_compute", $3); };
o_mode_check : MODE_CHECK { driver.option_num("mode_check.status", "true"); };
o_mode_check_neighbourhood_size : MODE_CHECK_NEIGHBOURHOOD_SIZE EQUAL signed_number_w_inf { driver.option_num("mode_check.neighbourhood_size", $3); };
o_mode_check_number_of_points : MODE_CHECK_NUMBER_OF_POINTS EQUAL INT_NUMBER { driver.option_num("mode_check.number_of_points", $3); };
o_mode_check_symmetric_plots : MODE_CHECK_SYMMETRIC_PLOTS EQUAL INT_NUMBER { driver.option_num("mode_check.symmetric_plots", $3); };
o_prior_trunc : PRIOR_TRUNC EQUAL non_negative_number { driver.option_num("prior_trunc", $3); };
o_mh_posterior_mode_estimation : MH_POSTERIOR_MODE_ESTIMATION { driver.option_num("mh_posterior_mode_estimation", "true"); };
o_mh_nblocks : MH_NBLOCKS EQUAL INT_NUMBER { driver.option_num("mh_nblck", $3); };
o_load_mh_file : LOAD_MH_FILE { driver.option_num("load_mh_file", "true"); };
o_load_results_after_load_mh : LOAD_RESULTS_AFTER_LOAD_MH { driver.option_num("load_results_after_load_mh", "true"); };
o_loglinear : LOGLINEAR { driver.option_num("loglinear", "true"); };
o_linear_approximation : LINEAR_APPROXIMATION { driver.option_num("linear_approximation", "true"); };
o_logdata : LOGDATA { driver.option_num("logdata", "true"); };
o_conditional_likelihood : CONDITIONAL_LIKELIHOOD {driver.option_num("conditional_likelihood.status", "true"); };
o_nodiagnostic : NODIAGNOSTIC { driver.option_num("nodiagnostic", "true"); };
o_bayesian_irf : BAYESIAN_IRF { driver.option_num("bayesian_irf", "true"); };
o_dsge_var : DSGE_VAR EQUAL non_negative_number
             { driver.option_num("dsge_var", $3); }
           | DSGE_VAR EQUAL INF_CONSTANT
             { driver.option_num("dsge_var", "Inf"); }
           | DSGE_VAR
             { driver.option_str("dsge_var", "NaN"); }
           ;
o_dsge_varlag : DSGE_VARLAG EQUAL INT_NUMBER { driver.option_num("dsge_varlag", $3); };
o_tex : TEX { driver.option_num("TeX", "true"); };
o_forecast : FORECAST EQUAL INT_NUMBER { driver.option_num("forecast", $3); };
o_smoother : SMOOTHER { driver.option_num("smoother", "true"); };
o_moments_varendo : MOMENTS_VARENDO { driver.option_num("moments_varendo", "true"); };
o_contemporaneous_correlation : CONTEMPORANEOUS_CORRELATION { driver.option_num("contemporaneous_correlation", "true"); };
o_filtered_vars : FILTERED_VARS { driver.option_num("filtered_vars", "true"); };
o_relative_irf : RELATIVE_IRF { driver.option_num("relative_irf", "true"); };
o_fast_kalman_filter : FAST_KALMAN_FILTER  { driver.option_num("fast_kalman_filter", "true"); };
o_kalman_algo : KALMAN_ALGO EQUAL INT_NUMBER { driver.option_num("kalman_algo", $3); };
o_kalman_tol : KALMAN_TOL EQUAL non_negative_number { driver.option_num("kalman_tol", $3); };
o_diffuse_kalman_tol : DIFFUSE_KALMAN_TOL EQUAL non_negative_number { driver.option_num("diffuse_kalman_tol", $3); };
o_schur_vec_tol : SCHUR_VEC_TOL EQUAL non_negative_number { driver.option_num("schur_vec_tol", $3); };
o_marginal_density : MARGINAL_DENSITY EQUAL LAPLACE
                     { driver.option_str("mc_marginal_density", "laplace"); }
                   | MARGINAL_DENSITY EQUAL MODIFIEDHARMONICMEAN
                     { driver.option_str("mc_marginal_density", "modifiedharmonicmean"); }
                   ;
o_print : PRINT { driver.option_num("noprint", "false"); };
o_noprint : NOPRINT { driver.option_num("noprint", "true"); };
o_xls_sheet : XLS_SHEET EQUAL symbol { driver.option_str("xls_sheet", $3); } // Kept for backward compatibility, but no longeer recommended (see #67)
            | XLS_SHEET EQUAL QUOTED_STRING { driver.option_str("xls_sheet", $3); }
            ;
o_xls_range : XLS_RANGE EQUAL range { driver.option_str("xls_range", $3); };
o_filter_step_ahead : FILTER_STEP_AHEAD EQUAL vec_int { driver.option_vec_int("filter_step_ahead", $3); };
o_taper_steps : TAPER_STEPS EQUAL vec_int { driver.option_vec_int("convergence.geweke.taper_steps", $3); };
o_geweke_interval : GEWEKE_INTERVAL EQUAL vec_value { driver.option_vec_value("convergence.geweke.geweke_interval",$3); };
o_raftery_lewis_diagnostics : RAFTERY_LEWIS_DIAGNOSTICS { driver.option_num("convergence.rafterylewis.indicator", "true"); };
o_raftery_lewis_qrs : RAFTERY_LEWIS_QRS EQUAL vec_value { driver.option_vec_value("convergence.rafterylewis.qrs",$3); };
o_brooks_gelman_plotrows: BROOKS_GELMAN_PLOTROWS EQUAL INT_NUMBER { driver.option_num("convergence.brooksgelman.plotrows", $3); };
o_constant : CONSTANT { driver.option_num("noconstant", "false"); };
o_noconstant : NOCONSTANT { driver.option_num("noconstant", "true"); };
o_mh_recover : MH_RECOVER { driver.option_num("mh_recover", "true"); };
o_mh_initialize_from_previous_mcmc : MH_INITIALIZE_FROM_PREVIOUS_MCMC { driver.option_num("mh_initialize_from_previous_mcmc.status", "true"); }
o_mh_initialize_from_previous_mcmc_directory : MH_INITIALIZE_FROM_PREVIOUS_MCMC_DIRECTORY EQUAL filename { driver.option_str("mh_initialize_from_previous_mcmc.directory", $3); };
o_mh_initialize_from_previous_mcmc_record : MH_INITIALIZE_FROM_PREVIOUS_MCMC_RECORD EQUAL filename { driver.option_str("mh_initialize_from_previous_mcmc.record", $3); };
o_mh_initialize_from_previous_mcmc_prior : MH_INITIALIZE_FROM_PREVIOUS_MCMC_PRIOR EQUAL filename { driver.option_str("mh_initialize_from_previous_mcmc.prior", $3); };
o_diffuse_filter: DIFFUSE_FILTER {driver.option_num("diffuse_filter", "true"); };
o_plot_priors: PLOT_PRIORS EQUAL INT_NUMBER {driver.option_num("plot_priors", $3); };
o_aim_solver: AIM_SOLVER {driver.option_num("aim_solver", "true"); };
o_partial_information : PARTIAL_INFORMATION {driver.option_num("partial_information", "true"); };
o_sub_draws: SUB_DRAWS EQUAL INT_NUMBER {driver.option_num("sub_draws",$3);};
o_planner_discount : PLANNER_DISCOUNT EQUAL expression { driver.set_planner_discount($3); };
o_planner_discount_latex_name : PLANNER_DISCOUNT_LATEX_NAME EQUAL TEX_NAME { driver.set_planner_discount_latex_name($3); };
o_sylvester : SYLVESTER EQUAL FIXED_POINT {driver.option_num("sylvester_fp", "true"); }
               | SYLVESTER EQUAL DEFAULT {driver.option_num("sylvester_fp", "false"); };
o_sylvester_fixed_point_tol : SYLVESTER_FIXED_POINT_TOL EQUAL non_negative_number {driver.option_num("sylvester_fixed_point_tol",$3);};
o_lyapunov : LYAPUNOV EQUAL FIXED_POINT {driver.option_num("lyapunov_fp", "true"); }
              | LYAPUNOV EQUAL DOUBLING {driver.option_num("lyapunov_db", "true"); }
              | LYAPUNOV EQUAL SQUARE_ROOT_SOLVER {driver.option_num("lyapunov_srs", "true"); }
              | LYAPUNOV EQUAL DEFAULT {driver.option_num("lyapunov_fp", "false"); driver.option_num("lyapunov_db", "false"); driver.option_num("lyapunov_srs", "false");};
o_lyapunov_complex_threshold : LYAPUNOV_COMPLEX_THRESHOLD EQUAL non_negative_number {driver.option_num("lyapunov_complex_threshold",$3);};
o_lyapunov_fixed_point_tol : LYAPUNOV_FIXED_POINT_TOL EQUAL non_negative_number {driver.option_num("lyapunov_fixed_point_tol",$3);};
o_lyapunov_doubling_tol : LYAPUNOV_DOUBLING_TOL EQUAL non_negative_number {driver.option_num("lyapunov_doubling_tol",$3);};
o_dr : DR EQUAL CYCLE_REDUCTION {driver.option_num("dr_cycle_reduction", "true"); }
       | DR EQUAL LOGARITHMIC_REDUCTION {driver.option_num("dr_logarithmic_reduction", "true"); }
       | DR EQUAL DEFAULT {driver.option_num("dr_cycle_reduction", "false"); driver.option_num("dr_logarithmic_reduction", "false");};
o_dr_cycle_reduction_tol : DR_CYCLE_REDUCTION_TOL EQUAL non_negative_number {driver.option_num("dr_cycle_reduction_tol",$3);};
o_dr_logarithmic_reduction_tol : DR_LOGARITHMIC_REDUCTION_TOL EQUAL non_negative_number {driver.option_num("dr_logarithmic_reduction_tol",$3);};
o_dr_logarithmic_reduction_maxiter : DR_LOGARITHMIC_REDUCTION_MAXITER EQUAL INT_NUMBER {driver.option_num("dr_logarithmic_reduction_maxiter",$3);};
o_psd_detail_plot : DETAIL_PLOT { driver.option_num("plot_shock_decomp.detail_plot", "true"); };
o_icd_detail_plot : DETAIL_PLOT { driver.option_num("initial_condition_decomp.detail_plot", "true"); };
o_psd_interactive : INTERACTIVE { driver.option_num("plot_shock_decomp.interactive", "true"); };
o_psd_screen_shocks : SCREEN_SHOCKS { driver.option_num("plot_shock_decomp.screen_shocks", "true"); };
o_psd_steadystate : STEADYSTATE { driver.option_num("plot_shock_decomp.steadystate", "true"); };
o_icd_steadystate : STEADYSTATE { driver.option_num("initial_condition_decomp.steadystate", "true"); };
o_icd_fig_name : FIG_NAME EQUAL filename { driver.option_str("initial_condition_decomp.fig_name", $3); };
o_psd_fig_name : FIG_NAME EQUAL filename { driver.option_str("plot_shock_decomp.fig_name", $3); };
o_psd_type : TYPE EQUAL QOQ
             { driver.option_str("plot_shock_decomp.type", "qoq"); }
           | TYPE EQUAL YOY
             { driver.option_str("plot_shock_decomp.type", "yoy"); }
           | TYPE EQUAL AOA
             { driver.option_str("plot_shock_decomp.type", "aoa"); }
           ;
o_icd_type : TYPE EQUAL QOQ
             { driver.option_str("initial_condition_decomp.type", "qoq"); }
           | TYPE EQUAL YOY
             { driver.option_str("initial_condition_decomp.type", "yoy"); }
           | TYPE EQUAL AOA
             { driver.option_str("initial_condition_decomp.type", "aoa"); }
           ;
o_icd_plot_init_date : PLOT_INIT_DATE EQUAL date_expr { driver.option_date("initial_condition_decomp.plot_init_date", $3); } ;
o_icd_plot_end_date : PLOT_END_DATE EQUAL date_expr { driver.option_date("initial_condition_decomp.plot_end_date", $3); } ;
o_psd_plot_init_date : PLOT_INIT_DATE EQUAL date_expr { driver.option_date("plot_shock_decomp.plot_init_date", $3); } ;
o_psd_plot_end_date : PLOT_END_DATE EQUAL date_expr { driver.option_date("plot_shock_decomp.plot_end_date", $3); } ;
o_icd_write_xls : WRITE_XLS { driver.option_num("initial_condition_decomp.write_xls", "true"); };
o_psd_write_xls : WRITE_XLS { driver.option_num("plot_shock_decomp.write_xls", "true"); };
o_psd_realtime : REALTIME EQUAL INT_NUMBER { driver.option_num("plot_shock_decomp.realtime", $3); };
o_psd_vintage : VINTAGE EQUAL INT_NUMBER { driver.option_num("plot_shock_decomp.vintage", $3); };
o_psd_diff : DIFF { driver.option_num("plot_shock_decomp.diff", "true"); };
o_icd_diff : DIFF { driver.option_num("initial_condition_decomp.diff", "true"); };
o_psd_flip : FLIP { driver.option_num("plot_shock_decomp.flip", "true"); };
o_icd_flip : FLIP { driver.option_num("initial_condition_decomp.flip", "true"); };
o_bvar_prior_tau : BVAR_PRIOR_TAU EQUAL signed_number { driver.option_num("bvar_prior_tau", $3); };
o_bvar_prior_decay : BVAR_PRIOR_DECAY EQUAL non_negative_number { driver.option_num("bvar_prior_decay", $3); };
o_bvar_prior_lambda : BVAR_PRIOR_LAMBDA EQUAL signed_number { driver.option_num("bvar_prior_lambda", $3); };
o_bvar_prior_mu : BVAR_PRIOR_MU EQUAL non_negative_number { driver.option_num("bvar_prior_mu", $3); };
o_bvar_prior_omega : BVAR_PRIOR_OMEGA EQUAL INT_NUMBER { driver.option_num("bvar_prior_omega", $3); };
o_bvar_prior_flat : BVAR_PRIOR_FLAT { driver.option_num("bvar_prior_flat", "true"); };
o_bvar_prior_train : BVAR_PRIOR_TRAIN EQUAL INT_NUMBER { driver.option_num("bvar_prior_train", $3); };
o_bvar_replic : BVAR_REPLIC EQUAL INT_NUMBER { driver.option_num("bvar_replic", $3); };
o_stderr_multiples : STDERR_MULTIPLES { driver.option_num("irf_opt.stderr_multiples", "true"); };
o_diagonal_only : DIAGONAL_ONLY { driver.option_num("irf_opt.diagonal_only", "true"); };
o_number_of_particles : NUMBER_OF_PARTICLES EQUAL INT_NUMBER { driver.option_num("particle.number_of_particles", $3); };
o_particle_filter_options : PARTICLE_FILTER_OPTIONS EQUAL '(' name_value_pair_with_boolean_list ')' { driver.option_str("particle.particle_filter_options", $4); }
o_resampling : RESAMPLING EQUAL SYSTEMATIC
              | RESAMPLING EQUAL NONE {driver.option_num("particle.resampling.status.systematic", "false"); driver.option_num("particle.resampling.status.none", "true"); }
              | RESAMPLING EQUAL GENERIC {driver.option_num("particle.resampling.status.systematic", "false"); driver.option_num("particle.resampling.status.generic", "true"); };
o_resampling_threshold : RESAMPLING_THRESHOLD EQUAL non_negative_number { driver.option_num("particle.resampling.threshold", $3); };
o_resampling_method : RESAMPLING_METHOD EQUAL KITAGAWA {driver.option_num("particle.resampling.method.kitagawa", "true"); driver.option_num("particle.resampling.method.smooth", "false"); driver.option_num("particle.resampling.smethod.stratified", "false"); }
              | RESAMPLING_METHOD EQUAL SMOOTH {driver.option_num("particle.resampling.method.kitagawa", "false"); driver.option_num("particle.resampling.method.smooth", "true"); driver.option_num("particle.resampling.smethod.stratified", "false"); }
              | RESAMPLING_METHOD EQUAL STRATIFIED {driver.option_num("particle.resampling.method.kitagawa", "false"); driver.option_num("particle.resampling.method.smooth", "false"); driver.option_num("particle.resampling.method.stratified", "true"); };
o_cpf_weights : CPF_WEIGHTS EQUAL AMISANOTRISTANI {driver.option_num("particle.cpf_weights_method.amisanotristani", "true"); driver.option_num("particle.cpf_weights_method.murrayjonesparslow", "false"); }
              | CPF_WEIGHTS EQUAL MURRAYJONESPARSLOW {driver.option_num("particle.cpf_weights_method.amisanotristani", "false"); driver.option_num("particle.cpf_weights_method.murrayjonesparslow", "true"); };
o_filter_algorithm : FILTER_ALGORITHM EQUAL symbol { driver.option_str("particle.filter_algorithm", $3); };
o_nonlinear_filter_initialization : NONLINEAR_FILTER_INITIALIZATION EQUAL INT_NUMBER { driver.option_num("particle.initialization", $3); };
o_proposal_approximation : PROPOSAL_APPROXIMATION EQUAL CUBATURE {driver.option_num("particle.proposal_approximation.cubature", "true"); driver.option_num("particle.proposal_approximation.unscented", "false"); driver.option_num("particle.proposal_approximation.montecarlo", "false");}
		| PROPOSAL_APPROXIMATION EQUAL UNSCENTED {driver.option_num("particle.proposal_approximation.cubature", "false"); driver.option_num("particle.proposal_approximation.unscented", "true"); driver.option_num("particle.proposal_approximation.montecarlo", "false");}
		| PROPOSAL_APPROXIMATION EQUAL MONTECARLO {driver.option_num("particle.proposal_approximation.cubature", "false"); driver.option_num("particle.proposal_approximation.unscented", "false"); driver.option_num("particle.proposal_approximation.montecarlo", "true");} ;
o_distribution_approximation : DISTRIBUTION_APPROXIMATION EQUAL CUBATURE {driver.option_num("particle.distribution_approximation.cubature", "true"); driver.option_num("particle.distribution_approximation.unscented", "false"); driver.option_num("particle.distribution_approximation.montecarlo", "false");}
		| DISTRIBUTION_APPROXIMATION EQUAL UNSCENTED {driver.option_num("particle.distribution_approximation.cubature", "false"); driver.option_num("particle.distribution_approximation.unscented", "true"); driver.option_num("particle.distribution_approximation.montecarlo", "false");}
		| DISTRIBUTION_APPROXIMATION EQUAL MONTECARLO {driver.option_num("particle.distribution_approximation.cubature", "false"); driver.option_num("particle.distribution_approximation.unscented", "false"); driver.option_num("particle.distribution_approximation.montecarlo", "true");} ;
o_gsa_identification : IDENTIFICATION EQUAL INT_NUMBER { driver.option_num("identification", $3); }; /*not in doc */
o_gsa_morris : MORRIS EQUAL INT_NUMBER { driver.option_num("morris", $3); };
o_gsa_stab : STAB  EQUAL INT_NUMBER { driver.option_num("stab", $3); };
o_gsa_redform : REDFORM  EQUAL INT_NUMBER { driver.option_num("redform", $3); };
o_gsa_pprior : PPRIOR  EQUAL INT_NUMBER { driver.option_num("pprior", $3); };
o_gsa_prior_range : PRIOR_RANGE  EQUAL INT_NUMBER { driver.option_num("prior_range", $3); };
o_gsa_ppost : PPOST  EQUAL INT_NUMBER { driver.option_num("ppost", $3); };
o_gsa_ilptau : ILPTAU EQUAL INT_NUMBER { driver.option_num("ilptau", $3); };
o_gsa_morris_nliv : MORRIS_NLIV EQUAL INT_NUMBER { driver.option_num("morris_nliv", $3); };
o_gsa_morris_ntra : MORRIS_NTRA EQUAL INT_NUMBER { driver.option_num("morris_ntra", $3); };
o_gsa_nsam : NSAM EQUAL INT_NUMBER { driver.option_num("Nsam", $3); }; /* not in doc ??*/
o_gsa_load_redform : LOAD_REDFORM EQUAL INT_NUMBER { driver.option_num("load_redform", $3); };
o_gsa_load_rmse : LOAD_RMSE EQUAL INT_NUMBER { driver.option_num("load_rmse", $3); };
o_gsa_load_stab : LOAD_STAB EQUAL INT_NUMBER { driver.option_num("load_stab", $3); };
o_gsa_alpha2_stab : ALPHA2_STAB EQUAL non_negative_number { driver.option_num("alpha2_stab", $3); };
o_gsa_logtrans_redform : LOGTRANS_REDFORM EQUAL INT_NUMBER { driver.option_num("logtrans_redform", $3); };
o_gsa_threshold_redform : THRESHOLD_REDFORM EQUAL vec_value_w_inf { driver.option_vec_value("threshold_redform",$3); };
o_gsa_ksstat_redform : KSSTAT_REDFORM EQUAL non_negative_number { driver.option_num("ksstat_redform", $3); };
o_gsa_alpha2_redform : ALPHA2_REDFORM EQUAL non_negative_number { driver.option_num("alpha2_redform", $3); };
o_gsa_namendo : NAMENDO EQUAL '(' symbol_list_or_wildcard ')' { driver.option_symbol_list("namendo", $4); };
o_gsa_namlagendo : NAMLAGENDO EQUAL '(' symbol_list_or_wildcard ')' { driver.option_symbol_list("namlagendo", $4); };
o_gsa_namexo : NAMEXO EQUAL '(' symbol_list_or_wildcard ')' { driver.option_symbol_list("namexo", $4); };
o_gsa_rmse : RMSE EQUAL INT_NUMBER { driver.option_num("rmse", $3); };
o_gsa_lik_only : LIK_ONLY EQUAL INT_NUMBER { driver.option_num("lik_only", $3); };
o_gsa_var_rmse : VAR_RMSE EQUAL '(' symbol_list_or_wildcard ')' { driver.option_symbol_list("var_rmse", $4); };
o_gsa_pfilt_rmse : PFILT_RMSE EQUAL non_negative_number { driver.option_num("pfilt_rmse", $3); };
o_gsa_istart_rmse : ISTART_RMSE EQUAL INT_NUMBER { driver.option_num("istart_rmse", $3); };
o_gsa_alpha_rmse : ALPHA_RMSE EQUAL non_negative_number { driver.option_num("alpha_rmse", $3); };
o_gsa_alpha2_rmse : ALPHA2_RMSE EQUAL non_negative_number { driver.option_num("alpha2_rmse", $3); };
o_gsa_sample_file : GSA_SAMPLE_FILE EQUAL INT_NUMBER
                    { driver.option_num("gsa_sample_file", $3); }
                  | GSA_SAMPLE_FILE EQUAL filename
                    { driver.option_str("gsa_sample_file", $3); }
                  ;
o_gsa_neighborhood_width : NEIGHBORHOOD_WIDTH EQUAL non_negative_number { driver.option_num("neighborhood_width", $3); };
o_gsa_pvalue_ks : PVALUE_KS EQUAL  non_negative_number { driver.option_num("pvalue_ks", $3); };
o_gsa_pvalue_corr : PVALUE_CORR EQUAL  non_negative_number { driver.option_num("pvalue_corr", $3); };
o_load_ident_files : LOAD_IDENT_FILES EQUAL INT_NUMBER { driver.option_num("load_ident_files", $3); }
o_useautocorr : USEAUTOCORR EQUAL INT_NUMBER { driver.option_num("useautocorr", $3); }
o_prior_mc : PRIOR_MC EQUAL INT_NUMBER { driver.option_num("prior_mc", $3); }
o_advanced : ADVANCED EQUAL signed_integer { driver.option_num("advanced", $3); }
o_max_dim_cova_group : MAX_DIM_COVA_GROUP EQUAL INT_NUMBER { driver.option_num("max_dim_cova_group", $3); }

o_homotopy_mode : HOMOTOPY_MODE EQUAL INT_NUMBER {driver.option_num("homotopy_mode",$3); };
o_homotopy_steps : HOMOTOPY_STEPS EQUAL INT_NUMBER {driver.option_num("homotopy_steps",$3); };
o_homotopy_force_continue: HOMOTOPY_FORCE_CONTINUE EQUAL INT_NUMBER { driver.option_num("homotopy_force_continue",$3); };
o_nocheck : NOCHECK {driver.option_num("steadystate.nocheck","true"); };

o_controlled_varexo : CONTROLLED_VAREXO EQUAL '(' symbol_list ')' { driver.option_symbol_list("controlled_varexo", $4); };
o_parameter_set : PARAMETER_SET EQUAL PRIOR_MODE
                  { driver.option_str("parameter_set", "prior_mode"); }
                | PARAMETER_SET EQUAL PRIOR_MEAN
                  { driver.option_str("parameter_set", "prior_mean"); }
                | PARAMETER_SET EQUAL POSTERIOR_MEAN
                  { driver.option_str("parameter_set", "posterior_mean"); }
                | PARAMETER_SET EQUAL POSTERIOR_MODE
                  { driver.option_str("parameter_set", "posterior_mode"); }
                | PARAMETER_SET EQUAL POSTERIOR_MEDIAN
                  { driver.option_str("parameter_set", "posterior_median"); }
                | PARAMETER_SET EQUAL MLE_MODE
                  { driver.option_str("parameter_set", "mle_mode"); }
                | PARAMETER_SET EQUAL CALIBRATION
                  { driver.option_str("parameter_set", "calibration"); }
                ;
o_nodecomposition : NODECOMPOSITION { driver.option_num("nodecomposition", "true"); };
o_spectral_density : SPECTRAL_DENSITY { driver.option_num("SpectralDensity.trigger", "true"); };
o_ms_drop : DROP EQUAL INT_NUMBER { driver.option_num("ms.drop", $3); };
o_ms_mh_replic : MH_REPLIC EQUAL INT_NUMBER { driver.option_num("ms.mh_replic", $3); };
o_freq : FREQ EQUAL INT_NUMBER
         { driver.option_num("ms.freq",$3); }
       | FREQ EQUAL MONTHLY
         { driver.option_num("ms.freq","12"); }
       | FREQ EQUAL QUARTERLY
         { driver.option_num("ms.freq","4"); }
       ;
o_initial_year : INITIAL_YEAR EQUAL INT_NUMBER {driver.option_num("ms.initial_year",$3); };
o_initial_subperiod : INITIAL_SUBPERIOD EQUAL INT_NUMBER {driver.option_num("ms.initial_subperiod",$3); };
o_final_year : FINAL_YEAR EQUAL INT_NUMBER {driver.option_num("ms.final_year",$3); };
o_final_subperiod : FINAL_SUBPERIOD EQUAL INT_NUMBER {driver.option_num("ms.final_subperiod",$3); };
o_data : DATA EQUAL filename { driver.option_str("ms.data", $3); };
o_vlist : VLIST EQUAL INT_NUMBER {driver.option_num("ms.vlist",$3); };
o_vlistlog : VLISTLOG EQUAL '(' symbol_list ')' {driver.option_symbol_list("ms.vlistlog", $4); };
o_vlistper : VLISTPER EQUAL INT_NUMBER {driver.option_num("ms.vlistper",$3); };
o_restriction_fname : RESTRICTION_FNAME EQUAL NAME
                      {
                        driver.warning("restriction_fname is now deprecated, and may be removed in a future version of Dynare. Use svar_identification instead.");
                        driver.option_str("ms.restriction_fname",$3);
                      }
                    | RESTRICTION_FNAME EQUAL UPPER_CHOLESKY
                      {
                        driver.warning("restriction_fname is now deprecated, and may be removed in a future version of Dynare. Use svar_identification instead.");
                        driver.option_str("ms.restriction_fname","upper_cholesky");
                      }
                    | RESTRICTION_FNAME EQUAL LOWER_CHOLESKY
                      {
                        driver.warning("restriction_fname is now deprecated, and may be removed in a future version of Dynare. Use svar_identification instead.");
                        driver.option_str("ms.restriction_fname","lower_cholesky");
                      }
                    ;
o_nlags : NLAGS EQUAL INT_NUMBER {driver.option_num("ms.nlags",$3); };
o_cross_restrictions : CROSS_RESTRICTIONS {driver.option_num("ms.cross_restrictions","true"); };
o_contemp_reduced_form : CONTEMP_REDUCED_FORM {driver.option_num("ms.contemp_reduced_form","true"); };
o_real_pseudo_forecast : REAL_PSEUDO_FORECAST EQUAL INT_NUMBER {driver.option_num("ms.real_pseudo_forecast",$3); };
o_no_bayesian_prior : NO_BAYESIAN_PRIOR {driver.option_num("ms.bayesian_prior","false"); };
o_dummy_obs : DUMMY_OBS EQUAL INT_NUMBER {driver.option_num("ms.dummy_obs",$3); };
o_nstates : NSTATES EQUAL INT_NUMBER {driver.option_num("ms.nstates",$3); };
o_indxscalesstates : INDXSCALESSTATES EQUAL INT_NUMBER {driver.option_num("ms.indxscalesstates",$3); };
o_alpha : ALPHA EQUAL non_negative_number {driver.option_num("ms.alpha",$3); };
o_beta : BETA EQUAL non_negative_number {driver.option_num("ms.beta",$3); };
o_gsig2_lmdm : GSIG2_LMDM EQUAL INT_NUMBER {driver.option_num("ms.gsig2_lmdm",$3); };
o_specification : SPECIFICATION EQUAL SIMS_ZHA
                  {driver.option_num("ms.specification","1"); }
                | SPECIFICATION EQUAL NONE
                  {driver.option_num("ms.specification","0"); }
                ;
o_q_diag : Q_DIAG EQUAL non_negative_number {driver.option_num("ms.q_diag",$3); };
o_flat_prior : FLAT_PRIOR EQUAL INT_NUMBER {driver.option_num("ms.flat_prior",$3); };
o_ncsk : NCSK EQUAL INT_NUMBER {driver.option_num("ms.ncsk",$3); };
o_nstd : NSTD EQUAL INT_NUMBER {driver.option_num("ms.nstd",$3); };
o_ninv : NINV EQUAL INT_NUMBER {driver.option_num("ms.ninv",$3); };
o_indxparr : INDXPARR EQUAL INT_NUMBER {driver.option_num("ms.indxparr",$3); };
o_indxovr : INDXOVR EQUAL INT_NUMBER {driver.option_num("ms.indxovr",$3); };
o_aband : ABAND EQUAL INT_NUMBER {driver.option_num("ms.aband",$3); };
o_indxap : INDXAP EQUAL INT_NUMBER {driver.option_num("ms.indxap",$3); };
o_apband : APBAND EQUAL INT_NUMBER {driver.option_num("ms.apband",$3); };
o_indximf : INDXIMF EQUAL INT_NUMBER {driver.option_num("ms.indximf",$3); };
o_indxfore : INDXFORE EQUAL INT_NUMBER {driver.option_num("ms.indxfore",$3); };
o_foreband : FOREBAND EQUAL INT_NUMBER {driver.option_num("ms.foreband",$3); };
o_indxgforhat : INDXGFOREHAT EQUAL INT_NUMBER {driver.option_num("ms.indxgforehat",$3); };
o_indxgimfhat : INDXGIMFHAT EQUAL INT_NUMBER {driver.option_num("ms.indxgimfhat",$3); };
o_indxestima : INDXESTIMA EQUAL INT_NUMBER {driver.option_num("ms.indxestima",$3); };
o_indxgdls : INDXGDLS EQUAL INT_NUMBER {driver.option_num("ms.indxgdls",$3); };
o_eq_ms : EQ_MS EQUAL INT_NUMBER {driver.option_num("ms.eq_ms",$3); };
o_cms : CMS EQUAL INT_NUMBER {driver.option_num("ms.cms",$3); };
o_ncms : NCMS EQUAL INT_NUMBER {driver.option_num("ms.ncms",$3); };
o_eq_cms : EQ_CMS EQUAL INT_NUMBER {driver.option_num("ms.eq_cms",$3); };
o_tlindx : TLINDX EQUAL INT_NUMBER {driver.option_num("ms.tlindx",$3); };
o_tlnumber : TLNUMBER EQUAL INT_NUMBER {driver.option_num("ms.tlnumber",$3); };
o_cnum : CNUM EQUAL INT_NUMBER {driver.option_num("ms.cnum",$3); };
o_k_order_solver : K_ORDER_SOLVER {driver.option_num("k_order_solver","true"); };
o_pruning : PRUNING { driver.option_num("pruning", "true"); };
o_chain : CHAIN EQUAL INT_NUMBER { driver.option_num("ms.chain",$3); };
o_restrictions : RESTRICTIONS EQUAL vec_of_vec_value
                 { driver.option_vec_of_vec_value("ms.restrictions",$3); }
               ;
o_duration : DURATION EQUAL non_negative_number
             { driver.option_num("ms.duration",$3); }
           | DURATION EQUAL vec_value_w_inf
             { driver.option_vec_value("ms.duration",$3); }
           ;
o_number_of_regimes : NUMBER_OF_REGIMES EQUAL INT_NUMBER { driver.option_num("ms.number_of_regimes",$3); };
o_number_of_lags : NUMBER_OF_LAGS EQUAL INT_NUMBER { driver.option_num("ms.number_of_lags",$3); };
o_parameters : PARAMETERS EQUAL '[' symbol_list ']' { driver.option_symbol_list("ms.parameters", $4); };
o_coefficients : COEFFICIENTS { driver.option_str("ms.coefficients","svar_coefficients"); };
o_variances : VARIANCES { driver.option_str("ms.variances","svar_variances"); };
o_equations : EQUATIONS EQUAL vec_int
              { driver.option_vec_int("ms.equations",$3); }
            | EQUATIONS EQUAL vec_int_number
              { driver.option_vec_int("ms.equations",$3); }
            ;
o_silent_optimizer : SILENT_OPTIMIZER { driver.option_num("silent_optimizer", "true"); };
o_instruments : INSTRUMENTS EQUAL '(' symbol_list ')' {driver.option_symbol_list("instruments", $4); };

o_ext_func_name : EXT_FUNC_NAME EQUAL namespace_qualified_filename { driver.external_function_option("name", $3); };
o_ext_func_nargs : EXT_FUNC_NARGS EQUAL INT_NUMBER { driver.external_function_option("nargs",$3); };
o_first_deriv_provided : FIRST_DERIV_PROVIDED EQUAL namespace_qualified_filename
                         { driver.external_function_option("first_deriv_provided", $3); }
                       | FIRST_DERIV_PROVIDED
                         { driver.external_function_option("first_deriv_provided", ""); }
                       ;
o_second_deriv_provided : SECOND_DERIV_PROVIDED EQUAL namespace_qualified_filename
                          { driver.external_function_option("second_deriv_provided", $3); }
                        | SECOND_DERIV_PROVIDED
                          { driver.external_function_option("second_deriv_provided", ""); }
                        ;
o_filter_covariance : FILTER_COVARIANCE
                        { driver.option_num("filter_covariance","true");}
                      ;
o_updated_covariance : UPDATED_COVARIANCE
                        { driver.option_num("updated_covariance","true");}

o_filter_decomposition : FILTER_DECOMPOSITION
                           { driver.option_num("filter_decomposition","true");}
                         ;
o_smoothed_state_uncertainty : SMOOTHED_STATE_UNCERTAINTY
                           { driver.option_num("smoothed_state_uncertainty","true");}
                         ;
o_smoother_redux : SMOOTHER_REDUX
                           { driver.option_num("smoother_redux","true");}

o_selected_variables_only : SELECTED_VARIABLES_ONLY
                           { driver.option_num("selected_variables_only","true");}
                          ;
o_cova_compute : COVA_COMPUTE EQUAL INT_NUMBER
                 { driver.option_num("cova_compute",$3);}
               ;
o_output_file_tag : OUTPUT_FILE_TAG EQUAL filename {driver.option_str("ms.output_file_tag", $3); };
o_file_tag : FILE_TAG EQUAL filename { driver.option_str("ms.file_tag", $3); };
o_no_create_init : NO_CREATE_INIT { driver.option_num("ms.create_init", "false"); };
o_simulation_file_tag : SIMULATION_FILE_TAG EQUAL filename { driver.option_str("ms.simulation_file_tag", $3); };
o_coefficients_prior_hyperparameters : COEFFICIENTS_PRIOR_HYPERPARAMETERS EQUAL vec_value
                                       { driver.option_vec_value("ms.coefficients_prior_hyperparameters",$3); };
o_convergence_starting_value : CONVERGENCE_STARTING_VALUE EQUAL non_negative_number
                               { driver.option_num("ms.convergence_starting_value",$3); };
o_convergence_ending_value : CONVERGENCE_ENDING_VALUE EQUAL non_negative_number
                             { driver.option_num("ms.convergence_ending_value",$3); };
o_convergence_increment_value : CONVERGENCE_INCREMENT_VALUE EQUAL non_negative_number
                                { driver.option_num("ms.convergence_increment_value",$3); };
o_max_iterations_starting_value : MAX_ITERATIONS_STARTING_VALUE EQUAL INT_NUMBER
                                  { driver.option_num("ms.max_iterations_starting_value",$3); };
o_max_iterations_increment_value : MAX_ITERATIONS_INCREMENT_VALUE EQUAL non_negative_number
                                   { driver.option_num("ms.max_iterations_increment_value",$3); };
o_max_block_iterations : MAX_BLOCK_ITERATIONS EQUAL INT_NUMBER
                         { driver.option_num("ms.max_block_iterations",$3); };
o_max_repeated_optimization_runs : MAX_REPEATED_OPTIMIZATION_RUNS EQUAL INT_NUMBER
                                   { driver.option_num("ms.max_repeated_optimization_runs",$3); };
o_function_convergence_criterion : FUNCTION_CONVERGENCE_CRITERION EQUAL non_negative_number
                                   { driver.option_num("ms.function_convergence_criterion",$3); };
o_parameter_convergence_criterion : PARAMETER_CONVERGENCE_CRITERION EQUAL non_negative_number
                                    { driver.option_num("ms.parameter_convergence_criterion",$3); };
o_number_of_large_perturbations : NUMBER_OF_LARGE_PERTURBATIONS EQUAL INT_NUMBER
                                  { driver.option_num("ms.number_of_large_perturbations",$3); };
o_number_of_small_perturbations : NUMBER_OF_SMALL_PERTURBATIONS EQUAL INT_NUMBER
                                  { driver.option_num("ms.number_of_small_perturbations",$3); };
o_number_of_posterior_draws_after_perturbation : NUMBER_OF_POSTERIOR_DRAWS_AFTER_PERTURBATION EQUAL INT_NUMBER
                                                 { driver.option_num("ms.number_of_posterior_draws_after_perturbation",$3); };
o_max_number_of_stages : MAX_NUMBER_OF_STAGES EQUAL INT_NUMBER
                         { driver.option_num("ms.max_number_of_stages",$3); };
o_random_function_convergence_criterion : RANDOM_FUNCTION_CONVERGENCE_CRITERION EQUAL non_negative_number
                                          { driver.option_num("ms.random_function_convergence_criterion",$3); };
o_random_parameter_convergence_criterion : RANDOM_PARAMETER_CONVERGENCE_CRITERION EQUAL non_negative_number
                                           { driver.option_num("ms.random_parameter_convergence_criterion",$3); };
o_thinning_factor : THINNING_FACTOR EQUAL INT_NUMBER { driver.option_num("ms.thinning_factor",$3); };
o_adaptive_mh_draws : ADAPTIVE_MH_DRAWS EQUAL INT_NUMBER { driver.option_num("ms.adaptive_mh_draws",$3); };
o_save_draws : SAVE_DRAWS { driver.option_num("ms.save_draws","true"); };
o_proposal_draws : PROPOSAL_DRAWS EQUAL INT_NUMBER { driver.option_num("ms.proposal_draws",$3); };
o_use_mean_center : USE_MEAN_CENTER { driver.option_num("ms.use_mean_center","true"); };
o_proposal_type : PROPOSAL_TYPE EQUAL INT_NUMBER { driver.option_num("ms.proposal_type",$3); }
o_proposal_lower_bound : PROPOSAL_LOWER_BOUND EQUAL signed_number { driver.option_num("ms.proposal_lower_bound",$3); }
o_proposal_upper_bound : PROPOSAL_UPPER_BOUND EQUAL signed_number { driver.option_num("ms.proposal_upper_bound",$3); }
o_parameter_uncertainty : PARAMETER_UNCERTAINTY { driver.option_num("ms.parameter_uncertainty","true"); };
o_horizon : HORIZON EQUAL INT_NUMBER { driver.option_num("ms.horizon",$3); };
o_filtered_probabilities : FILTERED_PROBABILITIES { driver.option_num("ms.filtered_probabilities","true"); };
o_real_time_smoothed : REAL_TIME_SMOOTHED { driver.option_num("ms.real_time_smoothed_probabilities","true"); };
o_no_error_bands : NO_ERROR_BANDS { driver.option_num("ms.error_bands","false"); };
o_error_band_percentiles : ERROR_BAND_PERCENTILES EQUAL vec_value { driver.option_vec_value("ms.percentiles",$3); };
o_shock_draws : SHOCK_DRAWS EQUAL INT_NUMBER { driver.option_num("ms.shock_draws",$3); };
o_shocks_per_parameter : SHOCKS_PER_PARAMETER EQUAL INT_NUMBER { driver.option_num("ms.shocks_per_parameter",$3); };
o_free_parameters : FREE_PARAMETERS EQUAL vec_value { driver.option_vec_value("ms.free_parameters",$3); };
o_median : MEDIAN { driver.option_num("ms.median","1"); }
         | MEDIAN EQUAL signed_number { driver.option_num("median", $3); };
o_regimes : REGIMES { driver.option_num("ms.regimes","true"); };
o_regime : REGIME EQUAL INT_NUMBER { driver.option_num("ms.regime",$3); };
o_data_obs_nbr : DATA_OBS_NBR EQUAL INT_NUMBER { driver.option_num("ms.forecast_data_obs",$3); };
o_discretionary_tol: DISCRETIONARY_TOL EQUAL non_negative_number { driver.option_num("discretionary_tol",$3); };
o_analytic_derivation : ANALYTIC_DERIVATION { driver.option_num("analytic_derivation", "1"); }
o_analytic_derivation_mode : ANALYTIC_DERIVATION_MODE EQUAL signed_number { driver.option_num("analytic_derivation_mode", $3); }
o_endogenous_prior : ENDOGENOUS_PRIOR { driver.option_num("endogenous_prior", "true"); }
o_use_univariate_filters_if_singularity_is_detected : USE_UNIVARIATE_FILTERS_IF_SINGULARITY_IS_DETECTED EQUAL INT_NUMBER { driver.option_num("use_univariate_filters_if_singularity_is_detected", $3); }
o_mcmc_jumping_covariance : MCMC_JUMPING_COVARIANCE EQUAL HESSIAN
                            { driver.option_str("MCMC_jumping_covariance", $3); }                                  | MCMC_JUMPING_COVARIANCE EQUAL PRIOR_VARIANCE
                            { driver.option_str("MCMC_jumping_covariance", $3); }
                          | MCMC_JUMPING_COVARIANCE EQUAL IDENTITY_MATRIX
                            { driver.option_str("MCMC_jumping_covariance", $3); }
                          | MCMC_JUMPING_COVARIANCE EQUAL filename
                            { driver.option_str("MCMC_jumping_covariance", $3); }
                          ;
o_rescale_prediction_error_covariance : RESCALE_PREDICTION_ERROR_COVARIANCE { driver.option_num("rescale_prediction_error_covariance", "true"); };
o_use_penalized_objective_for_hessian : USE_PENALIZED_OBJECTIVE_FOR_HESSIAN { driver.option_num("hessian.use_penalized_objective","true"); };
o_irf_plot_threshold : IRF_PLOT_THRESHOLD EQUAL non_negative_number { driver.option_num("impulse_responses.plot_threshold", $3); };
o_dr_display_tol : DR_DISPLAY_TOL EQUAL non_negative_number { driver.option_num("dr_display_tol", $3); };
o_consider_all_endogenous : CONSIDER_ALL_ENDOGENOUS { driver.option_str("endo_vars_for_moment_computations_in_estimation", "all_endogenous_variables"); };
o_consider_all_endogenous_and_auxiliary : CONSIDER_ALL_ENDOGENOUS_AND_AUXILIARY { driver.option_str("endo_vars_for_moment_computations_in_estimation", "all_endogenous_and_auxiliary_variables"); };
o_consider_only_observed : CONSIDER_ONLY_OBSERVED { driver.option_str("endo_vars_for_moment_computations_in_estimation", "only_observed_variables"); };
o_no_homotopy : NO_HOMOTOPY { driver.option_num("no_homotopy", "true"); };
o_endval_steady : ENDVAL_STEADY { driver.option_num("simul.endval_steady", "true"); }
o_homotopy_max_completion_share : HOMOTOPY_MAX_COMPLETION_SHARE EQUAL non_negative_number { driver.option_num("simul.homotopy_max_completion_share", $3); }
o_homotopy_min_step_size : HOMOTOPY_MIN_STEP_SIZE EQUAL non_negative_number { driver.option_num("simul.homotopy_min_step_size", $3); }
o_homotopy_initial_step_size : HOMOTOPY_INITIAL_STEP_SIZE EQUAL non_negative_number { driver.option_num("simul.homotopy_initial_step_size", $3); }
o_homotopy_step_size_increase_success_count : HOMOTOPY_STEP_SIZE_INCREASE_SUCCESS_COUNT EQUAL INT_NUMBER { driver.option_num("simul.homotopy_step_size_increase_success_count", $3); }
o_homotopy_linearization_fallback : HOMOTOPY_LINEARIZATION_FALLBACK { driver.option_num("simul.homotopy_linearization_fallback", "true"); }
o_homotopy_marginal_linearization_fallback : HOMOTOPY_MARGINAL_LINEARIZATION_FALLBACK { driver.option_num("simul.homotopy_marginal_linearization_fallback", "0.01"); }
                                           | HOMOTOPY_MARGINAL_LINEARIZATION_FALLBACK EQUAL non_negative_number { driver.option_num("simul.homotopy_marginal_linearization_fallback", $3); }

o_infile : INFILE EQUAL filename { driver.option_str("infile", $3); };
o_invars : INVARS EQUAL '(' symbol_list ')' { driver.option_symbol_list("invars", $4); };
o_period : PERIOD EQUAL INT_NUMBER { driver.option_num("period", $3); };
o_outfile : OUTFILE EQUAL filename { driver.option_str("outfile", $3); };
o_outvars : OUTVARS EQUAL '(' symbol_list ')' { driver.option_symbol_list("outvars", $4); };
o_lmmcp : LMMCP {driver.option_num("lmmcp.status", "true"); };
o_function : FUNCTION EQUAL filename { driver.option_str("function", $3); };
o_sampling_draws : SAMPLING_DRAWS EQUAL INT_NUMBER { driver.option_num("sampling_draws",$3); };
o_use_shock_groups : USE_SHOCK_GROUPS { driver.option_str("plot_shock_decomp.use_shock_groups","default"); }
                   | USE_SHOCK_GROUPS EQUAL symbol { driver.option_str("plot_shock_decomp.use_shock_groups", $3); }
                   ;
o_colormap : COLORMAP EQUAL symbol { driver.option_num("plot_shock_decomp.colormap",$3); };
o_icd_colormap : COLORMAP EQUAL symbol { driver.option_num("initial_condition_decomp.colormap",$3); };
o_no_init_estimation_check_first_obs : NO_INIT_ESTIMATION_CHECK_FIRST_OBS { driver.option_num("no_init_estimation_check_first_obs", "true"); };
o_heteroskedastic_filter : HETEROSKEDASTIC_FILTER { driver.option_num("heteroskedastic_filter", "true"); };
o_pfwee_constant_simulation_length : CONSTANT_SIMULATION_LENGTH { driver.option_num("pfwee.constant_simulation_length", "true"); };
o_fsolve_options : FSOLVE_OPTIONS EQUAL '(' name_value_pair_with_boolean_list ')' { driver.option_str("fsolve_options", $4); };

// Some options to "method_of_moments"
o_bartlett_kernel_lag : BARTLETT_KERNEL_LAG EQUAL INT_NUMBER { driver.option_num("mom.bartlett_kernel_lag", $3); };
o_weighting_matrix : WEIGHTING_MATRIX EQUAL vec_str { driver.option_vec_cellstr("mom.weighting_matrix", $3); }
o_weighting_matrix_scaling_factor : WEIGHTING_MATRIX_SCALING_FACTOR EQUAL non_negative_number { driver.option_num("mom.weighting_matrix_scaling_factor", $3); };
o_analytic_standard_errors : ANALYTIC_STANDARD_ERRORS { driver.option_num("mom.analytic_standard_errors", "true"); };
o_analytic_jacobian : ANALYTIC_JACOBIAN { driver.option_num("mom.analytic_jacobian", "true"); };
o_mom_method : MOM_METHOD EQUAL GMM
               { driver.option_str("mom.mom_method", "GMM"); }
             | MOM_METHOD EQUAL SMM
               { driver.option_str("mom.mom_method", "SMM"); }
             | MOM_METHOD EQUAL IRF_MATCHING
               { driver.option_str("mom.mom_method", "IRF_MATCHING"); };
o_simulation_method : SIMULATION_METHOD EQUAL STOCH_SIMUL
               { driver.option_str("mom.simulation_method", "STOCH_SIMUL"); };
o_penalized_estimator : PENALIZED_ESTIMATOR { driver.option_num("mom.penalized_estimator", "true"); };
o_verbose : VERBOSE { driver.option_num("mom.verbose", "true"); };
o_simulation_multiple : SIMULATION_MULTIPLE EQUAL INT_NUMBER { driver.option_num("mom.simulation_multiple", $3); };
o_mom_burnin : BURNIN EQUAL INT_NUMBER { driver.option_num("mom.burnin", $3); };
o_bounded_shock_support : BOUNDED_SHOCK_SUPPORT { driver.option_num("mom.bounded_shock_support", "true"); };
o_mom_seed : SEED EQUAL INT_NUMBER { driver.option_num("mom.seed", $3); };
o_additional_optimizer_steps : ADDITIONAL_OPTIMIZER_STEPS EQUAL vec_int { driver.option_vec_int("additional_optimizer_steps", $3); };
o_mom_se_tolx : SE_TOLX EQUAL non_negative_number { driver.option_num("mom.se_tolx", $3); };
o_irf_matching_file : IRF_MATCHING_FILE EQUAL filename { driver.option_str("mom.irf_matching_file.name", $3); };
o_add_tiny_number_to_cholesky : ADD_TINY_NUMBER_TO_CHOLESKY EQUAL non_negative_number { driver.option_num("add_tiny_number_to_cholesky", $3); };

o_analytical_girf : ANALYTICAL_GIRF { driver.option_num("irf_opt.analytical_GIRF", "true"); };
o_irf_in_percent : IRF_IN_PERCENT { driver.option_num("irf_opt.percent", "true"); };
o_emas_girf : EMAS_GIRF { driver.option_num("irf_opt.ergodic_mean_irf", "true"); };
o_emas_drop : EMAS_DROP EQUAL INT_NUMBER { driver.option_num("irf_opt.EM.drop", $3); };
o_emas_tolf : EMAS_TOLF EQUAL non_negative_number { driver.option_num("irf_opt.EM.tolf", $3); };
o_emas_max_iter : EMAS_MAX_ITER EQUAL INT_NUMBER { driver.option_num("irf_opt.EM.iter", $3); };
o_non_zero : NON_ZERO { driver.option_num("non_zero", "true"); };

// Some options to "identification"
o_no_identification_strength : NO_IDENTIFICATION_STRENGTH { driver.option_num("no_identification_strength", "true"); };
o_no_identification_reducedform : NO_IDENTIFICATION_REDUCEDFORM { driver.option_num("no_identification_reducedform", "true"); };
o_no_identification_moments : NO_IDENTIFICATION_MOMENTS { driver.option_num("no_identification_moments", "true"); };
o_no_identification_minimal : NO_IDENTIFICATION_MINIMAL { driver.option_num("no_identification_minimal", "true"); };
o_no_identification_spectrum : NO_IDENTIFICATION_SPECTRUM { driver.option_num("no_identification_spectrum", "true"); };
o_normalize_jacobians : NORMALIZE_JACOBIANS EQUAL INT_NUMBER { driver.option_num("normalize_jacobians", $3); };
o_grid_nbr : GRID_NBR EQUAL INT_NUMBER { driver.option_num("grid_nbr", $3); };
o_tol_rank : TOL_RANK EQUAL non_negative_number { driver.option_num("tol_rank", $3); };
o_tol_deriv : TOL_DERIV EQUAL non_negative_number { driver.option_num("tol_deriv", $3); };
o_tol_sv : TOL_SV EQUAL non_negative_number { driver.option_num("tol_sv", $3); };
o_checks_via_subsets : CHECKS_VIA_SUBSETS EQUAL INT_NUMBER { driver.option_num("checks_via_subsets", $3); };
o_max_dim_subsets_groups : MAX_DIM_SUBSETS_GROUPS EQUAL INT_NUMBER { driver.option_num("max_dim_subsets_groups", $3); };
o_block_static : BLOCK_STATIC { driver.option_num("block_static","true"); };
o_block_dynamic : BLOCK_DYNAMIC { driver.option_num("block_dynamic","true"); };
o_incidence : INCIDENCE { driver.option_num("incidence","true"); };

// Some options to "evaluate_planner_objective"
o_evaluate_planner_objective_periods : PERIODS EQUAL INT_NUMBER { driver.option_num("ramsey.periods", $3); };
o_evaluate_planner_objective_drop : DROP EQUAL INT_NUMBER { driver.option_num("ramsey.drop", $3); };

// Some options to "occbin_solver"
o_occbin_simul_maxit : SIMUL_MAXIT EQUAL INT_NUMBER { driver.option_num("simul.maxit", $3); };
o_occbin_simul_periods : SIMUL_PERIODS EQUAL INT_NUMBER { driver.option_num("simul.periods", $3); };
o_occbin_simul_curb_retrench : SIMUL_CURB_RETRENCH { driver.option_num("simul.curb_retrench", "true"); };
o_occbin_simul_check_ahead_periods : SIMUL_CHECK_AHEAD_PERIODS EQUAL INT_NUMBER { driver.option_num("simul.check_ahead_periods", $3); };
o_occbin_simul_max_check_ahead_periods : SIMUL_MAX_CHECK_AHEAD_PERIODS EQUAL INT_NUMBER { driver.option_num("simul.max_check_ahead_periods", $3); };
o_occbin_simul_reset_check_ahead_periods : SIMUL_RESET_CHECK_AHEAD_PERIODS { driver.option_num("simul.reset_check_ahead_periods_in_new_period", "true"); };
o_occbin_simul_debug : SIMUL_DEBUG { driver.option_num("simul.debug", "true"); };
o_occbin_simul_periodic_solution : SIMUL_PERIODIC_SOLUTION { driver.option_num("simul.periodic_solution", "true"); };

// Some options to "occbin_setup"
o_occbin_likelihood_inversion_filter : LIKELIHOOD_INVERSION_FILTER { driver.option_num("likelihood.inversion_filter", "true"); };
o_occbin_likelihood_piecewise_kalman_filter : LIKELIHOOD_PIECEWISE_KALMAN_FILTER { driver.option_num("likelihood.inversion_filter", "false"); };
o_occbin_likelihood_maxit : LIKELIHOOD_MAXIT EQUAL INT_NUMBER { driver.option_num("likelihood.maxit", $3); };
o_occbin_likelihood_periods : LIKELIHOOD_PERIODS EQUAL INT_NUMBER { driver.option_num("likelihood.periods", $3); };
o_occbin_likelihood_curb_retrench : LIKELIHOOD_CURB_RETRENCH { driver.option_num("likelihood.curb_retrench", "true"); };
o_occbin_likelihood_check_ahead_periods : LIKELIHOOD_CHECK_AHEAD_PERIODS EQUAL INT_NUMBER { driver.option_num("likelihood.check_ahead_periods", $3); };
o_occbin_likelihood_max_check_ahead_periods : LIKELIHOOD_MAX_CHECK_AHEAD_PERIODS EQUAL INT_NUMBER { driver.option_num("likelihood.max_check_ahead_periods", $3); };
o_occbin_likelihood_periodic_solution : LIKELIHOOD_PERIODIC_SOLUTION { driver.option_num("likelihood.periodic_solution", "true"); };
o_occbin_likelihood_max_kalman_iterations : LIKELIHOOD_MAX_KALMAN_ITERATIONS EQUAL INT_NUMBER { driver.option_num("likelihood.max_number_of_iterations", $3); };
o_occbin_smoother_inversion_filter : SMOOTHER_INVERSION_FILTER { driver.option_num("smoother.inversion_filter", "true"); };
o_occbin_smoother_piecewise_kalman_filter : SMOOTHER_PIECEWISE_KALMAN_FILTER { driver.option_num("smoother.inversion_filter", "false"); };
o_occbin_smoother_maxit : SMOOTHER_MAXIT EQUAL INT_NUMBER { driver.option_num("smoother.maxit", $3); };
o_occbin_smoother_periods : SMOOTHER_PERIODS EQUAL INT_NUMBER { driver.option_num("smoother.periods", $3); };
o_occbin_smoother_curb_retrench : SMOOTHER_CURB_RETRENCH { driver.option_num("smoother.curb_retrench", "true"); };
o_occbin_smoother_check_ahead_periods : SMOOTHER_CHECK_AHEAD_PERIODS EQUAL INT_NUMBER { driver.option_num("smoother.check_ahead_periods", $3); };
o_occbin_smoother_max_check_ahead_periods : SMOOTHER_MAX_CHECK_AHEAD_PERIODS EQUAL INT_NUMBER { driver.option_num("smoother.max_check_ahead_periods", $3); };
o_occbin_smoother_debug : SMOOTHER_DEBUG { driver.option_num("smoother.debug", "true"); };
o_occbin_smoother_periodic_solution : SMOOTHER_PERIODIC_SOLUTION { driver.option_num("smoother.periodic_solution", "true"); };
o_occbin_filter_use_relaxation : FILTER_USE_RELEXATION { driver.option_num("filter.use_relaxation", "true"); };

// Some options to "occbin_write_regimes"
o_occbin_write_regimes_periods : PERIODS EQUAL vec_int
                                 { driver.option_vec_int("write_regimes.periods", $3); };
                               | PERIODS EQUAL vec_int_number
                                 { driver.option_vec_int("write_regimes.periods", $3); }
o_occbin_write_regimes_filename : FILENAME EQUAL filename { driver.option_str("write_regimes.filename", $3); };
o_occbin_write_regimes_smoother : SMOOTHER { driver.option_str("write_regimes.type", "smoother"); };
o_occbin_write_regimes_simul : SIMUL { driver.option_str("write_regimes.type", "simul"); };

// Some options to "occbin_graph"
o_occbin_graph_noconstant : NOCONSTANT { driver.option_num("graph.steady_state", "false"); };

range : symbol ':' symbol
        { $$ = $1 + ':' + $3; }

integer_range : INT_NUMBER ':' INT_NUMBER
                { $$ = $1 + ':' + $3; }

integer_range_w_inf : INT_NUMBER ':' INT_NUMBER
                      { $$ = { $1, $3 }; }
                    | INT_NUMBER ':' INF_CONSTANT
                      { $$ = { $1, "Inf" }; }
                    ;

signed_integer_range : signed_integer ':' signed_integer
                       { $$ = $1 + ':' + $3; }
                     | MINUS '(' signed_integer ':' signed_integer ')'
                       { $$ = "-(" + $3 + ':' + $5 + ")"; };

vec_int_number : INT_NUMBER
                 { $$ = { stoi($1) }; }
               ;

vec_int_elem : vec_int_number
             | INT_NUMBER ':' INT_NUMBER
               {
                 $$ = {};
                 for (int i = stoi($1); i <= stoi($3); i++)
                   $$.push_back(i);
               }
             ;

vec_int_1 : '[' vec_int_elem
            { $$ = $2;}
          | '[' COMMA vec_int_elem
            { $$ = $3;}
          | vec_int_1 vec_int_elem
            {
              $$ = $1;
              $$.insert($$.end(), $2.begin(), $2.end());
            }
          | vec_int_1 COMMA vec_int_elem
            {
              $$ = $1;
              $$.insert($$.end(), $3.begin(), $3.end());
            }
          ;

vec_int : vec_int_1 ']'
        | vec_int_1 COMMA ']'
        ;

vec_str_1 : '[' QUOTED_STRING
            { $$ = { $2 }; }
          | '[' COMMA QUOTED_STRING
            { $$ = { $3 }; }
          | vec_str_1 QUOTED_STRING
            {
              $$ = $1;
              $$.push_back($2);
            }
          | vec_str_1 COMMA QUOTED_STRING
            {
              $$ = $1;
              $$.push_back($3);
            }
          ;

vec_str : vec_str_1 ']'
        | vec_str_1 COMMA ']'
        ;

vec_value_1 : '[' signed_number
              { $$ = { $2 }; }
            | '[' COMMA signed_number
              { $$ = { $3 }; }
            | vec_value_1 signed_number
              {
                $$ = $1;
                $$.push_back($2);
              }
            | vec_value_1 COMMA signed_number
              {
                $$ = $1;
                $$.push_back($3);
              }
            ;

vec_value : vec_value_1 ']'
          | vec_value_1 COMMA ']'
          ;

vec_value_w_inf_1 : signed_number_w_inf
                    { $$ = { $1 }; }
                  | vec_value_w_inf_1 signed_number_w_inf
                    {
                      $$ = $1;
                      $$.push_back($2);
                    }
                  | vec_value_w_inf_1 COMMA signed_number_w_inf
                    {
                      $$ = $1;
                      $$.push_back($3);
                    }
                  ;

vec_value_w_inf : '[' vec_value_w_inf_1 ']'
                  { $$ = $2; }
                ;

vec_of_vec_value_1 : vec_of_vec_value_1 COMMA vec_value
                     {
                       $$ = $1;
                       $$.push_back($3);
                     }
                   | vec_value
                     { $$ = { $1 }; }
                   ;

vec_of_vec_value : '[' vec_of_vec_value_1 ']'
                   { $$ = $2; }
                 ;

symbol : NAME
       | ALPHA
       | BETA
       | NINV
       | ABAND
       | CMS
       | NCMS
       | CNUM
       | GAMMA
       | INV_GAMMA
       | INV_GAMMA1
       | INV_GAMMA2
       | NORMAL
       | UNIFORM
       | EPS
       | PDF
       | FIG
       | NONE
       | DR
       | PRIOR
       | TRUE
       | FALSE
       | BIND
       | RELAX
       | ERROR_BIND
       | ERROR_RELAX
       | KIND
       | LL
       | DL
       | DD
       | ADD
       | MULTIPLY
       ;

%%

void
Dynare::parser::error(const Dynare::parser::location_type &l,
                      const string &m)
{
  driver.error(l, m);
}
