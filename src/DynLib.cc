#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;
namespace py = pybind11;

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

#include "ParsingDriver.hh"
#include "ExtendedPreprocessorTypes.hh"
#include "WarningConsolidation.hh"
#include "ModFile.hh"
#include "Exceptions.hh"

#include "DynLib.hh"

DynareModel::DynareModel(const string &modfile_string, int derivs_order, int params_derivs_order) {
    // Capture stderr
    ostringstream errss;
    auto cerr_original = cerr.rdbuf(errss.rdbuf());
    try{
        set_mod_file(modfile_string, derivs_order, params_derivs_order);   
        set_symbols();
        set_equations();
        set_calibration();
        set_exogenous();
        set_symbolic_derivatives();
    } catch (const PreprocessorException & ex){
        cerr.rdbuf(cerr_original);
        throw PreprocessorException(errss.str());
    } catch(...){
        cerr.rdbuf(cerr_original);
        throw;
    }
    // Stop capturing cerr
    cerr.rdbuf(cerr_original);
}


void DynareModel::set_mod_file(const string& modfile_string, int derivs_order, int params_derivs_order){
    stringstream modfile;
    modfile << modfile_string;
    cout.setstate(ios_base::failbit);
    // Do parsing and construct internal representation of mod file
    bool debug = false;
    bool no_warn = true;
    bool nostrict = true;
    WarningConsolidation warnings(no_warn);
    ParsingDriver p(warnings, nostrict);
    mod_file = p.parse(modfile, debug);    
    
    // Run checking pass
    bool stochastic = true;
    mod_file->checkPass(nostrict, stochastic);
    
    // Perform transformations on the model (creation of auxiliary vars and equations)
    bool compute_xrefs = false;
    bool transform_unary_ops = false;
    string exclude_eqs = "";
    string include_eqs = "";
    mod_file->transformPass(nostrict, stochastic, compute_xrefs,
                          transform_unary_ops, exclude_eqs, include_eqs);

    // Evaluate parameters initialization, initval and endval
    bool warn_uninit = false;
    mod_file->evalAllExpressions(warn_uninit);

    // Do computations (including derivatives)
    bool no_tmp_terms = true;
    // //! forces the preprocessor to compute derivative w.r.t. parameters
    // mod_file->mod_file_struct.identification_present = true;
    // mod_file->computingPass(no_tmp_terms, output_mode, params_derivs_order); 
    mod_file->static_model = static_cast<StaticModel>(mod_file->dynamic_model);
    mod_file->static_model.computingPass(
        derivs_order,
        params_derivs_order,
        mod_file->global_eval_context,
        no_tmp_terms,
        mod_file->block,
        mod_file->use_dll
    );
    mod_file->dynamic_model.computingPass(
        derivs_order,
        params_derivs_order,
        mod_file->global_eval_context,
        no_tmp_terms,
        mod_file->block,
        mod_file->use_dll
    );
    cout.clear();
}

void DynareModel::set_symbols(){
    const SymbolTable& table = mod_file->symbol_table;
    endogenous = vector<string>();
    for(int id: table.endo_ids){
        endogenous.push_back(table.getName(id));
    }
    exogenous = vector<string>();
    for(int id: table.exo_ids){
        exogenous.push_back(table.getName(id));
    }
    exogenous_det = vector<string>();
    for(int id: table.exo_det_ids){
        exogenous_det.push_back(table.getName(id));
    }
    parameters = vector<string>();
    for(int id: table.param_ids){
        parameters.push_back(table.getName(id));
    }
}

void DynareModel::set_equations(){
    const DynamicModel& dm = mod_file->dynamic_model;
    equations = vector<string>();
    for(expr_t eq : dm.equations){
        equations.push_back(eq->toString());
    }
}

void DynareModel::set_calibration(){
    const SymbolTable& table = mod_file->symbol_table;
    eval_context_t context = mod_file->global_eval_context;
    // Add steady state to eval context for uninitialized variables
    for(const auto& [vect, expr] : mod_file->steady_state_model.def_table){
        double val = expr->eval(context);
        for(int id : vect){
            if(context[id] == 0) context[id] = val;
        }
    }

    // Get calibration
    calibration = map<string,double>();
    for(const auto& [id,val] : context){
        calibration[table.getName(id)] = val;
    }
}


void DynareModel::set_exogenous(){
    const SymbolTable& table = mod_file->symbol_table;
    const eval_context_t& context = mod_file->global_eval_context;
    covariances = map<pair<string,string>, double>();
    trajectories = map<string, vector<tuple<int, int, double>>>();
    for(const auto &statement : mod_file->statements){
        const type_info& type = typeid(*statement);
        if(type == typeid(ShocksStatement)){
            ShocksStatement* shock = static_cast<ShocksStatement*>(statement.get());
            for(const auto& [id, expr] : shock->var_shocks){
                string s = table.getName(id);
                covariances[{s,s}] = expr->eval(context);
            }
            for(const auto& [id, expr] : shock->std_shocks){
                string s = table.getName(id);
                covariances[{s,s}] = pow(expr->eval(context),2);
            }
            for(const auto& [key, expr] : shock->covar_shocks){
                const auto& [id1,id2] = key;
                string s1 = table.getName(id1);
                string s2 = table.getName(id2);
                double covar = expr->eval(context);
                covariances[{s1,s2}] = covar;
            }
            for(const auto& [key, expr] : shock->corr_shocks){
                const auto& [id1,id2] = key;
                string s1 = table.getName(id1);
                string s2 = table.getName(id2);
                double corr = expr->eval(context);
                double std_1 = sqrt(covariances[{s1,s1}]);
                double std_2 = sqrt(covariances[{s2,s2}]);
                double covar = corr*std_1*std_2;
                covariances[{s1,s2}] = covar;
            }
        } else if(type == typeid(ShocksSurpriseStatement)){
            ShocksSurpriseStatement* shock = static_cast<ShocksSurpriseStatement*>(statement.get());
            if(shock->overwrite){
                covariances.clear();
            }
            for(const auto& [id, trajectory] : shock->surprise_shocks){
                string var = table.getName(id);
                for(const auto& [p1, p2, expr] : trajectory){
                    double val = expr->eval(context);
                    trajectories[var].push_back({p1, p2, val});
                }
            }
        }
    }
}

void DynareModel::set_symbolic_derivatives(){
    const DynamicModel& dm = mod_file->dynamic_model;
    // Get first order derivatives w.r.t variables
    symb_jacob_endo = vector<symb_jacobian_t>(3);
    symb_jacob_exo = symb_jacobian_t();
    symb_jacob_exo_det = symb_jacobian_t();
    for(const auto& [vect, expr] : dm.derivatives[1]){
        int eq = vect[0];
        for(int i = 1; i < vect.size(); i++){
            int derivID = vect[i];
            SymbolType st = dm.getTypeByDerivID(derivID);
            int tsid = dm.getTypeSpecificIDByDerivID(derivID);
            switch(st){
                case SymbolType::endogenous:
                {
                    int lag = dm.getLagByDerivID(derivID);
                    if(lag > 1 || lag < -1){
                        throw py::value_error{"Unsupported lag value"};
                    }
                    symb_jacob_endo[lag+1][{eq, tsid}] = expr;
                }
                break;
                case SymbolType::exogenous:
                    symb_jacob_exo[{eq,tsid}] = expr;
                    break;
                case SymbolType::exogenousDet:
                    symb_jacob_exo_det[{eq,tsid}] = expr;
                    break;
                default:
                    throw py::value_error{"Unknown symbol type"};
            }
        }
    }
    // Get first order derivatives wrt parameters
    symb_jacob_params = symb_jacobian_t();
    for (const auto &[indices, expr] : dm.params_derivatives.at({0,1})){
        int eq = indices[0];
        int param = dm.getTypeSpecificIDByDerivID(indices[1]);
        symb_jacob_params[{eq,param}] = expr;
    }
    // Get all derivatives w.r.t. variables
    symb_derivatives = dm.derivatives;

    // Set symbol_info table
    for(const auto& [id,lag] : dm.inv_deriv_id_table){
        SymbolType type = mod_file->symbol_table.getType(id);
        int tsid = mod_file->symbol_table.getTypeSpecificID(id);
        symbol_info.emplace_back(type, tsid, lag);
    }
}

double DynareModel::evaluate_with_lags(
    expr_t expression,
    const vector<vector<double>>& endo,
    const vector<double>& exo,
    const vector<double>& exo_det,
    const vector<double>& params
){
    switch(expression_type(expression)){
        case ExprNodeType::NumConstNode:
        {
            NumConstNode* expr = static_cast<NumConstNode*>(expression);
            return this->mod_file->num_constants.getDouble(expr->id);
        }
        case ExprNodeType::VariableNode:
        {
            VariableNode* expr = static_cast<VariableNode*>(expression);
            int lag = expr->lag;
            if(lag > 1 || lag < -1){
                throw py::value_error{"Unsupported lag value"};
            }
            int id = expr->symb_id;
            SymbolType type = mod_file->symbol_table.getType(id);
            int sid = mod_file->symbol_table.getTypeSpecificID(id);
            switch(type){
                case SymbolType::endogenous:
                    return endo[lag+1][sid];
                case SymbolType::exogenous:
                    return exo[sid];
                case SymbolType::exogenousDet:
                    return exo_det[sid];
                case SymbolType::parameter:
                    return params[sid];
                default:
                    throw py::value_error{"Unsupported variable type"};
            }
        }
        case ExprNodeType::UnaryOpNode:
        {
            UnaryOpNode* expr = static_cast<UnaryOpNode*>(expression);
            double arg = evaluate_with_lags(expr->arg, endo, exo, exo_det, params);
            return expr->eval_opcode(expr->op_code, arg);
        }
        case ExprNodeType::BinaryOpNode:
        {
            BinaryOpNode* expr = static_cast<BinaryOpNode*>(expression);
            double arg1 = evaluate_with_lags(expr->arg1, endo, exo, exo_det, params);
            double arg2 = evaluate_with_lags(expr->arg2, endo, exo, exo_det, params);
            BinaryOpcode opcode = expr->op_code;
            if(opcode == BinaryOpcode::equal){
                // special convention to make evaluation of residuals easier
                opcode = BinaryOpcode::minus;
            }
            return expr->eval_opcode(arg1, opcode, arg2, expr->powerDerivOrder);
        }
        case ExprNodeType::TrinaryOpNode:
        {
            TrinaryOpNode* expr = static_cast<TrinaryOpNode*>(expression);
            double arg1 = evaluate_with_lags(expr->arg1, endo, exo, exo_det, params);
            double arg2 = evaluate_with_lags(expr->arg2, endo, exo, exo_det, params);
            double arg3 = evaluate_with_lags(expr->arg3, endo, exo, exo_det, params);
            return expr->eval_opcode(arg1, expr->op_code, arg2, arg3);
        }
        default:
            throw py::type_error{"Unknown expression type"};
    }
}

double DynareModel::checked_evaluate_with_lags(
    expr_t expr,
    const vector<vector<double>>& endo,
    const vector<double>& exo,
    const vector<double>& exo_det,
    const vector<double>& params
){
    // Capture cerr
    ostringstream errss;
    auto cerr_original = cerr.rdbuf(errss.rdbuf());

    try{
    double res = evaluate_with_lags(expr,endo,exo,exo_det,params);
    // Stop cerr capture
    cerr.rdbuf(cerr_original);
    return res;
    }
    catch (const PreprocessorException & ex){
        cerr.rdbuf(cerr_original);
        throw PreprocessorException(errss.str());
    } catch(...){
        cerr.rdbuf(cerr_original);
        throw;
    }
}

vector<double> DynareModel::dynamic_function(
    vector<double> endo_future,
    vector<double> endo_present,
    vector<double> endo_past,
    vector<double> exo,
    vector<double> exo_det,
    vector<double> params
){
    vector<vector<double>> endo;
    endo.push_back(endo_past);
    endo.push_back(endo_present);
    endo.push_back(endo_future);
    vector<double> res;
    for(auto eq : this->mod_file->dynamic_model.equations){
        res.push_back(checked_evaluate_with_lags(eq, endo, exo, exo_det, params));
    }
    return res;
}

void DynareModel::eval_symb_jacob(
    const symb_jacobian_t& symb_jacob,
    jacobian_t& out,
    const vector<vector<double>>& endo,
    const vector<double>& exo,
    const vector<double>& exo_det,
    const vector<double>& params
){
    for(const auto& [key, expr] : symb_jacob){
        const auto& [eq,tsid] = key;
        out[{eq, tsid}] = checked_evaluate_with_lags(expr, endo, exo, exo_det, params);
    }
}

vector<jacobian_t> DynareModel::jacobians(
    vector<double> endo_future,
    vector<double> endo_present,
    vector<double> endo_past,
    vector<double> exo,
    vector<double> exo_det,
    vector<double> params
){
    vector<vector<double>> endo;
    endo.push_back(endo_past);
    endo.push_back(endo_present);
    endo.push_back(endo_future);

    vector<jacobian_t> res(6);
    eval_symb_jacob(this->symb_jacob_endo[2],res[0],endo, exo, exo_det, params);
    eval_symb_jacob(this->symb_jacob_endo[1],res[1],endo, exo, exo_det, params);
    eval_symb_jacob(this->symb_jacob_endo[0],res[2],endo, exo, exo_det, params);
    eval_symb_jacob(this->symb_jacob_exo,res[3],endo, exo, exo_det, params);
    eval_symb_jacob(this->symb_jacob_exo_det,res[4],endo, exo, exo_det, params);
    eval_symb_jacob(this->symb_jacob_params,res[5],endo, exo, exo_det, params);
    return res;
}

derivatives_t DynareModel::derivatives(
    vector<double> endo_future,
    vector<double> endo_present,
    vector<double> endo_past,
    vector<double> exo,
    vector<double> exo_det,
    vector<double> params
){
    vector<vector<double>> endo;
    endo.push_back(endo_past);
    endo.push_back(endo_present);
    endo.push_back(endo_future);

    derivatives_t res;

    for(const auto& fixed_order_derivs : symb_derivatives){
        vector<vector<int>> coords;
        vector<double> values;
        for(const auto& [vect, expr]: fixed_order_derivs){
            coords.push_back(vect);
            values.push_back(checked_evaluate_with_lags(expr, endo, exo, exo_det, params));
        }
        res.emplace_back(coords, values);
    }

    return res;
}