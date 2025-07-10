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

enum class ExprNodeType{
    NumConstNode,
    VariableNode,
    UnaryOpNode,
    BinaryOpNode,
    TrinaryOpNode
};

ExprNodeType expression_type(expr_t expression){
    const type_info& type = typeid(*expression);
    if(type == typeid(NumConstNode))
        return ExprNodeType::NumConstNode;
    else if(type == typeid(VariableNode))
        return ExprNodeType::VariableNode;
    else if(type == typeid(UnaryOpNode))
        return ExprNodeType::UnaryOpNode;
    else if(type == typeid(BinaryOpNode))
        return ExprNodeType::BinaryOpNode;
    else if(type == typeid(TrinaryOpNode))
        return ExprNodeType::TrinaryOpNode;
    else
        throw py::type_error{"Unknown expression type"};
}


class DynareModel{
    unique_ptr<ModFile> mod_file;
    double evaluate_with_lags(
        expr_t expression,
        // vector of size 3 containing endogenous variables at times t-1, t and t+1 respectively
        const vector<vector<double>>& endo,
        const vector<double>& exo,
        const vector<double>& exo_det,
        const vector<double>& params
    );

    public:
        DynareModel (const string& mod_string);
        vector<string> endogenous;
        vector<string> exogenous;
        vector<string> exogenous_det;
        vector<string> parameters;
        vector<string> equations;
        map<string,double> calibration;
        vector<double> dynamic_function(
            vector<double> endo_p1,
            vector<double> endo_0,
            vector<double> endo_m1,
            vector<double> exo,
            vector<double> exo_det,
            vector<double> params
        );
        map<pair<string,string>,double> covariances;
        map<string, vector<tuple<int, int, double>>> trajectories;
};

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

vector<double> DynareModel::dynamic_function(
    vector<double> endo_p1,
    vector<double> endo_0,
    vector<double> endo_m1,
    vector<double> exo,
    vector<double> exo_det,
    vector<double> params
){
    vector<vector<double>> endo;
    endo.push_back(endo_m1);
    endo.push_back(endo_0);
    endo.push_back(endo_p1);
    vector<double> res;
    for(auto eq : this->mod_file->dynamic_model.equations){
        res.push_back(evaluate_with_lags(eq, endo, exo, exo_det, params));
    }
    return res;
}

DynareModel::DynareModel(const string &modfile_string) {
    stringstream modfile;
    modfile << modfile_string;
    
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
    OutputType output_mode = OutputType::standard;
    int params_derivs_order = 1;
    mod_file->computingPass(no_tmp_terms, output_mode, params_derivs_order);    

    // Get symbols
    SymbolTable table = mod_file->symbol_table;
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

    // Get equations
    equations = vector<string>();
    for(auto eq : mod_file->dynamic_model.equations){
        equations.push_back(eq->toString());
    }

    // Get calibration
    calibration = map<string,double>();
    for(const auto& [id,val] : mod_file->global_eval_context){
        calibration[table.getName(id)] = val;
    }

    // Get exogenous variable definitions
    covariances = map<pair<string,string>, double>();
    trajectories = map<string, vector<tuple<int, int, double>>>();
    for(const auto &statement : mod_file->statements){
        const type_info& type = typeid(*statement);
        if(type == typeid(ShocksStatement)){
            ShocksStatement* shock = static_cast<ShocksStatement*>(statement.get());
            for(const auto& [id, expr] : shock->var_shocks){
                string s = table.getName(id);
                covariances[make_pair(s,s)] = expr->eval(mod_file->global_eval_context);
            }
            for(const auto& [id, expr] : shock->std_shocks){
                string s = table.getName(id);
                covariances[make_pair(s,s)] = pow(expr->eval(mod_file->global_eval_context),2);
            }
            for(const auto& [key, expr] : shock->covar_shocks){
                const auto& [id1,id2] = key;
                string s1 = table.getName(id1);
                string s2 = table.getName(id2);
                double covar = expr->eval(mod_file->global_eval_context);
                covariances[make_pair(s1,s2)] = covar;
            }
            for(const auto& [key, expr] : shock->corr_shocks){
                const auto& [id1,id2] = key;
                string s1 = table.getName(id1);
                string s2 = table.getName(id2);
                double corr = expr->eval(mod_file->global_eval_context);
                double std_1 = sqrt(covariances[make_pair(s1,s1)]);
                double std_2 = sqrt(covariances[make_pair(s2,s2)]);
                double covar = corr*std_1*std_2;
                covariances[make_pair(s1,s2)] = covar;
            }
        } else if(type == typeid(ShocksSurpriseStatement)){
            ShocksSurpriseStatement* shock = static_cast<ShocksSurpriseStatement*>(statement.get());
            if(shock->overwrite){
                covariances.clear();
            }
            for(const auto& [id, trajectory] : shock->surprise_shocks){
                string var = table.getName(id);
                for(const auto& [p1, p2, expr] : trajectory){
                    double val = expr->eval(mod_file->global_eval_context);
                    trajectories[var].push_back(make_tuple(p1, p2, val));
                }
            }
        }
    }
}




PYBIND11_MODULE(dynare_preprocessor, m) {
    m.doc() = "dynare preprocessor";
    py::class_<DynareModel>(m, "DynareModel")
    .def(py::init<const string &>())
    .def_readwrite("endogenous", &DynareModel::endogenous)
    .def_readwrite("exogenous", &DynareModel::exogenous)
    .def_readwrite("exogenous_det", &DynareModel::exogenous_det)
    .def_readwrite("parameters", &DynareModel::parameters)
    .def_readwrite("equations", &DynareModel::equations)
    .def_readwrite("calibration", &DynareModel::calibration)
    .def("dynamic_function", &DynareModel::dynamic_function)
    .def_readwrite("covariances", &DynareModel::covariances)
    .def_readwrite("trajectories", &DynareModel::trajectories);
}