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

// necessary to make json valid
#include <boost/algorithm/string.hpp>

#include "WarningConsolidation.hh"
#include "ParsingDriver.hh"
#include "ExtendedPreprocessorTypes.hh"
#include "ModFile.hh"
#include "Exceptions.hh"

#include "DynLib.hh"

DynareModel::DynareModel(
    const string &modfile_string,
    int derivs_order,
    int params_derivs_order
){
    set_mod_file(modfile_string, derivs_order, params_derivs_order);
    set_json_string();
    set_symbols();
    set_equations();
    set_context();
    set_exogenous();
    set_symbolic_derivatives();
}


void DynareModel::set_mod_file(const string& modfile_string, int derivs_order, int params_derivs_order){
    stringstream modfile;
    modfile << modfile_string;
    // Disable standard output
    cout.setstate(ios_base::failbit);
    // Do parsing and construct internal representation of mod file
    const bool nostrict = true;
    const bool nowarn = true;
    warnings = make_unique<WarningConsolidation>(nowarn);
    driver = make_unique<ParsingDriver>(*warnings,nostrict);
    const bool debug = false;
    mod_file = driver->parse(modfile, debug);    
    
    // Run checking pass
    const bool stochastic = true;
    mod_file->checkPass(nostrict, stochastic);
    
    // Perform transformations on the model (creation of auxiliary vars and equations)
    const bool compute_xrefs = false;
    const bool transform_unary_ops = false;
    string exclude_eqs = "";
    string include_eqs = "";
    mod_file->transformPass(nostrict, stochastic, compute_xrefs,
                          transform_unary_ops, exclude_eqs, include_eqs);

    // Evaluate parameters initialization, initval and endval
    const bool warn_uninit = false;
    mod_file->evalAllExpressions(warn_uninit);

    // Do computations (including derivatives)
    const bool no_tmp_terms = true;

    mod_file->mod_file_struct.order_option = derivs_order;

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

    for (auto& statement : mod_file->statements){
        statement->computingPass(mod_file->mod_file_struct);
    }
    
    // Those matrices can only be filled here, because we use derivatives
    mod_file->dynamic_model.fillVarModelTableMatrices();

    for (auto& hm : mod_file->heterogeneous_models){
        hm.computingPass(
            derivs_order,
            no_tmp_terms,
            mod_file->use_dll
        );
    }

    //Reenable standard output
    cout.clear();
}

void DynareModel::set_json_string(){
    const string basename = "model";
    JsonOutputPointType outputpoint = JsonOutputPointType::computingpass;
    JsonFileOutputType json_output_mode = JsonFileOutputType::standardout;
    bool onlyjson = false; // exits after json output if set to true
    // we capture output completely
    std::stringstream buffer;
    std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
    bool jsonderivsimple = true;
    mod_file->writeJsonOutput(basename, outputpoint, json_output_mode, onlyjson, jsonderivsimple);
    std::cout.rdbuf(old);
    json_string = buffer.str();
    // below needed otherwide output file would be invalid json (json 513: property expected)
    boost::replace_all(json_string , ", ,", ",");
    string prefix = "//-- BEGIN JSON --// \n";
    string suffix = "\n//-- END JSON --// \nJSON written after Computing step.\n";
    if(json_string.ends_with(suffix)){
        json_string.erase(json_string.length() - suffix.length(), suffix.length());
    }
    if(json_string.starts_with(prefix)){
        json_string.erase(0, prefix.length());
    }
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

void DynareModel::set_context(){
    const SymbolTable& table = mod_file->symbol_table;
    eval_context_t new_context = mod_file->global_eval_context;
    // Add steady state to eval context for uninitialized variables
    for(const auto& [vect, expr] : mod_file->steady_state_model.def_table){
        try{
            double val = expr->eval(new_context);
            for(int id : vect){
                if(new_context[id] == 0) new_context[id] = val;
            }
        }
        catch(const ExprNode::EvalExternalFunctionException& ex){
            throw UnsupportedFeatureException("External functions are not supported (yet).");
        }
        catch(const ExprNode::EvalException& ex){
            throw EvalException("Evaluation error in steady state");
        }
    }

    // Get evaluation context
    context = map<string,double>();
    for(const auto& [id,val] : new_context){
        context[table.getName(id)] = val;
    }
}


void DynareModel::set_exogenous(){
    const SymbolTable& table = mod_file->symbol_table;
    const eval_context_t& context = mod_file->global_eval_context;
    covariances = map<pair<string,string>, double>();
    trajectories = map<string, vector<tuple<int, int, double>>>();
    try{
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
                    for(const auto& [period_range, expr] : trajectory){
                        if(holds_alternative<pair<int,int>>(period_range)){
                            auto [p1, p2] = get<pair<int,int>>(period_range);
                            double val = expr->eval(context);
                            trajectories[var].emplace_back(p1, p2, val);
                        } else {
                            throw UnsupportedFeatureException("Date period ranges are not supported (yet)");
                        }
                    }
                }
            } else if(type == typeid(NativeStatement)){
                NativeStatement* native_statement = static_cast<NativeStatement*>(statement.get());
                throw UnsupportedFeatureException("Unsupported native statement: `" + native_statement->native_statement + "`");
            }
        }
    }
    catch(const ExprNode::EvalExternalFunctionException& ex){
        throw UnsupportedFeatureException("External functions are not supported (yet).");
    }
    catch(const ExprNode::EvalException& ex){
        throw EvalException("Evaluation error in steady state");
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
                        throw UnsupportedFeatureException("The only supported lag values are 1, 0 and -1");
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
                    throw PreprocessorException("Unknown symbol type");
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

    // Insert residuals into symbolic derivatives at order 0
    map<vector<int>, expr_t> residuals;
    for(int i = 0; i < dm.equations.size(); i++){
        vector<int> coordinate = {i};
        residuals[coordinate] = dm.equations[i];
    }
    symb_derivatives[0] = residuals;
    
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
    try{
        switch(expression_type(expression)){
            case ExprNodeType::NumConstNode:
            {
                NumConstNode* expr = static_cast<NumConstNode*>(expression);
                return mod_file->num_constants.getDouble(expr->id);
            }
            case ExprNodeType::VariableNode:
            {
                VariableNode* expr = static_cast<VariableNode*>(expression);
                int lag = expr->lag;
                if(lag > 1 || lag < -1){
                    throw PreprocessorException("Unsupported lag value");
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
                        throw PreprocessorException("Unsupported variable type");
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
                throw PreprocessorException("Unknown expression type");
        }
    }
    catch(const ExprNode::EvalExternalFunctionException& ex){
        throw UnsupportedFeatureException("External functions are not supported (yet).");
    }
    catch(const ExprNode::EvalException& ex){
        throw EvalException("Evaluation error in steady state");
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
    for(auto eq : mod_file->dynamic_model.equations){
        res.push_back(evaluate_with_lags(eq, endo, exo, exo_det, params));
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
        out[{eq, tsid}] = evaluate_with_lags(expr, endo, exo, exo_det, params);
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
    eval_symb_jacob(symb_jacob_endo[2],res[0],endo, exo, exo_det, params);
    eval_symb_jacob(symb_jacob_endo[1],res[1],endo, exo, exo_det, params);
    eval_symb_jacob(symb_jacob_endo[0],res[2],endo, exo, exo_det, params);
    eval_symb_jacob(symb_jacob_exo,res[3],endo, exo, exo_det, params);
    eval_symb_jacob(symb_jacob_exo_det,res[4],endo, exo, exo_det, params);
    eval_symb_jacob(symb_jacob_params,res[5],endo, exo, exo_det, params);
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
        vector<pair<vector<int>, double>> coords;
        for(const auto& [vect, expr]: fixed_order_derivs){
            coords.emplace_back(vect, evaluate_with_lags(expr, endo, exo, exo_det, params));
        }
        res.push_back(coords);
    }

    return res;
}


PYBIND11_MODULE(dynare_preprocessor, m) {
    m.doc() = "dynare preprocessor";
    auto PyExc_Preprocessor = py::register_exception<PreprocessorException>(m, "PreprocessorException", PyExc_RuntimeError);
    py::register_exception<ParserException>(m, "ParserException", PyExc_Preprocessor.ptr());
    py::register_exception<EvalException>(m, "EvalException", PyExc_Preprocessor.ptr());
    py::register_exception<UnsupportedFeatureException>(m, "UnsupportedFeatureException", PyExc_Preprocessor.ptr());
    py::enum_<SymbolType>(m, "SymbolType", "enum.Enum")
    .value("endogenous", SymbolType::endogenous)
    .value("exogenous", SymbolType::exogenous)
    .value("exogenous_det", SymbolType::exogenousDet)
    .value("parameter", SymbolType::parameter)
    .export_values();
    py::class_<DynareModel>(m, "DynareModel")
    .def(
        py::init<const string &, int, int>(),
        py::arg(),
        py::arg("derivs_order") = 1,
        py::arg("params_derivs_order") = 0
    )
    .def_readonly("json_string", &DynareModel::json_string)
    .def_readonly("endogenous", &DynareModel::endogenous)
    .def_readonly("exogenous", &DynareModel::exogenous)
    .def_readonly("exogenous_det", &DynareModel::exogenous_det)
    .def_readonly("parameters", &DynareModel::parameters)
    .def_readonly("equations", &DynareModel::equations)
    .def_readonly("context", &DynareModel::context)
    .def_readonly("covariances", &DynareModel::covariances)
    .def_readonly("trajectories", &DynareModel::trajectories)
    .def_readonly("symbol_info", &DynareModel::symbol_info)
    .def("dynamic_function", &DynareModel::dynamic_function)
    .def("jacobians", &DynareModel::jacobians)
    .def("derivatives", &DynareModel::derivatives);
}