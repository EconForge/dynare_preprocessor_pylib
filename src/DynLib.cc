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

enum class ExprNodeType{
    NumConstNode,
    VariableNode,
    UnaryOpNode,
    BinaryOpNode,
    TrinaryOpNode
};

// (equation id, type specific id) -> symbolic partial derivative
using symb_jacobian_t = map<pair<int,int>, expr_t>;

// (equation id, type specific id) -> evaluated partial derivative
using jacobian_t = map<pair<int,int>, double>;

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
    double checked_evaluate_with_lags(
        expr_t expr,
        const vector<vector<double>>& endo,
        const vector<double>& exo,
        const vector<double>& exo_det,
        const vector<double>& params
    );
    vector<symb_jacobian_t> symb_jacob_endo;
    symb_jacobian_t symb_jacob_exo;
    symb_jacobian_t symb_jacob_exo_det;
    symb_jacobian_t symb_jacob_params;
    
    void eval_symb_jacob(
        const symb_jacobian_t& symb_jacob,
        jacobian_t& out,
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
            vector<double> endo_future,
            vector<double> endo_present,
            vector<double> endo_past,
            vector<double> exo,
            vector<double> exo_det,
            vector<double> params
        );
        map<pair<string,string>,double> covariances;
        map<string, vector<tuple<int, int, double>>> trajectories;

        // returns vector of partial derivative matrices wrt endo_future, endo_present,
        // endo_past, exo, exo_det and params vectors respectively
        vector<jacobian_t> jacobians(
            vector<double> endo_future,
            vector<double> endo_present,
            vector<double> endo_past,
            vector<double> exo,
            vector<double> exo_det,
            vector<double> params
        );
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

DynareModel::DynareModel(const string &modfile_string) {
    // Capture stderr
    ostringstream errss;
    auto cerr_original = cerr.rdbuf(errss.rdbuf());
    
    stringstream modfile;
    modfile << modfile_string;
    
    try{
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
    //! forces the preprocessor to compute derivative w.r.t. parameters
    mod_file->mod_file_struct.identification_present = true;
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
    const DynamicModel& dm = mod_file->dynamic_model;
    equations = vector<string>();
    for(expr_t eq : dm.equations){
        equations.push_back(eq->toString());
    }

    // Add steady state to eval context for uninitialized variables
    eval_context_t context = mod_file->global_eval_context;
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

    // Get exogenous variable definitions
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
    // Get derivatives wrt variables
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
    // Get derivatives wrt parameters
    symb_jacob_params = symb_jacobian_t();
    for (const auto &[indices, expr] : dm.params_derivatives.at({0,1})){
        int eq = indices[0];
        int param = dm.getTypeSpecificIDByDerivID(indices[1]);
        symb_jacob_params[{eq,param}] = expr;
    }
    } catch (const PreprocessorException & ex){
        cerr.rdbuf(cerr_original);
        throw PreprocessorException(errss.str());
    }
    // Stop capturing cerr
    cerr.rdbuf(cerr_original);
}




PYBIND11_MODULE(dynare_preprocessor, m) {
    m.doc() = "dynare preprocessor";
    py::register_exception<PreprocessorException>(m, "PreprocessorException", PyExc_RuntimeError);
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
    .def_readwrite("trajectories", &DynareModel::trajectories)
    .def("jacobians", &DynareModel::jacobians)
    .doc() = R"(The DynareModel class is initialized by passing a mod file in string format.

**Fields:**
- `endogenous`, `exogenous`, `exogenous_det` and `parameters` (_list[str]_): lists of symbols
- `equations` (_list[str]_): list of model equations in string form
- `calibration` (_dict[str,float]_): dictionary of values indexed by the model symbols, filled using steady state then initval blocks, uninitialized values are set to 0
- `covariances` (_dict[tuple[str,str],float_): contains the variances and covariances of the exogenous variables declared in Shocks block
- `trajectories` (_dict[str,list[tuple[int,int,float]]]_): contains the trajectories of exogenous variables declared in surpise shocks block, encoded as lists of elements of the form (p1, p2, v) signifying that the variable in question takes on the value v in periods p1 to p2.

**Methods:**

- `dynamic_function`:
  - Input:
    - `endo_future` (_list[float]_): list of values of endogenous variables at time t+1, given in the same order as the `endogenous` field
    - `endo_present` (_list[float]_): list of values of endogenous variables at time t, given in the same order as the `endogenous` field
    - `endo_past` (_list[float]_): list of values of endogenous variables at time t-1, given in the same order as the `endogenous` field
    - `exo` (_list[float]_): list of values of exogenous variables at time t, given in the same order as the `exogenous` field
    - `exo_det` (_list[float]_): list of values of deterministic exogenous variables at time t, given in the same order as the `exogenous` field
    - `params` (_list[float]_): list of values of parameters, given in the same order as the `exogenous` field
  - Output: list of residuals at time t (_list[float]_)


- `jacobians`:
  - Input:
    - `endo_future` (_list[float]_): list of values of endogenous variables at time t+1, given in the same order as the `endogenous` field
    - `endo_present` (_list[float]_): list of values of endogenous variables at time t, given in the same order as the `endogenous` field
    - `endo_past` (_list[float]_): list of values of endogenous variables at time t-1, given in the same order as the `endogenous` field
    - `exo` (_list[float]_): list of values of exogenous variables at time t, given in the same order as the `exogenous` field
    - `exo_det` (_list[float]_): list of values of deterministic exogenous variables at time t, given in the same order as the `exogenous` field
    - `params` (_list[float]_): list of values of parameters, given in the same order as the `exogenous` field
  - Output: list of jacobians of residuals with regards to endo_future, endo_present, endo_past, exo, exo_det and params vectors respectively evaluated at time t, each jacobian has the type dict[tuple[int,int], float]. The first index is the equation number and the second is the symbol number. For example, `parameter_jacobian[(0,0)]` is supposed to be the partial derivative of the first residual w.r.t. the first parameter and it exists as long as the partial derivative in equation is not identically equal to zero.)";
}