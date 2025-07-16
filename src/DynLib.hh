#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// Currently unable to use this header because the used version of pybind11 is too outdated
// #include <pybind11/native_enum.h> 

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

//! Matrix in COO format
using coo_matrix_t = pair<vector<vector<int>>, vector<double>>;

/*! Index 0 is not used, index 1 contains first derivatives, ...
     For each derivation order, stores a matrix in COO form where coordinates are vectors of integer: the
     first integer is the equation index, the remaining ones are the derivation
     IDs of variables (in non-decreasing order, to avoid storing symmetric
     elements several times). Only non-zero derivatives are stored. */
using derivatives_t = vector<coo_matrix_t>;

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
    public:
        DynareModel (
            const string& mod_string,
            int derivs_order = 1,
            int params_derivs_order = 0
        );
        vector<string> endogenous;
        vector<string> exogenous;
        vector<string> exogenous_det;
        vector<string> parameters;
        vector<string> equations;
        map<string,double> calibration;
        map<pair<string,string>,double> covariances;
        map<string, vector<tuple<int, int, double>>> trajectories;

        // deriv ID -> (symbol type, index, lag)
        vector<tuple<SymbolType, int, int>> symbol_info;

        vector<double> dynamic_function(
            vector<double> endo_future,
            vector<double> endo_present,
            vector<double> endo_past,
            vector<double> exo,
            vector<double> exo_det,
            vector<double> params
        );

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

        /*! Returns a matrix in COO form for each derivation order (Index 0 is not used,
        index 1 contains first derivatives, ...) where coordinates are vectors of
        integer: the first integer is the equation index, the remaining ones are the
        derivation IDs of variables (in non-decreasing order, to avoid storing symmetric
        elements several times). Only non-zero derivatives are stored. */
        derivatives_t derivatives(
            vector<double> endo_future,
            vector<double> endo_present,
            vector<double> endo_past,
            vector<double> exo,
            vector<double> exo_det,
            vector<double> params
        );
    private:
        unique_ptr<ModFile> mod_file;
        void set_mod_file(const string& modfile_string, int derivs_order, int params_derivs_order);
        void set_symbols();
        void set_equations();
        void set_calibration();
        void set_exogenous();
        void set_symbolic_derivatives();
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

        vector<map<vector<int>, expr_t>> symb_derivatives;
        
        void eval_symb_jacob(
            const symb_jacobian_t& symb_jacob,
            jacobian_t& out,
            // vector of size 3 containing endogenous variables at times t-1, t and t+1 respectively
            const vector<vector<double>>& endo,
            const vector<double>& exo,
            const vector<double>& exo_det,
            const vector<double>& params
        );

};

PYBIND11_MODULE(dynare_preprocessor, m) {
    m.doc() = "dynare preprocessor";
    py::register_exception<PreprocessorException>(m, "PreprocessorException", PyExc_RuntimeError);
    py::register_exception<ParserException>(m, "ParserException", PyExc_RuntimeError);
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
    .def_readonly("endogenous", &DynareModel::endogenous)
    .def_readonly("exogenous", &DynareModel::exogenous)
    .def_readonly("exogenous_det", &DynareModel::exogenous_det)
    .def_readonly("parameters", &DynareModel::parameters)
    .def_readonly("equations", &DynareModel::equations)
    .def_readonly("calibration", &DynareModel::calibration)
    .def_readonly("covariances", &DynareModel::covariances)
    .def_readonly("trajectories", &DynareModel::trajectories)
    .def_readonly("symbol_info", &DynareModel::symbol_info)
    .def("dynamic_function", &DynareModel::dynamic_function)
    .def("jacobians", &DynareModel::jacobians)
    .def("derivatives", &DynareModel::derivatives)
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