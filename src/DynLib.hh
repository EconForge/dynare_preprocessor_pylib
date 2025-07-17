#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// Currently unable to use this header because the used version of pybind11 is not recent enough
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

class ParserException; // needed for pybind11

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
        map<string,double> context;
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
        //! String representing the preprocessor's json output
        string json_string;
    private:
        unique_ptr<ModFile> mod_file;
        void set_mod_file(const string& modfile_string, int derivs_order, int params_derivs_order);
        void set_json_string();
        void set_symbols();
        void set_equations();
        void set_context();
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