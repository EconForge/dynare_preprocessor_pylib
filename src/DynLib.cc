#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <thread>
#include <algorithm>
#include <filesystem>

#include <cstdlib>

#include <unistd.h>
#include <assert.h>

#include "ParsingDriver.hh"
#include "ExtendedPreprocessorTypes.hh"
#include "ConfigFile.hh"
#include "ModFile.hh"

#include <boost/algorithm/string.hpp>

std::string get_json(unique_ptr<ModFile> mod_file, JsonOutputPointType json){
    const string basename = "model";
    JsonFileOutputType json_output_mode = JsonFileOutputType::standardout;
    bool onlyjson = false; // hangs if set to true
    // # we capture output completely
    std::stringstream buffer;
    std::streambuf *old = std::cout.rdbuf(buffer.rdbuf()); // make cout's buffer point to buffer's and keep pointer to original
    std::cout.clear(); // unsilence standard output
    if(json == JsonOutputPointType::computingpass){
        bool jsonderivsimple = true;
        mod_file->writeJsonOutput(basename, json, json_output_mode, onlyjson, jsonderivsimple);
    } else{
        mod_file->writeJsonOutput(basename, json, json_output_mode, onlyjson);
    }
    buffer << std::flush;
    std::cout.rdbuf(old); // point back to original cout buffer (necessary for destructor)
    std::string output = buffer.str();
    // below needed otherwide output file would be invalid json (json 513: property expected)
    boost::replace_all(output , ", ,", ",");
    return output;
}

std::string preprocess(const std::string &modfile_string, int mode) {
    
    assert(mode >= 0 && mode <= 4 && "mode must be between 0 and 4 inclusive");
    // Allowed values for mode:
    // 0 -> no json (useless here)
    // 1 -> json generated after parsing
    // 2 -> json generated after checking
    // 3 -> json generated after transforming
    // 4 -> json generated after computing
    
    std::cout.setstate(std::ios_base::failbit); // silence standard output

    JsonOutputPointType json = static_cast<JsonOutputPointType>(mode);
    

    if(json == JsonOutputPointType::nojson) return "{}";

    stringstream modfile;
    modfile << modfile_string;
    
    // Do parsing and construct internal representation of mod file
    bool debug = false;
    bool no_warn = true;
    bool nostrict = true;
    WarningConsolidation warnings(no_warn);
    ParsingDriver p(warnings, nostrict);
    unique_ptr<ModFile> mod_file = p.parse(modfile, debug);    
    
    if(json == JsonOutputPointType::parsing){
        return get_json(std::move(mod_file), json);
    }
    // Run checking pass
    bool stochastic = true;
    mod_file->checkPass(nostrict, stochastic);
    if(json == JsonOutputPointType::checkpass){
        return get_json(std::move(mod_file), json);
    }


    // Perform transformations on the model (creation of auxiliary vars and equations)
    bool compute_xrefs = false;
    bool transform_unary_ops = false;
    std::string exclude_eqs = "";
    std::string include_eqs = "";
    mod_file->transformPass(nostrict, stochastic, compute_xrefs,
                          transform_unary_ops, exclude_eqs, include_eqs);

    if(json == JsonOutputPointType::transformpass){
        return get_json(std::move(mod_file), json);
    }
    
    // Evaluate parameters initialization, initval, endval and pounds
    bool warn_uninit = false;
    mod_file->evalAllExpressions(warn_uninit);

    // Do computations (including derivatives)
    bool no_tmp_terms = true;
    OutputType output_mode = OutputType::standard;
    int params_derivs_order = 1;
    mod_file->computingPass(no_tmp_terms, output_mode, params_derivs_order);
    return get_json(std::move(mod_file), json);
    
}


#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(dynare_preprocessor, m) {

    m.doc() = "dynare preprocessor";
    m.def("preprocess", &preprocess, "preprocess mod file using Dynare preprocessor");

}