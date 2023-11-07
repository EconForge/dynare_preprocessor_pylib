

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

#include "ParsingDriver.hh"
#include "ExtendedPreprocessorTypes.hh"
#include "ConfigFile.hh"
#include "ModFile.hh"


std::string preprocess(const std::string &modfile_string) {

    // takes the model file as a sttring
    // returns a json file with result of the preprocessing (only the parse step)

    // # we capture output completely
    std::stringstream buffer;
    std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());

    const string basename = "model";

    stringstream modfile;
    modfile << modfile_string;

    bool debug = false;
    bool no_warn = true;
    bool nostrict = true;
    
    WarningConsolidation warnings(no_warn);
    ParsingDriver p(warnings, nostrict);

    unique_ptr<ModFile> mod_file = p.parse(modfile, debug);

    JsonOutputPointType json{JsonOutputPointType::nojson};
    JsonFileOutputType json_output_mode{JsonFileOutputType::file};
  
    json = JsonOutputPointType::parsing;
    json_output_mode = JsonFileOutputType::standardout;
  
    bool onlyjson = true;

    mod_file->writeJsonOutput(basename, json, json_output_mode, true);

    std::string output = buffer.str();

    return output;

}


#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(dynare_preprocessor, m) {

    m.doc() = "Dynare preprocessor";
    m.def("preprocess", &preprocess, "Preprocess dynare model.");

}
