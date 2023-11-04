

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


// #include<bits/stdc++.h> 

// using namespace std; 

// extern "C" {

int just_try(int a) {

  return a+1;
}


std::string preprocess(const std::string &modfile_string) {

//   dup2(STDOUT_FILENO, STDERR_FILENO);


    // # we capture output completely
    std::stringstream buffer;
    std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
    // std::streambuf * old = std::cerr.rdbuf(buffer.rdbuf());


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
  
    bool onlyjson = false; // remark: why would writeJsonOutput ever decide to exit ?

    mod_file->writeJsonOutput(basename, json, json_output_mode, false);

    std::string output = buffer.str(); // text will now contain "Bla\n"

    return output;

}


#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(dynare_preprocessor, m) {

    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("preprocess", &preprocess, "Another one");

}

    // m.def("preprocess", &preprocess, "Another one");
    // m.def("preprocess",
    //     [](const std::string &s) {
    //         cout << "utf-8 is icing on the cake.\n";
    //         cout << s;
    //     }
    // );


//     // m.def("add", &just_try, "A function that adds two numbers");
//     // m.def("notmain", &notmain, "Another one");
//     // m.def("utf8_test", [](const std::string &s) {
//     //     cout << "utf-8 is icing on the cake.\n";
//     //     cout << s;
//     // }
// );