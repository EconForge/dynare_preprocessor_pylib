

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
#include "Configuration.hh"
#include "ModFile.hh"


int just_try(int a) {

  return a+1;
}

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(dynare_preprocessor, m) {

    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("just_try", &just_try, "AI Augmented Dynare Preprocessor");

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