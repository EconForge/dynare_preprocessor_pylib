#ifndef _EXCEPTIONS_HH
#define _EXCEPTIONS_HH

#include <sstream>
#include <exception>
#include "WarningConsolidation.hh"

using namespace std;

inline ostringstream err_msg;

class PreprocessorException : public exception {
  public:
    string message;
    PreprocessorException(): message("Unknown preprocessor exception") {}
    PreprocessorException(string msg): message(msg) {
      err_msg.str("");
      string prefix = "ERROR: ";
      if(message.starts_with(prefix)){
        message.erase(0, prefix.length());
      }
    }
    const char* what() const noexcept {
        return message.c_str();
    }
};

// defined in ParsingDriver.hh because dynare namespace is inaccessible from this file
// inherits from PreprocessorException
class ParserException;

class UnsupportedFeatureException : public PreprocessorException {
  public:
    UnsupportedFeatureException(string msg){
      message = msg;
    }
};

class EvalException : public PreprocessorException {
  public:
    EvalException(string msg){
      message = msg;
    }
};

#endif