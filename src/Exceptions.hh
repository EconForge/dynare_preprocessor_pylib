#ifndef _EXCEPTIONS_HH
#define _EXCEPTIONS_HH

#include <sstream>
#include <exception>


using namespace std;

inline ostringstream err_msg;

class PreprocessorException : public exception {
  string message;
  public:
    PreprocessorException(): message("Unknown preprocessor exception") {}
    PreprocessorException(string msg): message(msg) {
      err_msg.str("");
    }
    const char* what() const noexcept {
        return message.c_str();
    }
};

#endif