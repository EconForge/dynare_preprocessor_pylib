#ifndef _EXCEPTIONS_HH
#define _EXCEPTIONS_HH

#include <exception>
using namespace std;

class PreprocessorException : public exception {
  string message;
  public:
    PreprocessorException(): message("preprocessor exception") {}
    PreprocessorException(string msg): message(msg) {}
    const char* what() const noexcept {
        return message.c_str();
    }
};

#endif