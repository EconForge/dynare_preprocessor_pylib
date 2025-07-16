#ifndef _EXCEPTIONS_HH
#define _EXCEPTIONS_HH

#include <sstream>
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

class ParserException : public exception {
  string message;
  public:
    ParserException(int line_begin, int col_begin, int line_end, int col_end, const string& m)
      {
        ostringstream err_msg;
        err_msg << m << " at line " << line_begin << ", col " << col_begin;
        if(line_begin == line_end){
          if (col_begin < col_end -1){
            err_msg << " - " << col_end;
          }
        } else{
          err_msg << " - line " << line_end << ", col " << col_end;
        }
        message = err_msg.str();
      }
    const char* what() const noexcept {
      return message.c_str();
    }
};

#endif