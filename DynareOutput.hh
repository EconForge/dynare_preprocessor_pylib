#ifndef _DYNARE_OUTPUT_HH
#define _DYNARE_OUTPUT_HH

enum OutputType
  {
    none,                             // outputs files for Matlab/Octave processing
    dynamic,                          // outputs <fname>_dynamic.cc and related files
    first,                            // outputs <fname>_first_derivatives and related files 
    second,                           // outputs <fname>_first_derivatives, <fname>_second_derivatives.cc and related files 
    third,                            // outputs <fname>_first_derivatives, <fname>_second_derivatives.cc, <fname>_third_derivatives.cc  and related files 
  };
#endif
