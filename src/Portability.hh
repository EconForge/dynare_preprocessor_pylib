#ifndef PORTABILITY_HH
#define PORTABILITY_HH

#ifdef _MSC_VER
#ifndef __builtin_unreachable
#define __builtin_unreachable() __assume(0)
#endif
#endif

#endif

