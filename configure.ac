dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2009-2018 Dynare Team
dnl
dnl This file is part of Dynare.
dnl
dnl Dynare is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl Dynare is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.62])
AC_INIT([dynare-preprocessor], [4.6-unstable])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([src/DynareMain.cc])
AM_INIT_AUTOMAKE([1.11 -Wall -Wno-portability foreign no-dist-gzip dist-xz tar-pax])

AC_PROG_CC
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_14

AC_CANONICAL_HOST
case ${host_os} in
  *mingw32*)
    # On Windows, we don't want dynamic linking
    AM_LDFLAGS="-static"
    AC_SUBST([AM_LDFLAGS])
    ;;
  *cygwin*)
    AC_MSG_WARN([You are compiling for the Cygwin target. This means that the preprocessor will])
    AC_MSG_WARN([not run from MATLAB unless you add the Cygwin DLL to the path.])
    AC_MSG_WARN([This is probably not what you want. Consider using a MinGW cross-compiler.])
    ;;
esac

# Use C++ for testing headers
AC_LANG([C++])

AM_CXXFLAGS="-Wall -Wno-parentheses"
AC_SUBST([AM_CXXFLAGS])

# If default 'ar' is not available, try to find one with a host prefix (see ticket #145)
AC_CHECK_PROGS([AR], [ar ${host_alias}-ar])

AC_PROG_RANLIB
AM_PROG_AR

AM_PROG_LEX
# Hack to get lex include dir, ticket #575
AC_PATH_PROG([LEXPATH], [$LEX])
AC_SUBST([LEXINC], [$(eval "echo $LEXPATH | sed 's|\(.*\)$LEX$|\1../include|'")])

AC_CHECK_PROG([YACC], [bison], [bison])
if test -z "$YACC"; then
  unset YACC # AM_MISSING_PROG needs an unset variable: an empty variable won't do
  AM_MISSING_PROG([YACC], [bison])
fi

# We need 1.36 because of unordered_{set,hash} used by Dynare++
AX_BOOST_BASE([1.36], [], [AC_MSG_ERROR([Can't find Boost >= 1.36])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM

CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp], [], [AC_MSG_ERROR([Can't find Boost Graph Library])])
AC_CHECK_HEADERS([boost/algorithm/string/trim.hpp], [], [AC_MSG_ERROR([Can't find Boost String Library])])
AC_CHECK_HEADERS([boost/algorithm/string/split.hpp], [], [AC_MSG_ERROR([Can't find Boost String Library])])
CPPFLAGS="$CPPFLAGS_SAVED"

# Don't use deprecated hash structures
AC_DEFINE([BOOST_NO_HASH], [], [Don't use deprecated STL hash structures])

AC_CHECK_PROG([DOXYGEN], [doxygen], [doxygen])
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

AC_CHECK_PROG([PDFLATEX], [pdflatex], [pdflatex])
AM_CONDITIONAL([HAVE_PDFLATEX], [test -n "$PDFLATEX"])

AC_CHECK_PROG([BIBTEX], [bibtex], [bibtex])
AM_CONDITIONAL([HAVE_BIBTEX], [test -n "$BIBTEX"])

if test -n "$PDFLATEX" -a -n "$BIBTEX"; then
  AX_LATEX_CLASS([beamer], [ax_latex_have_beamer])
fi
AM_CONDITIONAL([HAVE_BEAMER], [test "$ax_latex_have_beamer" = yes])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/macro/Makefile
                 doc/Makefile
                 doc/preprocessor/Makefile
                 doc/macroprocessor/Makefile
])

if test -n "$DOXYGEN"; then
  BUILD_DYNARE_PREPROC_DOC="yes"
else
  BUILD_DYNARE_PREPROC_DOC="no (missing doxygen)"
fi


if test -n "$PDFLATEX" -a "$ax_latex_have_beamer" = yes; then
  BUILD_BEAMER_DOC="yes"
else
  BUILD_BEAMER_DOC="no (missing one of: pdflatex, beamer)"
fi

AC_MSG_NOTICE([

The Dynare Preprocessor is now configured to build...

Binaries (with "make")
 Dynare preprocessor:                         yes

PDF documentation (with "make pdf"):
 Preprocessor & Macroprocessor presentations: $BUILD_BEAMER_DOC

HTML documentation (with "make html"):
 Dynare preprocessor developer doc:           $BUILD_DYNARE_PREPROC_DOC

])

AC_OUTPUT
