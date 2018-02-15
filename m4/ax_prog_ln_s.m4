dnl Copyright (C) 2009-2010 Dynare Team
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

AC_DEFUN([AX_PROG_LN_S],
[
AC_REQUIRE([AC_CANONICAL_HOST])
AC_REQUIRE([AC_PROG_LN_S])
case ${host_os} in
  *cygwin*|*mingw32*)
    LN_S="cp -p" # Symbolic links are not understood by MATLAB under Windows
    ;;
esac
])
