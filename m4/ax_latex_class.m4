# ===========================================================================
#         http://www.nongnu.org/autoconf-archive/ax_latex_class.html
# ===========================================================================
#
# OBSOLETE MACRO
#
#   Deprecated because of licensing issues. The Lesser GPL imposes licensing
#   restrictions on the generated configure script unless it is augmented
#   with an Autoconf Exception clause.
#
# SYNOPSIS
#
#   AX_LATEX_CLASS(CLASSNAME,VARIABLETOSET[,ACTION-IF-FOUND[,ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macros test is class CLASSNAME exists and work and set
#   VARIABLETOSET to yes or no If ACTION-IF-FOUND (and ACTION-IF-NOT-FOUND)
#   are set, do the correct action
#
# LICENSE
#
#   Copyright (c) 2008 Boretti Mathieu <boretti@eig.unige.ch>
#   Copyright (c) 2009 Dynare Team
#
#   This library is free software; you can redistribute it and/or modify it
#   under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation; either version 2.1 of the License, or (at
#   your option) any later version.
#
#   This library is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
#   General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with this library. If not, see <http://www.gnu.org/licenses/>.

AC_DEFUN([AX_LATEX_CLASS],[
AC_CACHE_CHECK([for usability of class $1],[ac_cv_latex_class_]translit($1,[-],[_]),[
AX_LATEX_TEST([\documentclass{$1}
\begin{document}
\end{document}],[ac_cv_latex_class_]translit($1,[-],[_]))
])
$2=$[ac_cv_latex_class_]translit($1,[-],[_]) ; export $2;
AC_SUBST($2)
ifelse($#,2,[],$#,3,[
    if test "[$]$2" = "yes" ;
    then
        $3
    fi
],$#,4,[
    ifelse($3,[],[
        if test "[$]$2" = "no" ;
        then
            $4
        fi
    ],[
        if test "[$]$2" = "yes" ;
        then
            $3
        else
            $4
        fi
    ])
])

])
