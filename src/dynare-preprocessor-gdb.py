# GDB pretty-printer for ExprNode class hierarchy

# Copyright © 2022-2023 Dynare Team
#
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

class ExprNodePrinter:
    '''Pretty-prints an ExprNode value'''
    def __init__(self, val):
        self.val = val

    def to_string(self):
        # Call the toString() method on the pointer.
        # We must use the raw pretty-printer for the pointer itself, otherwise
        # we enter an infinite loop.
        # We retrieve a C string, because with C++ strings the gdb pretty-printer
        # insists on keeping the quotes around it.
        r = gdb.parse_and_eval("((ExprNode *) " + self.val.format_string(raw = True) + ")->toString().c_str()")
        typestr = "(" + str(self.val.type) + ") ";
        # Add dynamic type information between brackets, if different from static type
        if str(self.val.type) != str(self.val.dynamic_type):
            typestr += "[" + str(self.val.dynamic_type) + "] "
        return typestr + r.string()

class ExprNodePrinterControl(gdb.printing.PrettyPrinter):
    '''Determines whether a value can be pretty printed with ExprNodePrinter. To be directly registered within the GDB API.'''
    def __init__(self):
        # The name below will appear in “info pretty-printer”, and can be used with “enable/disable pretty-printer”
        super().__init__('ExprNode')

    def __call__(self, val):
        # Check if the value is a subtype of ExprNode *.
        # Doing a dynamic_cast on a non-pointer type triggers an exception, so we first check
        # whether it’s a pointer (after resolving for typedefs, such as “expr_t”).
        if val.type.strip_typedefs().code == gdb.TYPE_CODE_PTR and val.dynamic_cast(gdb.lookup_type('ExprNode').pointer()) != 0:
            return ExprNodePrinter(val)

# Register the pretty printer
gdb.pretty_printers.append(ExprNodePrinterControl())
