#################################################################################
#
# (c) Copyright 2010 William Stein
#
#  This file is part of PSAGE
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################


import function_field

def FunctionField(X, names=None):
    """
    Return the function field defined by X.

    INPUT:

        - `X` -- a field; return the function field in one variable over X.

        - ``names`` -- name of variable as a string
    
    EXAMPLES::

        sage: FunctionField(QQ,'alpha')
        Rational function field in alpha over Rational Field
        sage: K.<alpha> = FunctionField(GF(7)); K
        Rational function field in alpha over Finite Field of size 7
    """
    return function_field.RationalFunctionField(X, names=names)
    
