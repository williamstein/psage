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
r"""
Category of function fields 
"""


from sage.categories.category import Category
from sage.misc.cachefunc import cached_method
from sage.categories.basic import Fields
#from sage.rings.field import is_Field

class FunctionFields(Category):
    r"""
    The category of function fields.

    EXAMPLES:

    We create the category of function fields::

        sage: C = FunctionFields()
        sage: C
        Category of function fields

    TESTS::

        sage: TestSuite(FunctionFields()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: FunctionFields().super_categories()
            [Category of fields]
        """
        return[Fields()]

    def __contains__(self, x):
        r"""
        Returns True if ``x`` is a function field.

        EXAMPLES::

        """
        import sage.rings.all
        return sage.rings.all.is_FunctionField(x)

    def _call_(self, x):
        r"""
        Constructs an object in this category from the data in ``x``,
        or throws a TypeError.

        EXAMPLES::

            sage: C = FunctionFields()

        """
        try:
            return x.function_field()
        except AttributeError:
            raise  TypeError, "unable to canonically associate a function field to %s"%x


    class ParentMethods:
        pass

    class ElementMethods:
        pass
