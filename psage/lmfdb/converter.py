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


"""
The converter module defines code for converting between Sage objects
and MongoDB documents.

DBConverter is a class whose methods are for converting native Sage
objects into MongoDB documents (JSON-like dictionaries), and
conversely.  In addition to the DBConverter, this module defines two
objects: db_converter and to_db.  The first, db_converter, is simply
an instance of the DBConverter class.  The second, to_db is the method
of the DBConverter class that converts a native Sage object to a
document.  Conversion in the other direction, from document to Sage
object, must be done by explicitly calling a method of db_converter.
This is because no type codes are stored in the database document, in
order to not tie the database too tightly to Sage.
"""

from builtins import object
class DBConverter(object):
    def from_dirichlet_character(self, e):
        zeta_order = int(e.parent().base_ring().zeta_order())
        return {'modulus':int(e.modulus()),
                'order':int(e.order()),
                'even':e.is_even(),
                'element':[int(a) for a in e.element()],
                'zeta_order':zeta_order}

    def to_dirichlet_character(self, character):
        if character['order'] == 1:
            from sage.all import trivial_character
            return trivial_character(character['modulus'])
        from sage.all import DirichletGroup, CyclotomicField, QQ
        from sage.modular.dirichlet import DirichletCharacter
        zeta_order = character['zeta_order']
        R = QQ if zeta_order == 2 else CyclotomicField(zeta_order)
        G = DirichletGroup(character['modulus'], R, zeta_order=zeta_order)
        v = G.an_element().element().parent()(character['element'])
        return DirichletCharacter(G, v)
    
    def to_db(self, x):
        from sage.modular.dirichlet import DirichletCharacter
        from sage.all import Integer
        import six
        if six.PY2:
            integer_types = (Integer, int, long)
        else:
            integer_types = (Integer, int)
        if isinstance(x, DirichletCharacter):
            return self.from_dirichlet_character(x)
        elif isinstance(x, integer_types):
            return int(x)
        elif x is None:
            return x
        raise NotImplementedError

db_converter = DBConverter()
to_db = db_converter.to_db
        
