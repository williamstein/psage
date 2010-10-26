/*********************************************************************

 (c) Copyright 2006-2010 Salman Baig and Chris Hall

 This file is part of ELLFF

 ELLFF is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 ELLFF is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*********************************************************************/

class zz_pEXRat {
 private:

 public:
  zz_pEXRat();
  zz_pEXRat(const zz_PEXRat&);
  void operator=(const zz_PEXRat&);

};

// elliptic surface over F_q

class ellsurf_pEInfoT {
 private:
  ellsurf_pEInfoT();
  ellsurf_pEInfoT(const ellsurf_pEInfoT&);
  void operator=(const ellsurf_pEInfoT&);

 public: 
  class jTable {
  public:
    long ref_count;

    jTable(long q);
    ~jTable();

    zz_pEX *table; // needs to be rational functions
  
};

class ellfiber_pE {
 private:

 public:

};
