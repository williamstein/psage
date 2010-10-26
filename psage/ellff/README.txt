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

ELLFF (for elliptic curve L-function over function fields) computes
L-functions for non-constant elliptic curves over F_q(t) with
gcd(q,6)=1. It is a C++ library built on top of Victor Shoup's NTL
library and hooks into Sage and PSage via Cython.

Its files are as follows:

    * ell_surface.cpp/h: Contain the C++ class for elliptic curves
      over a function field, i.e. elliptic surfaces over a finite
      constant field, along with basic associated functions and
      context switching (a la NTL).

    * ell.cpp/h: Contain the C++ class for elliptic curves over finite
      fields and the requisite functionality needed to count
      points. These curves are assumed to have an elliptic surface
      context.

    * ellff.cpp: Auto-generated C++ code from ellff.pyx.

    * ellff.pyx: Cython interface for ELLFF. In particular, contains
      functions to create an elliptic curve over a funtion field as a
      Sage object, compute its L-function, form pullbacks and
      quadratic twists, and interface with the table of Euler factors
      (i.e. traces of Frobenius).

    * euler.cpp/h: C++ code to compute tables of Euler factors for a
      curve, along with a pullback or quadratic twist of the
      curve. This file also computes the coefficients used to compute
      the L-function.

    * helper.cpp/h: Contains various helper routines for working over
      a polynomial ring over a finite field: order elements of F_p[x],
      F_q, F_q[x]; convert elements of F_p[x], F_q, F_q[x] to unsigned
      longs; enumerate elements in F_p[x], F_q, F_q[x]; evaluate
      elements of F_q[x] at element of F_q; and exponentiation in
      F_q^*

   * interface.h: <SHB: This should be deprecated I think>

    * jacobi.cpp/h: Computes the (abstract) Jacobi sum for elements in
      F_q^*

    * lzz_pEExtra.cpp/h: Contains various helper routines for working
      over F_q: table of multiplicative inverses and square roots;
      evaluation of the Legendre character; the action of Frobenius;
      and context-switching

    * lzz_pEratX.cpp/h: Contains code for elements in F_q(t) as
      quotients of elements in F_q[t], building on top of NTL's
      lzz_pEX code.

