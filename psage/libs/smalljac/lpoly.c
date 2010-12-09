#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "smalljac.h"

/*
    Copyright 2007 Andrew V. Sutherland

    This file is part of smalljac.

    smalljac is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    smalljac is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with smalljac.  If not, see <http://www.gnu.org/licenses/>.
*/

// simple command-line program to compute the L-polynomial of a J(C/F_p)

int main (int argc, char *argv[])
{
	unsigned long p;
	long a[3];
	int sts, a1_only;
	
	if ( argc < 3 ) { puts ("lpoly curve q [a1-only-flag]");   puts ("   \"lpoly [1,2,3,4,5] 65537\"");    puts("    \"lpoly x^5+3x^3-2x+1 10000019\""); return 0; }

	p = atoll(argv[2]);
	a1_only = ( argc > 3 ? atoi(argv[3]) : 0 );

	sts = smalljac_Lpoly (a, argv[1], p, a1_only);
	if ( sts <= 0 ) { 
		switch (sts) {
		case 0: printf ("Specified curve has bad reduction at %lu\n", p);  break;
		case SMALLJAC_PARSE_ERROR: printf ("Unable to parse curve string: %s\n", argv[3]); break;
		case SMALLJAC_UNSUPPORTED_CURVE: puts ("Specified curve not supported - should be degree 3, 5, or 7\n");  break;
		case SMALLJAC_SINGULAR_CURVE: puts ("Specified curve is singular over Q\n");  break;
		case SMALLJAC_INVALID_PP: puts ("q must be a prime power");  break;
		case SMALLJAC_WRONG_GENUS: printf ("Specified curve has the wrong genus, compiled for genus %d\n", SMALLJAC_GENUS);  break;
		default: printf ("smalljac_Lpoly returned error %d\n", sts);
		}
		return 0;
	}
	if ( a1_only ) {
		printf ("a1 coefficient of L_p(T) is %ld\n", a[0]);
	} else {
		switch ( sts ) {
		case 1: printf ("L_q(T) = qT^2 + %ldT + 1\n", a[0]);  break;
		case 2: printf ("L_q(T) = q^2T^4 + %ldqT^3 + %ldT^2 + %ldT + 1\n", a[0], a[1], a[0]);  break;
		case 3: printf ("L_q(T) = q^3T^6 + %ldq^2T^5 + %ldqT^4 + %ldT^3 + %ldT^2 + %ldT + 1\n", a[0], a[1], a[2], a[1], a[0]);  break;
		default: printf ("Lpoly returned unexpected value %d\n", sts);
		}
	}
}
