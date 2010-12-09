#include <stdio.h>
#include <stdlib.h>
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
	unsigned long p, N;
	long m[8];
	int i, n;
	
	if ( argc < 3 ) { puts ("jgroup curve p");   puts ("   \"jgroup [1,2,3,4,5] 65537\"\n   \"jgroup x^5+3x^2-1 1000003\"");    return 0; }

	p = atol(argv[2]);

	n = smalljac_group (m, argv[1], p, SMALLJAC_SGROUP);
	if ( n <= 0 ) { 
		switch (n) {
		case 0: printf ("Specified curve has bad reduction at %lu\n", p);  break;
		case SMALLJAC_PARSE_ERROR: printf ("Unable to parse curve string: %s\n", argv[3]); break;
		case SMALLJAC_UNSUPPORTED_CURVE: puts ("Specified curve not supported - should be degree 3, 5, or 7\n");  break;
		case SMALLJAC_SINGULAR_CURVE: puts ("Specified curve is singular over Q\n");  break;
		case SMALLJAC_INVALID_PP: puts ("p must be prime");  break;
		case SMALLJAC_WRONG_GENUS: printf ("Specified curve has the wrong genus, compiled for genus %d\n", SMALLJAC_GENUS);  break;
		default: printf ("smalljac_group returned error %d\n", n);
		}
		return 0;
	}
	printf ("J(C/F_p) is isomorphic to ");
	for ( i = 0 ; i < n ; i++ ) if ( i ) printf (" x Z_%lu", m[i]);  else printf ("Z_%lu", m[i]);
	puts ("");
	if ( n > 1 ) {
		for ( i = 1, N = m[0] ; i < n ; i++ ) N *= m[i];
		printf ("Group order is %lu\n", N);
	}
}
