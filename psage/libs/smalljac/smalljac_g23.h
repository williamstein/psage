#ifndef _SMALLJAC_G23_INCLUDE_
#define _SMALLJAC_G23_INCLUDE_

#include "jac.h"

unsigned long smalljac_charpoly_gops;

int smalljac_genus2_charpoly (long a[3], long p, double sqrtp, long P1, long PN1);
int smalljac_genus3_charpoly (long a[3], long p, double sqrtp, long P1, long PN1, unsigned long pts);
int smalljac_genus2_charpoly_from_P1 (long a[2], long P1, long Min, long Max, curve_t c[1]);
int smalljac_genus2_charpoly_from_Pmodp (long a[2], curve_t c[1]);
int smalljac_genus3_charpoly_from_P1 (long o[3], long P1, long pts, long Min, long Max, curve_t c[1]);
int smalljac_genus3_charpoly_from_Pmodp (long a[3], curve_t c[1]);

#endif
