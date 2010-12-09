#ifndef _ASM_INCLUDE_
#define _ASM_INCLUDE_

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


// This code could all be optimized - it is just a first stab
static inline unsigned long _asm_highbit (unsigned long x) { asm ("bsrq %0, %0" : "=r" (x) : "0" (x)); return x; }

#define _asm_div_q_q(q,r,x,y)				asm ("divq %4" :"=a"(q) ,"=d"(r) : "0"(x), "1"(r), "rm"(y))
#define _asm_mult_1_1(z1,z0,x0,y0)		asm ("mulq %3" :"=a"(z0) ,"=d"(z1) : "a"(x0), "rm"(y0))
#define _asm_mult_2_2_1(z1,z0,x1,x0,y0)	asm ("mulq %3" :"=a"(z0) ,"=d"(z1) : "a"(x0), "rm"(y0));(z1)+=(y0)*(x1)
#define _asm_addto_2_2(z1,z0,x1,x0)		asm ("addq %3,%0;adcq %5,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1):"cc")
#define _asm_addto_2_1(z1,z0,x0)			asm ("addq %3,%0;adcq $0,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1):"cc")
#define _asm_addto_3_3(z2,z1,z0,x2,x1,x0)	asm ("addq %4,%0;adcq %6,%1;adcq %8,%2":"=r"(z0),"=r"(z1),"=r"(z2): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1), "2"(z2), "rim"(x2) :"cc")
#define _asm_addto_3_2(z2,z1,z0,x1,x0)		asm ("addq %4,%0;adcq %6,%1;adcq 0,%2":"=r"(z0),"=r"(z1),"=r"(z2): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1), "2"(z2) :"cc")
#define _asm_subfrom_2_2(z1,z0,x1,x0)		asm ("subq %3,%0;sbbq %5,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1):"cc")
// increment needs to propogate the carry - performance not critical here anyway
#define _asm_inc_2(z1,z0)				asm ("addq $1,%0;adcq $0,%1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")
//#define _asm_shiftl_2(z1,z0)				asm ("shlq %0,1;rclq %1,1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")
//#define _asm_shiftr_2(z1,z0)				asm ("shrq %1,1;rcrq %0,1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")


#define _asm_mult_3_2_1(z2,z1,z0,x1,x0,y0)	{ register unsigned long __u; \
									   _asm_mult_1_1 (__u,z0,x0,y0); \
									   _asm_mult_1_1 (z2,z1,x1,y0); \
									   _asm_addto_2_1 (z2,z1,__u); }

// This function assumes that x[1] and y[1] < 2^31
static inline unsigned long _asm_mult_3_2_2 (unsigned long z[3], unsigned long x[2], unsigned long y[2])
{
	register unsigned long U, V, R0,R1, R2;
	
	R1 = 0;
	_asm_mult_1_1(R0,z[0],x[0],y[0]);
	_asm_mult_1_1(U,V,x[0],y[1]);
	_asm_addto_2_2(R1,R0,U,V);
	_asm_mult_1_1(U,V,x[1],y[0]);
	_asm_addto_2_2(R1,R0,U,V);
	z[1] = R0;
	z[2] = x[1]*y[1]+R1;
}

// This function assumes that x[1] and y[1] < 2^31
static inline unsigned long _asm_mult_3_2_2r (unsigned long z2, unsigned long z1, unsigned long z0, unsigned long x[2], unsigned long y[2])
{
	register unsigned long U, V, R0,R1, R2;
	
	R1 = 0;
	_asm_mult_1_1(R0,z0,x[0],y[0]);
	_asm_mult_1_1(U,V,x[0],y[1]);
	_asm_addto_2_2(R1,R0,U,V);
	_asm_mult_1_1(U,V,x[1],y[0]);
	_asm_addto_2_2(R1,R0,U,V);
	z1 = R0;
	z2 = x[1]*y[1]+R1;
}


// This function assumes that x[1] < 2^31.  For no obvious reason, this is slower than multiplying?!
static inline unsigned long _asm_square_3_2 (unsigned long z[3], unsigned long x[2])
{
	register unsigned long U, V, R0,R1, R2;
	
	R1 = 0;
	_asm_mult_1_1(R0,z[0],x[0],x[0]);
	_asm_mult_1_1(U,V,x[0],x[1]);
	_asm_addto_2_2(R1,R0,U,V);
	_asm_addto_2_2(R1,R0,U,V);
	z[1] = R0;
	z[2] = x[1]*x[1]+R1;
}
#endif
