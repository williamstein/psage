#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include "gmp.h"
#include "cstd.h"

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

unsigned long mem_count;
unsigned long mem_bytes;

void *mem_alloc (unsigned long bytes)
{
	void *ptr;

	if ( bytes < 4 ) { err_printf ("mem_alloc to small - ptr error?!\n");  exit (0); }
	if ( bytes == 4 ) { dbg_printf ("malloc 4 warning\n"); }
	ptr = malloc (bytes);
	if ( ! ptr ) { err_printf ("Fatal error, attempted memory allocation of %d bytes failed.\n");  exit (0); }
	dbg_printf ("Allocated %lu bytes at %x\n", bytes, ptr);
	memset (ptr, 0, bytes);
	mem_count++;
	mem_bytes += bytes;
	return ptr;
}

void mem_free (void *ptr)
	{ dbg_printf ("freed %x\n", ptr);  free (ptr); }
