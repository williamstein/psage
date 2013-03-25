# -*- coding=utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""

Parallel routines for computing Maass forms.

"""

cdef int compute_V_cplx_dp_sym_par(double complex **V,
                           int N1,
                           double *Xm,
                           double ***Xpb,
                           double ***Ypb,
                           double complex ***Cvec,
                           double complex *cusp_evs,
                           double *alphas,
                           int **Mv,int **Qv,double *Qfak,
                           int *symmetric_cusps,
                           double R,double Y,
                           int nc, int ncols,
                           int cuspidal,
                           int verbose,
                           int ncpus=?,
                           int is_trivial=?)
