/*
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Copyright (C) 2026 Martin Köhler <koehlerm(AT)mpi-magdeburg.mpg.de>
 *
 * This file is part of qrupdate-ng.
 *
 * qrupdate is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "qrupdate.h"
#include "qrupdate_f77.h"

/**
  \brief Updates a row-pivoted LU factorization after rank-1 modification.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_clup1up( qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t *l, qrupdate_fortran_int_t ldl, qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t *p, qrupdate_fortran_complex_float_t *u, qrupdate_fortran_complex_float_t *v, qrupdate_fortran_complex_float_t *w)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  CLUP1UP updates a row-pivoted LU factorization after rank-1
  modification.  Given an m-by-k lower-triangular matrix L with
  unit diagonal, a k-by-n upper-trapezoidal matrix R, and a
  permutation vector p, where k = min(m,n), CLUP1UP
  updates L -> L1, R -> R1 and p -> p1 so that L1 is again
  lower unit triangular, R1 upper trapezoidal, p1 a permutation,
  and P1'*L1*R1 = P'*L*R + u*v', where v' denotes the conjugate
  transpose of v and P is the permutation matrix corresponding
  to p.
  \endverbatim

  \param[in] m
  \verbatim
           m is INTEGER
           The number of rows of the matrix L.  m >= 0.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The number of columns of the matrix R.  n >= 0.
  \endverbatim

  \param[in,out] l
  \verbatim
           L is COMPLEX array, dimension (ldl,k)
           On entry, the unit lower triangular matrix L.  On exit,
           the updated unit lower triangular matrix L1.
  \endverbatim

  \param[in] ldl
  \verbatim
           ldl is INTEGER
           The leading dimension of the array L.  ldl >= m.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is COMPLEX array, dimension (ldr,n)
           On entry, the upper trapezoidal m-by-n matrix R.
           On exit, the updated upper trapezoidal matrix R1.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of the array R.  ldr >= k,
           where k = min(m,n).
  \endverbatim

  \param[in] p
  \verbatim
           p is INTEGER array, dimension (m)
           The permutation vector representing the row pivoting.
           On exit, p is updated to reflect the new pivoting.
  \endverbatim

  \param[in] u
  \verbatim
           u is COMPLEX array, dimension (m)
           The left m-vector defining the rank-1 modification.
  \endverbatim

  \param[in] v
  \verbatim
           v is COMPLEX array, dimension (n)
           The right n-vector defining the rank-1 modification.
  \endverbatim

  \param[out] w
  \verbatim
           w is COMPLEX array, dimension (m)
           Workspace vector used during the update computation.
  \endverbatim

  \ingroup c_ludecomp
 */

QRUPDATE_EXPORT void qrupdate_clup1up(
    qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
    qrupdate_fortran_complex_float_t *l, qrupdate_fortran_int_t ldl,
    qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t ldr,
    qrupdate_fortran_int_t *p, qrupdate_fortran_complex_float_t *u,
    qrupdate_fortran_complex_float_t *v, qrupdate_fortran_complex_float_t *w) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(clup1up, CLUP1UP)(&m, &n, l, &ldl, r, &ldr, p, u, v,
                                            w);
}