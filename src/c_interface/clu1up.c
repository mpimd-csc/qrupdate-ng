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
  \brief Updates an LU factorization after a rank-1 modification.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_clu1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t *l, qrupdate_fortran_int_t ldl, qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_float_t *u, qrupdate_fortran_complex_float_t *v)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  CLU1UP updates an LU factorization after rank-1 modification.
  Given an m-by-k lower-triangular matrix L with unit diagonal and
  a k-by-n upper-trapezoidal matrix R, where k = min(m,n), this
  CLU1UP updates L -> L1 and R -> R1 so that L1 is again
  lower unit triangular, R1 upper trapezoidal, and
  L1*R1 = L*R + u*v', where v' denotes the conjugate transpose of
  v.

  The update is performed using the Bennett algorithm with
  column-major access, which processes the leading k-by-k block
  first and then finishes the trailing part of R if needed.
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

  \param[in,out] u
  \verbatim
           u is COMPLEX array, dimension (m)
           On entry, the left m-vector defining the rank-1
           modification.  On exit, if k < m, u is destroyed;
           otherwise, u contains the updated vector.
  \endverbatim

  \param[in,out] v
  \verbatim
           v is COMPLEX array, dimension (n)
           On entry, the right n-vector defining the rank-1
           modification.  On exit, v is destroyed.
  \endverbatim

  \ingroup c_ludecomp
 */

QRUPDATE_EXPORT void
qrupdate_clu1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
                qrupdate_fortran_complex_float_t *l, qrupdate_fortran_int_t ldl,
                qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t ldr,
                qrupdate_fortran_complex_float_t *u,
                qrupdate_fortran_complex_float_t *v) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(clu1up, CLU1UP)(&m, &n, l, &ldl, r, &ldr, u, v);
}