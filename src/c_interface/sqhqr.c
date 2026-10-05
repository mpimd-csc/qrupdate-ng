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
  \brief Reduces an upper Hessenberg matrix to upper trapezoidal form.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_sqhqr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t *c, qrupdate_fortran_float_t *s)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  SQHQR reduces an m-by-n upper Hessenberg matrix R to upper
  trapezoidal form.  Given an m-by-n upper Hessenberg matrix R,
  SQHQR applies min(m-1,n) Givens rotations from the
  left to eliminate the subdiagonal elements, producing an upper
  trapezoidal matrix.

  On exit, c contains the cosine parts and s contains the sine
  parts of the Givens rotations used in the reduction.
  \endverbatim

  \param[in] m
  \verbatim
           m is INTEGER
           The number of rows of the matrix R.  m >= 0.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The number of columns of the matrix R.  n >= 0.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is REAL array, dimension (ldr,n)
           On entry, the upper Hessenberg matrix R.  On exit, the
           updated upper trapezoidal matrix.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of the array R.  ldr >= m.
  \endverbatim

  \param[out] c
  \verbatim
           c is REAL array, dimension (min(m-1,n))
           On exit, the cosine parts of the Givens rotations used
           to reduce R to upper trapezoidal form.
  \endverbatim

  \param[out] s
  \verbatim
           s is REAL array, dimension (min(m-1,n))
           On exit, the sine parts of the Givens rotations used
           to reduce R to upper trapezoidal form.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void
qrupdate_sqhqr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
               qrupdate_fortran_float_t *r, qrupdate_fortran_int_t ldr,
               qrupdate_fortran_float_t *c, qrupdate_fortran_float_t *s) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(sqhqr, SQHQR)(&m, &n, r, &ldr, c, s);
}