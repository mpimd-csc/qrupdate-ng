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
  \brief Converts an upper trapezoidal matrix to upper Hessenberg form.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_cqrqh(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t *c, qrupdate_fortran_complex_float_t *s)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  CQRQH brings an m-by-n upper trapezoidal matrix R into upper
  Hessenberg form.  Given an m-by-n upper trapezoidal matrix R,
  CQRQH applies min(m-1,n) inverse Givens rotations
  from the right to introduce subdiagonal elements, producing an
  upper Hessenberg matrix.

  On exit, c contains the cosine parts and s contains the sine
  parts of the Givens rotations used in the transformation.
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
           R is COMPLEX array, dimension (ldr,n)
           On entry, the upper trapezoidal matrix R.  On exit, the
           upper Hessenberg matrix.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of the array R.  ldr >= m.
  \endverbatim

  \param[in] c
  \verbatim
           c is REAL array, dimension (min(m-1,n))
           The cosine parts of the Givens rotations.
  \endverbatim

  \param[in] s
  \verbatim
           s is COMPLEX array, dimension (min(m-1,n))
           The sine parts of the Givens rotations.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void qrupdate_cqrqh(qrupdate_fortran_int_t m,
                                    qrupdate_fortran_int_t n,
                                    qrupdate_fortran_complex_float_t *r,
                                    qrupdate_fortran_int_t ldr,
                                    qrupdate_fortran_float_t *c,
                                    qrupdate_fortran_complex_float_t *s) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(cqrqh, CQRQH)(&m, &n, r, &ldr, c, s);
}