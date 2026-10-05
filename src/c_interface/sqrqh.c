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
    QRUPDATE_EXPORT void qrupdate_sqrqh(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t *c, qrupdate_fortran_float_t *s)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  SQRQH brings an upper trapezoidal matrix R into upper Hessenberg form
  using min(m-1,n) Givens rotations. (real version)
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
           R is REAL array, dimension (ldr,*)
           On entry, the upper Hessenberg matrix R.  On exit,
           the updated upper trapezoidal matrix.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of R.  ldr >= m.
  \endverbatim

  \param[in] c
  \verbatim
           c is REAL array, dimension (*)
           The rotation cosines.  Must contain at least
           min(m-1,n) elements.
  \endverbatim

  \param[in] s
  \verbatim
           s is REAL array, dimension (*)
           The rotation sines.  Must contain at least
           min(m-1,n) elements.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void
qrupdate_sqrqh(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
               qrupdate_fortran_float_t *r, qrupdate_fortran_int_t ldr,
               qrupdate_fortran_float_t *c, qrupdate_fortran_float_t *s) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(sqrqh, SQRQH)(&m, &n, r, &ldr, c, s);
}