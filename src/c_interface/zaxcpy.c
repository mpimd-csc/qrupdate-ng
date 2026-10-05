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
  \brief Performs scaled conjugate vector addition.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_zaxcpy(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t a, qrupdate_fortran_complex_double_t *x, qrupdate_fortran_int_t incx, qrupdate_fortran_complex_double_t *y, qrupdate_fortran_int_t incy)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  ZAXCPY performs the operation y := y + a * conjg(x), where a is a
  double complex scalar, x is a vector of length n, conjg(x) denotes
  the element-wise complex conjugate of x, and y is a vector of the
  same length.  On entry, y contains the existing values; on exit, y
  is overwritten with the result.  This is the double precision
  complex analogue of the BLAS zaxpy, with the x argument conjugated
  before scaling.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The number of elements in vectors x and y.  If n <= 0,
           the subroutine returns immediately without modification.
  \endverbatim

  \param[in] a
  \verbatim
           a is COMPLEX*16
           The double complex scalar used to scale the conjugated
           vector conjg(x) before accumulation into y.
  \endverbatim

  \param[in] x
  \verbatim
           x is COMPLEX*16 array, dimension (*)
           The vector of length n whose complex conjugate is scaled
           by a and added to y.  x is not modified.
  \endverbatim

  \param[in] incx
  \verbatim
           incx is INTEGER
           The stride (increment) for elements of x.  If incx > 0,
           elements are accessed starting from x(1); if incx < 0,
           elements are accessed starting from
           x(1 + (-n+1)*incx).  A value of 1 accesses
           contiguous elements.
  \endverbatim

  \param[in,out] y
  \verbatim
           y is COMPLEX*16 array, dimension (*)
           On entry, the vector y of length n.  On exit, y is
           overwritten with y + a * conjg(x).
  \endverbatim

  \param[in] incy
  \verbatim
           incy is INTEGER
           The stride (increment) for elements of y.  If incy > 0,
           elements are accessed starting from y(1); if incy < 0,
           elements are accessed starting from
           y(1 + (-n+1)*incy).  A value of 1 accesses
           contiguous elements.
  \endverbatim
  \ingroup aux
 */

QRUPDATE_EXPORT void qrupdate_zaxcpy(qrupdate_fortran_int_t n,
                                     qrupdate_fortran_complex_double_t a,
                                     qrupdate_fortran_complex_double_t *x,
                                     qrupdate_fortran_int_t incx,
                                     qrupdate_fortran_complex_double_t *y,
                                     qrupdate_fortran_int_t incy) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(zaxcpy, ZAXCPY)(&n, &a, x, &incx, y, &incy);
}