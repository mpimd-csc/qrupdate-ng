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
  \brief Updates a QR factorization after a rank-1 modification.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_zqr1up( qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t *q, qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_double_t *u, qrupdate_fortran_complex_double_t *v, qrupdate_fortran_complex_double_t *w, qrupdate_fortran_double_t *rw)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  ZQR1UP updates a QR factorization after rank-1 modification i.e.,
  given a m-by-k unitary Q and m-by-n upper trapezoidal R, an m-vector
  u and n-vector v, ZQR1UP updates Q -> Q1 and R -> R1 so that
  Q1*R1 = Q*R + u*v', and Q1 is again unitary and R1 upper trapezoidal.
  (complex version)
  \endverbatim

  \param[in] m
  \verbatim
           m is INTEGER
           The number of rows of the matrix Q.  m >= 0.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The number of columns of the matrix R.  n >= 0.
  \endverbatim

  \param[in] k
  \verbatim
           k is INTEGER
           The number of columns of Q, and rows of R.  Must be
           either k = m (full Q) or k = n < m (economical form).
  \endverbatim

  \param[in,out] q
  \verbatim
           Q is COMPLEX*16 array, dimension (ldq,*)
           On entry, the unitary m-by-k matrix Q.  On exit,
           the updated matrix Q1.
  \endverbatim

  \param[in] ldq
  \verbatim
           ldq is INTEGER
           The leading dimension of Q.  ldq >= m.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is COMPLEX*16 array, dimension (ldr,*)
           On entry, the upper trapezoidal m-by-n matrix R.  On
           exit, the updated matrix R1.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of R.  ldr >= k.
  \endverbatim

  \param[in,out] u
  \verbatim
           u is COMPLEX*16 array, dimension (*)
           On entry, the left m-vector.  On exit, if k < m,
           u is destroyed.
  \endverbatim

  \param[in,out] v
  \verbatim
           v is COMPLEX*16 array, dimension (*)
           On entry, the right n-vector.  On exit, v is
           destroyed.
  \endverbatim

  \param[out] w
  \verbatim
           w is COMPLEX*16 array, dimension (*)
           A workspace vector of size k.
  \endverbatim

  \param[out] rw
  \verbatim
           rw is DOUBLE PRECISION array, dimension (*)
           A real workspace vector of size k.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void qrupdate_zqr1up(
    qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
    qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t *q,
    qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_double_t *u,
    qrupdate_fortran_complex_double_t *v, qrupdate_fortran_complex_double_t *w,
    qrupdate_fortran_double_t *rw) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(zqr1up, ZQR1UP)(&m, &n, &k, q, &ldq, r, &ldr, u, v, w,
                                          rw);
}