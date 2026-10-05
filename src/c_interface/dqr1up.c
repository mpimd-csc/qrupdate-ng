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
    QRUPDATE_EXPORT void qrupdate_dqr1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_double_t *q, qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_double_t *u, qrupdate_fortran_double_t *v, qrupdate_fortran_double_t *w)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  DQR1UP updates a QR factorization after rank-1 modification i.e.,
  given a m-by-k orthogonal Q and m-by-n upper trapezoidal R, an
  m-vector u and n-vector v, DQR1UP updates Q -> Q1 and R ->
  R1 so that Q1*R1 = Q*R + u*v', and Q1 is again orthonormal and R1
  upper trapezoidal. (real version)
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
           Q is DOUBLE PRECISION array, dimension (ldq,*)
           On entry, the orthogonal m-by-k matrix Q.  On exit,
           the updated matrix Q1.
  \endverbatim

  \param[in] ldq
  \verbatim
           ldq is INTEGER
           The leading dimension of Q.  ldq >= m.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is DOUBLE PRECISION array, dimension (ldr,*)
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
           u is DOUBLE PRECISION array, dimension (*)
           On entry, the left m-vector.  On exit, if k < m,
           u is destroyed.
  \endverbatim

  \param[in,out] v
  \verbatim
           v is DOUBLE PRECISION array, dimension (*)
           On entry, the right n-vector.  On exit, v is
           destroyed.
  \endverbatim

  \param[out] w
  \verbatim
           w is DOUBLE PRECISION array, dimension (*)
           A workspace vector of size 2*k.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void
qrupdate_dqr1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
                qrupdate_fortran_int_t k, qrupdate_fortran_double_t *q,
                qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t *r,
                qrupdate_fortran_int_t ldr, qrupdate_fortran_double_t *u,
                qrupdate_fortran_double_t *v, qrupdate_fortran_double_t *w) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(dqr1up, DQR1UP)(&m, &n, &k, q, &ldq, r, &ldr, u, v,
                                          w);
}