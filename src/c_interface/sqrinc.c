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
  \brief Updates a QR factorization after inserting a new column.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_sqrinc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_float_t *q, qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_float_t *x, qrupdate_fortran_float_t *w)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  SQRINC updates a QR factorization after inserting a new column. i.e.,
  given an m-by-k orthogonal matrix Q, an m-by-n upper trapezoidal
  matrix R and index j in the range 1:n+1, SQRINC updates the
  matrix Q -> Q1 and R -> R1 so that Q1 is again orthogonal, R1 upper
  trapezoidal, and Q1*R1 = [A(:,1:j-1); x; A(:,j:n)], where A = Q*R.
  (real version)
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
           either k = m (full Q) or k = n <= m (economical form,
           basis dimension will increase).
  \endverbatim

  \param[in,out] q
  \verbatim
           Q is REAL array, dimension (ldq,*)
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
           R is REAL array, dimension (ldr,*)
           On entry, the original matrix R.  On exit, the
           updated matrix R1.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of R.  ldr >= min(m,n+1).
  \endverbatim

  \param[in] j
  \verbatim
           j is INTEGER
           The position of the new column in R1.  1 <= j <= n+1.
  \endverbatim

  \param[in] x
  \verbatim
           x is REAL array, dimension (*)
           The column being inserted.
  \endverbatim

  \param[out] w
  \verbatim
           w is REAL array, dimension (*)
           A workspace vector of size k.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void
qrupdate_sqrinc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
                qrupdate_fortran_int_t k, qrupdate_fortran_float_t *q,
                qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t *r,
                qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
                qrupdate_fortran_float_t *x, qrupdate_fortran_float_t *w) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(sqrinc, SQRINC)(&m, &n, &k, q, &ldq, r, &ldr, &j, x,
                                          w);
}