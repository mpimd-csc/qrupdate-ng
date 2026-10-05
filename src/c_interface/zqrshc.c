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
  \brief Updates a QR factorization after a circular shift of columns.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_zqrshc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t *q, qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i, qrupdate_fortran_int_t j, qrupdate_fortran_complex_double_t *w, qrupdate_fortran_double_t *rw)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  ZQRSHC updates a QR factorization after circular shift of columns.
  i.e., given an m-by-k unitary matrix Q, an k-by-n upper trapezoidal
  matrix R and index j in the range 1:n+1, ZQRSHC updates the
  matrix Q -> Q1 and R -> R1 so that Q1 is again unitary, R1 upper
  trapezoidal, and Q1*R1 = A(:,p), where A = Q*R and p is the
  permutation [1:i-1,shift(i:j,-1),j+1:n] if i < j or
  [1:j-1,shift(j:i,+1),i+1:n] if j < i. (complex version)
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
           The number of columns of Q1, and rows of R1.  Must be
           either k = m (full Q) or k = n <= m (economical form).
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
           On entry, the original matrix R.  On exit, the
           updated matrix R1.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of R.  ldr >= k.
  \endverbatim

  \param[in] i
  \verbatim
           i is INTEGER
           The first index determining the range (see above).
  \endverbatim

  \param[in] j
  \verbatim
           j is INTEGER
           The second index determining the range (see above).
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

QRUPDATE_EXPORT void
qrupdate_zqrshc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
                qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t *q,
                qrupdate_fortran_int_t ldq,
                qrupdate_fortran_complex_double_t *r,
                qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
                qrupdate_fortran_int_t j, qrupdate_fortran_complex_double_t *w,
                qrupdate_fortran_double_t *rw) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(zqrshc, ZQRSHC)(&m, &n, &k, q, &ldq, r, &ldr, &i, &j,
                                          w, rw);
}