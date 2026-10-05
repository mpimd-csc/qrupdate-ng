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
  \brief Updates a QR factorization after deleting a row.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_dqrder(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t *q, qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_double_t *w)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  DQRDER updates a QR factorization after deleting a row. i.e., given
  an m-by-m orthogonal matrix Q, an m-by-n upper trapezoidal matrix R
  and index j in the range 1:m, DQRDER updates Q ->Q1 and an R
  -> R1 so that Q1 is again orthogonal, R1 upper trapezoidal, and Q1*R1
  = [A(1:j-1,:); A(j+1:m,:)], where A = Q*R. (real version)
  \endverbatim

  \param[in] m
  \verbatim
           m is INTEGER
           The number of rows of the matrix Q.  m >= 1.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The number of columns of the matrix R.  n >= 0.
  \endverbatim

  \param[in,out] q
  \verbatim
           Q is DOUBLE PRECISION array, dimension (ldq,*)
           On entry, the orthogonal matrix Q.  On exit, the
           updated matrix Q1.
  \endverbatim

  \param[in] ldq
  \verbatim
           ldq is INTEGER
           The leading dimension of Q.  ldq >= m.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is DOUBLE PRECISION array, dimension (ldr,*)
           On entry, the original matrix R.  On exit, the
           updated matrix R1.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of R.  ldr >= m.
  \endverbatim

  \param[in] j
  \verbatim
           j is INTEGER
           The position of the deleted row.  1 <= j <= m.
  \endverbatim

  \param[out] w
  \verbatim
           w is DOUBLE PRECISION array, dimension (*)
           A workspace vector of size 2*m.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void
qrupdate_dqrder(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
                qrupdate_fortran_double_t *q, qrupdate_fortran_int_t ldq,
                qrupdate_fortran_double_t *r, qrupdate_fortran_int_t ldr,
                qrupdate_fortran_int_t j, qrupdate_fortran_double_t *w) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(dqrder, DQRDER)(&m, &n, q, &ldq, r, &ldr, &j, w);
}