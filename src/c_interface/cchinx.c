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
  \brief Updates a Cholesky factorization after inserting a row and column.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_cchinx(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_complex_float_t *u, qrupdate_fortran_float_t *rw, qrupdate_fortran_int_t *info)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  CCHINX updates the Cholesky factorization of a hermitian
  positive definite matrix A after inserting a row and column.
  Given an upper triangular matrix R that is a Cholesky factor of
  A, i.e., A = R'*R, where R' denotes the conjugate transpose of
  R, CCHINX updates R -> R1 so that R1'*R1 = A1, where
  A1(jj,jj) = A, A1(j,:) = u', A1(:,j) = u, and
  jj = [1:j-1, j+1:n+1].

  On exit, u is destroyed and R is extended by one row and column.
  The insertion is performed by first solving R'*u = v, checking
  positive definiteness, and then retriangularizing.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The order of matrix R.  n >= 0.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is COMPLEX array, dimension (ldr,n+1)
           On entry, the upper triangular matrix R, the Cholesky
           factor of A.  On exit, the updated upper triangular
           matrix R1, the Cholesky factor of A1.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of the array R.  ldr >= n+1.
  \endverbatim

  \param[in] j
  \verbatim
           j is INTEGER
           The position of the inserted row and column.
           1 <= j <= n+1.
  \endverbatim

  \param[in,out] u
  \verbatim
           u is COMPLEX array, dimension (n+1)
           On entry, the vector defining the inserted row/column.
           On exit, u is destroyed.
  \endverbatim

  \param[out] rw
  \verbatim
           rw is REAL array, dimension (n)
           Workspace vector used to store rotation cosines during
           the retriangularization.
  \endverbatim

  \param[out] info
  \verbatim
           info is INTEGER
           = 0:  successful exit
           = 1:  the update would violate positive-definiteness
           = 2:  R is singular
           = 3:  the diagonal element of u is not real
  \endverbatim

  \ingroup c_choldecomp
 */

QRUPDATE_EXPORT void
qrupdate_cchinx(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t *r,
                qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
                qrupdate_fortran_complex_float_t *u,
                qrupdate_fortran_float_t *rw, qrupdate_fortran_int_t *info) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(cchinx, CCHINX)(&n, r, &ldr, &j, u, rw, info);
}