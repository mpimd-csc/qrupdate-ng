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
  \brief Updates a Cholesky factorization after deleting a row and column.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_dchdex(qrupdate_fortran_int_t n, qrupdate_fortran_double_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_double_t *w)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  DCHDEX updates the Cholesky factorization of a symmetric
  positive definite matrix A after deleting a row/column.
  Given an upper triangular matrix R that is a Cholesky
  factor of A, i.e., A = R.'*R, where R.' denotes the
  transpose of R, DCHDEX updates R -> R1 so that
  R1.'*R1 = A(jj,jj), where jj = [1:j-1, j+1:n+1].
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The order of matrix R.  n >= 0.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is DOUBLE PRECISION array, dimension (ldr,n)
           On entry, the upper triangular matrix R, the Cholesky
           factor of A.  On exit, the updated upper triangular
           matrix R1, the Cholesky factor of A(jj,jj).
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of the array R.  ldr >= n.
  \endverbatim

  \param[in] j
  \verbatim
           j is INTEGER
           The position of the deleted row/column.
  \endverbatim

  \param[out] w
  \verbatim
           w is DOUBLE PRECISION array, dimension (n)
           A workspace vector.
  \endverbatim

  \ingroup c_choldecomp
 */

QRUPDATE_EXPORT void qrupdate_dchdex(qrupdate_fortran_int_t n,
                                     qrupdate_fortran_double_t *r,
                                     qrupdate_fortran_int_t ldr,
                                     qrupdate_fortran_int_t j,
                                     qrupdate_fortran_double_t *w) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(dchdex, DCHDEX)(&n, r, &ldr, &j, w);
}