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
  \brief Updates a Cholesky factorization after a rank-1 modification.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_sch1up(qrupdate_fortran_int_t n, qrupdate_fortran_float_t *r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t *u, qrupdate_fortran_float_t *w)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  SCH1UP updates the Cholesky factorization of a symmetric
  positive definite matrix A after a rank-1 modification.  Given an
  upper triangular matrix R that is a Cholesky factor of A, i.e.,
  A = R.'*R, where R.' denotes the transpose of R, SCH1UP
  updates R -> R1 so that R1.'*R1 = A + u*u.', where u is a given
  vector.

  The update is performed by applying a sequence of Givens rotations
  to restore the upper triangular structure of R.  On exit, u
  contains the rotation sines and w contains the rotation cosines
  used in the transformation.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The order of matrix R.  n >= 0.
  \endverbatim

  \param[in,out] r
  \verbatim
           R is REAL array, dimension (ldr,n)
           On entry, the upper triangular matrix R, the Cholesky
           factor of A.  On exit, the updated upper triangular
           matrix R1, the Cholesky factor of A + u*u.'.
  \endverbatim

  \param[in] ldr
  \verbatim
           ldr is INTEGER
           The leading dimension of the array R.  ldr >= n.
  \endverbatim

  \param[in,out] u
  \verbatim
           u is REAL array, dimension (n)
           On entry, the vector determining the rank-1 update.
           On exit, u contains the rotation sines used to
           transform R to R1.
  \endverbatim

  \param[out] w
  \verbatim
           w is REAL array, dimension (n)
           On exit, w contains the cosine parts of the Givens
           rotations used to transform R to R1.
  \endverbatim

  \ingroup c_choldecomp
 */

QRUPDATE_EXPORT void qrupdate_sch1up(qrupdate_fortran_int_t n,
                                     qrupdate_fortran_float_t *r,
                                     qrupdate_fortran_int_t ldr,
                                     qrupdate_fortran_float_t *u,
                                     qrupdate_fortran_float_t *w) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(sch1up, SCH1UP)(&n, r, &ldr, u, w);
}