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
  \brief Generates a unit vector orthogonal to the column space of a unitary matrix.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_zgqvec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t *q, qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t *u)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  ZGQVEC generates a vector u in the orthogonal complement of the
  column space of a unitary matrix Q.  Given an m-by-n unitary
  matrix Q with n < m, ZGQVEC generates a vector u of
  length m such that Q'*u = 0 and norm(u) = 1, where Q' denotes
  the conjugate transpose of Q.

  The algorithm projects canonical unit vectors onto the orthogonal
  complement of Q's column space until a nonzero result is found.
  If n = 0, the first canonical unit vector is returned.
  \endverbatim

  \param[in] m
  \verbatim
           m is INTEGER
           The number of rows of the matrix Q.  m >= 0.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The number of columns of the matrix Q.  n >= 0 and
           n < m.
  \endverbatim

  \param[in] q
  \verbatim
           Q is COMPLEX*16 array, dimension (ldq,n)
           The unitary m-by-n matrix Q.
  \endverbatim

  \param[in] ldq
  \verbatim
           ldq is INTEGER
           The leading dimension of the array Q.  ldq >= m.
  \endverbatim

  \param[out] u
  \verbatim
           u is COMPLEX*16 array, dimension (m)
           The generated vector such that Q'*u = 0 and norm(u) = 1.
  \endverbatim

  \ingroup c_qrdecomp
 */

QRUPDATE_EXPORT void qrupdate_zgqvec(qrupdate_fortran_int_t m,
                                     qrupdate_fortran_int_t n,
                                     qrupdate_fortran_complex_double_t *q,
                                     qrupdate_fortran_int_t ldq,
                                     qrupdate_fortran_complex_double_t *u) {

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(zgqvec, ZGQVEC)(&m, &n, q, &ldq, u);
}