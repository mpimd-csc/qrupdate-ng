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
  \brief Applies a sequence of Givens rotations from the right to a matrix.

  \par C Interface:
  ==============
  \verbatim
    QRUPDATE_EXPORT void qrupdate_sqrot(char *dir, qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t *q, qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t *c, qrupdate_fortran_float_t *s)
  \endverbatim

  \par Purpose:
  =============
  \verbatim

  SQROT applies a sequence of Givens rotations from the right
  side to an m-by-n matrix Q.  Given a direction indicator
  dir, the rotation cosine and sine vectors c and s, SQROT
  applies the rotations to Q, updating it in place.  If dir
  is 'F' (forward), rotations are applied from the first to
  the last; if dir is 'B' (backward), from the last to the
  first.
  \endverbatim

  \param[in] dir
  \verbatim
           dir is CHARACTER
           If 'B' or 'b', rotations are applied backwards
           (from the last to the first).  If 'F' or 'f',
           rotations are applied forwards (from the first to
           the last).
  \endverbatim

  \param[in] m
  \verbatim
           m is INTEGER
           The number of rows of matrix Q.  m >= 0.
  \endverbatim

  \param[in] n
  \verbatim
           n is INTEGER
           The number of columns of the matrix Q.  n >= 0.
  \endverbatim

  \param[in,out] q
  \verbatim
           Q is REAL array, dimension (ldq,*)
           On entry, the matrix Q.  On exit, the updated
           matrix Q1.
  \endverbatim

  \param[in] ldq
  \verbatim
           ldq is INTEGER
           The leading dimension of Q.  ldq >= m.
  \endverbatim

  \param[in] c
  \verbatim
           c is REAL array, dimension (*)
           The rotation cosines.  Must contain at least
           n-1 elements.
  \endverbatim

  \param[in] s
  \verbatim
           s is REAL array, dimension (*)
           The rotation sines.  Must contain at least
           n-1 elements.
  \endverbatim

  \ingroup givens
 */

QRUPDATE_EXPORT void
qrupdate_sqrot(char *dir, qrupdate_fortran_int_t m, qrupdate_fortran_int_t n,
               qrupdate_fortran_float_t *q, qrupdate_fortran_int_t ldq,
               qrupdate_fortran_float_t *c, qrupdate_fortran_float_t *s) {

  char _dir[2] = {dir[0], 0};

  // Call QRUPDATE Fortran
  QRUPDATE_FORTRAN_GLOBAL(sqrot, SQROT)(_dir, &m, &n, q, &ldq, c, s, 1);
}