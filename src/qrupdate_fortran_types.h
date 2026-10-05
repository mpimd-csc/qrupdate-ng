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


#ifndef QRUPDATE_FORTRAN_TYPES_H
#define QRUPDATE_FORTRAN_TYPES_H
#include <stdint.h>
#include <stdlib.h>
#include <complex.h>

#include "qrupdate_config.h"

#ifndef QRUPDATE_CHARLEN_T
#define QRUPDATE_CHARLEN_T

#if defined(__INTEL_LLVM_COMPILER) || defined(__ICC)
/* Intel Compiler (oneAPI and classic) */
typedef size_t qrupdate_fortran_charlen_t;
#elif defined (__PGI) || defined(__NVCOMPILER)
typedef int qrupdate_fortran_charlen_t;
#elif defined(__aocc__)
/* CLANG/FLANG as AMD AOCC */
typedef size_t qrupdate_fortran_charlen_t;
#elif defined(__clang__)
/* CLANG/FLANG */
typedef size_t qrupdate_fortran_charlen_t;
#elif __GNUC__ > 7
/* GNU 8.x and newer */
typedef size_t qrupdate_fortran_charlen_t;
#else
/* GNU 4.x - 7.x */
typedef int qrupdate_fortran_charlen_t;
#endif
#endif

/* Define the LOGICAL dataype */
#ifndef qrupdate_fortran_logical_t
#if __GNUC__ >= 5 && !defined (__clang__)
#ifdef QRUPDATE_INTEGER8
#define qrupdate_fortran_logical_t int_fast64_t
#else
#define qrupdate_fortran_logical_t int_least32_t
#endif
#else
#ifdef QRUPDATE_INTEGER8
#define qrupdate_fortran_logical_t int64_t
#else
#define qrupdate_fortran_logical_t int
#endif
#endif
#endif

/* Define the BLAS integer */
#ifndef qrupdate_fortran_int_t
#ifdef QRUPDATE_INTEGER8
#define qrupdate_fortran_int_t int64_t
#else
#define qrupdate_fortran_int_t int32_t
#endif
#endif

/* Define the DOUBLE type */
#ifndef qrupdate_fortran_double_t
#define qrupdate_fortran_double_t double
#endif

/* Define the FLOAT type */
#ifndef qrupdate_fortran_float_t
#define qrupdate_fortran_float_t float
#endif


/* Define the complex*16 */
#ifndef qrupdate_fortran_complex_double_t
#define qrupdate_fortran_complex_double_t double complex
#endif


/* Define the complex*8 */
#ifndef qrupdate_fortran_complex_float_t
#define qrupdate_fortran_complex_float_t float complex
#endif

#endif

