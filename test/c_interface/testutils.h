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
 * (at your option) any later version.
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
#ifndef QRUPDATE_TESTUTILS_H
#define QRUPDATE_TESTUTILS_H

#include "qrupdate_fortran_mangle.h"

extern struct { int passed, failed; } stats_;
extern qrupdate_fortran_int_t xrand_[4];

/* random matrix/vector generators */
extern void QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float *x, qrupdate_fortran_int_t *ldx);
extern void QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double *x, qrupdate_fortran_int_t *ldx);
extern void QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float _Complex *x, qrupdate_fortran_int_t *ldx);
extern void QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double _Complex *x, qrupdate_fortran_int_t *ldx);

/* BLAS copy */
extern void QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(qrupdate_fortran_int_t *n, float *x, qrupdate_fortran_int_t *incx, float *y, qrupdate_fortran_int_t *incy);
extern void QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(qrupdate_fortran_int_t *n, double *x, qrupdate_fortran_int_t *incx, double *y, qrupdate_fortran_int_t *incy);
extern void QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(qrupdate_fortran_int_t *n, float _Complex *x, qrupdate_fortran_int_t *incx, float _Complex *y, qrupdate_fortran_int_t *incy);
extern void QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(qrupdate_fortran_int_t *n, double _Complex *x, qrupdate_fortran_int_t *incx, double _Complex *y, qrupdate_fortran_int_t *incy);

/* Cholesky generators/checkers */
extern void QRUPDATE_FORTRAN_GLOBAL(schgen,SCHGEN)(qrupdate_fortran_int_t *n, float *a, qrupdate_fortran_int_t *lda, float *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(dchgen,DCHGEN)(qrupdate_fortran_int_t *n, double *a, qrupdate_fortran_int_t *lda, double *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(cchgen,CCHGEN)(qrupdate_fortran_int_t *n, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(zchgen,ZCHGEN)(qrupdate_fortran_int_t *n, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *r, qrupdate_fortran_int_t *ldr);

extern void QRUPDATE_FORTRAN_GLOBAL(schchk,SCHCHK)(qrupdate_fortran_int_t *n, float *a, qrupdate_fortran_int_t *lda, float *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(dchchk,DCHCHK)(qrupdate_fortran_int_t *n, double *a, qrupdate_fortran_int_t *lda, double *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(cchchk,CCHCHK)(qrupdate_fortran_int_t *n, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(zchchk,ZCHCHK)(qrupdate_fortran_int_t *n, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *r, qrupdate_fortran_int_t *ldr);

/* Cholesky rank-1 update (Fortran internal) */
extern void QRUPDATE_FORTRAN_GLOBAL(sch1up,SCH1UP)(qrupdate_fortran_int_t *n, float *r, qrupdate_fortran_int_t *ldr, float *u, float *w);
extern void QRUPDATE_FORTRAN_GLOBAL(dch1up,DCH1UP)(qrupdate_fortran_int_t *n, double *r, qrupdate_fortran_int_t *ldr, double *u, double *w);
extern void QRUPDATE_FORTRAN_GLOBAL(cch1up,CCH1UP)(qrupdate_fortran_int_t *n, float _Complex *r, qrupdate_fortran_int_t *ldr, float _Complex *u, float *w);
extern void QRUPDATE_FORTRAN_GLOBAL(zch1up,ZCH1UP)(qrupdate_fortran_int_t *n, double _Complex *r, qrupdate_fortran_int_t *ldr, double _Complex *u, double *w);

/* Cholesky symmetric delete */
extern void QRUPDATE_FORTRAN_GLOBAL(schdex,SCHDEX)(qrupdate_fortran_int_t *n, float *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j, float *w);
extern void QRUPDATE_FORTRAN_GLOBAL(dchdex,DCHDEX)(qrupdate_fortran_int_t *n, double *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j, double *w);
extern void QRUPDATE_FORTRAN_GLOBAL(cchdex,CCHDEX)(qrupdate_fortran_int_t *n, float _Complex *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j, float *rw);
extern void QRUPDATE_FORTRAN_GLOBAL(zchdex,ZCHDEX)(qrupdate_fortran_int_t *n, double _Complex *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j, double *rw);

/* QR generators/checkers */
extern void QRUPDATE_FORTRAN_GLOBAL(sqrgen,SQRGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float *a, qrupdate_fortran_int_t *lda, float *q, qrupdate_fortran_int_t *ldq, float *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(dqrgen,DQRGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double *a, qrupdate_fortran_int_t *lda, double *q, qrupdate_fortran_int_t *ldq, double *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(cqrgen,CQRGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *q, qrupdate_fortran_int_t *ldq, float _Complex *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(zqrgen,ZQRGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *q, qrupdate_fortran_int_t *ldq, double _Complex *r, qrupdate_fortran_int_t *ldr);

extern void QRUPDATE_FORTRAN_GLOBAL(sqrchk,SQRCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, qrupdate_fortran_int_t *k, float *a, qrupdate_fortran_int_t *lda, float *q, qrupdate_fortran_int_t *ldq, float *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(dqrchk,DQRCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, qrupdate_fortran_int_t *k, double *a, qrupdate_fortran_int_t *lda, double *q, qrupdate_fortran_int_t *ldq, double *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(cqrchk,CQRCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, qrupdate_fortran_int_t *k, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *q, qrupdate_fortran_int_t *ldq, float _Complex *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(zqrchk,ZQRCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, qrupdate_fortran_int_t *k, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *q, qrupdate_fortran_int_t *ldq, double _Complex *r, qrupdate_fortran_int_t *ldr);

/* LU generators/checkers */
extern void QRUPDATE_FORTRAN_GLOBAL(slugen,SLUGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float *a, qrupdate_fortran_int_t *lda, float *l, qrupdate_fortran_int_t *ldl, float *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(dlugen,DLUGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double *a, qrupdate_fortran_int_t *lda, double *l, qrupdate_fortran_int_t *ldl, double *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(clugen,CLUGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *l, qrupdate_fortran_int_t *ldl, float _Complex *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(zlugen,ZLUGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *l, qrupdate_fortran_int_t *ldl, double _Complex *r, qrupdate_fortran_int_t *ldr);

extern void QRUPDATE_FORTRAN_GLOBAL(sluchk,SLUCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float *a, qrupdate_fortran_int_t *lda, float *l, qrupdate_fortran_int_t *ldl, float *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(dluchk,DLUCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double *a, qrupdate_fortran_int_t *lda, double *l, qrupdate_fortran_int_t *ldl, double *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(cluchk,CLUCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *l, qrupdate_fortran_int_t *ldl, float _Complex *r, qrupdate_fortran_int_t *ldr);
extern void QRUPDATE_FORTRAN_GLOBAL(zluchk,ZLUCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *l, qrupdate_fortran_int_t *ldl, double _Complex *r, qrupdate_fortran_int_t *ldr);

/* pivoted LU generators/checkers */
extern void QRUPDATE_FORTRAN_GLOBAL(slupgen,SLUPGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float *a, qrupdate_fortran_int_t *lda, float *l, qrupdate_fortran_int_t *ldl, float *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);
extern void QRUPDATE_FORTRAN_GLOBAL(dlupgen,DLUPGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double *a, qrupdate_fortran_int_t *lda, double *l, qrupdate_fortran_int_t *ldl, double *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);
extern void QRUPDATE_FORTRAN_GLOBAL(clupgen,CLUPGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *l, qrupdate_fortran_int_t *ldl, float _Complex *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);
extern void QRUPDATE_FORTRAN_GLOBAL(zlupgen,ZLUPGEN)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *l, qrupdate_fortran_int_t *ldl, double _Complex *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);

extern void QRUPDATE_FORTRAN_GLOBAL(slupchk,SLUPCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float *a, qrupdate_fortran_int_t *lda, float *l, qrupdate_fortran_int_t *ldl, float *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);
extern void QRUPDATE_FORTRAN_GLOBAL(dlupchk,DLUPCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double *a, qrupdate_fortran_int_t *lda, double *l, qrupdate_fortran_int_t *ldl, double *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);
extern void QRUPDATE_FORTRAN_GLOBAL(clupchk,CLUPCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, float _Complex *a, qrupdate_fortran_int_t *lda, float _Complex *l, qrupdate_fortran_int_t *ldl, float _Complex *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);
extern void QRUPDATE_FORTRAN_GLOBAL(zlupchk,ZLUPCHK)(qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n, double _Complex *a, qrupdate_fortran_int_t *lda, double _Complex *l, qrupdate_fortran_int_t *ldl, double _Complex *r, qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *p);

/* LAPACK machine precision */
extern float QRUPDATE_FORTRAN_GLOBAL(slamch,SLAMCH)(char *cmach);
extern double QRUPDATE_FORTRAN_GLOBAL(dlamch,DLAMCH)(char *cmach);

/* statistics */
extern void QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)(void);

#endif
