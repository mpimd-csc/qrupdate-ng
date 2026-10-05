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

#ifndef QRUPDATE_F77_H

#ifdef __cplusplus
extern "C" {
#endif

#include "qrupdate_config.h"
#include "qrupdate_fortran_mangle.h"
#include "qrupdate_fortran_types.h"

void QRUPDATE_FORTRAN_GLOBAL(caxcpy, CAXCPY)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_float_t *a,
    qrupdate_fortran_complex_float_t *x, qrupdate_fortran_int_t *incx,
    qrupdate_fortran_complex_float_t *y, qrupdate_fortran_int_t *incy);
void QRUPDATE_FORTRAN_GLOBAL(cch1dn, CCH1DN)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_complex_float_t *u,
    qrupdate_fortran_float_t *rw, qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(cch1up,
                             CCH1UP)(qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_float_t *r,
                                     qrupdate_fortran_int_t *ldr,
                                     qrupdate_fortran_complex_float_t *u,
                                     qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(cchdex,
                             CCHDEX)(qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_float_t *r,
                                     qrupdate_fortran_int_t *ldr,
                                     qrupdate_fortran_int_t *j,
                                     qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cchinx, CCHINX)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_complex_float_t *u, qrupdate_fortran_float_t *rw,
    qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(cchshx, CCHSHX)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *i,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_float_t *w,
    qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cgqvec,
                             CGQVEC)(qrupdate_fortran_int_t *m,
                                     qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_float_t *q,
                                     qrupdate_fortran_int_t *ldq,
                                     qrupdate_fortran_complex_float_t *u);
void QRUPDATE_FORTRAN_GLOBAL(clu1up, CLU1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_float_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_complex_float_t *u, qrupdate_fortran_complex_float_t *v);
void QRUPDATE_FORTRAN_GLOBAL(clup1up, CLUP1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_float_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *p, qrupdate_fortran_complex_float_t *u,
    qrupdate_fortran_complex_float_t *v, qrupdate_fortran_complex_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(cqhqr, CQHQR)(qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_complex_float_t *r,
                                           qrupdate_fortran_int_t *ldr,
                                           qrupdate_fortran_float_t *c,
                                           qrupdate_fortran_complex_float_t *s);
void QRUPDATE_FORTRAN_GLOBAL(cqr1up, CQR1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_complex_float_t *u,
    qrupdate_fortran_complex_float_t *v, qrupdate_fortran_complex_float_t *w,
    qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cqrdec, CQRDEC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cqrder, CQRDER)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_float_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_float_t *w,
    qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cqrinc, CQRINC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_complex_float_t *x, qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cqrinr, CQRINR)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_float_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_complex_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_float_t *x,
    qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cqrot, CQROT)(char *dir, qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_complex_float_t *q,
                                           qrupdate_fortran_int_t *ldq,
                                           qrupdate_fortran_float_t *c,
                                           qrupdate_fortran_complex_float_t *s,
                                           qrupdate_fortran_charlen_t len0);
void QRUPDATE_FORTRAN_GLOBAL(cqrqh, CQRQH)(qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_complex_float_t *r,
                                           qrupdate_fortran_int_t *ldr,
                                           qrupdate_fortran_float_t *c,
                                           qrupdate_fortran_complex_float_t *s);
void QRUPDATE_FORTRAN_GLOBAL(cqrshc, CQRSHC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *i,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_float_t *w,
    qrupdate_fortran_float_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(cqrtv1,
                             CQRTV1)(qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_float_t *u,
                                     qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dch1dn, DCH1DN)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_double_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_double_t *u,
                                             qrupdate_fortran_double_t *w,
                                             qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(dch1up, DCH1UP)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_double_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_double_t *u,
                                             qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dchdex, DCHDEX)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_double_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_int_t *j,
                                             qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dchinx, DCHINX)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_double_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_int_t *j,
                                             qrupdate_fortran_double_t *u,
                                             qrupdate_fortran_double_t *w,
                                             qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(dchshx, DCHSHX)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_double_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_int_t *i,
                                             qrupdate_fortran_int_t *j,
                                             qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dgqvec, DGQVEC)(qrupdate_fortran_int_t *m,
                                             qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_double_t *q,
                                             qrupdate_fortran_int_t *ldq,
                                             qrupdate_fortran_double_t *u);
void QRUPDATE_FORTRAN_GLOBAL(dlu1up, DLU1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_double_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_double_t *u, qrupdate_fortran_double_t *v);
void QRUPDATE_FORTRAN_GLOBAL(dlup1up, DLUP1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_double_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *p, qrupdate_fortran_double_t *u,
    qrupdate_fortran_double_t *v, qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dqhqr, DQHQR)(qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_double_t *r,
                                           qrupdate_fortran_int_t *ldr,
                                           qrupdate_fortran_double_t *c,
                                           qrupdate_fortran_double_t *s);
void QRUPDATE_FORTRAN_GLOBAL(dqr1up, DQR1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_double_t *u,
    qrupdate_fortran_double_t *v, qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dqrdec, DQRDEC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dqrder, DQRDER)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_double_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dqrinc, DQRINC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_double_t *x, qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dqrinr, DQRINR)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_double_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_double_t *x,
    qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dqrot, DQROT)(char *dir, qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_double_t *q,
                                           qrupdate_fortran_int_t *ldq,
                                           qrupdate_fortran_double_t *c,
                                           qrupdate_fortran_double_t *s,
                                           qrupdate_fortran_charlen_t len0);
void QRUPDATE_FORTRAN_GLOBAL(dqrqh, DQRQH)(qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_double_t *r,
                                           qrupdate_fortran_int_t *ldr,
                                           qrupdate_fortran_double_t *c,
                                           qrupdate_fortran_double_t *s);
void QRUPDATE_FORTRAN_GLOBAL(dqrshc, DQRSHC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *i,
    qrupdate_fortran_int_t *j, qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(dqrtv1, DQRTV1)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_double_t *u,
                                             qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sch1dn, SCH1DN)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_float_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_float_t *u,
                                             qrupdate_fortran_float_t *w,
                                             qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(sch1up, SCH1UP)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_float_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_float_t *u,
                                             qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(schdex, SCHDEX)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_float_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_int_t *j,
                                             qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(schinx, SCHINX)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_float_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_int_t *j,
                                             qrupdate_fortran_float_t *u,
                                             qrupdate_fortran_float_t *w,
                                             qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(schshx, SCHSHX)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_float_t *r,
                                             qrupdate_fortran_int_t *ldr,
                                             qrupdate_fortran_int_t *i,
                                             qrupdate_fortran_int_t *j,
                                             qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sgqvec, SGQVEC)(qrupdate_fortran_int_t *m,
                                             qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_float_t *q,
                                             qrupdate_fortran_int_t *ldq,
                                             qrupdate_fortran_float_t *u);
void QRUPDATE_FORTRAN_GLOBAL(slu1up, SLU1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_float_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_float_t *u, qrupdate_fortran_float_t *v);
void QRUPDATE_FORTRAN_GLOBAL(slup1up, SLUP1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_float_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *p, qrupdate_fortran_float_t *u,
    qrupdate_fortran_float_t *v, qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sqhqr, SQHQR)(qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_float_t *r,
                                           qrupdate_fortran_int_t *ldr,
                                           qrupdate_fortran_float_t *c,
                                           qrupdate_fortran_float_t *s);
void QRUPDATE_FORTRAN_GLOBAL(sqr1up, SQR1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_float_t *u,
    qrupdate_fortran_float_t *v, qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sqrdec, SQRDEC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sqrder, SQRDER)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_float_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sqrinc, SQRINC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_float_t *x, qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sqrinr, SQRINR)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_float_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_float_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_float_t *x,
    qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sqrot, SQROT)(char *dir, qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_float_t *q,
                                           qrupdate_fortran_int_t *ldq,
                                           qrupdate_fortran_float_t *c,
                                           qrupdate_fortran_float_t *s,
                                           qrupdate_fortran_charlen_t len0);
void QRUPDATE_FORTRAN_GLOBAL(sqrqh, SQRQH)(qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_float_t *r,
                                           qrupdate_fortran_int_t *ldr,
                                           qrupdate_fortran_float_t *c,
                                           qrupdate_fortran_float_t *s);
void QRUPDATE_FORTRAN_GLOBAL(sqrshc, SQRSHC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_float_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_float_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *i,
    qrupdate_fortran_int_t *j, qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(sqrtv1, SQRTV1)(qrupdate_fortran_int_t *n,
                                             qrupdate_fortran_float_t *u,
                                             qrupdate_fortran_float_t *w);
void QRUPDATE_FORTRAN_GLOBAL(zaxcpy, ZAXCPY)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_double_t *a,
    qrupdate_fortran_complex_double_t *x, qrupdate_fortran_int_t *incx,
    qrupdate_fortran_complex_double_t *y, qrupdate_fortran_int_t *incy);
void QRUPDATE_FORTRAN_GLOBAL(zch1dn, ZCH1DN)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_complex_double_t *u,
    qrupdate_fortran_double_t *rw, qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(zch1up,
                             ZCH1UP)(qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_double_t *r,
                                     qrupdate_fortran_int_t *ldr,
                                     qrupdate_fortran_complex_double_t *u,
                                     qrupdate_fortran_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(zchdex,
                             ZCHDEX)(qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_double_t *r,
                                     qrupdate_fortran_int_t *ldr,
                                     qrupdate_fortran_int_t *j,
                                     qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zchinx, ZCHINX)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_complex_double_t *u, qrupdate_fortran_double_t *rw,
    qrupdate_fortran_int_t *info);
void QRUPDATE_FORTRAN_GLOBAL(zchshx, ZCHSHX)(
    qrupdate_fortran_int_t *n, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *i,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_double_t *w,
    qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zgqvec,
                             ZGQVEC)(qrupdate_fortran_int_t *m,
                                     qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_double_t *q,
                                     qrupdate_fortran_int_t *ldq,
                                     qrupdate_fortran_complex_double_t *u);
void QRUPDATE_FORTRAN_GLOBAL(zlu1up, ZLU1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_double_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_complex_double_t *u, qrupdate_fortran_complex_double_t *v);
void QRUPDATE_FORTRAN_GLOBAL(zlup1up, ZLUP1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_double_t *l, qrupdate_fortran_int_t *ldl,
    qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *p, qrupdate_fortran_complex_double_t *u,
    qrupdate_fortran_complex_double_t *v, qrupdate_fortran_complex_double_t *w);
void QRUPDATE_FORTRAN_GLOBAL(zqhqr, ZQHQR)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_double_t *c, qrupdate_fortran_complex_double_t *s);
void QRUPDATE_FORTRAN_GLOBAL(zqr1up, ZQR1UP)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_complex_double_t *u,
    qrupdate_fortran_complex_double_t *v, qrupdate_fortran_complex_double_t *w,
    qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zqrdec, ZQRDEC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zqrder, ZQRDER)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_double_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_double_t *w,
    qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zqrinc, ZQRINC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *j,
    qrupdate_fortran_complex_double_t *x, qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zqrinr, ZQRINR)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_double_t *q, qrupdate_fortran_int_t *ldq,
    qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_double_t *x,
    qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zqrot, ZQROT)(char *dir, qrupdate_fortran_int_t *m,
                                           qrupdate_fortran_int_t *n,
                                           qrupdate_fortran_complex_double_t *q,
                                           qrupdate_fortran_int_t *ldq,
                                           qrupdate_fortran_double_t *c,
                                           qrupdate_fortran_complex_double_t *s,
                                           qrupdate_fortran_charlen_t len0);
void QRUPDATE_FORTRAN_GLOBAL(zqrqh, ZQRQH)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_complex_double_t *r, qrupdate_fortran_int_t *ldr,
    qrupdate_fortran_double_t *c, qrupdate_fortran_complex_double_t *s);
void QRUPDATE_FORTRAN_GLOBAL(zqrshc, ZQRSHC)(
    qrupdate_fortran_int_t *m, qrupdate_fortran_int_t *n,
    qrupdate_fortran_int_t *k, qrupdate_fortran_complex_double_t *q,
    qrupdate_fortran_int_t *ldq, qrupdate_fortran_complex_double_t *r,
    qrupdate_fortran_int_t *ldr, qrupdate_fortran_int_t *i,
    qrupdate_fortran_int_t *j, qrupdate_fortran_complex_double_t *w,
    qrupdate_fortran_double_t *rw);
void QRUPDATE_FORTRAN_GLOBAL(zqrtv1,
                             ZQRTV1)(qrupdate_fortran_int_t *n,
                                     qrupdate_fortran_complex_double_t *u,
                                     qrupdate_fortran_double_t *w);

#ifdef __cplusplus
}
#endif

#endif