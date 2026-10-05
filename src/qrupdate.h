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

#ifndef QRUPDATE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "qrupdate_export.h"
#include "qrupdate_config.h"
#include "qrupdate_fortran_mangle.h"
#include "qrupdate_fortran_types.h"
#include "qrupdate_error.h"

    QRUPDATE_EXPORT void qrupdate_caxcpy(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t a, qrupdate_fortran_complex_float_t* x, qrupdate_fortran_int_t incx,
            qrupdate_fortran_complex_float_t* y, qrupdate_fortran_int_t incy);

    QRUPDATE_EXPORT void qrupdate_cch1dn(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_float_t* u,
            qrupdate_fortran_float_t* rw, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_cch1up(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_float_t* u,
            qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_cchdex(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cchinx(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_complex_float_t* u, qrupdate_fortran_float_t* rw, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_cchshx(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_complex_float_t* w, qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cgqvec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_complex_float_t* u);

    QRUPDATE_EXPORT void qrupdate_clu1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_float_t* u, qrupdate_fortran_complex_float_t* v);

    QRUPDATE_EXPORT void qrupdate_clup1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t* p, qrupdate_fortran_complex_float_t* u,
            qrupdate_fortran_complex_float_t* v, qrupdate_fortran_complex_float_t* w);

    QRUPDATE_EXPORT void qrupdate_cqhqr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_float_t* c, qrupdate_fortran_complex_float_t* s);

    QRUPDATE_EXPORT void qrupdate_cqr1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_float_t* u,
            qrupdate_fortran_complex_float_t* v, qrupdate_fortran_complex_float_t* w, qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cqrdec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cqrder(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_complex_float_t* w,
            qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cqrinc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_complex_float_t* x, qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cqrinr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_complex_float_t* x,
            qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cqrot(char* dir, qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t* c, qrupdate_fortran_complex_float_t* s);

    QRUPDATE_EXPORT void qrupdate_cqrqh(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_float_t* c, qrupdate_fortran_complex_float_t* s);

    QRUPDATE_EXPORT void qrupdate_cqrshc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_complex_float_t* w, qrupdate_fortran_float_t* rw);

    QRUPDATE_EXPORT void qrupdate_cqrtv1(qrupdate_fortran_int_t n, qrupdate_fortran_complex_float_t* u, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_dch1dn(qrupdate_fortran_int_t n, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_double_t* u,
            qrupdate_fortran_double_t* w, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_dch1up(qrupdate_fortran_int_t n, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_double_t* u,
            qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dchdex(qrupdate_fortran_int_t n, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dchinx(qrupdate_fortran_int_t n, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_double_t* u, qrupdate_fortran_double_t* w, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_dchshx(qrupdate_fortran_int_t n, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dgqvec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_double_t* u);

    QRUPDATE_EXPORT void qrupdate_dlu1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_double_t* u, qrupdate_fortran_double_t* v);

    QRUPDATE_EXPORT void qrupdate_dlup1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t* p, qrupdate_fortran_double_t* u,
            qrupdate_fortran_double_t* v, qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dqhqr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_double_t* c, qrupdate_fortran_double_t* s);

    QRUPDATE_EXPORT void qrupdate_dqr1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_double_t* u,
            qrupdate_fortran_double_t* v, qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dqrdec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dqrder(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dqrinc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_double_t* x, qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dqrinr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_double_t* x,
            qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dqrot(char* dir, qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t* c, qrupdate_fortran_double_t* s);

    QRUPDATE_EXPORT void qrupdate_dqrqh(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_double_t* c, qrupdate_fortran_double_t* s);

    QRUPDATE_EXPORT void qrupdate_dqrshc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_dqrtv1(qrupdate_fortran_int_t n, qrupdate_fortran_double_t* u, qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_sch1dn(qrupdate_fortran_int_t n, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t* u,
            qrupdate_fortran_float_t* w, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_sch1up(qrupdate_fortran_int_t n, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t* u,
            qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_schdex(qrupdate_fortran_int_t n, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_schinx(qrupdate_fortran_int_t n, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_float_t* u, qrupdate_fortran_float_t* w, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_schshx(qrupdate_fortran_int_t n, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sgqvec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_float_t* u);

    QRUPDATE_EXPORT void qrupdate_slu1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t* u, qrupdate_fortran_float_t* v);

    QRUPDATE_EXPORT void qrupdate_slup1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t* p, qrupdate_fortran_float_t* u,
            qrupdate_fortran_float_t* v, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sqhqr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_float_t* c, qrupdate_fortran_float_t* s);

    QRUPDATE_EXPORT void qrupdate_sqr1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_float_t* u,
            qrupdate_fortran_float_t* v, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sqrdec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sqrder(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sqrinc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_float_t* x, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sqrinr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_float_t* x,
            qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sqrot(char* dir, qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t* c, qrupdate_fortran_float_t* s);

    QRUPDATE_EXPORT void qrupdate_sqrqh(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_float_t* c, qrupdate_fortran_float_t* s);

    QRUPDATE_EXPORT void qrupdate_sqrshc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_float_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_float_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_sqrtv1(qrupdate_fortran_int_t n, qrupdate_fortran_float_t* u, qrupdate_fortran_float_t* w);

    QRUPDATE_EXPORT void qrupdate_zaxcpy(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t a, qrupdate_fortran_complex_double_t* x, qrupdate_fortran_int_t incx,
            qrupdate_fortran_complex_double_t* y, qrupdate_fortran_int_t incy);

    QRUPDATE_EXPORT void qrupdate_zch1dn(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_double_t* u,
            qrupdate_fortran_double_t* rw, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_zch1up(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_double_t* u,
            qrupdate_fortran_double_t* w);

    QRUPDATE_EXPORT void qrupdate_zchdex(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zchinx(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_complex_double_t* u, qrupdate_fortran_double_t* rw, qrupdate_fortran_int_t* info);

    QRUPDATE_EXPORT void qrupdate_zchshx(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_complex_double_t* w, qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zgqvec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_complex_double_t* u);

    QRUPDATE_EXPORT void qrupdate_zlu1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_double_t* u, qrupdate_fortran_complex_double_t* v);

    QRUPDATE_EXPORT void qrupdate_zlup1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* l, qrupdate_fortran_int_t ldl,
            qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t* p, qrupdate_fortran_complex_double_t* u,
            qrupdate_fortran_complex_double_t* v, qrupdate_fortran_complex_double_t* w);

    QRUPDATE_EXPORT void qrupdate_zqhqr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_double_t* c, qrupdate_fortran_complex_double_t* s);

    QRUPDATE_EXPORT void qrupdate_zqr1up(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_complex_double_t* u,
            qrupdate_fortran_complex_double_t* v, qrupdate_fortran_complex_double_t* w, qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zqrdec(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zqrder(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_complex_double_t* w,
            qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zqrinc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j,
            qrupdate_fortran_complex_double_t* x, qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zqrinr(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* q, qrupdate_fortran_int_t ldq,
            qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t j, qrupdate_fortran_complex_double_t* x,
            qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zqrot(char* dir, qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_double_t* c, qrupdate_fortran_complex_double_t* s);

    QRUPDATE_EXPORT void qrupdate_zqrqh(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr,
            qrupdate_fortran_double_t* c, qrupdate_fortran_complex_double_t* s);

    QRUPDATE_EXPORT void qrupdate_zqrshc(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t k, qrupdate_fortran_complex_double_t* q,
            qrupdate_fortran_int_t ldq, qrupdate_fortran_complex_double_t* r, qrupdate_fortran_int_t ldr, qrupdate_fortran_int_t i,
            qrupdate_fortran_int_t j, qrupdate_fortran_complex_double_t* w, qrupdate_fortran_double_t* rw);

    QRUPDATE_EXPORT void qrupdate_zqrtv1(qrupdate_fortran_int_t n, qrupdate_fortran_complex_double_t* u, qrupdate_fortran_double_t* w);


#ifdef __cplusplus
}
#endif

#endif
