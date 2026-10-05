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
#include <stdio.h>
#include <stdlib.h>
#include "qrupdate.h"

#include "testutils.h"

static void stest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    qrupdate_fortran_int_t k = m < n ? m : n;
    float *A = (float *)malloc(m * n * sizeof(float));
    float *L = (float *)malloc(m * k * sizeof(float));
    float *R = (float *)malloc(k * n * sizeof(float));
    float *u = (float *)malloc(m * sizeof(float));
    float *v = (float *)malloc(n * sizeof(float));
    float *w = (float *)malloc(m * sizeof(float));
    qrupdate_fortran_int_t *p = (qrupdate_fortran_int_t *)malloc(m * sizeof(qrupdate_fortran_int_t));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(slupgen,SLUPGEN)(&m, &n, A, &m, L, &m, R, &k, p);
    /* sger(m,n,1e0,u,1,v,1,A,m) */
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0f * u[i] * v[j];
    qrupdate_slup1up(m, n, L, m, R, k, p, u, v, w);
    QRUPDATE_FORTRAN_GLOBAL(slupchk,SLUPCHK)(&m, &n, A, &m, L, &m, R, &k, p);

    free(A); free(L); free(R); free(u); free(v); free(w); free(p);
}

static void dtest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    qrupdate_fortran_int_t k = m < n ? m : n;
    double *A = (double *)malloc(m * n * sizeof(double));
    double *L = (double *)malloc(m * k * sizeof(double));
    double *R = (double *)malloc(k * n * sizeof(double));
    double *u = (double *)malloc(m * sizeof(double));
    double *v = (double *)malloc(n * sizeof(double));
    double *w = (double *)malloc(m * sizeof(double));
    qrupdate_fortran_int_t *p = (qrupdate_fortran_int_t *)malloc(m * sizeof(qrupdate_fortran_int_t));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(dlupgen,DLUPGEN)(&m, &n, A, &m, L, &m, R, &k, p);
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0 * u[i] * v[j];
    qrupdate_dlup1up(m, n, L, m, R, k, p, u, v, w);
    QRUPDATE_FORTRAN_GLOBAL(dlupchk,DLUPCHK)(&m, &n, A, &m, L, &m, R, &k, p);

    free(A); free(L); free(R); free(u); free(v); free(w); free(p);
}

static void ctest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    qrupdate_fortran_int_t k = m < n ? m : n;
    float _Complex *A = (float _Complex *)malloc(m * n * sizeof(float _Complex));
    float _Complex *L = (float _Complex *)malloc(m * k * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(k * n * sizeof(float _Complex));
    float _Complex *u = (float _Complex *)malloc(m * sizeof(float _Complex));
    float _Complex *v = (float _Complex *)malloc(n * sizeof(float _Complex));
    float _Complex *w = (float _Complex *)malloc(m * sizeof(float _Complex));
    qrupdate_fortran_int_t *p = (qrupdate_fortran_int_t *)malloc(m * sizeof(qrupdate_fortran_int_t));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(clupgen,CLUPGEN)(&m, &n, A, &m, L, &m, R, &k, p);
    /* cgeru(m,n,(1e0,0e0),u,1,v,1,A,m) */
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0f * u[i] * v[j];
    qrupdate_clup1up(m, n, L, m, R, k, p, u, v, w);
    QRUPDATE_FORTRAN_GLOBAL(clupchk,CLUPCHK)(&m, &n, A, &m, L, &m, R, &k, p);

    free(A); free(L); free(R); free(u); free(v); free(w); free(p);
}

static void ztest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    qrupdate_fortran_int_t k = m < n ? m : n;
    double _Complex *A = (double _Complex *)malloc(m * n * sizeof(double _Complex));
    double _Complex *L = (double _Complex *)malloc(m * k * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(k * n * sizeof(double _Complex));
    double _Complex *u = (double _Complex *)malloc(m * sizeof(double _Complex));
    double _Complex *v = (double _Complex *)malloc(n * sizeof(double _Complex));
    double _Complex *w = (double _Complex *)malloc(m * sizeof(double _Complex));
    qrupdate_fortran_int_t *p = (qrupdate_fortran_int_t *)malloc(m * sizeof(qrupdate_fortran_int_t));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(zlupgen,ZLUPGEN)(&m, &n, A, &m, L, &m, R, &k, p);
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0 * u[i] * v[j];
    qrupdate_zlup1up(m, n, L, m, R, k, p, u, v, w);
    QRUPDATE_FORTRAN_GLOBAL(zlupchk,ZLUPCHK)(&m, &n, A, &m, L, &m, R, &k, p);

    free(A); free(L); free(R); free(u); free(v); free(w); free(p);
}

int main(void)
{
    qrupdate_fortran_int_t m, n, one = 1;

    printf("\n");
    printf("testing pivoted LU rank-1 update routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    m = 60; n = 40;
    printf("slup1up test (rows > columns):\n");
    stest(m, n);
    printf("dlup1up test (rows > columns):\n");
    dtest(m, n);
    printf("clup1up test (rows > columns):\n");
    ctest(m, n);
    printf("zlup1up test (rows > columns):\n");
    ztest(m, n);

    m = 40; n = 60;
    printf("slup1up test (rows < columns):\n");
    stest(m, n);
    printf("dlup1up test (rows < columns):\n");
    dtest(m, n);
    printf("clup1up test (rows < columns):\n");
    ctest(m, n);
    printf("zlup1up test (rows < columns):\n");
    ztest(m, n);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
