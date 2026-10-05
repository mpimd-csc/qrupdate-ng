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

static void stest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    float *A = (float *)malloc(m * maxmn * sizeof(float));
    float *Q = (float *)malloc(m * m * sizeof(float));
    float *R = (float *)malloc(m * n * sizeof(float));
    float *u = (float *)malloc(m * sizeof(float));
    float *v = (float *)malloc(n * sizeof(float));
    float *wrk = (float *)malloc(2 * m * sizeof(float));
    qrupdate_fortran_int_t k, i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(sqrgen,SQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    /* sger(m,n,1e0,u,1,v,1,A,m) */
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0f * u[i] * v[j];
    k = ec ? n : m;
    qrupdate_sqr1up(m, n, k, Q, m, R, m, u, v, wrk);
    QRUPDATE_FORTRAN_GLOBAL(sqrchk,SQRCHK)(&m, &n, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(u); free(v); free(wrk);
}

static void dtest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    double *A = (double *)malloc(m * maxmn * sizeof(double));
    double *Q = (double *)malloc(m * m * sizeof(double));
    double *R = (double *)malloc(m * n * sizeof(double));
    double *u = (double *)malloc(m * sizeof(double));
    double *v = (double *)malloc(n * sizeof(double));
    double *wrk = (double *)malloc(2 * m * sizeof(double));
    qrupdate_fortran_int_t k, i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(dqrgen,DQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0 * u[i] * v[j];
    k = ec ? n : m;
    qrupdate_dqr1up(m, n, k, Q, m, R, m, u, v, wrk);
    QRUPDATE_FORTRAN_GLOBAL(dqrchk,DQRCHK)(&m, &n, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(u); free(v); free(wrk);
}

static void ctest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    float _Complex *A = (float _Complex *)malloc(m * maxmn * sizeof(float _Complex));
    float _Complex *Q = (float _Complex *)malloc(m * m * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(m * n * sizeof(float _Complex));
    float _Complex *u = (float _Complex *)malloc(m * sizeof(float _Complex));
    float _Complex *v = (float _Complex *)malloc(n * sizeof(float _Complex));
    float _Complex *wrk = (float _Complex *)malloc(m * sizeof(float _Complex));
    float *rwrk = (float *)malloc(m * sizeof(float));
    qrupdate_fortran_int_t k, i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(cqrgen,CQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    /* cgerc(m,n,(1e0,0e0),u,1,v,1,A,m) */
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0f * u[i] * __builtin_conjf(v[j]);
    k = ec ? n : m;
    qrupdate_cqr1up(m, n, k, Q, m, R, m, u, v, wrk, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(cqrchk,CQRCHK)(&m, &n, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(u); free(v); free(wrk); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    double _Complex *A = (double _Complex *)malloc(m * maxmn * sizeof(double _Complex));
    double _Complex *Q = (double _Complex *)malloc(m * m * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(m * n * sizeof(double _Complex));
    double _Complex *u = (double _Complex *)malloc(m * sizeof(double _Complex));
    double _Complex *v = (double _Complex *)malloc(n * sizeof(double _Complex));
    double _Complex *wrk = (double _Complex *)malloc(m * sizeof(double _Complex));
    double *rwrk = (double *)malloc(m * sizeof(double));
    qrupdate_fortran_int_t k, i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &one, v, &n);
    QRUPDATE_FORTRAN_GLOBAL(zqrgen,ZQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            A[i + j * m] += 1.0 * u[i] * __builtin_conj(v[j]);
    k = ec ? n : m;
    qrupdate_zqr1up(m, n, k, Q, m, R, m, u, v, wrk, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(zqrchk,ZQRCHK)(&m, &n, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(u); free(v); free(wrk); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t m, n, one = 1;

    printf("\n");
    printf("testing QR rank-1 update routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    m = 60; n = 40;
    printf("sqr1up test (full factorization):\n");
    stest(m, n, 0);
    printf("dqr1up test (full factorization):\n");
    dtest(m, n, 0);
    printf("cqr1up test (full factorization):\n");
    ctest(m, n, 0);
    printf("zqr1up test (full factorization):\n");
    ztest(m, n, 0);

    printf("sqr1up test (economized factorization):\n");
    stest(m, n, 1);
    printf("dqr1up test (economized factorization):\n");
    dtest(m, n, 1);
    printf("cqr1up test (economized factorization):\n");
    ctest(m, n, 1);
    printf("zqr1up test (economized factorization):\n");
    ztest(m, n, 1);

    m = 40; n = 60;
    printf("sqr1up test (rows < columns):\n");
    stest(m, n, 0);
    printf("dqr1up test (rows < columns):\n");
    dtest(m, n, 0);
    printf("cqr1up test (rows < columns):\n");
    ctest(m, n, 0);
    printf("zqr1up test (rows < columns):\n");
    ztest(m, n, 0);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
