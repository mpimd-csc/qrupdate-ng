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

static void stest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn1 = (m > n + 1) ? m : n + 1;
    float *A = (float *)malloc(m * maxmn1 * sizeof(float));
    float *Q = (float *)malloc(m * m * sizeof(float));
    float *R = (float *)malloc(m * (n + 1) * sizeof(float));
    float *u = (float *)malloc(m * sizeof(float));
    float *wrk = (float *)malloc(m * sizeof(float));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(sqrgen,SQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = n; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(&m, &A[(i - 1) * m], &one, &A[i * m], &one);
    QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(&m, u, &one, &A[(j - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_sqrinc(m, n, k, Q, m, R, m, j, u, wrk);
    if (ec) k = n + 1;
    {
        qrupdate_fortran_int_t np1 = n + 1;
        QRUPDATE_FORTRAN_GLOBAL(sqrchk,SQRCHK)(&m, &np1, &k, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(u); free(wrk);
}

static void dtest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn1 = (m > n + 1) ? m : n + 1;
    double *A = (double *)malloc(m * maxmn1 * sizeof(double));
    double *Q = (double *)malloc(m * m * sizeof(double));
    double *R = (double *)malloc(m * (n + 1) * sizeof(double));
    double *u = (double *)malloc(m * sizeof(double));
    double *wrk = (double *)malloc(m * sizeof(double));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(dqrgen,DQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = n; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(&m, &A[(i - 1) * m], &one, &A[i * m], &one);
    QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(&m, u, &one, &A[(j - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_dqrinc(m, n, k, Q, m, R, m, j, u, wrk);
    if (ec) k = n + 1;
    {
        qrupdate_fortran_int_t np1 = n + 1;
        QRUPDATE_FORTRAN_GLOBAL(dqrchk,DQRCHK)(&m, &np1, &k, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(u); free(wrk);
}

static void ctest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn1 = (m > n + 1) ? m : n + 1;
    float _Complex *A = (float _Complex *)malloc(m * maxmn1 * sizeof(float _Complex));
    float _Complex *Q = (float _Complex *)malloc(m * m * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(m * (n + 1) * sizeof(float _Complex));
    float _Complex *u = (float _Complex *)malloc(m * sizeof(float _Complex));
    float *rwrk = (float *)malloc(m * sizeof(float));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(cqrgen,CQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = n; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(&m, &A[(i - 1) * m], &one, &A[i * m], &one);
    QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(&m, u, &one, &A[(j - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_cqrinc(m, n, k, Q, m, R, m, j, u, rwrk);
    if (ec) k = n + 1;
    {
        qrupdate_fortran_int_t np1 = n + 1;
        QRUPDATE_FORTRAN_GLOBAL(cqrchk,CQRCHK)(&m, &np1, &k, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(u); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn1 = (m > n + 1) ? m : n + 1;
    double _Complex *A = (double _Complex *)malloc(m * maxmn1 * sizeof(double _Complex));
    double _Complex *Q = (double _Complex *)malloc(m * m * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(m * (n + 1) * sizeof(double _Complex));
    double _Complex *u = (double _Complex *)malloc(m * sizeof(double _Complex));
    double *rwrk = (double *)malloc(m * sizeof(double));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &one, u, &m);
    QRUPDATE_FORTRAN_GLOBAL(zqrgen,ZQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = n; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(&m, &A[(i - 1) * m], &one, &A[i * m], &one);
    QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(&m, u, &one, &A[(j - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_zqrinc(m, n, k, Q, m, R, m, j, u, rwrk);
    if (ec) k = n + 1;
    {
        qrupdate_fortran_int_t np1 = n + 1;
        QRUPDATE_FORTRAN_GLOBAL(zqrchk,ZQRCHK)(&m, &np1, &k, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(u); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t m = 60, n = 40, j = 28;

    printf("\n");
    printf("testing QR column insert routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("sqrinc test (full factorization):\n");
    stest(m, n, j, 0);
    printf("dqrinc test (full factorization):\n");
    dtest(m, n, j, 0);
    printf("cqrinc test (full factorization):\n");
    ctest(m, n, j, 0);
    printf("zqrinc test (full factorization):\n");
    ztest(m, n, j, 0);

    printf("sqrinc test (economized factorization):\n");
    stest(m, n, j, 1);
    printf("dqrinc test (economized factorization):\n");
    dtest(m, n, j, 1);
    printf("cqrinc test (economized factorization):\n");
    ctest(m, n, j, 1);
    printf("zqrinc test (economized factorization):\n");
    ztest(m, n, j, 1);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
