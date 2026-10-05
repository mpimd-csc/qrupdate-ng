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
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    float *A = (float *)malloc(m * maxmn * sizeof(float));
    float *Q = (float *)malloc(m * m * sizeof(float));
    float *R = (float *)malloc(m * n * sizeof(float));
    float *wrk = (float *)malloc(m * sizeof(float));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(sqrgen,SQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < n; i++)
        QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(&m, &A[i * m], &one, &A[(i - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_sqrdec(m, n, k, Q, m, R, m, j, wrk);
    if (ec) k = n + 1;
    QRUPDATE_FORTRAN_GLOBAL(sqrchk,SQRCHK)(&(qrupdate_fortran_int_t){m}, &(qrupdate_fortran_int_t){n - 1}, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(wrk);
}

static void dtest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    double *A = (double *)malloc(m * maxmn * sizeof(double));
    double *Q = (double *)malloc(m * m * sizeof(double));
    double *R = (double *)malloc(m * n * sizeof(double));
    double *wrk = (double *)malloc(m * sizeof(double));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(dqrgen,DQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < n; i++)
        QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(&m, &A[i * m], &one, &A[(i - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_dqrdec(m, n, k, Q, m, R, m, j, wrk);
    if (ec) k = n + 1;
    QRUPDATE_FORTRAN_GLOBAL(dqrchk,DQRCHK)(&(qrupdate_fortran_int_t){m}, &(qrupdate_fortran_int_t){n - 1}, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(wrk);
}

static void ctest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    float _Complex *A = (float _Complex *)malloc(m * maxmn * sizeof(float _Complex));
    float _Complex *Q = (float _Complex *)malloc(m * m * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(m * n * sizeof(float _Complex));
    float *rwrk = (float *)malloc(m * sizeof(float));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(cqrgen,CQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < n; i++)
        QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(&m, &A[i * m], &one, &A[(i - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_cqrdec(m, n, k, Q, m, R, m, j, rwrk);
    if (ec) k = n + 1;
    QRUPDATE_FORTRAN_GLOBAL(cqrchk,CQRCHK)(&(qrupdate_fortran_int_t){m}, &(qrupdate_fortran_int_t){n - 1}, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j, qrupdate_fortran_int_t ec)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    double _Complex *A = (double _Complex *)malloc(m * maxmn * sizeof(double _Complex));
    double _Complex *Q = (double _Complex *)malloc(m * m * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(m * n * sizeof(double _Complex));
    double *rwrk = (double *)malloc(m * sizeof(double));
    qrupdate_fortran_int_t k, i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(zqrgen,ZQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < n; i++)
        QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(&m, &A[i * m], &one, &A[(i - 1) * m], &one);
    k = ec ? n : m;
    qrupdate_zqrdec(m, n, k, Q, m, R, m, j, rwrk);
    if (ec) k = n + 1;
    QRUPDATE_FORTRAN_GLOBAL(zqrchk,ZQRCHK)(&(qrupdate_fortran_int_t){m}, &(qrupdate_fortran_int_t){n - 1}, &k, A, &m, Q, &m, R, &m);

    free(A); free(Q); free(R); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t m = 60, n = 40, j = 12;

    printf("\n");
    printf("testing QR column delete routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("sqrdec test (full factorization):\n");
    stest(m, n, j, 0);
    printf("dqrdec test (full factorization):\n");
    dtest(m, n, j, 0);
    printf("cqrdec test (full factorization):\n");
    ctest(m, n, j, 0);
    printf("zqrdec test (full factorization):\n");
    ztest(m, n, j, 0);

    printf("sqrdec test (economized factorization):\n");
    stest(m, n, j, 1);
    printf("dqrdec test (economized factorization):\n");
    dtest(m, n, j, 1);
    printf("cqrdec test (economized factorization):\n");
    ctest(m, n, j, 1);
    printf("zqrdec test (economized factorization):\n");
    ztest(m, n, j, 1);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
