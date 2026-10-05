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

static void stest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    qrupdate_fortran_int_t mp1 = m + 1;
    qrupdate_fortran_int_t maxmn = mp1 > n ? mp1 : n;
    float *A = (float *)malloc(mp1 * maxmn * sizeof(float));
    float *Q = (float *)malloc(mp1 * mp1 * sizeof(float));
    float *R = (float *)malloc(mp1 * n * sizeof(float));
    float *u = (float *)malloc(n * sizeof(float));
    float *wrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &n, A, &mp1);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(sqrgen,SQRGEN)(&m, &n, A, &mp1, Q, &mp1, R, &mp1);
    for (i = m; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(&n, &A[i - 1], &mp1, &A[i], &mp1);
    QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(&n, u, &one, &A[j - 1], &mp1);
    qrupdate_sqrinr(m, n, Q, mp1, R, mp1, j, u, wrk);
    QRUPDATE_FORTRAN_GLOBAL(sqrchk,SQRCHK)(&mp1, &n, &mp1, A, &mp1, Q, &mp1, R, &mp1);

    free(A); free(Q); free(R); free(u); free(wrk);
}

static void dtest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    qrupdate_fortran_int_t mp1 = m + 1;
    qrupdate_fortran_int_t maxmn = mp1 > n ? mp1 : n;
    double *A = (double *)malloc(mp1 * maxmn * sizeof(double));
    double *Q = (double *)malloc(mp1 * mp1 * sizeof(double));
    double *R = (double *)malloc(mp1 * n * sizeof(double));
    double *u = (double *)malloc(n * sizeof(double));
    double *wrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &n, A, &mp1);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(dqrgen,DQRGEN)(&m, &n, A, &mp1, Q, &mp1, R, &mp1);
    for (i = m; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(&n, &A[i - 1], &mp1, &A[i], &mp1);
    QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(&n, u, &one, &A[j - 1], &mp1);
    qrupdate_dqrinr(m, n, Q, mp1, R, mp1, j, u, wrk);
    QRUPDATE_FORTRAN_GLOBAL(dqrchk,DQRCHK)(&mp1, &n, &mp1, A, &mp1, Q, &mp1, R, &mp1);

    free(A); free(Q); free(R); free(u); free(wrk);
}

static void ctest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    qrupdate_fortran_int_t mp1 = m + 1;
    qrupdate_fortran_int_t maxmn = mp1 > n ? mp1 : n;
    float _Complex *A = (float _Complex *)malloc(mp1 * maxmn * sizeof(float _Complex));
    float _Complex *Q = (float _Complex *)malloc(mp1 * mp1 * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(mp1 * n * sizeof(float _Complex));
    float _Complex *u = (float _Complex *)malloc(n * sizeof(float _Complex));
    float *rwrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &n, A, &mp1);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(cqrgen,CQRGEN)(&m, &n, A, &mp1, Q, &mp1, R, &mp1);
    for (i = m; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(&n, &A[i - 1], &mp1, &A[i], &mp1);
    QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(&n, u, &one, &A[j - 1], &mp1);
    qrupdate_cqrinr(m, n, Q, mp1, R, mp1, j, u, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(cqrchk,CQRCHK)(&mp1, &n, &mp1, A, &mp1, Q, &mp1, R, &mp1);

    free(A); free(Q); free(R); free(u); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    qrupdate_fortran_int_t mp1 = m + 1;
    qrupdate_fortran_int_t maxmn = mp1 > n ? mp1 : n;
    double _Complex *A = (double _Complex *)malloc(mp1 * maxmn * sizeof(double _Complex));
    double _Complex *Q = (double _Complex *)malloc(mp1 * mp1 * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(mp1 * n * sizeof(double _Complex));
    double _Complex *u = (double _Complex *)malloc(n * sizeof(double _Complex));
    double *rwrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &n, A, &mp1);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(zqrgen,ZQRGEN)(&m, &n, A, &mp1, Q, &mp1, R, &mp1);
    for (i = m; i >= j; i--)
        QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(&n, &A[i - 1], &mp1, &A[i], &mp1);
    QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(&n, u, &one, &A[j - 1], &mp1);
    qrupdate_zqrinr(m, n, Q, mp1, R, mp1, j, u, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(zqrchk,ZQRCHK)(&mp1, &n, &mp1, A, &mp1, Q, &mp1, R, &mp1);

    free(A); free(Q); free(R); free(u); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t m = 60, n = 40, j = 30;

    printf("\n");
    printf("testing QR row insert routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("sqrinr test (full factorization):\n");
    stest(m, n, j);
    printf("dqrinr test (full factorization):\n");
    dtest(m, n, j);
    printf("cqrinr test (full factorization):\n");
    ctest(m, n, j);
    printf("zqrinr test (full factorization):\n");
    ztest(m, n, j);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
