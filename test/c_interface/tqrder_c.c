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
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    float *A = (float *)malloc(m * maxmn * sizeof(float));
    float *Q = (float *)malloc(m * m * sizeof(float));
    float *R = (float *)malloc(m * n * sizeof(float));
    float *wrk = (float *)malloc(2 * m * sizeof(float));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(sqrgen,SQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < m; i++)
        QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(&n, &A[i], &m, &A[i - 1], &m);
    qrupdate_sqrder(m, n, Q, m, R, m, j, wrk);
    {
        qrupdate_fortran_int_t mm1 = m - 1, nm1 = m - 1;
        QRUPDATE_FORTRAN_GLOBAL(sqrchk,SQRCHK)(&mm1, &n, &nm1, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(wrk);
}

static void dtest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    double *A = (double *)malloc(m * maxmn * sizeof(double));
    double *Q = (double *)malloc(m * m * sizeof(double));
    double *R = (double *)malloc(m * n * sizeof(double));
    double *wrk = (double *)malloc(2 * m * sizeof(double));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(dqrgen,DQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < m; i++)
        QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(&n, &A[i], &m, &A[i - 1], &m);
    qrupdate_dqrder(m, n, Q, m, R, m, j, wrk);
    {
        qrupdate_fortran_int_t mm1 = m - 1, nm1 = m - 1;
        QRUPDATE_FORTRAN_GLOBAL(dqrchk,DQRCHK)(&mm1, &n, &nm1, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(wrk);
}

static void ctest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    float _Complex *A = (float _Complex *)malloc(m * maxmn * sizeof(float _Complex));
    float _Complex *Q = (float _Complex *)malloc(m * m * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(m * n * sizeof(float _Complex));
    float _Complex *wrk = (float _Complex *)malloc(m * sizeof(float _Complex));
    float *rwrk = (float *)malloc(m * sizeof(float));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(cqrgen,CQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < m; i++)
        QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(&n, &A[i], &m, &A[i - 1], &m);
    qrupdate_cqrder(m, n, Q, m, R, m, j, wrk, rwrk);
    {
        qrupdate_fortran_int_t mm1 = m - 1, nm1 = m - 1;
        QRUPDATE_FORTRAN_GLOBAL(cqrchk,CQRCHK)(&mm1, &n, &nm1, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(wrk); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    qrupdate_fortran_int_t maxmn = m > n ? m : n;
    double _Complex *A = (double _Complex *)malloc(m * maxmn * sizeof(double _Complex));
    double _Complex *Q = (double _Complex *)malloc(m * m * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(m * n * sizeof(double _Complex));
    double _Complex *wrk = (double _Complex *)malloc(m * sizeof(double _Complex));
    double *rwrk = (double *)malloc(m * sizeof(double));
    qrupdate_fortran_int_t i, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(zqrgen,ZQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    for (i = j; i < m; i++)
        QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(&n, &A[i], &m, &A[i - 1], &m);
    qrupdate_zqrder(m, n, Q, m, R, m, j, wrk, rwrk);
    {
        qrupdate_fortran_int_t mm1 = m - 1, nm1 = m - 1;
        QRUPDATE_FORTRAN_GLOBAL(zqrchk,ZQRCHK)(&mm1, &n, &nm1, A, &m, Q, &m, R, &m);
    }

    free(A); free(Q); free(R); free(wrk); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t m = 60, n = 40, j = 30;

    printf("\n");
    printf("testing QR row delete routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("sqrder test (full factorization):\n");
    stest(m, n, j);
    printf("dqrder test (full factorization):\n");
    dtest(m, n, j);
    printf("cqrder test (full factorization):\n");
    ctest(m, n, j);
    printf("zqrder test (full factorization):\n");
    ztest(m, n, j);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
