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

static void stest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    float *A = (float *)malloc(n * n * sizeof(float));
    float *R = (float *)malloc(n * n * sizeof(float));
    float *wrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t i, k;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(schgen,SCHGEN)(&n, A, &n, R, &n);
    for (k = j; k < n; k++) {
        for (i = 0; i < j - 1; i++)
            A[i + (k - 1) * n] = A[i + k * n];
        for (i = j - 1; i < k; i++)
            A[i + (k - 1) * n] = A[(i + 1) + k * n];
    }
    qrupdate_schdex(n, R, n, j, wrk);
    {
        qrupdate_fortran_int_t nm1 = n - 1;
        QRUPDATE_FORTRAN_GLOBAL(schchk,SCHCHK)(&nm1, A, &n, R, &n);
    }

    free(A); free(R); free(wrk);
}

static void dtest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    double *A = (double *)malloc(n * n * sizeof(double));
    double *R = (double *)malloc(n * n * sizeof(double));
    double *wrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t i, k;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(dchgen,DCHGEN)(&n, A, &n, R, &n);
    for (k = j; k < n; k++) {
        for (i = 0; i < j - 1; i++)
            A[i + (k - 1) * n] = A[i + k * n];
        for (i = j - 1; i < k; i++)
            A[i + (k - 1) * n] = A[(i + 1) + k * n];
    }
    qrupdate_dchdex(n, R, n, j, wrk);
    {
        qrupdate_fortran_int_t nm1 = n - 1;
        QRUPDATE_FORTRAN_GLOBAL(dchchk,DCHCHK)(&nm1, A, &n, R, &n);
    }

    free(A); free(R); free(wrk);
}

static void ctest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    float _Complex *A = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float *rwrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t i, k;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(cchgen,CCHGEN)(&n, A, &n, R, &n);
    for (k = j; k < n; k++) {
        for (i = 0; i < j - 1; i++)
            A[i + (k - 1) * n] = A[i + k * n];
        for (i = j - 1; i < k; i++)
            A[i + (k - 1) * n] = A[(i + 1) + k * n];
    }
    qrupdate_cchdex(n, R, n, j, rwrk);
    {
        qrupdate_fortran_int_t nm1 = n - 1;
        QRUPDATE_FORTRAN_GLOBAL(cchchk,CCHCHK)(&nm1, A, &n, R, &n);
    }

    free(A); free(R); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t j)
{
    double _Complex *A = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double *rwrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t i, k;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(zchgen,ZCHGEN)(&n, A, &n, R, &n);
    for (k = j; k < n; k++) {
        for (i = 0; i < j - 1; i++)
            A[i + (k - 1) * n] = A[i + k * n];
        for (i = j - 1; i < k; i++)
            A[i + (k - 1) * n] = A[(i + 1) + k * n];
    }
    qrupdate_zchdex(n, R, n, j, rwrk);
    {
        qrupdate_fortran_int_t nm1 = n - 1;
        QRUPDATE_FORTRAN_GLOBAL(zchchk,ZCHCHK)(&nm1, A, &n, R, &n);
    }

    free(A); free(R); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t n = 50, j = 15;

    printf("\n");
    printf("testing Cholesky symmetric delete routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("schdex test:\n");
    stest(n, j);
    printf("dchdex test:\n");
    dtest(n, j);
    printf("cchdex test:\n");
    ctest(n, j);
    printf("zchdex test:\n");
    ztest(n, j);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
