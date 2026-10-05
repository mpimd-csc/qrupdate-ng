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

/* Common blocks shared with Fortran utils */
#include "testutils.h"

static void stest(qrupdate_fortran_int_t n)
{
    float *A = (float *)malloc(n * n * sizeof(float));
    float *R = (float *)malloc(n * n * sizeof(float));
    float *u = (float *)malloc(n * sizeof(float));
    float *wrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(schgen,SCHGEN)(&n, A, &n, R, &n);
    /* ssyr('U',n,1e0,u,1,A,n) */
    for (j = 0; j < n; j++)
        for (i = 0; i <= j; i++)
            A[i + j * n] += 1.0f * u[i] * u[j];
    qrupdate_sch1up(n, R, n, u, wrk);
    QRUPDATE_FORTRAN_GLOBAL(schchk,SCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(wrk);
}

static void dtest(qrupdate_fortran_int_t n)
{
    double *A = (double *)malloc(n * n * sizeof(double));
    double *R = (double *)malloc(n * n * sizeof(double));
    double *u = (double *)malloc(n * sizeof(double));
    double *wrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(dchgen,DCHGEN)(&n, A, &n, R, &n);
    /* dsyr('U',n,1d0,u,1,A,n) */
    for (j = 0; j < n; j++)
        for (i = 0; i <= j; i++)
            A[i + j * n] += 1.0 * u[i] * u[j];
    qrupdate_dch1up(n, R, n, u, wrk);
    QRUPDATE_FORTRAN_GLOBAL(dchchk,DCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(wrk);
}

static void ctest(qrupdate_fortran_int_t n)
{
    float _Complex *A = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float _Complex *u = (float _Complex *)malloc(n * sizeof(float _Complex));
    float *rwrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(cchgen,CCHGEN)(&n, A, &n, R, &n);
    /* cher('U',n,1e0,u,1,A,n) */
    for (j = 0; j < n; j++)
        for (i = 0; i <= j; i++)
            A[i + j * n] += 1.0f * u[i] * __builtin_conjf(u[j]);
    qrupdate_cch1up(n, R, n, u, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(cchchk,CCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t n)
{
    double _Complex *A = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double _Complex *u = (double _Complex *)malloc(n * sizeof(double _Complex));
    double *rwrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t i, j, one = 1;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(zchgen,ZCHGEN)(&n, A, &n, R, &n);
    /* zher('U',n,1d0,u,1,A,n) */
    for (j = 0; j < n; j++)
        for (i = 0; i <= j; i++)
            A[i + j * n] += 1.0 * u[i] * __builtin_conj(u[j]);
    qrupdate_zch1up(n, R, n, u, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(zchchk,ZCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t n = 50;

    printf("\n");
    printf("testing Cholesky rank-1 update routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("sch1up test:\n");
    stest(n);
    printf("dch1up test:\n");
    dtest(n);
    printf("cch1up test:\n");
    ctest(n);
    printf("zch1up test:\n");
    ztest(n);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
