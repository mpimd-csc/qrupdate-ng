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
#include <string.h>
#include "qrupdate.h"

#include "testutils.h"

static void stest(qrupdate_fortran_int_t n)
{
    float *A = (float *)malloc(n * n * sizeof(float));
    float *R = (float *)malloc(n * n * sizeof(float));
    float *u = (float *)malloc(n * sizeof(float));
    float *wrk = (float *)malloc(2 * n * sizeof(float));
    qrupdate_fortran_int_t one = 1, info;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(scopy,SCOPY)(&n, u, &one, wrk, &one);
    QRUPDATE_FORTRAN_GLOBAL(schgen,SCHGEN)(&n, A, &n, R, &n);
    QRUPDATE_FORTRAN_GLOBAL(sch1up,SCH1UP)(&n, R, &n, u, wrk + n);
    qrupdate_sch1dn(n, R, n, wrk, wrk + n, &info);
    QRUPDATE_FORTRAN_GLOBAL(schchk,SCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(wrk);
}

static void dtest(qrupdate_fortran_int_t n)
{
    double *A = (double *)malloc(n * n * sizeof(double));
    double *R = (double *)malloc(n * n * sizeof(double));
    double *u = (double *)malloc(n * sizeof(double));
    double *wrk = (double *)malloc(2 * n * sizeof(double));
    qrupdate_fortran_int_t one = 1, info;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(dcopy,DCOPY)(&n, u, &one, wrk, &one);
    QRUPDATE_FORTRAN_GLOBAL(dchgen,DCHGEN)(&n, A, &n, R, &n);
    QRUPDATE_FORTRAN_GLOBAL(dch1up,DCH1UP)(&n, R, &n, u, wrk + n);
    qrupdate_dch1dn(n, R, n, wrk, wrk + n, &info);
    QRUPDATE_FORTRAN_GLOBAL(dchchk,DCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(wrk);
}

static void ctest(qrupdate_fortran_int_t n)
{
    float _Complex *A = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float _Complex *u = (float _Complex *)malloc(n * sizeof(float _Complex));
    float _Complex *wrk = (float _Complex *)malloc(n * sizeof(float _Complex));
    float *rwrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t one = 1, info;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(ccopy,CCOPY)(&n, u, &one, wrk, &one);
    QRUPDATE_FORTRAN_GLOBAL(cchgen,CCHGEN)(&n, A, &n, R, &n);
    QRUPDATE_FORTRAN_GLOBAL(cch1up,CCH1UP)(&n, R, &n, u, rwrk);
    qrupdate_cch1dn(n, R, n, wrk, rwrk, &info);
    QRUPDATE_FORTRAN_GLOBAL(cchchk,CCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(wrk); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t n)
{
    double _Complex *A = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double _Complex *u = (double _Complex *)malloc(n * sizeof(double _Complex));
    double _Complex *wrk = (double _Complex *)malloc(n * sizeof(double _Complex));
    double *rwrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t one = 1, info;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &one, u, &n);
    QRUPDATE_FORTRAN_GLOBAL(zcopy,ZCOPY)(&n, u, &one, wrk, &one);
    QRUPDATE_FORTRAN_GLOBAL(zchgen,ZCHGEN)(&n, A, &n, R, &n);
    QRUPDATE_FORTRAN_GLOBAL(zch1up,ZCH1UP)(&n, R, &n, u, rwrk);
    qrupdate_zch1dn(n, R, n, wrk, rwrk, &info);
    QRUPDATE_FORTRAN_GLOBAL(zchchk,ZCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(u); free(wrk); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t n = 50;

    printf("\n");
    printf("testing Cholesky rank-1 downdate routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("sch1dn test:\n");
    stest(n);
    printf("dch1dn test:\n");
    dtest(n);
    printf("cch1dn test:\n");
    ctest(n);
    printf("zch1dn test:\n");
    ztest(n);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
