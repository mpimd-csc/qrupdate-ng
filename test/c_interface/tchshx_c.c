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

static void swap_float(float *a, float *b) { float t = *a; *a = *b; *b = t; }
static void swap_double(double *a, double *b) { double t = *a; *a = *b; *b = t; }
static void swap_cfloat(float _Complex *a, float _Complex *b) { float _Complex t = *a; *a = *b; *b = t; }
static void swap_cdouble(double _Complex *a, double _Complex *b) { double _Complex t = *a; *a = *b; *b = t; }

static void stest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t si, qrupdate_fortran_int_t sj)
{
    float *A = (float *)malloc(n * n * sizeof(float));
    float *R = (float *)malloc(n * n * sizeof(float));
    float *wrk = (float *)malloc(2 * n * sizeof(float));
    qrupdate_fortran_int_t k;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(schgen,SCHGEN)(&n, A, &n, R, &n);
    if (si < sj) {
        for (k = si; k < sj; k++) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_float(&A[i + (k - 1) * n], &A[i + k * n]);
            for (i = 0; i < n; i++) swap_float(&A[(k - 1) + i * n], &A[k + i * n]);
        }
    } else if (si > sj) {
        for (k = si; k > sj; k--) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_float(&A[i + (k - 1) * n], &A[i + (k - 2) * n]);
            for (i = 0; i < n; i++) swap_float(&A[(k - 1) + i * n], &A[(k - 2) + i * n]);
        }
    }
    qrupdate_schshx(n, R, n, si, sj, wrk);
    QRUPDATE_FORTRAN_GLOBAL(schchk,SCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(wrk);
}

static void dtest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t si, qrupdate_fortran_int_t sj)
{
    double *A = (double *)malloc(n * n * sizeof(double));
    double *R = (double *)malloc(n * n * sizeof(double));
    double *wrk = (double *)malloc(2 * n * sizeof(double));
    qrupdate_fortran_int_t k;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(dchgen,DCHGEN)(&n, A, &n, R, &n);
    if (si < sj) {
        for (k = si; k < sj; k++) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_double(&A[i + (k - 1) * n], &A[i + k * n]);
            for (i = 0; i < n; i++) swap_double(&A[(k - 1) + i * n], &A[k + i * n]);
        }
    } else if (si > sj) {
        for (k = si; k > sj; k--) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_double(&A[i + (k - 1) * n], &A[i + (k - 2) * n]);
            for (i = 0; i < n; i++) swap_double(&A[(k - 1) + i * n], &A[(k - 2) + i * n]);
        }
    }
    qrupdate_dchshx(n, R, n, si, sj, wrk);
    QRUPDATE_FORTRAN_GLOBAL(dchchk,DCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(wrk);
}

static void ctest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t si, qrupdate_fortran_int_t sj)
{
    float _Complex *A = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(n * n * sizeof(float _Complex));
    float _Complex *wrk = (float _Complex *)malloc(n * sizeof(float _Complex));
    float *rwrk = (float *)malloc(n * sizeof(float));
    qrupdate_fortran_int_t k;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(cchgen,CCHGEN)(&n, A, &n, R, &n);
    if (si < sj) {
        for (k = si; k < sj; k++) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_cfloat(&A[i + (k - 1) * n], &A[i + k * n]);
            for (i = 0; i < n; i++) swap_cfloat(&A[(k - 1) + i * n], &A[k + i * n]);
        }
    } else if (si > sj) {
        for (k = si; k > sj; k--) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_cfloat(&A[i + (k - 1) * n], &A[i + (k - 2) * n]);
            for (i = 0; i < n; i++) swap_cfloat(&A[(k - 1) + i * n], &A[(k - 2) + i * n]);
        }
    }
    qrupdate_cchshx(n, R, n, si, sj, wrk, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(cchchk,CCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(wrk); free(rwrk);
}

static void ztest(qrupdate_fortran_int_t n, qrupdate_fortran_int_t si, qrupdate_fortran_int_t sj)
{
    double _Complex *A = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(n * n * sizeof(double _Complex));
    double _Complex *wrk = (double _Complex *)malloc(n * sizeof(double _Complex));
    double *rwrk = (double *)malloc(n * sizeof(double));
    qrupdate_fortran_int_t k;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&n, &n, A, &n);
    QRUPDATE_FORTRAN_GLOBAL(zchgen,ZCHGEN)(&n, A, &n, R, &n);
    if (si < sj) {
        for (k = si; k < sj; k++) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_cdouble(&A[i + (k - 1) * n], &A[i + k * n]);
            for (i = 0; i < n; i++) swap_cdouble(&A[(k - 1) + i * n], &A[k + i * n]);
        }
    } else if (si > sj) {
        for (k = si; k > sj; k--) {
            qrupdate_fortran_int_t i;
            for (i = 0; i < n; i++) swap_cdouble(&A[i + (k - 1) * n], &A[i + (k - 2) * n]);
            for (i = 0; i < n; i++) swap_cdouble(&A[(k - 1) + i * n], &A[(k - 2) + i * n]);
        }
    }
    qrupdate_zchshx(n, R, n, si, sj, wrk, rwrk);
    QRUPDATE_FORTRAN_GLOBAL(zchchk,ZCHCHK)(&n, A, &n, R, &n);

    free(A); free(R); free(wrk); free(rwrk);
}

int main(void)
{
    qrupdate_fortran_int_t n = 50;

    printf("\n");
    printf("testing QR column shift routines.\n");
    printf("All residual errors are expected to be small.\n");
    printf("\n");

    printf("schshx test (left shift):\n");
    stest(n, 20, 40);
    printf("dchshx test (left shift):\n");
    dtest(n, 20, 40);
    printf("cchshx test (left shift):\n");
    ctest(n, 20, 40);
    printf("zchshx test (left shift):\n");
    ztest(n, 20, 40);

    printf("schshx test (right shift):\n");
    stest(n, 40, 20);
    printf("dchshx test (right shift):\n");
    dtest(n, 40, 20);
    printf("cchshx test (right shift):\n");
    ctest(n, 40, 20);
    printf("zchshx test (right shift):\n");
    ztest(n, 40, 20);

    QRUPDATE_FORTRAN_GLOBAL(pstats,PSTATS)();
    return 0;
}
