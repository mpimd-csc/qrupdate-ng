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
#include <math.h>
#include "qrupdate.h"

#include "testutils.h"

static void stest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    float *A = (float *)malloc(m * m * sizeof(float));
    float *Q = (float *)malloc(m * m * sizeof(float));
    float *R = (float *)malloc(m * m * sizeof(float));
    float *u = (float *)malloc(m * sizeof(float));
    float tol, err, prod;
    qrupdate_fortran_int_t i, j;

    QRUPDATE_FORTRAN_GLOBAL(srandg,SRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(sqrgen,SQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    qrupdate_sgqvec(m, n, Q, m, u);
    tol = 500.0f * QRUPDATE_FORTRAN_GLOBAL(slamch,SLAMCH)("p");

    if (n == 0) {
        /* check u(1) = 1 */
        err = fabsf(u[0] - 1.0f);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   s (n=0, ldq=%d) u(1)=1  : %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        /* check u(2:m) = 0 */
        if (m > 1) {
            err = 0.0f;
            for (i = 1; i < m; i++) if (fabsf(u[i]) > err) err = fabsf(u[i]);
        } else {
            err = 0.0f;
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   s (n=0, ldq=%d) u(2:m)=0: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    } else {
        /* check Q'*u = 0 */
        err = 0.0f;
        for (i = 0; i < n; i++) {
            prod = 0.0f;
            for (j = 0; j < m; j++) prod += Q[j + i * m] * u[j];
            if (fabsf(prod) > err) err = fabsf(prod);
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   s (n>0, ldq=%d) max|Q'u|: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        /* check norm(u) = 1 */
        err = 0.0f;
        for (i = 0; i < m; i++) err += u[i] * u[i];
        err = fabsf(sqrtf(err) - 1.0f);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   s (n>0, ldq=%d) norm(u) : %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    }

    free(A); free(Q); free(R); free(u);
}

static void dtest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    double *A = (double *)malloc(m * m * sizeof(double));
    double *Q = (double *)malloc(m * m * sizeof(double));
    double *R = (double *)malloc(m * m * sizeof(double));
    double *u = (double *)malloc(m * sizeof(double));
    double tol, err, prod;
    qrupdate_fortran_int_t i, j;

    QRUPDATE_FORTRAN_GLOBAL(drandg,DRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(dqrgen,DQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    qrupdate_dgqvec(m, n, Q, m, u);
    tol = 500.0 * QRUPDATE_FORTRAN_GLOBAL(dlamch,DLAMCH)("p");

    if (n == 0) {
        err = fabs(u[0] - 1.0);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   d (n=0, ldq=%d) u(1)=1  : %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        if (m > 1) {
            err = 0.0;
            for (i = 1; i < m; i++) if (fabs(u[i]) > err) err = fabs(u[i]);
        } else {
            err = 0.0;
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   d (n=0, ldq=%d) u(2:m)=0: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    } else {
        err = 0.0;
        for (i = 0; i < n; i++) {
            prod = 0.0;
            for (j = 0; j < m; j++) prod += Q[j + i * m] * u[j];
            if (fabs(prod) > err) err = fabs(prod);
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   d (n>0, ldq=%d) max|Q'u|: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        err = 0.0;
        for (i = 0; i < m; i++) err += u[i] * u[i];
        err = fabs(sqrt(err) - 1.0);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   d (n>0, ldq=%d) norm(u) : %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    }

    free(A); free(Q); free(R); free(u);
}

static void ctest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    float _Complex *A = (float _Complex *)malloc(m * m * sizeof(float _Complex));
    float _Complex *Q = (float _Complex *)malloc(m * m * sizeof(float _Complex));
    float _Complex *R = (float _Complex *)malloc(m * m * sizeof(float _Complex));
    float _Complex *u = (float _Complex *)malloc(m * sizeof(float _Complex));
    float tol, err, rr;
    float _Complex rc;
    qrupdate_fortran_int_t i, j;

    QRUPDATE_FORTRAN_GLOBAL(crandg,CRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(cqrgen,CQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    qrupdate_cgqvec(m, n, Q, m, u);
    tol = 500.0f * QRUPDATE_FORTRAN_GLOBAL(slamch,SLAMCH)("p");

    if (n == 0) {
        err = fabsf(sqrtf(crealf(u[0]) * crealf(u[0]) + cimagf(u[0]) * cimagf(u[0])) - 1.0f);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   c (n=0, ldq=%d) |u(1)|=1: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        if (m > 1) {
            err = 0.0f;
            for (i = 1; i < m; i++) {
                rr = sqrtf(crealf(u[i]) * crealf(u[i]) + cimagf(u[i]) * cimagf(u[i]));
                if (rr > err) err = rr;
            }
        } else {
            err = 0.0f;
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   c (n=0, ldq=%d) u(2:m)=0: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    } else {
        err = 0.0f;
        for (i = 0; i < n; i++) {
            rc = 0.0f;
            for (j = 0; j < m; j++) rc += __builtin_conjf(Q[j + i * m]) * u[j];
            rr = sqrtf(crealf(rc) * crealf(rc) + cimagf(rc) * cimagf(rc));
            if (rr > err) err = rr;
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   c (n>0, ldq=%d) max|Q'u|: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        err = 0.0f;
        for (i = 0; i < m; i++) err += crealf(u[i]) * crealf(u[i]) + cimagf(u[i]) * cimagf(u[i]);
        err = fabsf(sqrtf(err) - 1.0f);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   c (n>0, ldq=%d) norm(u) : %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    }

    free(A); free(Q); free(R); free(u);
}

static void ztest(qrupdate_fortran_int_t m, qrupdate_fortran_int_t n)
{
    double _Complex *A = (double _Complex *)malloc(m * m * sizeof(double _Complex));
    double _Complex *Q = (double _Complex *)malloc(m * m * sizeof(double _Complex));
    double _Complex *R = (double _Complex *)malloc(m * m * sizeof(double _Complex));
    double _Complex *u = (double _Complex *)malloc(m * sizeof(double _Complex));
    double tol, err, rr;
    double _Complex rc;
    qrupdate_fortran_int_t i, j;

    QRUPDATE_FORTRAN_GLOBAL(zrandg,ZRANDG)(&m, &n, A, &m);
    QRUPDATE_FORTRAN_GLOBAL(zqrgen,ZQRGEN)(&m, &n, A, &m, Q, &m, R, &m);
    qrupdate_zgqvec(m, n, Q, m, u);
    tol = 500.0 * QRUPDATE_FORTRAN_GLOBAL(dlamch,DLAMCH)("p");

    if (n == 0) {
        err = fabs(sqrt(creal(u[0]) * creal(u[0]) + cimag(u[0]) * cimag(u[0])) - 1.0);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   z (n=0, ldq=%d) |u(1)|=1: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        if (m > 1) {
            err = 0.0;
            for (i = 1; i < m; i++) {
                rr = sqrt(creal(u[i]) * creal(u[i]) + cimag(u[i]) * cimag(u[i]));
                if (rr > err) err = rr;
            }
        } else {
            err = 0.0;
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   z (n=0, ldq=%d) u(2:m)=0: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    } else {
        err = 0.0;
        for (i = 0; i < n; i++) {
            rc = 0.0;
            for (j = 0; j < m; j++) rc += __builtin_conj(Q[j + i * m]) * u[j];
            rr = sqrt(creal(rc) * creal(rc) + cimag(rc) * cimag(rc));
            if (rr > err) err = rr;
        }
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   z (n>0, ldq=%d) max|Q'u|: %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
        err = 0.0;
        for (i = 0; i < m; i++) err += creal(u[i]) * creal(u[i]) + cimag(u[i]) * cimag(u[i]);
        err = fabs(sqrt(err) - 1.0);
        if (err < tol) stats_.passed++; else stats_.failed++;
        printf("   z (n>0, ldq=%d) norm(u) : %12.4E %s\n", (int)m, err, err < tol ? "PASS" : "FAIL");
    }

    free(A); free(Q); free(R); free(u);
}

int main(void)
{
    qrupdate_fortran_int_t m, n, i;

    printf("\n");
    printf("Testing gqvec routines.\n");
    printf("\n");
    printf("----------------------------------------------------------------------\n");

    stats_.passed = 0;
    stats_.failed = 0;

    for (i = 0; i < 6; i++) {
        switch (i) {
        case 0: m = 10; n = 0; break;
        case 1: m = 10; n = 1; break;
        case 2: m = 10; n = m - 1; break;
        case 3: m = 10; n = m / 2; break;
        case 4: m = 10; n = m / 2 - 1; break;
        case 5: m = 100; n = 50; break;
        default: m = 10; n = 0; break;
        }
        printf("Test case (m,n) = (%d, %d):\n", (int)m, (int)n);
        stest(m, n);
        dtest(m, n);
        ctest(m, n);
        ztest(m, n);
    }

    printf("----------------------------------------------------------------------\n");
    printf(" total: PASSED %6d FAILED %6d\n", stats_.passed, stats_.failed);
    printf("\n");

    if (stats_.failed != 0) return 1;
    return 0;
}
