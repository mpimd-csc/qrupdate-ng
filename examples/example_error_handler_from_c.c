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
 *  (at your option) any later version.
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
#include "qrupdate_error.h"

typedef struct {
    int aux1, aux2;
} error_data_t;

void c_error_handler(const char * srname, int info, void * aux)
{
    error_data_t * ed = (error_data_t *) aux;
    printf("I got an error from %s. (CODE = %d ), The address of the aux-parameter is 0x%lX and (aux1 = %d, aux2 = %d)\n",
            srname, (int) info, (unsigned long) aux, ed->aux1, ed->aux2);
    return;
}

int main(int argc, char *argv[])
{
    qrupdate_fortran_int_t info = 0;
    qrupdate_error_handler_t old_handler;
    error_data_t ed;

    ed.aux1 = 4711;
    ed.aux2 = 1337;

    old_handler = qrupdate_get_error();
    printf("Initial error handler %p\n", old_handler);

    qrupdate_set_error((void*) c_error_handler);
    qrupdate_set_error_data(&ed);
    qrupdate_cch1dn(-1, NULL, 1, NULL, NULL, &info);

    old_handler = qrupdate_get_error();
    printf("error handler before reset: %p\n", (void*) old_handler);


    printf("reset the error handler.\n");

    qrupdate_set_error(NULL);

    qrupdate_cch1dn(-1, NULL, 1, NULL, NULL, &info);

    printf("restore the handler from var.\n");
    qrupdate_set_error(old_handler);
    qrupdate_cch1dn(-1, NULL, 1, NULL, NULL, &info);

    return 0;
}

