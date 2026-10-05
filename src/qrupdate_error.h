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

#ifndef QRUPDATE_ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "qrupdate_export.h"
#include "qrupdate_config.h"
#include "qrupdate_fortran_mangle.h"
#include "qrupdate_fortran_types.h"

    /**
     * @brief Type definition for C error handler function pointers.
     * \ingroup error
     *
     * This typedef defines the signature for C error handler functions
     * used with qrupdate_set_error. The handler is called when an error
     * occurs in qrupdate routines.
     *
     * @param srname The name of the routine that encountered the error.
     *               Passed as a null-terminated C string.
     * @param info   The error code.
     * @param aux    Optional pointer to auxiliary data (NULL if not set).
     */
    typedef void (*qrupdate_error_handler_t)(const char*, int, void*);

    /**
     * @copydoc qrupdate_error::qrupdate_set_error_c
     * \ingroup error
     * \sa qrupdate_error::qrupdate_set_error_c
     */
    QRUPDATE_EXPORT void qrupdate_set_error(qrupdate_error_handler_t p_handler);

    /**
     * @copydoc qrupdate_error::qrupdate_set_error_data_c
     * \ingroup error
     * \sa qrupdate_error::qrupdate_set_error_data_c
     * */
    QRUPDATE_EXPORT void qrupdate_set_error_data(void *p_aux);

    /**
     * @copydoc qrupdate_error::qrupdate_get_error_c
     * \ingroup error
     * \sa qrupdate_error::qrupdate_get_error_c
     */
    QRUPDATE_EXPORT qrupdate_error_handler_t qrupdate_get_error(void);
#ifdef __cplusplus
}
#endif

#endif
