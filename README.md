qrupdate-ng: A Library for Fast Updating of QR and Cholesky Decompositions
==========================================================================

* License: GPL-3.0-or-later (opensource.org/licenses/GPL-3.0)
* Version: 1.1.3
* Authors: M. Koehler (0000-0003-2338-9904), C. Himpe (0000-0003-2194-6754)
* Summary: A FORTRAN library for rank-1 matrix decomposition updates

## Supported Matrix Decompositions

### [QR Decomposition](https://en.wikipedia.org/wiki/QR_decomposition)

* QR rank-1 update (`qr1up`)

    Updates the QR factorization after an additive rank-1 update to the original matrix (`A = A + u*v'`).
    Works for full & economized factorization.

* QR column insert (`qrinc`)

    Updates the QR factorization after an inserting a column to the original matrix.
    Works for full & economized factorization.

* QR column delete (`qrdec`)

    Updates the QR factorization after an deleting a column from the original matrix.
    Works for full & economized factorization.

* QR column shift (`qrshc`)

    Updates the QR factorization after a circular shift of columns in  the original matrix.
    Works for full & economized factorization.

* QR row insert (`qrinr`)

    Updates the QR factorization after an inserting a row to the original matrix.
    Works for full factorization only.

* QR row delete (`qrder`)

    Updates the QR factorization after an deleting a row from the original matrix.
    Works for full factorization only.

### [Cholesky Decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition)

* Cholesky rank-1 update (`ch1up`)

    Updates the Cholesky factorization after positive rank-1 update `A = A + u*u'`.

* Cholesky rank-1 downdate (`ch1dn`)

    Updates the Cholesky factorization after positive rank-1 downdate `A = A - u*u'`.

* Cholesky symmetric insert (`chinx`)

    Updates the Cholesky factorization after a symmetric column/row insertion.

* Cholesky symmetric insert (`chdex`)

    Updates the Cholesky factorization after a symmetric column/row deletion.

* Cholesky symmetric shift (`chshx`)

    Updates the Cholesky factorization after a symmetric column/row left/right circular shift.

### [LU Decomposition](https://en.wikipedia.org/wiki/LU_decomposition)

* LU rank-1 update (`lu1up`)

    Updates the LU factorization after a rank-1 update (`A = A + u*v.'`).
    No pivoting available. Faster than `lup1up`, but less stable.

* LU pivoted rank-1 update (`lup1up`)

    Updates a row-pivoted LUP factorization after a rank-1 update (`A = A + u*v.'`).
    Also updates the row permutation matrix. Slower than `lu1up`, but more stable.

### Based on

qrupdate-ng is a modernized and maintained fork of:

[qrupdate](https://sourceforge.net/p/qrupdate) (Version 1.1.2)

For a list of original authors see [CHANGELOG](CHANGELOG) file.

### License

See the [LICENSE](LICENSE) file.

### Installation

See the [INSTALL](INSTALL.md) file for installation details.
Consult individual routines for documentation.
