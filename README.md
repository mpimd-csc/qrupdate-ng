qrupdate-ng
===========

A Library for Fast Updating of QR, Cholesky, and LU  Decompositions

Version: 1.1.90 (2026-05-28)

* Authors: M. Koehler (0000-0003-2338-9904), C. Himpe (0000-0003-2194-6754)
* Original Author: Jaroslav Hajek (VZLU Prague)

## Supported Matrix Decompositions

The library supports the follwoing matrix decompositions

### QR Decomposition

See [QR Decompositions](\ref qrdecomp) for details.

This following operations are supported:

* QR rank-1 update (`qr1up`)

    Updates the QR factorization after an additive rank-1 update to the original matrix.
    Works for full and economized factorization.

* QR column insert (`qrinc`)

    Updates the QR factorization after an inserting a column to the original matrix.
    Works for full and economized factorization.

* QR column delete (`qrdec`)

    Updates the QR factorization after an deleting a column from the original matrix.
    Works for full and economized factorization.

* QR column shift (`qrshc`)

    Updates the QR factorization after a circular shift of columns in  the original matrix.
    Works for full and economized factorization.

* QR row insert (`qrinr`)

    Updates the QR factorization after an inserting a row to the original matrix.
    Works for full factorization only.

* QR row delete (`qrder`)

    Updates the QR factorization after an deleting a row from the original matrix.
    Works for full factorization only.

### Cholesky Decomposition

See [Cholesky Decompositions](\ref choldecomp) for details.

This following operations are supported:

* Cholesky rank-1 update (`ch1up`)

    Updates the Cholesky factorization after positive rank-1 update.

* Cholesky rank-1 downdate (`ch1dn`)

    Updates the Cholesky factorization after positive rank-1 downdate.

* Cholesky symmetric insert (`chinx`)

    Updates the Cholesky factorization after a symmetric column/row insertion.

* Cholesky symmetric insert (`chdex`)

    Updates the Cholesky factorization after a symmetric column/row deletion.

* Cholesky symmetric shift (`chshx`)

    Updates the Cholesky factorization after a symmetric column/row left/right circular shift.

### LU Decomposition

See [LU Decompositions](\ref ludecomp) for details.

This following operations are supported:

* LU rank-1 update (`lu1up`)

    Updates the LU factorization after a rank-1 update.
    No pivoting available. Faster than `lup1up`, but less stable.

* LU pivoted rank-1 update (`lup1up`)

    Updates a row-pivoted LUP factorization after a rank-1 update.
    Also updates the row permutation matrix. Slower than `lu1up`, but more stable.

## Installation

See the [INSTALL](INSTALL.md) file for installation details.
Consult individual routines for documentation.

### Fortran Package Manager

If you are using the Fortran Package Manager [FPM](https://fpm.fortran-lang.org/)
you can add `qrupdate-ng` as dependency by including the following in your
`fpm.toml`:
```toml
[dependencies]
qrupdate-ng = { git = "https://github.com/mpimd-csc/qrupdate-ng.git", branch = "main" }
```


## Based on

qrupdate-ng is a modernized and maintained fork of:

[qrupdate](https://sourceforge.net/p/qrupdate) (Version 1.1.2)

For a list of original authors see [CHANGELOG](CHANGELOG.md) file.

## License

License: GPL-3.0-or-later (opensource.org/licenses/GPL-3.0)
See `LICENSE` for details.


