ESSL-LAPACK 
===========

IBM's ESSL library does only provide a small subset of the LAPACK library's
routines. Even simple functions like DLACPY are not available if one uses the
ESSL library. For this reason this project checks for missing LAPACK routines in
a provided library and creates a small static library containg the missing
functions. 


