LAPACK-MISSING 
==============
Copyright 2020 by Martin Koehler 


This library checks at compile-time if all required LAPACK routines are
available and creates a small library containing the missing ones. This is
required by projects that use IBM's ESSL library since this library does not
contain the full set of LAPACK routines. 

The routines are taken from LAPACK without any modification and thus uses the
same license. 


The cmake logic is licensed under the same terms as. 

