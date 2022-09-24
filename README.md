### eigenevectors_of_real_symmetric_matrices_with_eispack
Fortran 90/95 adaptation of a couple of EISPACK Fortran 77 routines

This micro library computes the eigenvalues and the eigenvectors of
real symmetric matrices.

These routines have been taken from https://netlib.org/eispack/ on Sept., 2022
and minimally modified to adapt their syntax from Fortran 77 to Fortran 90.
However, there is no guarantee that this operation did not introduce bugs which
were not present in the original version.
Please, use the code exclusively at your own risk!!!

Copyright: since the original version did not contain a Copyright note and it
was developed by people working for public US institutions, it is reasonable
to assume that this is PUBLIC DOMAIN code.
The accompanying example code is released as PUBLIC DOMAIN, too.

Each subroutine or function contains the description of its usage.

To compile the library with gfortran (tested with gcc version 9.4.0):
gfortran -c eispack_symm_matrix.f90

To compile the example program:
gfortran test_eispack.f90 eispack.a -o test

To execute the example program:
./test
