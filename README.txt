### eigenevectors_of_general_and_real_symmetric_matrices_with_eispack

Fortran 90/95 adaptation of a couple of EISPACK Fortran 77 routines

This (micro) library computes the eigenvalues and the eigenvectors of
general and real symmetric matrices.

These routines have been taken from https://netlib.org/eispack/ on Sept., 2022
(for the real symmetric mat. solver) and on July, 2023 (for the general solver)
and minimally modified to adapt their syntax from Fortran 77 to Fortran 90.
However, there is no guarantee that this operation did not introduce bugs which
were not present in the original version.
Please, use the code exclusively at your own risk!!!

Copyright: since the original version did not contain a Copyright note and it
was developed by people working for public US institutions, it is reasonable
to assume that this is PUBLIC DOMAIN code.
The accompanying example code is released as PUBLIC DOMAIN, too.

Each subroutine or function contains the description of its usage.

To compile the (micro) library code with gfortran (tested with gcc version 9.4.0):
gfortran -c eispack_real_symm_matrix.f90

To create the (micro) library:
ar r eispack_real_symm_matrix.a eispack_real_symm_matrix.o

To compile the example program:
gfortran test_eispack.f90 eispack_real_symm_matrix.a -o test

To execute the example program:
./test
