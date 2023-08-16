### eigenvectors_of_general_and_real_symmetric_matrices_with_eispack

Fortran 90/95 adaptation of a few selected EISPACK Fortran 77 routines

This tiny library computes the eigenvalues and the eigenvectors of
real symmetric matrices and it solves the generalized eigenvalues problem
for A . x = lambda B . x, where A and B are generic matrices, x are the
eigenvectors and lambda the eigenvalues.

These routines have been taken from https://netlib.org/eispack/ on Sept., 2022
(for the real symmetric mat. solver) and on July, 2023 (for the general solver)
and minimally modified to adapt their syntax from Fortran 77 to Fortran 90.
However, there is no guarantee that this operation did not introduce bugs which
were not present in the original version.
Please, use the code exclusively at your own risk!!!

**Copyright**: since the original version did not contain a Copyright note and it
was developed by people working for public US institutions, it is reasonable
to assume that **this is PUBLIC DOMAIN code**.
The accompanying example code is released as PUBLIC DOMAIN, too.

Each subroutine or function contains the description of its usage.

#### Compilation and usage ####
 
To compile the library code with gfortran (tested with gcc version 11.3.0)
simply type:

```
make
```

You will get a static library file `eispack.a` and an executable `test.exe`
with the sample test program.

To clean the directory with the sources type:

```
make clean
```

To compile a f90 program named "my_program.f90" tpye:

```
gfortran my_program.f90 -o my_program.exe eispack.a
```

#### About the test/example code

In the first two tests we use the subroutine rs to solve the eigenvalue problem:

$ A . x = \lambda x $

where A is a real symmetric matric, x is an eigenvector and \lambda its eigenvalue.
 
In the other tests we use the subroutine rgg to solve the generalized eigevalue problem:

$ A . x = \lambda B . x $

where A and B are generic matrices, x is an eigenvector and \lambda its eigenvalue.
In tests 3 and 4 B is simply the identity matrix, so we effectively solving a standard
eigenvalue problem, but with generic matrices, which can have also complex eigenvalues
and eigenvectors. In the last test B is a generic matrix and we really solve the generalized problem.
 
The advantage of the rs solver is that it is a bit faster and simpler to handle than the rgg solver
(if A is a symmetric real matrix, the eigenvalues and the eigenvector components are real numbers).
The rgg solver is more generic, but a bit slower and one must be prepared to handle different types of outputs
(the output fields are different depending on whether the results are real or complex number,
see the documentation of the subroutine in the comments of the code).

The source code of test_eispack.f90 provides a few examples of how to use the library, while the execution
of the executable serves also as a test, because we verify that the solutions lambda and x satisfy

$ A . x - \lambda x = 0 $ and $ A . x - \lambda B . x = 0 $

by selecting, for each eigenvector, the biggest absolute value of the LHS of the expressions above and
displying it. If the library works fine, the results must be very small (~ < 10^-14) numbers.
Please, note that, in case of wrong results, no warning in printed, one has always to check the output.
Please, also note that currently no example of complex eigenvalue/eigenvector is given.
