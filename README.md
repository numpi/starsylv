# StarSylv

StarSylv implements a solver for Sylvester-like equations as described 
in the following paper: 

> <em>Unique solution of systems of generalized Sylvester equations: 
> an algorithmic approach</em>, F. De Terán, B. Iannazzo, F. Poloni, 
> L. Robol, to be submitted.

The algorithm solves a system of <code>R</code> T-Sylvester equations, 
with <code>N x N</code> coefficients. Only the last equation has a transposition,
and all are required to have triangular coefficients (upper or lower 
depending on which side of the unkonwns they multiply). This form is not 
restrictive, and can be obtained by performing a preliminary re-ordering 
of the equations and a periodic Schur decomposition. 
The details can be found in the above paper, and the code for the 
preliminary reduction is not included in this repository.

The code in the <code>fortran</code> subdirectory is the one used 
in the paper to run the tests. It can be compiled and run entering 
in the subdirectory and running 
```
 make ; ./starsylv N R
 ```
where <code>N</code> and <code>R</code> are two integers that specify 
the dimension of the coefficients and the number of equations to be used 
in the test. The program generates a random example of the given size, 
and returns the time needed to solve it and the backward error on the linear 
system equivalent to the matrix equation.
