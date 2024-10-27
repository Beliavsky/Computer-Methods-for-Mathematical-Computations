# Computer Methods for Mathematical Computations
Code from the book by Forsythe, Malcolm, and Moler, translated to Fortran 90 by [Ralph Carmichael](https://www.pdas.com/fmm.html). Compile with `gfortran fmm.f90 testfunctions.f90 samples.f90`. Also compiles and runs with the ifx, ifort, and g95 compilers. The original routines in fixed format are in `original.f`. Below is a list of procedures, taken from [here](https://www.pdas.com/fmm2.html).

**Name** Description<br>
**Decomp** LU-decomposition of a square matrix<br>
**Solve**	Solves a system of linear equations. Use after Decomp.<br>
**FMMspline**	Fit a cubic spline to data. FMM end conditions.<br>
**NaturalSpline**	Same as FMMspline, but with zero second derivatives at endpoints.<br>
**Seval**	Evaluate a cubic spline at a given point.<br>
**Seval3** Evaluate a cubic spline at a given point. Returns value of spline plus 1st, 2nd, and 3rd derivatives.<br>
**Quanc8** Numerical integration of a function.<br>
**Rkf45**	Solves a system of ordinary differential equations as an initial value problem.<br>
**Zeroin** Find a zero of a function.<br>
**BrentZero**	Same as Zeroin, but with additional dummy arguments.<br>
**Fmin** Find the minimum of a function.<br>
**BrentMin** Same as Fmin, but with additional dummy arguments.<br>
**SVD**	Singular Value Decomposition of a matrix.<br>
