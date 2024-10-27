

COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS                   /fmm/readme.txt

The files for this program are in the directory fmm 

  readme.txt      general description
  fmm.f90         the source code in modern Fortran
  original.src    the source code as published with the book
  samples.f90     the source code for the sample problems in the book
  samples.out     the output expected from executing samples.exe

To compile the sample case program, use the command
   gfortran samples.f90 -o samples.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.


This collection of Fortran procedures for mathematical computation is based 
on the procedures from the book "Computer Methods for Mathematical 
Computations", by George E. Forsythe, Michael A. Malcolm, and Cleve B. Moler. 
Prentice-Hall, 1977.
The book has been a classic textbook for numerical analysis since its publication
in 1977. You will need a copy of the book to make good use of the procedures,
although the subroutine headers are quite informative and together with the
samples, should provide a minimal set of instructions. The PDAS web site has 
additional information at the page  https://www.pdas.com/fmm.html
