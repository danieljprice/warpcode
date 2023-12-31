# warpcode

## About

This is a `ring code' for solving the diffusion equation for warped accretion discs
in astrophysics as first derived by Pringle (1992) and with the additional
precession term from Ogilvie (1999). It is modified from an
original code by Jim Pringle, passed down through Giuseppe Lodato
and rewritten by Daniel Price. The code has been extensively rewritten
in Fortran 90 and modernised.

Non-linear coefficients for large amplitude warps are computed
via a routine kindly provided by Gordon Ogilvie. This bit can be 
quite slow so this loop has been parallelised using openMP. For
linear warps the diffusion coefficients are taken from Ogilvie (1999)
as described in Lodato & Price (2010).

## License
This code is Public Domain and free to use, copy and modify.
I offer no warranties.

## Install and run
To compile it you will need a Fortran compiler (gfortran), then:

make
./disc

The default setup is taken from one of the runs in Lodato & Price (2010)

## See also
https://github.com/danieljprice/wavecode


Daniel Price, Dec. 2013
