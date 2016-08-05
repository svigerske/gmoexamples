# gmoexamples
some simple example on using GAMS gmo

To build, create a symlink "gams" that points to the GAMS system directory.
Makefile assumes it can find the GAMS API files in gams/apifiles/C/api.

## loadgms

These are help routines to read a problem instance from a .gms file
that has exactly one Solve statement into a GMO in-memory representation.

## addrow

An executable that given a .gms file with a nonlinear model (e.g., from
MINLPLib2), loads the model into memory (GMO) and adds a nonlinear row
to it.
The instance is printed (using GAMS/ConvertD) to files before.gms and
after.gms before and after the row has been added and the diff is displayed.
