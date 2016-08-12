# gmoexamples
Simple examples on using the GAMS GMO API.

To build, create a symlink "gams" that points to the GAMS system directory.
Makefile assumes it can find the GAMS API files in gams/apifiles/C/api.

## loadcntr
A simple main() program that shows how to load a compiled GAMS model
from a GAMS control file (225? directory) into a GMO object.

Building with make should create an executable ```./loadcntr```.

Execute

    rm -r 225*
    gamslib trnsport
    gams trnsport keep=1
    ./loadcntr 225a/gamscntr.dat

This executes GAMS for the trnsport model from the GAMS model library
and keeps the control file.
Then it calls the loadcntr executable on the generated control file.

This example may serve as a very basic starting point for setting up an
own solver in GAMS.
To install a solver executable (gmsxx_ux.out) in GAMS, do the following:

1. Create a batch file gmsxx_us.run with the line

        ${5}gmsxx_ux.out $4

2. Rename your executable to ```gmsxx_ux.out```

3. Copy ```gmsxx_us.run``` and ```gmsxx_ux.out``` to the GAMS system
   directory and make them executable.

4. Edit ```gmscmpun.txt``` in the GAMS system directory by adding

        MYSOLVER 2111 5 XX 1 0 1 MIP QCP MIQCP RMIQCP NLP DNLP RMINLP MINLP CNS
        gmsxx_us.run
        gmsxx_ux.out

5. Enjoy:
   You can now call your solver in GAMS via ```<MODELTYPE>=MYSOLVER```.

In steps 1-4, replace xx by your own abbreviation.
In step 4, replace MYSOLVER by your own solver name.
In step 4, edit the list of support model types as you see fit.



## loadgms

These are help routines to read a problem instance from a .gms file
that has exactly one Solve statement into a GMO in-memory representation.

## addrow

An executable that given a .gms file with a nonlinear model (e.g., from
MINLPLib2), loads the model into memory (GMO) and adds a nonlinear row
to it.
The instance is printed (using GAMS/ConvertD) to files before.gms and
after.gms before and after the row has been added and the diff is displayed.
