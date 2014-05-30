truthstudy
==========

truthstudy is a program which reads in a root file and extracts particle information from the Particle branch of an Event. The root file, Event number, and number of particles to extract information from are supplied by the user.

A version of truthstudy to be used with root's C++ interpreter, CINT, as well as a compiled version are offered here. To run the interpretted version, cd to its directory and run "root truthstudy.cpp". To compile the compiled version, cd to its directory and run "make".

To compile/run truthstudy, you must source root, install Delphes, and set the environment variable DELPHESPATH to point to your install of Delphes. (export DELPHESPATH=/path/to/Delphes/)

After make succeeds, you can run ./truthstudy. Here is a typical output from it:

Enter the full path to a root file: /Research/outputs/back-large-triple-102.root
There are 50000 Events.
Enter the number of an Event you would like to read: 2
There are 293 particles.
Enter the number of particles you would like information on: 5
Particle: 293
------------------------------------------------------------------------------------------------------------------------------------
   #  PID  Status  M1  M2  D1  D2  Charge    Mass        E        Px        Py        Pz        PT       Eta       Phi     Rapidity 
------------------------------------------------------------------------------------------------------------------------------------
   0 2212       3  -1  -1  -1  -1       1  9.38e-01  4.00e+03  0.00e+00  0.00e+00  4.00e+03  0.00e+00  1.00e+03  0.00e+00  1.00e+03
   1 2212       3  -1  -1  -1  -1       1  9.38e-01  4.00e+03  0.00e+00  0.00e+00 -4.00e+03  0.00e+00 -1.00e+03  0.00e+00 -1.00e+03
   2   21       3   0  -1  -1  -1       0  0.00e+00  1.16e+03 -2.83e-01  3.74e+00  1.16e+03  3.75e+00  6.43e+00  1.65e+00  6.43e+00
   3   21       3   1  -1  -1  -1       0  0.00e+00  3.01e+02 -6.99e-01  3.95e-01 -3.01e+02  8.03e-01 -6.62e+00  2.63e+00 -6.62e+00
   4   21       3   2  -1  -1  -1       0  0.00e+00  7.35e+01 -1.71e+01  3.80e+00  7.13e+01  1.76e+01  2.11e+00  2.92e+00  2.11e+00
------------------------------------------------------------------------------------------------------------------------------------

You will be asked for a path to a root file - you need to specify a root file produced by Delphes. You can then choose an event and how many particles from event are going to be shown. The table shown should be rather obvious. The M1 M2 values are the #s of the particles which produced the particl in question. In the table above, we started with #0 and #1. #2 was produced from #0, #3 was produced from #1, #4 was produced from #3 and so on...
