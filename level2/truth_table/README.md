truth_table
===========

truth_table is a program which reads in a root file and extracts particle information from the Particle branch of an Event. The root file, Event number, and number of particles to extract information from are supplied by the user.

To compile and run truth_table, run "make" and then "./truth_table". A usage statement will be printed showing the necessary command line input.

To compile/run truth_table, you must source root, install Delphes, and set the environment variable DELPHESPATH to point to your install of Delphes. (export DELPHESPATH=/path/to/Delphes/)

Here is a typical output from truth_table:

Enter the full path to a root file: /Research/outputs/back-large-triple-102.root
There are 50000 Events.
Enter the number of an Event you would like to read: 2
There are 293 particles.
Enter the number of particles you would like information on: 5
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

You can also supply the program with an Event number and the number of particles to read without the command line interaction, for example:

[nikola@corneria truth_table]$ ./truth_table ~/Research/ATLAS/Analysis/data/tag_1_pythia_events.root 10 6
------------------------------------------------------------------------------------------------------------------------------------
   #  PID  Status  M1  M2  D1  D2  Charge    Mass        E        Px        Py        Pz        PT       Eta       Phi     Rapidity 
------------------------------------------------------------------------------------------------------------------------------------
   0 2212       3  -1  -1  -1  -1       1  9.38e-01  7.00e+03  0.00e+00  0.00e+00  7.00e+03  0.00e+00  1.00e+03  0.00e+00  1.00e+03
   1 2212       3  -1  -1  -1  -1       1  9.38e-01  7.00e+03  0.00e+00  0.00e+00 -7.00e+03  0.00e+00 -1.00e+03  0.00e+00 -1.00e+03
   2   21       3   0  -1  -1  -1       0  0.00e+00  2.92e+02 -3.07e-01  1.07e-01  2.92e+02  3.25e-01  7.49e+00  2.81e+00  7.49e+00
   3   21       3   1  -1  -1  -1       0  0.00e+00  4.64e+02 -9.28e-01  1.39e+00 -4.64e+02  1.67e+00 -6.32e+00  2.16e+00 -6.32e+00
   4   21       3   2  -1  -1  -1       0  0.00e+00  1.11e+02  1.62e+01 -4.87e+00  1.10e+02  1.69e+01  2.57e+00 -2.91e-01  2.57e+00
   5   21       3   3  -1  -1  -1       0  0.00e+00  2.04e+02  1.79e+00  8.52e+00 -2.03e+02  8.71e+00 -3.84e+00  1.36e+00 -3.84e+00
------------------------------------------------------------------------------------------------------------------------------------

