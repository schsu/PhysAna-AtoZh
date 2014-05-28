truthstudy
==========

truthstudy is a program which reads in a root file and extracts particle information from the Particle branch of an Event. The root file, Event number, and number of particles to extract information from are supplied by the user.

A version of truthstudy to be used with root's C++ interpreter, CINT, as well as a compiled version are offered here. To run the interpretted version, cd to its directory and run "root truthstudy.cpp". To compile the compiled version, cd to its directory and run "make".

To compile/run truthstudy, you must source root, install Delphes, and set the environment variable DELPHESPATH to point to your install of Delphes. (export DELPHESPATH=/path/to/Delphes/)
