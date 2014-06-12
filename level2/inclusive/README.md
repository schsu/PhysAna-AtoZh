inclusive
=========

inclusive reads a ROOT file with Particle, Electron, and Muon branches, and outputs a series of plots and other information. To compile and run the code, simply run "make". Then run the program with "./inclusive inputFile".

It outputs the following histograms of the Electron and Muon branch kinematics:

electron pt
electron eta
electron phi
muon pt
muon eta
muon phi

It also outputs pt histograms of the electrons and muons in the Particle branch. These are truth electron pt and truth muon pt histograms.

It also outputs histograms of the Electron and Muon multiplicity.

It also compares the number of electrons and muons in the Particle branch to the number of electrons and muons in the Electron and Muon branches. The electron and muon efficiencies are then calculated and output to stdout.

To study the electron and muon efficiency as a function of pt, the electron pt and muon pt histograms are divided by the truth electron pt and truth muon pt histograms.
