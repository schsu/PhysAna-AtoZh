inclusive
=========

inclusive takes a ROOT file with Particle, Electron, Muon, and Jet branches as input, creates the directory "inclusive" in the directory where the program is executed, and outputs histograms in the directory "inclusive". To compile and run the code, simply run "make". Then run the program with "./inclusive inputFile".

It outputs the following histograms of the Electron, Muon, and Jet branch kinematics:

electron pt
electron E
electron eta
electron phi
muon pt
muon E
muon eta
muon phi
jet pt
jet E
jet eta
jet phi

It also outputs pt histograms of the electrons and muons in the Particle branch. These are truth electron pt and truth muon pt histograms.

It also outputs histograms of the Electron and Muon multiplicity.
The multiplicity of an electron is the number of electrons per event. Since in the process we are studying, our final state particles are 2 leptons and 2 jets, we expect to see many events with either 2 electrons and 0 muons or 2 muons and 0 electrons. Thus, multiplicity plots offer a good sanity check.

It also compares the number of electrons and muons in the Particle branch to the number of electrons and muons in the Electron and Muon branches. The electron and muon efficiencies are then calculated and output to stdout.
Not all particles are detected by the ATLAS detector, so it is important to do efficiency studies with simulated data where we can compare truth particle information with the final output of our detector simulation (Delphes).

To study the electron and muon efficiency as a function of pt, the electron pt and muon pt histograms are divided by the truth electron pt and truth muon pt histograms.
One would expect electrons and muons of low pt to be detected less often than electrons and muons of high pt. This is evident from the cylindrical geometry of the ATLAS detector - particles of high pt are likely to hit/pass through the walls of the cylindar (being detected), while particles of low pt are likely to hit/pass through the end caps (not being detected). Comparing the detection efficiency of particles which do go through the cylindar walls is another important measure which is not studied here at this time.

The studies done on electron and muon efficiency well motivate efficiency, multiplicity, and kinematic studies after cutting out events with low pt and high eta. These types of studies are carried out by the program cutflow. Understanding the program cutflow should be the next step in this analysis tutorial.
