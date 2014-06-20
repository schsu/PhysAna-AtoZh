cut_study
=========

cut_study takes a ROOT file with Electron, Muon, Jet4, and Particle branches as input, and creates a series of directories and plots as output. The following is a list of the directories and plots it outputs.

cut_study/electron_channel/dielectron_mass.eps
cut_study/electron_channel/dielectron_pt.eps
cut_study/electron_channel/ec_A_mass.eps
cut_study/electron_channel/ec_A_pt.eps
cut_study/electron_channel/ec_A_rapidity.eps
...

The plots in cut_study/electron_channel and cut_study/muon_channel occur after 3 cuts, called ll, jj, and bb. The meanings of these cuts are explained below.

ll
====

Events which pass ll have at exactly 2 leptons, 1 positive and 1 negative. These leptons must have pt greater than some minimum pt and |eta| less than some maximum eta. If the 2 leptons are electrons, the Event is considered as part of the electron channel, and its analysis is carried out separately from Events for which the 2 leptons are muons (this would be the muon channel).

jj
====

Events which pass ll also pass cut1. They also have exactly 2 jets with pt greater than some minimum pt and |eta| less than some maximum eta.

bb
====

Events which pass jj also pass cut2. Their 2 jets must also be b-tagged.

cut_flow
========

In addition to the final plots of physical quantities of interest, such as the reconstructed mass of A, various kinematic plots for each channel for each level of cut are output under cut_study/electron_channel/cut_flow_kinematics and cut_study/muon_channel/cut_flow_kinematics. The efficiencies of electron and muon detection for each cut level are also output by the program to stdout. Comparing the efficiencies and kinematics for each cut level is critical in determining useful cuts.
