***********************************************************************
*                                                                     *
*            Instructions on How to Use the TCas2HDM Model            *
*                                                                     *
***********************************************************************

I.    Introduction

The 2HDM4TC model is intended to be used as a model for strong dynamics in the framework of a two Higgs doublet model.  In certain regions of parameter space, a two Higgs doublet model closely resembles a low energy effective theory of strongly interacting electroweak symmetry breaking.  Some of the more exotic multiple electroweak boson and top quark final states can be simulated easily here.  It was for predominantly this reason that the model was created.

II.   Basic Use

Go into the 2HDM4TCCalc/bin folder and modify TCCalc.txt to contain the parameters you desire.  After modifying the file, run the program by typing,

./2HDM4TCCalc TCCalc.txt param_card.dat

The new param_card is ready to use with the model.  A by hand modification of param_card is discouraged because this will generically not have the correct width.  Major contributions for the widths (including some 3-body decays) are included in the calculator.  Additionally, warnings about lagrangian parameters growing too large are included.

After using the calculator, the model may be used as normal.  The only advantage of using the calculator over a program such as BRIDGE (by Patrick Meade and Matt Reece) is speed.  For a detailed analysis, the width computation by BRIDGE is recommended.

III.  Program Information

a.  GGF vertices - The program includes GGF vertices for production of the scalar resonances.  These are not perturbative and will work to any energy.  Additionally, these vertices can be made momentum dependent, which improves accuracy when dealing with processes involving broad Higgs widths (see instructions in couplings.f for more details).  Inital state radiation is not supported.  Generically, adding ISR to the diagram will make it incorrect.

More info... Coming Soon!!! ...maybe?

Please email questions and comments to Jared Evans (jaevans@ucdavis.edu)
