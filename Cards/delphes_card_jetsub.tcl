#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger
  Calorimeter
  EFlowMerger

  PhotonEfficiency
  PhotonIsolation
  
  ElectronEfficiency
  
  MuonEfficiency
  
  MissingET
  
  GenJetFinder
  
  FastJetFinder4
  FastJetFinder6
  FastJetFinder10
  
  JetEnergyScale4
  JetEnergyScale6
  JetEnergyScale10
  
  BTagging4
  BTagging6
  BTagging10
  TauTagging4
  TauTagging6
  TauTagging10
  
  UniqueObjectFinder4
  UniqueObjectFinder6
  UniqueObjectFinder10

  ScalarHT4
  ScalarHT6
  ScalarHT10
  
  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.15
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.51

  # magnetic field
  set Bz 2.0
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.60) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.85) + \
                         (abs(eta) > 2.5)                                                  * (0.00)}
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
              
  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.73) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2) * (0.95) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.99) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.50) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.83) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.90) + \
                         (abs(eta) > 2.5)                                                  * (0.00)}
}

##############################
# Muon tracking efficiency
##############################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
              
  # tracking efficiency formula for muons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.99) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.98) + \
                         (abs(eta) > 2.5)                                                  * (0.00)}
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.01) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.03) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.03) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.02) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05)}
}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  set ResolutionFormula {                  (abs(eta) <= 2.5) * (energy > 0.1   && energy <= 2.5e1) * (energy*0.015) + \
                                           (abs(eta) <= 2.5) * (energy > 2.5e1)                    * sqrt(energy^2*0.005^2 + energy*0.05^2 + 0.25^2) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0)                                       * sqrt(energy^2*0.005^2 + energy*0.05^2 + 0.25^2) + \
                         (abs(eta) > 3.0 && abs(eta) <= 5.0)                                       * sqrt(energy^2*0.107^2 + energy*2.08^2)}

}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for muons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.03) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.03) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.03) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05)}
}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}

#############
# Calorimeter
#############

module Calorimeter Calorimeter {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray towers
  set PhotonOutputArray photons

  set EFlowTrackOutputArray eflowTracks
  set EFlowPhotonOutputArray eflowPhotons
  set EFlowNeutralHadronOutputArray eflowNeutralHadrons

  set pi [expr {acos(-1)}]
  
  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 10 degrees towers
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} { 
    add PhiBins [expr {$i * $pi/18.0}]
  }
  foreach eta {-3.2 -2.5 -2.4 -2.3 -2.2 -2.1 -2 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.3 -1.2 -1.1 -1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 3.3} { 
    add EtaPhiBins $eta $PhiBins
  }
 
  # 20 degrees towers
  set PhiBins {}
  for {set i -9} {$i <= 9} {incr i} { 
    add PhiBins [expr {$i * $pi/9.0}]
  }  
  foreach eta {-4.9 -4.7 -4.5 -4.3 -4.1 -3.9 -3.7 -3.5 -3.3 -3 -2.8 -2.6 2.8 3 3.2 3.3 3.5 3.7 3.9 4.1 4.3 4.5 4.7 4.9} { 
    add EtaPhiBins $eta $PhiBins
  }

  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {0.0 1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0 0.0}
  add EnergyFraction {22} {1.0 0.0}
  add EnergyFraction {111} {1.0 0.0}
  # energy fractions for muon and neutrinos
  add EnergyFraction {12} {0.0 0.0}
  add EnergyFraction {13} {0.0 0.0}
  add EnergyFraction {14} {0.0 0.0}
  add EnergyFraction {16} {0.0 0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3 0.7}
  add EnergyFraction {3122} {0.3 0.7}
  
  # set ECalResolutionFormula {resolution formula as a function of eta and energy}
  # http://arxiv.org/pdf/physics/0608012v1 jinst8_08_s08003
  # http://villaolmo.mib.infn.it/ICATPP9th_2005/Calorimetry/Schram.p.pdf
  set ECalResolutionFormula {                  (abs(eta) <= 3.2) * sqrt(energy^2*0.0017^2 + energy*0.101^2) + \
                             (abs(eta) > 3.2 && abs(eta) <= 4.9) * sqrt(energy^2*0.0350^2 + energy*2.085^2)}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  # http://arxiv.org/pdf/hep-ex/0004009v1
  # http://villaolmo.mib.infn.it/ICATPP9th_2005/Calorimetry/Schram.p.pdf
  set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) + \
                             (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.0500^2 + energy*0.706^2) + \
                             (abs(eta) > 3.2 && abs(eta) <= 4.9) * sqrt(energy^2*0.9420^2 + energy*0.075^2)}
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray Calorimeter/eflowTracks
  add InputArray Calorimeter/eflowPhotons
  add InputArray Calorimeter/eflowNeutralHadrons
  set OutputArray eflow
}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray Calorimeter/photons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.85) + \
                         (abs(eta) > 2.5)                                   * (0.00)}
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowMerger/eflow
  
  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.1
}

#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronEnergySmearing/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.85) + \
                         (abs(eta) > 2.5)                                   * (0.00)}
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowMerger/eflow
  
  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.1
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray MuonMomentumSmearing/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.7) * (pt > 10.0)  * (0.85) + \
                         (abs(eta) > 2.7)                                   * (0.00)}
}

################
# Muon isolation
################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set IsolationInputArray EFlowMerger/eflow
  
  set OutputArray muons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.1
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT4 {
# add InputArray InputArray
  add InputArray UniqueObjectFinder4/jets4
  add InputArray UniqueObjectFinder4/electrons
  add InputArray UniqueObjectFinder4/photons
  add InputArray MuonEfficiency/muons
  set EnergyOutputArray energy
}

module Merger ScalarHT6 {
# add InputArray InputArray
  add InputArray UniqueObjectFinder6/jets6
  add InputArray UniqueObjectFinder6/electrons
  add InputArray UniqueObjectFinder6/photons
  add InputArray MuonEfficiency/muons
  set EnergyOutputArray energy
}

module Merger ScalarHT10 {
# add InputArray InputArray
  add InputArray UniqueObjectFinder10/jets10
  add InputArray UniqueObjectFinder10/electrons
  add InputArray UniqueObjectFinder10/photons
  add InputArray MuonEfficiency/muons
  set EnergyOutputArray energy
}

#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray Delphes/stableParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.6

  set JetPTMin 20.0
}

############
# Jet finder
############

module FastJetFinder FastJetFinder4 {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow
  
  set OutputArray jets4

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 N-subjettiness, 8 N-jettiness 
  set JetAlgorithm 6
  set ParameterR 0.4

  set ConeRadius 0.5
  set SeedThreshold 1.0
  set ConeAreaFraction 1.0
  set AdjacencyCut 2
  set OverlapThreshold 0.75

  set MaxIterations 100
  set MaxPairSize 2
  set Iratch 1
  
  set ComputeNsubjettiness 1
  set Beta  2.0
  # AxisMode - 1 wta kt axes, 2 optimised wta kt axes, 3 kt axes, 4 optimized kt axes 
  set AxisMode 3
  
  set JetPTMin 20.0
}

module FastJetFinder FastJetFinder6 {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow
  
  set OutputArray jets6

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.6

  set ConeRadius 0.5
  set SeedThreshold 1.0
  set ConeAreaFraction 1.0
  set AdjacencyCut 2
  set OverlapThreshold 0.75

  set MaxIterations 100
  set MaxPairSize 2
  set Iratch 1
  
  set ComputeNsubjettiness 1
  set Beta  2.0
  # AxisMode - 1 wta kt axes, 2 optimised wta kt axes, 3 kt axes, 4 optimized kt axes 
  set AxisMode 3
  
  set JetPTMin 20.0
}

module FastJetFinder FastJetFinder10 {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow
  
  set OutputArray jets10

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.0

  set ConeRadius 0.5
  set SeedThreshold 1.0
  set ConeAreaFraction 1.0
  set AdjacencyCut 2
  set OverlapThreshold 0.75

  set MaxIterations 100
  set MaxPairSize 2
  set Iratch 1
  
  set ComputeNsubjettiness 1
  set Beta  2.0
  # AxisMode - 1 wta kt axes, 2 optimised wta kt axes, 3 kt axes, 4 optimized kt axes 
  set AxisMode 3
  
  set JetPTMin 20.0
}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale4 {
  set InputArray FastJetFinder4/jets4
  set OutputArray jets4

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale6 {
  set InputArray FastJetFinder6/jets6
  set OutputArray jets6

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale10 {
  set InputArray FastJetFinder10/jets10
  set OutputArray jets10

 # scale formula for jets
  set ScaleFormula {1.00}
}


###########
# b-tagging
###########

module BTagging BTagging4 {
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale4/jets4

  set DeltaR 0.3

  set PartonPTMin 1.0

  set PartonEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {0.1}
  # efficiency formula for b-jets
  add EfficiencyFormula {5} {0.4}
}

module BTagging BTagging6 {
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale6/jets6

  set DeltaR 0.3

  set PartonPTMin 1.0

  set PartonEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {0.1}
  # efficiency formula for b-jets
  add EfficiencyFormula {5} {0.4}
}

module BTagging BTagging10 {
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale10/jets10

  set DeltaR 0.3

  set PartonPTMin 1.0

  set PartonEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {0.1}
  # efficiency formula for b-jets
  add EfficiencyFormula {5} {0.4}
}

module TauTagging TauTagging4 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale4/jets4

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging6 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale6/jets6

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging10 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale10/jets10

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder4 {
# earlier arrays take precedence over later ones
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronEfficiency/electrons electrons
  add InputArray JetEnergyScale4/jets4 jets4
}

module UniqueObjectFinder UniqueObjectFinder6 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronEfficiency/electrons electrons
  add InputArray JetEnergyScale6/jets6 jets6
}

module UniqueObjectFinder UniqueObjectFinder10 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronEfficiency/electrons electrons
  add InputArray JetEnergyScale10/jets10 jets10
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
  
  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower
  
  add Branch Calorimeter/eflowTracks EFlowTrack Track
  add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
  add Branch Calorimeter/eflowNeutralHadrons EFlowNeutralHadron Tower
  
  add Branch GenJetFinder/jets GenJet Jet
  
  add Branch UniqueObjectFinder4/jets4 Jet4 Jet
  add Branch UniqueObjectFinder6/jets6 Jet6 Jet
  add Branch UniqueObjectFinder10/jets10 Jet10 Jet
  
  add Branch UniqueObjectFinder4/electrons Electron Electron
  add Branch UniqueObjectFinder4/photons Photon Photon
  add Branch MuonEfficiency/muons Muon Muon
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT4/energy ScalarHT ScalarHT
}

