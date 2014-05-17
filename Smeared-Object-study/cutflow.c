#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TImage.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>


void cutflow(char* inputFile){


    //Create chain of root tree
    TChain* tc = new TChain("Delphes");
    tc->Add(inputFile);
    ExRootTreeReader* tr = new ExRootTreeReader(tc);
    TH1F *h1 = new TH1F("Cut flow","Cut flow",101,0,100);

    //Create object of class ExRootTreeReader
    Long64_t numberofEntries = tr->GetEntries();

    //Get pointer to branches
    TClonesArray* branchElectron = tr->UseBranch("Electron");
    TClonesArray* branchJet4 = tr->UseBranch("Jet4");
    TClonesArray* branchJet6 = tr->UseBranch("Jet6");
    TClonesArray* branchJet10 = tr->UseBranch("Jet10");
    TClonesArray* branchMuon = tr->UseBranch("Muon");
    //TClonesArray* branchParticle = tr->UseBranch("Particle");


    double electronMinPt = 25.0;
    double electronMaxEta = 2.4;
    double muonMinPt = 5.0;
    double muonMaxEta = 2.5;
    double jetMinPt = 25.0;

    double jetMaxEta = 2;
    double dileptonMaxmass = 100;
    double dileptonMinmass = 80;
    double dijetMaxmass = 165;
    double dijetMinmass = 85;




  for(Int_t entry = 0;entry < numberofEntries; entry++) {


    //Load selected branches
    tr->ReadEntry(entry);
    h1->Fill(10);
    bool positive = false;
    bool negative = false;
    int leadingleptonAt = 0;
    int subleadingleptonAt = 0;

    std::vector<Electron> goodElectrons;
    std::vector<Muon> goodMuons;
    std::vector<Jet> goodJets6;
    std::vector<Jet> goodJets4;
    std::vector<Jet> goodJets10;

    //have a pair of lepton pt > minPT , Eta < MaxEta
    for(int i=0; i<branchElectron->GetEntries();++i){
            Electron *electron = (Electron*)branchElectron->At(i);
            if(electron->PT > electronMinPt && fabs(electron->Eta) < electronMaxEta ){
                //goodElectrons.push_back(*electron);
                if(electron->Charge == 1){
                    positive = true ;
                }
                else if(electron->Charge == -1){
                    negative = true ;
                }
            }

    }

    for(int i=0;i<branchMuon->GetEntries();++i){
        Muon *muon = (Muon *)branchMuon->At(i);
        if(muon->PT > muonMinPt && fabs(muon->Eta) < muonMaxEta){
            goodMuons.push_back(*muon);
            if(muon->Charge == 1){
                positive = true;
            }
            else if(muon->Charge == -1) {
                negative = true;
            }
        }
    }
    int leptonNumber = 0;
    double dileptonMass = 0;
    if (goodElectrons.size() == 2 && goodMuons.size() == 0) {
    	dileptonMass = (goodElectrons[0].P4() + goodElectrons[1].P4()).M();    
	leptonNumber = 2 ;
	//dilepton.push_back(goodElectrons[0]);
        //dilepton.push_back(goodElectrons[1]);
    }
    else if (goodElectrons.size() == 0 && goodMuons.size() == 2) {
        //dilepton.push_back(goodMuons[0]);
        //dilepton.push_back(goodMuons[1]);
        leptonNumber = 2 ;
        dileptonMass = (goodMuons[0].P4() + goodMuons[1].P4()).M();
    }

    if( positive == false || negative == false|| leptonNumber != 2){    continue;    }


    h1->Fill(20);

    // minMass< mll < maxMass
    //double dileptonMass = (dilepton[0]->P4() + dilepton[1]->P4()).M();
    if(dileptonMass < dileptonMinmass || dileptonMass > dileptonMaxmass ){  continue;  }
    h1->Fill(30);

    // at least two jet4 (jet pT>25GeV, |eta|<2)
    double numerofgoodJet4 = 0;
    for(int i=0;i<branchJet4->GetEntries();++i){
            Jet* jet4=(Jet*) branchJet4->At(i);
            if(jet4->PT > jetMinPt && fabs(jet4->Eta) < jetMaxEta){
                numerofgoodJet4++;
            }
    }
    h1->Fill(40);

    // at least two b-jet4 (jet pT>25GeV, |eta|<2)

    float electronJetCloseness = 0.1;
    for (int i = 0; i < branchJet4->GetEntries(); ++i) {
      Jet* jet =(Jet *) branchJet4->At(i);
      if (jet->PT > jetMinPt && fabs(jet->Eta) < jetMaxEta) {
      // electrons are detected as jets too
      // we want to skip these, so we will compare each jet to the previously found electrons
      // and skip it if it is one of them
          bool electronJetMatch = false;
          for (size_t ii = 0; ii < goodElectrons.size(); ++ii) {
             if (goodElectrons[ii].P4().DeltaR(jet->P4()) < electronJetCloseness){
                 electronJetMatch = true;
             }
          }
          if (!electronJetMatch) {
             goodJets4.push_back(*jet);
          }
      }
    }
    
    for (int i = 0; i < branchJet6->GetEntries(); ++i) {
      Jet* jet = (Jet *) branchJet6->At(i);
      if (jet->PT > jetMinPt && fabs(jet->Eta) < jetMaxEta) {
        bool electronJetMatch = false;
          for (size_t ii = 0; ii < goodElectrons.size(); ++ii) {
             if (goodElectrons[ii].P4().DeltaR(jet->P4()) < electronJetCloseness){
                 electronJetMatch = true;
             }
          }
          if (!electronJetMatch) {
             goodJets6.push_back(*jet);
          }
      }
    }

    if(goodJets4.size() < 2 ||  goodJets6.size() < 2){  continue;    }
    h1->Fill(50);

    // 85<m(jj)<165 GeV
    double dijet4mass = (goodJets4[0].P4()+goodJets4[1].P4()).M();
    double dijet6mass = (goodJets6[0].P4()+goodJets6[1].P4()).M();
    if((dijet4mass > dijetMaxmass || dijet4mass < dijetMinmass) && (dijet6mass > dijetMaxmass || dijet6mass < dijetMinmass ) ){
        continue;
    }
    h1->Fill(60);

  

  }
  
  h1->Draw();
  //TAxis *x = h->GetXaxis()();
  //float t = h1->GetBinContent(10.1);
  //std::cout << t << std::endl ;
  //std::cout << h1->GetBinContent(10) << std::end;

}

