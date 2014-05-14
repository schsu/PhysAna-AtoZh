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



void SelectGoodJets(std::vector<BaseParticle*>& goodJets, TClonesArray* branchJet, const std::vector<BaseParticle*>& goodElectrons) {
  float electronJetCloseness = 0.1;
  for (int i = 0; i < branchJet->GetEntries(); ++i) {
    BaseParticle* jet = GetJetFromBranch(i,branchJet);
    if (jet->PT() > jetMinPt && fabs(jet->Eta()) < jetMaxEta) {
      // electrons are detected as jets too
      // we want to skip these, so we will compare each jet to the previously found electrons
      // and skip it if it is one of them
      bool electronJetMatch = false;
      for (size_t ii = 0; ii < goodElectrons.size(); ++ii) {
        if (goodElectrons[ii]->P4().DeltaR(jet->P4()) < electronJetCloseness){
          electronJetMatch = true;
        }
      }

      if (!electronJetMatch) {
        goodJets.push_back(jet);
      }
    }
  }
}



void cutflow(char* inpuFile){

    //Load share library
    gSystem->Load("/root/root/lib/libDelphes.so");


    //Create chain of root tree
    TChain* tc = new TChain("Delphes");
    tc->Add(inputFile);
    ExRootTreeReader* tr = new ExRootTreeReader(tc);
    TH1F *h1 = new TH1F("Cut flow","Cut flow",100,0,100);

    //Create object of class ExRootTreeReader
    Long64_t numberofEntries = tr->GetEntries();

    //Get pointer to branches
    TClonesArray* branchElectron = tr->UseBranch("Electron");
    TClonesArray* branchJet4 = tr->UseBranch("Jet4");
    TClonesArray* branchJet6 = tr->UseBranch("Jet6");
    TClonesArray* branchJet10 = tr->UseBranch("Jet10");
    TClonesArray* branchMuon = tr->UseBranch("Muon");
    TClonesArray* branchParticle = tr->UseBranch("Particle");


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




  for(Int_t entry;entry < numberofEntries; entry++) {


    //Load selected branches
    tr->ReadEntry(entry);
    Fill("Cut flow", 10);
    bool positive = false;
    bool negative = false;
    int leadingleptonAt = 0;
    int subleadingleptonAt = 0;

    std::vector<BaseParticle*> goodElectrons;
    std::vector<BaseParticle*> goodMuons;
    std::vector<BaseParticle*> Jets4;
    std::vector<BaseParticle*> goodJets6;
    std::vector<BaseParticle*> goodJets4;
    std::vector<BaseParticle*> goodJets10;
    std::vector<BaseParticle*> dilepton;


    //have a pair of lepton pt > minPT , Eta < MaxEta
    for(int i;branchElectron->GetEntries();++i){

            BaseParticle *e1ectron;
            electron = GetElectron(i);

            if(electron->PT() > electronMinPt && fabs(electron->Eta()) < electronMaxEta ){
                goodElectrons.push_back(eletron);
                if(electron->Charge() == 1){
                    positive = true ;
                }
                else if(electron->Charge() == -1){
                    negative = true;
                }
            }

    }
    for(int i;branchMuon->GetEntries();++i){

        BaseParticle *muon;
        muon = GetMuon(i);
        if(muon->PT > muonMinPt && fabs(muon->Eta()) < muonMaxEta){
            goodMuons.push_back(muon);
            if(muon->Charge() == 1){
                positive = true;
            }
            else if(muon->Charge() == -1) {
                negative = true;
            }
        }
    }

    if (goodElectrons.size() == 2 && goodMuons.size() == 0) {
        dilepton.push_back(goodElectrons[0]);
        dilepton.push_back(goodElectrons[1]);
    }
    else if (goodElectrons.size() == 0 && goodMuons.size() == 2) {
        dilepton.push_back(goodMuons[0]);
        dilepton.push_back(goodMuons[1]);
    }

    if((positive == false || negative == false)&& dilepton.size()!=2){  break;    }

    Fill("Cut flow",20);

    // minMass< mll < maxMass
    double dileptonMass = (dilepton[0]->P4() + dilepton[1]->P4()).M();
    if(dileptonMass < dileptonMinmass || dileptonMass > dileptonMaxmass ){  break;  }
    Fill("Cut flow",30);

    // at least two jet4 (jet pT>25GeV, |eta|<2)
    double numerofgoodJet4 = 0;
    for(int i;i<branchJet4->GetEntries();++i){
            Jet* jet4=GetJet4(i);
            if(jet4->PT() > jetMinPt && fabs(jet4->Eta()) < jetMaxEta){
                numerofgoodJet4++;
            }
    }
    Fill("Cut flow",40);

    // at least two b-jet4 (jet pT>25GeV, |eta|<2)
    SelectGoodJets(goodJets4, branchJet4, goodElectrons);
    SelectGoodJets(goodJets6, branchJet6, goodElectrons);
    if(goodJets4.size() < 2 ||goodJets6.size() < 2){  break;    }
    Fill("Cut flow",50);

    // 85<m(jj)<165 GeV
    double dijet4mass = (goodJets4[0]->P4()+goodJets4[1]->P4()).M();
    double dijet6mass = (goodJets6[0]->P4()+goodJets6[1]->P4()).M();
    if((dijet4mass > dijetMaxmass || dijet4mass < dijetMinmass) && (dijet6mass > dijetMaxmass || dijet6mass < dijetMinmass  ){
        break;
    }
    Fill("Cut flow",60);

    h1->Draw();



  }










}
