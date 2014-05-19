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
    TH1F *h1 = new TH1F("Cut flow ,r","Cut flow ,r",101,0,100);
    TH1F *h2 = new TH1F("Cut flow ,b","Cut flow ,b",111,0,110);

    //Create object of class ExRootTreeReader
    Long64_t numberofEntries = tr->GetEntries();

    //Get pointer to branches
    TClonesArray* branchElectron = tr->UseBranch("Electron");
    TClonesArray* branchJet4 = tr->UseBranch("Jet4");
    TClonesArray* branchJet8 = tr->UseBranch("Jet8");
    TClonesArray* branchJet12 = tr->UseBranch("Jet12");
    TClonesArray* branchMuon = tr->UseBranch("Muon");
    //TClonesArray* branchParticle = tr->UseBranch("Particle");


    double electronMinPt = 5.0;
    double electronMaxEta = 2.5;
    double muonMinPt = 5.0;
    double muonMaxEta = 2.5;
    double jetMinPt = 25.0;
    double jet8MinPt = 315.0;
    double jet12MinPt = 210.0;


    double jetMaxEta = 2.0;
    double dileptonMaxmass = 100;
    double dileptonMinmass = 80;
    double dijetMaxmass = 141;
    double dijetMinmass = 111;
    double boostjetMaxmass = 146;
    double boostjetMinmass = 106;




  for(Int_t entry = 0;entry < numberofEntries; entry++) {


    //Load selected branches
    tr->ReadEntry(entry);

    //all the events
    h1->Fill(10);
    h2->Fill(10);
    bool positive = false;
    bool negative = false;
    

    std::vector<Electron> goodElectrons;
    std::vector<Muon> goodMuons;
    std::vector<Jet> goodJets8;
    std::vector<Jet> goodJets4;
    std::vector<Jet> goodJets12;

    //events have a pair of lepton pt > minPT , |Eta| < MaxEta
    for(int i=0; i < branchElectron->GetEntries() ;++i){
            Electron *electron = (Electron*)branchElectron->At(i);
            if(electron->PT > electronMinPt && fabs(electron->Eta) < electronMaxEta ){
                goodElectrons.push_back(*electron);
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
    //std::cout << goodElectrons.size() << "  "  <<  goodMuons.size() ;
    if (goodElectrons.size() == 2 && goodMuons.size() == 0) {
    	dileptonMass = (goodElectrons[0].P4() + goodElectrons[1].P4()).M();
	leptonNumber = 2 ;

    }
    else if (goodElectrons.size() == 0 && goodMuons.size() == 2) {
        leptonNumber = 2 ;
        dileptonMass = (goodMuons[0].P4() + goodMuons[1].P4()).M();
    }

    if( positive == false || negative == false|| leptonNumber != 2){    continue;    }
    h1->Fill(20);
    h2->Fill(20);


    // minMass< m(ll) < maxMass
    if(dileptonMass < dileptonMinmass || dileptonMass > dileptonMaxmass ){  continue;  }
    h1->Fill(30);
    h2->Fill(30);

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

    bool resolved = false ;
    // at least two jet4 (jet pT>25GeV, |eta|<2)
    //resolved channel
    if(goodJets4.size() >= 2 ){
        h1->Fill(40);
        std::vector<Jet> btagJet4;
        for(size_t i=0 ; i < goodJets4.size() ; i++){
            if(goodJets4[i].BTag == 1){
                btagJet4.push_back(goodJets4[i]);
            }
        }
	//std::cout << " batg4 " << btagJet4.size() ;

        // at least two b-jet4 (jet pT>25GeV, |eta|<2)
        if( btagJet4.size() >=2 ){
            h1->Fill(50);
            // 111GeV < m(jj) <141 GeV
            double dijet4mass = (btagJet4[0].P4()+btagJet4[1].P4()).M();
            if((dijet4mass > dijetMaxmass || dijet4mass < dijetMinmass)){
            h1->Fill(60);
	    resolved = true;
            }
	    
        }

    }

    //merge channel
    if(resolved == false ){
        h2->Fill(40);


        //select goodJet8 (jet pT>315GeV, |eta|<2)
        for (int i = 0; i < branchJet8->GetEntries(); ++i) {
            Jet* jet = (Jet *) branchJet8->At(i);
            if (jet->PT > jet8MinPt && fabs(jet->Eta) < jetMaxEta) {
                bool electronJetMatch = false;
                for (size_t ii = 0; ii < goodElectrons.size(); ++ii) {
                    if (goodElectrons[ii].P4().DeltaR(jet->P4()) < electronJetCloseness){
                    electronJetMatch = true;
                    }
                }
                if (!electronJetMatch) {
                    goodJets8.push_back(*jet);
                }
            }
        }
        //select goodJet12 (jet pT>210GeV, |eta|<2)
        for (int i = 0; i < branchJet12->GetEntries(); ++i) {
            Jet* jet = (Jet *) branchJet12->At(i);
            if (jet->PT > jet12MinPt && fabs(jet->Eta) < jetMaxEta) {
                bool electronJetMatch = false;
                for (size_t ii = 0; ii < goodElectrons.size(); ++ii) {
                    if (goodElectrons[ii].P4().DeltaR(jet->P4()) < electronJetCloseness){
                    electronJetMatch = true;
                    }
                }
                if (!electronJetMatch) {
                    goodJets12.push_back(*jet);
                }
            }
        }

        //at least 1 Jet8
        if(goodJets8.size() >= 1){
            h2->Fill(50);
            std::vector<Jet> btagJet8;
            for(size_t i=0 ; i < goodJets8.size() ; i++){
                if(goodJets8[i].BTag == 1){
                    btagJet8.push_back(goodJets8[i]);
                }
            }
            //std::cout << "  batg8 " << btagJet8.size() ;
            //at least 1 b-tagged jet8
            if(btagJet8.size() >= 1 ){
                h2->Fill(60);
                //if(btagJet8.size() > 1){std::cout << "error8" ;}

                 //boost jet8 mass  in 106GeV 146GeV
                double jet8Mass = btagJet8[0].P4().M();
                if( jet8Mass < boostjetMaxmass && jet8Mass > boostjetMinmass ){
                    h2->Fill(70);
                }

            }

        }

        //at lest 1 jet 12
        if(goodJets12.size() >= 1){
            h2->Fill(80);
            std::vector<Jet> btagJet12;
            for(size_t i=0 ; i < goodJets12.size() ; i++){
                if(goodJets12[i].BTag == 1){
                    btagJet12.push_back(goodJets12[i]);
                }
            }
            //std::cout << "  batg12 " << btagJet12.size();

            //at least 1 b-tagged jet12
            if(btagJet12.size() >= 1){
                h2->Fill(90);
		//if(btagJet12.size() > 1){std::cout<< "error12" ;}
                //boost jet12 mass in 106GeV 146GeV
                double jet12Mass = btagJet12[0].P4().M();
                if( jet12Mass < boostjetMaxmass && jet12Mass > boostjetMinmass ){
                    h2->Fill(100);
                }

            }

        }


    }
	
  //std::cout << std::endl ;

  }

  std::cout << "Flow r " << std::endl;
  int nbins = h1->GetNbinsX();
  for (int i = 1; i <= nbins; ++i) {
    double items = h1->GetBinContent(i);
    if (items > 0.0) {
      std::cout << "n=" << i << "  " <<  items << " ";
    }
  }

  std::cout << std::endl;

  std::cout << "Flow b " << std::endl;
  int nbinsb = h2->GetNbinsX();
  for (int i = 1; i <= nbinsb; ++i) {
    double items = h2->GetBinContent(i);
    if (items > 0.0) {
      std::cout << "n=" << i <<"  " << items << " ";
    }
  }

  std::cout << std::endl;

}
