#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>

#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

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

using namespace std;

int main (int argc, char * argv[])
{
	// command line input
	if (argc != 2)
	{
		fprintf(stderr, "ERROR - wrong number of arguments\n");
		fprintf(stderr, "usage: %s inputFile\n", argv[0]);
		return 1;
	}

	char inputFile[1024];
	sprintf(inputFile, argv[1]);

	// system check if the input file exists
	if (access(inputFile, F_OK) != 0)
	{
		if (errno == ENOENT)
		{
			cout << inputFile << " does not exist" << endl;
		}
		else if (errno == EACCES)
		{
			cout << inputFile << " is not accessible" << endl;
		}
		return 1;
	}

	// create object of class TChain
	TChain * tc = new TChain("Delphes");
	tc->Add(inputFile);

	// create object of class ExRootTreeReader
	ExRootTreeReader * tr = new ExRootTreeReader(tc);
	Long64_t numberofEntries = tr->GetEntries();

	// create histograms to contain cut flow information
	TH1F * h1 = new TH1F("cut flow, r", "cut flow, r", 11, 0, 10);
	TH1F * h2 = new TH1F("cut flow, b", "cut flow, b", 12, 0, 11);

	// get pointer to branches
	TClonesArray * branchElectron = tr->UseBranch("Electron");
	TClonesArray * branchJet4 = tr->UseBranch("Jet4");
	TClonesArray * branchJet8 = tr->UseBranch("Jet6");
	TClonesArray * branchJet12 = tr->UseBranch("Jet10");
	TClonesArray * branchMuon = tr->UseBranch("Muon");

	// some cut parameters
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

	for(Int_t entry = 0; entry < numberofEntries; entry++)
	{
		//Load selected branches
		tr->ReadEntry(entry);

		// add all the events to the 1st bin of both channels
		h1->Fill(1);
		h2->Fill(1);

		bool positive = false;
		bool negative = false;

		std::vector<Electron> goodElectrons;
		std::vector<Muon> goodMuons;
		std::vector<Jet> goodJets8;
		std::vector<Jet> goodJets4;
		std::vector<Jet> goodJets12;

		// cut events which don't have a pair of leptons with pt > minPT, |Eta| < MaxEta
		for (int i = 0; i < branchElectron->GetEntries(); i++)
		{
			Electron * electron = (Electron*) branchElectron->At(i);
			if (electron->PT > electronMinPt && fabs(electron->Eta) < electronMaxEta)
			{
				goodElectrons.push_back(*electron);
				if (electron->Charge == 1)
				{
					positive = true;
				}
				else if (electron->Charge == -1)
				{
					negative = true;
				}
			}
		}

		for (int i = 0; i < branchMuon->GetEntries(); i++)
		{
			Muon * muon = (Muon *) branchMuon->At(i);
			if (muon->PT > muonMinPt && fabs(muon->Eta) < muonMaxEta)
			{
				goodMuons.push_back(*muon);
				if (muon->Charge == 1)
				{
					positive = true;
				}
				else if (muon->Charge == -1)
				{
					negative = true;
				}
			}
		}

		// more cut parameters
		int leptonNumber = 0;
		double dileptonMass = 0;

		if (goodElectrons.size() == 2 && goodMuons.size() == 0)
		{
			dileptonMass = (goodElectrons[0].P4() + goodElectrons[1].P4()).M();
			leptonNumber = 2;
		}
		else if (goodElectrons.size() == 0 && goodMuons.size() == 2)
		{
			leptonNumber = 2;
			dileptonMass = (goodMuons[0].P4() + goodMuons[1].P4()).M();
		}

		if (positive == false || negative == false || leptonNumber != 2)
		{
			continue;
		}

		// if we have exactly 2 leptons, 1 positive and 1 negative, fill the 2nd bin of both channels
		h1->Fill(2);
		h2->Fill(2);

		// if the dilepton mass is in the range minMass < m(ll) < maxMass
		if (dileptonMass < dileptonMinmass || dileptonMass > dileptonMaxmass)
		{
			continue;
		}
		// fill the 3rd bin of both channels
		h1->Fill(3);
		h2->Fill(3);

		// another cut parameter
		float electronJetCloseness = 0.1;

		for (int i = 0; i < branchJet4->GetEntries(); i++)
		{
			Jet * jet = (Jet *) branchJet4->At(i);
			if (jet->PT > jetMinPt && fabs(jet->Eta) < jetMaxEta)
			{
				// electrons are detected as jets too
				// we want to skip these, so we will compare each jet to the previously found electrons
				// and skip it if it is one of them

				bool electronJetMatch = false;
				for (size_t ii = 0; ii < goodElectrons.size(); ii++)
				{
					if (goodElectrons[ii].P4().DeltaR(jet->P4()) < electronJetCloseness)
					{
						electronJetMatch = true;
					}
				}
				if (!electronJetMatch)
				{
					goodJets4.push_back(*jet);
				}
			}
		}

		bool resolved = false;
		// cut events which don't have at least two jet4 (jet pT > 25GeV, |eta| < 2)
		// resolved channel
		if (goodJets4.size() >= 2)
		{
			// add events to 4th bin of resolved the channel
			h1->Fill(4);
			std::vector<Jet> btagJet4;
			for (size_t i = 0; i < goodJets4.size(); i++)
			{
				if (goodJets4[i].BTag == 1)
				{
					btagJet4.push_back(goodJets4[i]);
				}
			}

			// cut events whose jet4's aren't b-tagged
			if (btagJet4.size() >= 2)
			{
				// add these events to the 5th bin of the resolved channel
				h1->Fill(5);
				// cut events whose dijet4mass is not within the range 111GeV < m(jj) < 141GeV
				double dijet4mass = (btagJet4[0].P4() + btagJet4[1].P4()).M();
				if (dijet4mass > dijetMaxmass || dijet4mass < dijetMinmass)
				{
					// add these events to the 6th bin of the resolved channel
					h1->Fill(6);
					resolved = true;	// this seems to be when the channel is determined to be resolved - wouldn't that make bin 4 and 5 of the resolved channel cut flow histogram not actually represent a resolved channel?
				}
			}
		}

		// if the jets were not resolved, start adding events to the merged channel
		// merged channel
		if (resolved == false)
		{
			// add events to the 4th bin of the merged channel
			h2->Fill(4);

			// select goodJet8 (jet pT > 315GeV, |eta| < 2)
			for (int i = 0; i < branchJet8->GetEntries(); i++)
			{
				Jet * jet = (Jet *) branchJet8->At(i);
				if (jet->PT > jet8MinPt && fabs(jet->Eta) < jetMaxEta)
				{
					bool electronJetMatch = false;
					for (size_t ii = 0; ii < goodElectrons.size(); ii++)
					{
						if (goodElectrons[ii].P4().DeltaR(jet->P4()) < electronJetCloseness)
						{
							electronJetMatch = true;
						}
					}
					if (!electronJetMatch)
					{
						goodJets8.push_back(*jet);
					}
				}
			}

			// select goodJet12 (jet pT > 210GeV, |eta| < 2)
			for (int i = 0; i < branchJet12->GetEntries(); i++)
			{
				Jet * jet = (Jet *) branchJet12->At(i);
				if (jet->PT > jet12MinPt && fabs(jet->Eta) < jetMaxEta)
				{
					bool electronJetMatch = false;
					for (size_t ii = 0; ii < goodElectrons.size(); ++ii)
					{
						if (goodElectrons[ii].P4().DeltaR(jet->P4()) < electronJetCloseness)
						{
							electronJetMatch = true;
						}
					}
					if (!electronJetMatch)
					{
						goodJets12.push_back(*jet);
					}
				}
			}

			// at least 1 Jet8
			if(goodJets8.size() >= 1)
			{
				// add these events to the 5th bin of the merged channel
				h2->Fill(5);
				std::vector<Jet> btagJet8;
				for (size_t i = 0; i < goodJets8.size(); i++)
				{
					if (goodJets8[i].BTag == 1)
					{
						btagJet8.push_back(goodJets8[i]);
					}
				}
				// at least 1 b-tagged jet8
				if (btagJet8.size() >= 1)
				{
					// add these events to the 6th bin of the merged channel
					h2->Fill(6);

					// boost jet8 mass in the range 106GeV -> 146GeV
					double jet8Mass = btagJet8[0].P4().M();
					if (jet8Mass < boostjetMaxmass && jet8Mass > boostjetMinmass)
					{
						// add these events to the 7th bin of the merged channel
						h2->Fill(7);
					}
				}
			}

			// at least 1 jet12
			if (goodJets12.size() >= 1)
			{
				// add these events to the 8th bin of the merged channel
				h2->Fill(8);
				std::vector<Jet> btagJet12;
				for (size_t i = 0; i < goodJets12.size(); i++)
				{
					if (goodJets12[i].BTag == 1)
					{
						btagJet12.push_back(goodJets12[i]);
					}
				}

				// at least 1 b-tagged jet12
				if (btagJet12.size() >= 1)
				{
					// add these events to the 9th bin of the merged channel
					h2->Fill(9);
					// boost jet12 mass in the range 106GeV -> 146GeV
					double jet12Mass = btagJet12[0].P4().M();
					if (jet12Mass < boostjetMaxmass && jet12Mass > boostjetMinmass)
					{
						// add these events to the 10th bin of the merged channel
						h2->Fill(10);
					}
				}
			}
		}
	}

	cout << "Flow r " << std::endl;
	int nbins = h1->GetNbinsX();
	for (int i = 1; i <= nbins; i++)
	{
		double items = 0.0;
		items = h1->GetBinContent(i);
		cout << "cut" << i << "\t" <<  items << endl;
	}

	cout << endl;

	cout << "Flow b " << endl;
	int nbinsb = h2->GetNbinsX();

	for (int i = 1; i <= nbinsb; i++)
	{
		double items = 0.0;
		items = h2->GetBinContent(i);
		cout << "cut" << i << "\t" <<  items << endl;
	}

	cout << endl;

	return 0;
}
