// C++ headers
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>

// C headers
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

// ROOT headers
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TImage.h"

// Delphes headers (requires that Delphes' library path be added to your LD_LIBRARY_PATH)
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

	// get pointer to electron, muon, jet, and particle branches
	TClonesArray * branchElectron = tr->UseBranch("Electron");
	TClonesArray * branchMuon = tr->UseBranch("Muon");
	TClonesArray * branchJet = tr->UseBranch("Jet");
	TClonesArray * branchParticle = tr->UseBranch("Particle");

	// create the histograms
	TH1F * h_dielectron_mass = new TH1F("dielectron mass", "dielectron mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_dielectron_pt = new TH1F("dielectron pt", "dielectron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_electron_deltar = new TH1F("electron delta r", "electron delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_dimuon_mass = new TH1F("dimuon mass", "dimuon mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_dimuon_pt = new TH1F("dimuon pt", "dimuon pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_muon_deltar = new TH1F("muon delta r", "muon delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_dijet_mass = new TH1F("dijet mass", "dijet mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_dijet_pt = new TH1F("dijet pt", "dijet pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_jet_deltar = new TH1F("jet delta r", "jet delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_A_mass = new TH1F("A mass", "A mass; mass (GeV, 100 bins); count", 100, 100, 500);
	TH1F * h_A_pt = new TH1F("A pt", "A pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_A_rapidity = new TH1F("A rapidity", "A rapidity; rapitidy (GeV, 100 bins); count", 100, -3.4, 3.4);

	// pt/eta cut parameters
	double electronMinPt = 20.0;
	double electronMaxEta = 2.5;
	double muonMinPt = 20.0;
	double muonMaxEta = 2.5;
	double jetMinPt = 20.0;
	double jetMaxEta = 2.5;

	// counters for electron and muon efficiencies
	int64_t truth_electrons = 0;
	int64_t truth_muons = 0;

	int64_t detected_electrons = 0;
	int64_t detected_muons = 0;

	// cut 1 - 2 leptons, 1 positive, 1 negative, in pt/eta range
	int64_t cut1_electrons = 0;
	int64_t cut1_muons = 0;

	// cut 2 - 2 jets in pt/eta range
	int64_t cut2_electrons = 0;
	int64_t cut2_muons = 0;

	// cut 3 - 2 b-tagged jets
	int64_t cut3_electrons = 0;
	int64_t cut3_muons = 0;

	// loop over each event
	for(Int_t entry = 0; entry < numberofEntries; entry++)
	{
		// load selected branches
		tr->ReadEntry(entry);

		// count all detected electrons and muons from the electron and muon branches
		detected_electrons += branchElectron->GetEntries();
		detected_muons += branchMuon->GetEntries();

		// count truth electrons and muons
		for (int i = 0; i < branchParticle->GetEntries(); i++)
		{
			GenParticle * particle = (GenParticle*) branchParticle->At(i);

			if (particle->Status == 3)
			{
				if (particle->PID == 11 || particle->PID == -11)
				{
					truth_electrons += 1;
				}
				if (particle->PID == 13 || particle->PID == -13)
				{
					truth_muons += 1;
				}
			}
		}

		// good electrons, muons, and jets are electrons, muons, and jets which pass the pt/eta cut
		vector<Electron> goodElectrons;
		vector<Muon> goodMuons;
		vector<Jet> goodJets;

		// good b-jets are jets which pass the pt/eta cut and are b-tagged
		vector<Jet> goodBJets;

		// bools to keep track of whether or not both a positive and a negative lepton pass the pt/eta cut
		bool havePositiveLepton = false;
		bool haveNegativeLepton = false;

		// apply the pt/eta cut to electrons
		for (int i = 0; i < branchElectron->GetEntries(); i++)
		{
			Electron * electron = (Electron*) branchElectron->At(i);

			if (electron->PT > electronMinPt && fabs(electron->Eta) < electronMaxEta)
			{
				goodElectrons.push_back(*electron);

				if (electron->Charge == 1)
				{
					havePositiveLepton = true;
				}
				else if (electron->Charge == -1)
				{
					haveNegativeLepton = true;
				}
			}
		}

		// apply the pt/eta cut to muons
		for (int i = 0; i < branchMuon->GetEntries(); i++)
		{
			Muon * muon = (Muon*) branchMuon->At(i);

			if (muon->PT > muonMinPt && fabs(muon->Eta) < muonMaxEta)
			{
				goodMuons.push_back(*muon);

				if (muon->Charge == 1)
				{
					havePositiveLepton = true;
				}
				else if (muon->Charge == -1)
				{
					haveNegativeLepton = true;
				}
			}
		}

		// apply the pt/eta cut to jets
		for (int i = 0; i < branchJet->GetEntries(); i++)
		{
			Jet * jet = (Jet*) branchJet->At(i);

			if (jet->PT > jetMinPt && fabs(jet->Eta) < jetMaxEta)
			{
				goodJets.push_back(*jet);
			}
		}

		// apply b-tag cut
		for (int i = 0; i < goodJets.size(); i++)
		{
			if (goodJets[i].BTag == 1)
			{
				goodBJets.push_back(goodJets[i]);
			}
		}

		// if there is exactly 1 positive and 1 negative lepton (which passed the pt/eta cut)
		if (havePositiveLepton && haveNegativeLepton)
		{
			// electron channel
			// 2 electrons
			if (goodElectrons.size() == 2 && goodMuons.size() == 0)
			{
				cut1_electrons += 2;

				// 2 jets
				if (goodJets.size() == 2)
				{
					cut2_electrons += 2;

					// 2 b-jets
					if (goodBJets.size() == 2)
					{
						cut3_electrons += 2;

						h_dielectron_mass->Fill((goodElectrons[0].P4() + goodElectrons[1].P4()).M());
						h_dielectron_pt->Fill((goodElectrons[0].P4() + goodElectrons[1].P4()).Pt());
						h_electron_deltar->Fill(goodElectrons[0].P4().DeltaR(goodElectrons[1].P4()));

						h_dijet_mass->Fill((goodBJets[0].P4() + goodBJets[1].P4()).M());
						h_dijet_pt->Fill((goodBJets[0].P4() + goodBJets[1].P4()).Pt());
						h_jet_deltar->Fill(goodBJets[0].P4().DeltaR(goodBJets[1].P4()));

						h_A_mass->Fill((goodBJets[0].P4() + goodBJets[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).M());
						h_A_pt->Fill((goodBJets[0].P4() + goodBJets[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Pt());
						h_A_rapidity->Fill((goodBJets[0].P4() + goodBJets[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Rapidity());
					}
				}
			}

			// muon channel
			// 2 muons
			else if (goodElectrons.size() == 0 && goodMuons.size() == 2)
			{
				cut1_muons += 2;

				// 2 jets
				if (goodJets.size() == 2)
				{
					cut2_muons += 2;

					// 2 b-jets
					if (goodBJets.size() == 2)
					{
						cut3_muons += 2;

						h_dimuon_mass->Fill((goodMuons[0].P4() + goodMuons[1].P4()).M());
						h_dimuon_pt->Fill((goodMuons[0].P4() + goodMuons[1].P4()).Pt());
						h_muon_deltar->Fill(goodMuons[0].P4().DeltaR(goodMuons[1].P4()));

						h_dijet_mass->Fill((goodBJets[0].P4() + goodBJets[1].P4()).M());
						h_dijet_pt->Fill((goodBJets[0].P4() + goodBJets[1].P4()).Pt());
						h_jet_deltar->Fill(goodBJets[0].P4().DeltaR(goodBJets[1].P4()));

						h_A_mass->Fill((goodBJets[0].P4() + goodBJets[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).M());
						h_A_pt->Fill((goodBJets[0].P4() + goodBJets[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Pt());
						h_A_rapidity->Fill((goodBJets[0].P4() + goodBJets[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Rapidity());
					}
				}
			}
		}
	}

	// make the directory for the plots
	mkdir("cut_plots", 0777);

	// output the histograms
	TCanvas * c1 = new TCanvas("c1", "c1", 640, 480);

	// mass
	h_dielectron_mass->Draw();
	c1->SaveAs("cut_plots/dielectron_mass.eps");

	h_dimuon_mass->Draw();
	c1->SaveAs("cut_plots/dimuon_mass.eps");

	h_dijet_mass->Draw();
	c1->SaveAs("cut_plots/dijet_mass.eps");

	// pt
	h_dielectron_pt->Draw();
	c1->SaveAs("cut_plots/dielectron_pt.eps");

	h_dimuon_pt->Draw();
	c1->SaveAs("cut_plots/dimuon_pt.eps");

	h_dijet_pt->Draw();
	c1->SaveAs("cut_plots/dijet_pt.eps");

	// A
	h_A_mass->Draw();
	c1->SaveAs("cut_plots/A_mass.eps");

	h_A_pt->Draw();
	c1->SaveAs("cut_plots/A_pt.eps");

	h_A_rapidity->Draw();
	c1->SaveAs("cut_plots/A_rapidity.eps");

	// deltar
	h_electron_deltar->Draw();
	c1->SaveAs("cut_plots/electron_deltar.eps");

	h_muon_deltar->Draw();
	c1->SaveAs("cut_plots/muon_deltar.eps");

	h_jet_deltar->Draw();
	c1->SaveAs("jet_deltar.eps");

	// output efficiency information
	cout << endl;
	cout << "detected_electron_efficiency" << "\t" << (float) (detected_electrons) / (float) (truth_electrons) << endl;
	cout << "cut1_electron_efficiency" << "\t" << (float) (cut1_electrons) / (float) (truth_electrons) << endl;
	cout << "cut2_electron_efficiency" << "\t" << (float) (cut2_electrons) / (float) (truth_electrons) << endl;
	cout << "cut3_electron_efficiency" << "\t" << (float) (cut3_electrons) / (float) (truth_electrons) << endl;
	cout << endl;
	cout << "detected_muon_efficiency" << "\t" << (float) (detected_muons) / (float) (truth_muons) << endl;
	cout << "cut1_muon_efficiency" << "\t\t" << (float) (cut1_muons) / (float) (truth_muons) << endl;
	cout << "cut2_muon_efficiency" << "\t\t" << (float) (cut2_muons) / (float) (truth_muons) << endl;
	cout << "cut3_muon_efficiency" << "\t\t" << (float) (cut3_muons) / (float) (truth_muons) << endl;
	cout << endl;

	return 0;
}
