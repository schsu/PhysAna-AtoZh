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

	// get pointer to electron, muon, jet4, and particle branches
	TClonesArray * branchElectron = tr->UseBranch("Electron");
	TClonesArray * branchMuon = tr->UseBranch("Muon");
	TClonesArray * branchJet4 = tr->UseBranch("Jet4");
	TClonesArray * branchParticle = tr->UseBranch("Particle");

	// create the histograms
	// cut flow kinematics
	TH1F * h_ll_dielectron_pt = new TH1F("ll dielectron pt", "ll dielectron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_jj_dielectron_pt = new TH1F("jj dielectron pt", "jj dielectron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_bb_dielectron_pt = new TH1F("bb dielectron pt", "bb dielectron pt; pt (GeV, 100 bins); count", 100, 0, 300);

	// electrons channel
	TH1F * h_dielectron_mass = new TH1F("dielectron mass", "dielectron mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_dielectron_pt = new TH1F("dielectron pt", "dielectron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_electron_deltar = new TH1F("electron delta r", "electron delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_ec_dijet4_mass = new TH1F("ec dijet4 mass", "ec dijet4 mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_ec_dijet4_pt = new TH1F("ec dijet4 pt", "ec dijet4 pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_ec_jet4_deltar = new TH1F("ec jet4 delta r", "ec jet4 delta r; delta r (100 bins); count", 100, 0, 5);
	TH2F * h_ec_dijet4_pt_vs_deltar = new TH2F("ec dijet4 pt vs delta r", "ec dijet4 pt vs delta r; dijet4 pt (GeV); delta r", 100, 0, 300, 100, 0, 5);

	TH1F * h_ec_A_mass = new TH1F("ec A mass", "ec A mass; mass (GeV, 100 bins); count", 100, 100, 500);
	TH1F * h_ec_A_pt = new TH1F("ec A pt", "ec A pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_ec_A_rapidity = new TH1F("ec A rapidity", "ec A rapidity; rapitidy (GeV, 100 bins); count", 100, -3.4, 3.4);

	// muon channel
	TH1F * h_dimuon_mass = new TH1F("dimuon mass", "dimuon mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_dimuon_pt = new TH1F("dimuon pt", "dimuon pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_muon_deltar = new TH1F("muon delta r", "muon delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_mc_dijet4_mass = new TH1F("mc dijet4 mass", "mc dijet4 mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mc_dijet4_pt = new TH1F("mc dijet4 pt", "mc dijet4 pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mc_jet4_deltar = new TH1F("mc jet4 delta r", "mc jet4 delta r; delta r (100 bins); count", 100, 0, 5);
	TH2F * h_mc_dijet4_pt_vs_deltar = new TH2F("mc dijet4 pt vs delta r", "mc dijet4 pt vs delta r; dijet4 pt (GeV); delta r", 100, 0, 300, 100, 0, 5);

	TH1F * h_mc_A_mass = new TH1F("mc A mass", "mc A mass; mass (GeV, 100 bins); count", 100, 100, 500);
	TH1F * h_mc_A_pt = new TH1F("mc A pt", "mc A pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mc_A_rapidity = new TH1F("mc A rapidity", "mc A rapidity; rapitidy (GeV, 100 bins); count", 100, -3.4, 3.4);

	// pt/eta cut parameters
	double electronMinPt = 20.0;
	double electronMaxEta = 2.5;
	double muonMinPt = 20.0;
	double muonMaxEta = 2.5;
	double jet4MinPt = 20.0;
	double jet4MaxEta = 2.5;

	// efficiency study for muons
	double es_minpt = 10.0;
	double es_mineta = 1.5;
	double es_maxeta = 2.5;
	int64_t es_muons = 0;
	int64_t es_truth_muons = 0;

	// counters for electron and muon efficiencies
	int64_t truth_electrons = 0;
	int64_t truth_muons = 0;

	int64_t detected_electrons = 0;
	int64_t detected_muons = 0;

	// ll cut - 2 leptons, 1 positive, 1 negative, in pt/eta range
	int64_t ll_electrons = 0;
	int64_t ll_muons = 0;

	// jj cut - 2 jets in pt/eta range
	int64_t jj_electrons = 0;
	int64_t jj_muons = 0;

	// bb cut - 2 b-tagged jets
	int64_t bb_electrons = 0;
	int64_t bb_muons = 0;

	// loop over each event
	for(Int_t entry = 0; entry < numberofEntries; entry++)
	{
		// load selected branches
		tr->ReadEntry(entry);

		// count all detected electrons and muons from the electron and muon branches
		detected_electrons += branchElectron->GetEntries();
		detected_muons += branchMuon->GetEntries();

		// efficiency study for muons
		for (int i = 0; i < branchMuon->GetEntries(); i++)
		{
			Muon * muon = (Muon*) branchMuon->At(i);

			if (muon->PT > es_minpt && fabs(muon->Eta) < es_maxeta && fabs(muon->Eta) > es_mineta)
			{
				es_muons += 1;
			}
		}
		for (int i = 0; i < branchParticle->GetEntries(); i++)
		{
			GenParticle * particle = (GenParticle*) branchParticle->At(i);

			if (particle->Status == 1)
			{
				if (particle->PT > es_minpt && fabs(particle->Eta) < es_maxeta && fabs(particle->Eta) > es_mineta)
				{
					if (particle->PID == 13 || particle->PID == -13)
					{
						es_truth_muons += 1;
					}
				}
			}
		}

		// count truth electrons and muons
		for (int i = 0; i < branchParticle->GetEntries(); i++)
		{
			GenParticle * particle = (GenParticle*) branchParticle->At(i);

			if (particle->Status == 1)
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
		vector<Jet> goodJets4;

		// good b-jets are jets which pass the pt/eta cut and are b-tagged
		vector<Jet> goodBJets4;

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

		// apply the pt/eta cut to jets4
		for (int i = 0; i < branchJet4->GetEntries(); i++)
		{
			Jet * jet4 = (Jet*) branchJet4->At(i);

			if (jet4->PT > jet4MinPt && fabs(jet4->Eta) < jet4MaxEta)
			{
				goodJets4.push_back(*jet4);
			}
		}

		// apply b-tag cut
		for (int i = 0; i < goodJets4.size(); i++)
		{
			if (goodJets4[i].BTag == 1)
			{
				goodBJets4.push_back(goodJets4[i]);
			}
		}

		// if there is at least 1 positive and 1 negative lepton (which passed the pt/eta cut)
		if (havePositiveLepton && haveNegativeLepton)
		{
			// electron channel
			// 2 electrons
			if (goodElectrons.size() == 2 && goodMuons.size() == 0)
			{
				ll_electrons += 2;
				h_ll_dielectron_pt->Fill((goodElectrons[0].P4() + goodElectrons[1].P4()).Pt());

				// 2 jets
				if (goodJets4.size() == 2)
				{
					jj_electrons += 2;
					h_jj_dielectron_pt->Fill((goodElectrons[0].P4() + goodElectrons[1].P4()).Pt());

					// 2 b-jets
					if (goodBJets4.size() == 2)
					{
						bb_electrons += 2;
						h_bb_dielectron_pt->Fill((goodElectrons[0].P4() + goodElectrons[1].P4()).Pt());

						h_dielectron_mass->Fill((goodElectrons[0].P4() + goodElectrons[1].P4()).M());
						h_dielectron_pt->Fill((goodElectrons[0].P4() + goodElectrons[1].P4()).Pt());
						h_electron_deltar->Fill(goodElectrons[0].P4().DeltaR(goodElectrons[1].P4()));

						h_ec_dijet4_mass->Fill((goodBJets4[0].P4() + goodBJets4[1].P4()).M());
						h_ec_dijet4_pt->Fill((goodBJets4[0].P4() + goodBJets4[1].P4()).Pt());
						h_ec_jet4_deltar->Fill(goodBJets4[0].P4().DeltaR(goodBJets4[1].P4()));
						h_ec_dijet4_pt_vs_deltar->Fill((goodBJets4[0].P4() + goodBJets4[1].P4()).Pt(), goodBJets4[0].P4().DeltaR(goodBJets4[1].P4()));

						h_ec_A_mass->Fill((goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).M());
						h_ec_A_pt->Fill((goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Pt());
						h_ec_A_rapidity->Fill((goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Rapidity());
					}
				}
			}

			// muon channel
			// 2 muons
			else if (goodElectrons.size() == 0 && goodMuons.size() == 2)
			{
				ll_muons += 2;

				// 2 jets
				if (goodJets4.size() == 2)
				{
					jj_muons += 2;

					// 2 b-jets
					if (goodBJets4.size() == 2)
					{
						bb_muons += 2;

						h_dimuon_mass->Fill((goodMuons[0].P4() + goodMuons[1].P4()).M());
						h_dimuon_pt->Fill((goodMuons[0].P4() + goodMuons[1].P4()).Pt());
						h_muon_deltar->Fill(goodMuons[0].P4().DeltaR(goodMuons[1].P4()));

						h_mc_dijet4_mass->Fill((goodBJets4[0].P4() + goodBJets4[1].P4()).M());
						h_mc_dijet4_pt->Fill((goodBJets4[0].P4() + goodBJets4[1].P4()).Pt());
						h_mc_jet4_deltar->Fill(goodBJets4[0].P4().DeltaR(goodBJets4[1].P4()));
						h_mc_dijet4_pt_vs_deltar->Fill((goodBJets4[0].P4() + goodBJets4[1].P4()).Pt(), goodBJets4[0].P4().DeltaR(goodBJets4[1].P4()));

						h_mc_A_mass->Fill((goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).M());
						h_mc_A_pt->Fill((goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Pt());
						h_mc_A_rapidity->Fill((goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Rapidity());
					}
				}
			}
		}
	}

	// make the directories for the plots
	mkdir("cut_plots", 0777);
	mkdir("cut_plots/electron_channel", 0777);
	mkdir("cut_plots/muon_channel", 0777);
	mkdir("cut_plots/electron_channel/cut_flow_kinematics", 0777);

	// output the histograms
	TCanvas * c1 = new TCanvas("c1", "c1", 640, 480);

	// electron channel
	// mass
	h_dielectron_mass->Draw();
	c1->SaveAs("cut_plots/electron_channel/dielectron_mass.eps");

	h_ec_dijet4_mass->Draw();
	c1->SaveAs("cut_plots/electron_channel/ec_dijet4_mass.eps");

	h_ec_A_mass->Draw();
	c1->SaveAs("cut_plots/electron_channel/ec_A_mass.eps");

	// pt
	h_dielectron_pt->Draw();
	c1->SaveAs("cut_plots/electron_channel/dielectron_pt.eps");

	h_ec_dijet4_pt->Draw();
	c1->SaveAs("cut_plots/electron_channel/ec_dijet4_pt.eps");

	h_ec_A_pt->Draw();
	c1->SaveAs("cut_plots/electron_channel/ec_A_pt.eps");

	// rapidity
	h_ec_A_rapidity->Draw();
	c1->SaveAs("cut_plots/electron_channel/ec_A_rapidity.eps");

	// deltar
	h_electron_deltar->Draw();
	c1->SaveAs("cut_plots/electron_channel/electron_deltar.eps");

	h_ec_jet4_deltar->Draw();
	c1->SaveAs("cut_plots/electron_channel/ec_jet4_deltar.eps");

	h_ec_dijet4_pt_vs_deltar->Draw();
	c1->SaveAs("cut_plots/electron_channel/ec_dijet4_pt_vs_deltar.eps");

	// muon channel
	// mass
	h_dimuon_mass->Draw();
	c1->SaveAs("cut_plots/muon_channel/dimuon_mass.eps");

	h_mc_dijet4_mass->Draw();
	c1->SaveAs("cut_plots/muon_channel/mc_dijet4_mass.eps");

	h_mc_A_mass->Draw();
	c1->SaveAs("cut_plots/muon_channel/mc_A_mass.eps");

	// pt
	h_dimuon_pt->Draw();
	c1->SaveAs("cut_plots/muon_channel/dimuon_pt.eps");

	h_mc_dijet4_pt->Draw();
	c1->SaveAs("cut_plots/muon_channel/mc_dijet4_pt.eps");

	h_mc_A_pt->Draw();
	c1->SaveAs("cut_plots/muon_channel/mc_A_pt.eps");

	// rapidity
	h_mc_A_rapidity->Draw();
	c1->SaveAs("cut_plots/muon_channel/mc_A_rapidity.eps");

	// deltar
	h_muon_deltar->Draw();
	c1->SaveAs("cut_plots/muon_channel/muon_deltar.eps");

	h_mc_jet4_deltar->Draw();
	c1->SaveAs("cut_plots/muon_channel/mc_jet4_deltar.eps");

	h_mc_dijet4_pt_vs_deltar->Draw();
	c1->SaveAs("cut_plots/muon_channel/mc_dijet4_pt_vs_deltar.eps");

	// cut flow kinematics
	h_ll_dielectron_pt->Draw();
	c1->SaveAs("cut_plots/electron_channel/cut_flow_kinematics/ll_dielectron_pt.eps");
	h_jj_dielectron_pt->Draw();
	c1->SaveAs("cut_plots/electron_channel/cut_flow_kinematics/jj_dielectron_pt.eps");
	h_bb_dielectron_pt->Draw();
	c1->SaveAs("cut_plots/electron_channel/cut_flow_kinematics/bb_dielectron_pt.eps");

	// output efficiency information
	cout << endl;
	cout << "detected_electron_efficiency" << "\t" << (float) (detected_electrons) / (float) (truth_electrons) << endl;
	cout << "ll_electron_efficiency" << "\t" << (float) (ll_electrons) / (float) (truth_electrons) << endl;
	cout << "jj_electron_efficiency" << "\t" << (float) (jj_electrons) / (float) (truth_electrons) << endl;
	cout << "bb_electron_efficiency" << "\t" << (float) (bb_electrons) / (float) (truth_electrons) << endl;
	cout << endl;
	cout << "detected_muon_efficiency" << "\t" << (float) (detected_muons) / (float) (truth_muons) << endl;
	cout << "ll_muon_efficiency" << "\t\t" << (float) (ll_muons) / (float) (truth_muons) << endl;
	cout << "jj_muon_efficiency" << "\t\t" << (float) (jj_muons) / (float) (truth_muons) << endl;
	cout << "bb_muon_efficiency" << "\t\t" << (float) (bb_muons) / (float) (truth_muons) << endl;
	cout << endl;

	// efficiency study for muons
	cout << "detected muons in range (abs(eta) > 1.5 && abs(eta) <= 2.7) * (pt > 10.0)" << "\t\t\t" << es_muons << endl;
	cout << "truth muons in range (abs(eta) > 1.5 && abs(eta) <= 2.7) * (pt > 10.0)" << "\t\t\t\t" << es_truth_muons << endl;
	cout << "detection efficiency for muons in range (abs(eta) > 1.5 && abs(eta) <= 2.7) * (pt > 10.0)" << "\t" << (float) es_muons / (float) es_truth_muons << endl;
	cout << endl;

	return 0;
}
