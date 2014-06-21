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
//	TClonesArray * branchParticle = tr->UseBranch("Particle");

/*
	// create the histograms
	// electron channel

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
*/

	// pt/eta cut parameters
	double electronMinPt = 5.0;
	double electronMaxEta = 2.5;
	double muonMinPt = 5.0;
	double muonMaxEta = 2.5;
	double jet4MinPt = 25.0;
	double jet4MaxEta = 2.0;

	// mass cut parameters
	double llMinMass = 80;
	double llMaxMass = 100;
	double bbMinMass = 111;
	double bbMaxMass = 141;

	// counters for event efficiencies
	int total_events = 0;

	int ee_events = 0;
	int mee_events = 0;
	int eejj_events = 0;
	int meejj_events = 0;
	int eebb_events = 0;
	int meebb_events = 0;
	int meembb_events = 0;

	int mumu_events = 0;
	int mmumu_events = 0;
	int mumujj_events = 0;
	int mmumujj_events = 0;
	int mumubb_events = 0;
	int mmumubb_events = 0;
	int mmumumbb_events = 0;

	// loop over each event
	for (Int_t entry = 0; entry < numberofEntries; entry++)
	{
		// load selected branches
		tr->ReadEntry(entry);

		// add to the total number of events
		total_events += 1;

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
			// ee cut
			if (goodElectrons.size() == 2 && goodMuons.size() == 0)
			{
				ee_events += 1;

				double mee = (goodElectrons[0].P4() + goodElectrons[1].P4()).M();
//				double ptee = (goodElectrons[0].P4() + goodElectrons[1].P4()).Pt();
//				double dree = goodElectrons[0].P4().DeltaR(goodElectrons[1].P4());

				// mee cut
				if (mee < llMaxMass && mee > llMinMass)
				{
					mee_events += 1;
				}

				// eejj cut
				if (goodJets4.size() >= 2)
				{
					eejj_events += 1;

					// meejj cut
					if (mee < llMaxMass && mee > llMinMass)
					{
						meejj_events += 1;
					}

					// eebb cut
					if (goodBJets4.size() >= 2)
					{
						eebb_events += 1;

						// meebb cut
						if (mee < llMaxMass && mee > llMinMass)
						{
							meebb_events += 1;
						}

						double mbb = (goodBJets4[0].P4() + goodBJets4[1].P4()).M();
//						double ptbb = (goodBJets4[0].P4() + goodBJets4[1].P4()).Pt();
//						double drbb = goodBJets4[0].P4().DeltaR(goodBJets4[1].P4());

//						double mA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).M();
//						double ptA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Pt();
//						double rA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Rapidity();

						// meembb cut
						if (mee < llMaxMass && mee > llMinMass && mbb < bbMaxMass && mbb > bbMinMass)
						{
							meembb_events += 1;

//							h_dielectron_mass->Fill(mee);
//							h_dielectron_pt->Fill(ptee);
//							h_electron_deltar->Fill(dree);

//							h_ec_dijet4_mass->Fill(mbb);
//							h_ec_dijet4_pt->Fill(ptbb);
//							h_ec_jet4_deltar->Fill(drbb);
//							h_ec_dijet4_pt_vs_deltar->Fill(ptbb, drbb);

//							h_ec_A_mass->Fill(mA);
//							h_ec_A_pt->Fill(ptA);
//							h_ec_A_rapidity->Fill(rA);
						}
					}
				}
			}

			// muon channel
			// mumu cut
			else if (goodMuons.size() == 2 && goodElectrons.size() == 0)
			{
				mumu_events += 1;

				double mmumu = (goodMuons[0].P4() + goodMuons[1].P4()).M();
//				double ptmumu = (goodMuons[0].P4() + goodMuons[1].P4()).Pt();
//				double drmumu = goodMuons[0].P4().DeltaR(goodMuons[1].P4());

				// mmumu cut
				if (mmumu < llMaxMass && mmumu > llMinMass)
				{
					mmumu_events += 1;
				}

				// mumujj cut
				if (goodJets4.size() >= 2)
				{
					mumujj_events += 1;

					// mmumujj cut
					if (mmumu < llMaxMass && mmumu > llMinMass)
					{
						mmumujj_events += 1;
					}

					// mumubb cut
					if (goodBJets4.size() >= 2)
					{
						mumubb_events += 1;

						// mmumubb cut
						if (mmumu < llMaxMass && mmumu > llMinMass)
						{
							mmumubb_events += 1;
						}

						double mbb = (goodBJets4[0].P4() + goodBJets4[1].P4()).M();

//						double ptbb = (goodBJets4[0].P4() + goodBJets4[1].P4()).Pt();
//						double drbb = goodBJets4[0].P4().DeltaR(goodBJets4[1].P4());

//						double mA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).M();
//						double ptA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Pt();
//						double rA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Rapidity();

						// mmumumbb cut
						if (mmumu < llMaxMass && mmumu > llMinMass && mbb < bbMaxMass && mbb > bbMinMass)
						{
							mmumumbb_events += 1;

//							h_dimuon_mass->Fill(mmumu);
//							h_dimuon_pt->Fill(ptmumu);
//							h_muon_deltar->Fill(drmumu);

//							h_mc_dijet4_mass->Fill(mbb);
//							h_mc_dijet4_pt->Fill(ptbb);
//							h_mc_jet4_deltar->Fill(drbb);
//							h_mc_dijet4_pt_vs_deltar->Fill(ptbb, drbb);

//							h_mc_A_mass->Fill(mA);
//							h_mc_A_pt->Fill(ptA);
//							h_mc_A_rapidity->Fill(rA);
						}
					}
				}
			}
		}
	}

/*
	// make the directories for the plots
	mkdir("cut_plots", 0777);
	mkdir("cut_plots/electron_channel", 0777);
	mkdir("cut_plots/muon_channel", 0777);

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
*/

	// output event efficiency information
	cout << endl;
	cout << "total_events" << "\t\t" << total_events << endl;
	cout << "ee_events" << "\t\t" << ee_events << endl;
	cout << "mee_events" << "\t\t" << mee_events << endl;
	cout << "eejj_events" << "\t\t" << eejj_events << endl;
	cout << "meejj_events" << "\t\t" << meejj_events << endl;
	cout << "eebb_events" << "\t\t" << eebb_events << endl;
	cout << "meebb_events" << "\t\t" << meebb_events << endl;
	cout << "meembb_events" << "\t\t" << meembb_events << endl;
	cout << endl;
	cout << "ee_event_efficiency" << "\t" << (float) ee_events / (float) total_events << endl;
	cout << "mee_event_efficiency" << "\t" << (float) mee_events / (float) total_events << endl;
	cout << "eejj_event_efficiency" << "\t" << (float) eejj_events / (float) total_events << endl;
	cout << "meejj_events_efficiency" << "\t" << (float) meejj_events / (float) total_events << endl;
	cout << "eebb_event_efficiency" << "\t" << (float) eebb_events / (float) total_events << endl;
	cout << "meebb_event_efficiecy" << "\t" << (float) meebb_events / (float) total_events << endl;
	cout << "meembb_event_efficiency" << "\t" << (float) meembb_events / (float) total_events << endl;
	cout << endl;
	cout << endl;
	cout << "total_events" << "\t\t" << total_events << endl;
	cout << "mumu_events" << "\t\t" << mumu_events << endl;
	cout << "mmumu_events" << "\t\t" << mmumu_events << endl;
	cout << "mumujj_events" << "\t\t" << mumujj_events << endl;
	cout << "mmumujj_events" << "\t\t" << mmumujj_events << endl;
	cout << "mumubb_events" << "\t\t" << mumubb_events << endl;
	cout << "mmumubb_events" << "\t\t" << mmumubb_events << endl;
	cout << "mmumumbb_events" << "\t\t" << mmumumbb_events << endl;
	cout << endl;
	cout << "mumu_event_efficiency" << "\t\t" << (float) mumu_events / (float) total_events << endl;
	cout << "mmumu_event_efficiency" << "\t\t" << (float) mmumu_events / (float) total_events << endl;
	cout << "mumujj_event_efficiency" << "\t\t" << (float) mumujj_events / (float) total_events << endl;
	cout << "mmumujj_events_efficiency" << "\t" << (float) mmumujj_events / (float) total_events << endl;
	cout << "mumubb_event_efficiency" << "\t\t" << (float) mumubb_events / (float) total_events << endl;
	cout << "mmumubb_event_efficiecy" << "\t\t" << (float) mmumubb_events / (float) total_events << endl;
	cout << "mmumumbb_event_efficiency" << "\t" << (float) mmumumbb_events / (float) total_events << endl;
	cout << endl;
	int ll_events = ee_events + mumu_events;
	int mll_events = mee_events + mmumu_events;
	int lljj_events = eejj_events + mumujj_events;
	int mlljj_events = meejj_events + mmumujj_events;
	int llbb_events = eebb_events + mumubb_events;
	int mllbb_events = meebb_events + mmumubb_events;
	int mllmbb_events = meembb_events + mmumumbb_events;
	cout << endl;
	cout << "total_events" << "\t\t" << total_events << endl;
	cout << "ll_events" << "\t\t" << ll_events << endl;
	cout << "mll_events" << "\t\t" << mll_events << endl;
	cout << "lljj_events" << "\t\t" << lljj_events << endl;
	cout << "mlljj_events" << "\t\t" << mlljj_events << endl;
	cout << "llbb_events" << "\t\t" << llbb_events << endl;
	cout << "mllbb_events" << "\t\t" << mllbb_events << endl;
	cout << "mllmbb_events" << "\t\t" << mllmbb_events << endl;
	cout << endl;
	cout << "ll_event_efficiency" << "\t" << (float) ll_events / (float) total_events << endl;
	cout << "mll_event_efficiency" << "\t" << (float) mll_events / (float) total_events << endl;
	cout << "lljj_event_efficiency" << "\t" << (float) lljj_events / (float) total_events << endl;
	cout << "mlljj_events_efficiency" << "\t" << (float) mlljj_events / (float) total_events << endl;
	cout << "llbb_event_efficiency" << "\t" << (float) llbb_events / (float) total_events << endl;
	cout << "mllbb_event_efficiecy" << "\t" << (float) mllbb_events / (float) total_events << endl;
	cout << "mllmbb_event_efficiency" << "\t" << (float) mllmbb_events / (float) total_events << endl;
	cout << endl;

	return 0;
}
