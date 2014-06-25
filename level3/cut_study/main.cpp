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
	TClonesArray * branchJet8 = tr->UseBranch("Jet8");

	// create the histograms
	// electron resolved channel
	TH1F * h_res_dielectron_mass = new TH1F("res dielectron mass", "res dielectron mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_res_dielectron_pt = new TH1F("res dielectron pt", "res dielectron pt; pt (GeV, 100 bins); count", 100, 0, 3000);
	TH1F * h_res_electron_deltar = new TH1F("res electron delta r", "res electron delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_ec_res_dijet4_mass = new TH1F("ec res dijet4 mass", "ec res dijet4 mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_ec_res_dijet4_pt = new TH1F("ec res dijet4 pt", "ec res dijet4 pt; pt (GeV, 100 bins); count", 100, 0, 3000);
	TH1F * h_ec_res_jet4_deltar = new TH1F("ec res jet4 delta r", "ec res jet4 delta r; delta r (100 bins); count", 100, 0, 5);
	TH2F * h_ec_res_dijet4_pt_vs_deltar = new TH2F("ec res dijet4 pt vs delta r", "ec res dijet4 pt vs delta r; dijet4 pt (GeV); delta r", 100, 0, 300, 100, 0, 5);

	TH1F * h_ec_res_A_mass = new TH1F("ec res A mass", "ec res A mass; mass (GeV, 100 bins); count", 100, 100, 3000);
	TH1F * h_ec_res_A_pt = new TH1F("ec res A pt", "ec res A pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_ec_res_A_rapidity = new TH1F("ec res A rapidity", "ec res A rapidity; rapitidy (GeV, 100 bins); count", 100, -3.4, 3.4);

	// electron merged channel
	TH1F * h_mer_dielectron_mass = new TH1F("mer dielectron mass", "mer dielectron mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mer_dielectron_pt = new TH1F("mer dielectron pt", "mer dielectron pt; pt (GeV, 100 bins); count", 100, 0, 3000);
	TH1F * h_mer_electron_deltar = new TH1F("mer electron delta r", "mer electron delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_ec_mer_Jet8_mass = new TH1F("ec mer Jet8 mass", "ec mer Jet8 mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_ec_mer_Jet8_pt = new TH1F("ec mer Jet8 pt", "ec mer Jet8 pt; pt (GeV, 100 bins); count", 100, 0, 3000);

	TH1F * h_ec_mer_A_mass = new TH1F("ec mer A mass", "ec mer A mass; mass (GeV, 100 bins); count", 100, 100, 3000);
	TH1F * h_ec_mer_A_pt = new TH1F("ec mer A pt", "ec mer A pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_ec_mer_A_rapidity = new TH1F("ec mer A rapidity", "ec mer A rapidity; rapitidy (GeV, 100 bins); count", 100, -3.4, 3.4);

	// muon resolved channel
	TH1F * h_res_dimuon_mass = new TH1F("res dimuon mass", "res dimuon mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_res_dimuon_pt = new TH1F("res dimuon pt", "res dimuon pt; pt (GeV, 100 bins); count", 100, 0, 3000);
	TH1F * h_res_muon_deltar = new TH1F("res muon delta r", "res muon delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_mc_res_dijet4_mass = new TH1F("mc res dijet4 mass", "mc res dijet4 mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mc_res_dijet4_pt = new TH1F("mc res dijet4 pt", "mc res dijet4 pt; pt (GeV, 100 bins); count", 100, 0, 3000);
	TH1F * h_mc_res_jet4_deltar = new TH1F("mc res jet4 delta r", "mc res jet4 delta r; delta r (100 bins); count", 100, 0, 5);
	TH2F * h_mc_res_dijet4_pt_vs_deltar = new TH2F("mc res dijet4 pt vs delta r", "mc res dijet4 pt vs delta r; dijet4 pt (GeV); delta r", 100, 0, 300, 100, 0, 5);

	TH1F * h_mc_res_A_mass = new TH1F("mc res A mass", "mc res A mass; mass (GeV, 100 bins); count", 100, 100, 3000);
	TH1F * h_mc_res_A_pt = new TH1F("mc res A pt", "mc res A pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mc_res_A_rapidity = new TH1F("mc res A rapidity", "mc res A rapidity; rapitidy (GeV, 100 bins); count", 100, -3.4, 3.4);

	// muon merged channel
	TH1F * h_mer_dimuon_mass = new TH1F("mer dimuon mass", "mer dimuon mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mer_dimuon_pt = new TH1F("mer dimuon pt", "mer dimuon pt; pt (GeV, 100 bins); count", 100, 0, 3000);
	TH1F * h_mer_muon_deltar = new TH1F("mer muon delta r", "mer muon delta r; delta r (100 bins); count", 100, 0, 5);

	TH1F * h_mc_mer_Jet8_mass = new TH1F("mc mer Jet8 mass", "mc mer Jet8 mass; mass (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mc_mer_Jet8_pt = new TH1F("mc mer Jet8 pt", "mc mer Jet8 pt; pt (GeV, 100 bins); count", 100, 0, 3000);

	TH1F * h_mc_mer_A_mass = new TH1F("mc mer A mass", "mc mer A mass; mass (GeV, 100 bins); count", 100, 100, 3000);
	TH1F * h_mc_mer_A_pt = new TH1F("mc mer A pt", "mc mer A pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_mc_mer_A_rapidity = new TH1F("mc mer A rapidity", "mc mer A rapidity; rapitidy (GeV, 100 bins); count", 100, -3.4, 3.4);

	// pt/eta cut parameters
	double electronMinPt = 20.0;
	double electronMaxEta = 2.5;
	double muonMinPt = 20.0;
	double muonMaxEta = 2.5;
	double jet4MinPt = 30.0;
	double jet4MaxEta = 2.5;
	double jet8MinPt = 350.0;
	double jet8MaxEta = 2.5;

	// mass cut parameters
	double llMinMass = 0;
	double llMaxMass = 10000;
	double bbMinMass = 0;
	double bbMaxMass = 10000;
	double BMinMass = 0;
	double BMaxMass = 10000;

	// counters for event efficiencies
	int total_events = 0;

	int ee_events = 0;
	int total_truth_ee_events = 0;
	int truth_ee_events = 0;
	int mee_events = 0;
	int eejj_events = 0;
	int meejj_events = 0;
	int eebb_events = 0;
	int meebb_events = 0;
	int meembb_events = 0;

	int eeJ_events = 0;
	int meeJ_events = 0;
	int eeB_events = 0;
	int meeB_events = 0;
	int meemB_events = 0;

	int mumu_events = 0;
	int total_truth_mumu_events = 0;
	int truth_mumu_events = 0;
	int mmumu_events = 0;
	int mumujj_events = 0;
	int mmumujj_events = 0;
	int mumubb_events = 0;
	int mmumubb_events = 0;
	int mmumumbb_events = 0;

	int mumuJ_events = 0;
	int mmumuJ_events = 0;
	int mumuB_events = 0;
	int mmumuB_events = 0;
	int mmumumB_events = 0;

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
		vector<Jet> goodJets8;

		// good b-jets are jets which pass the pt/eta cut and are b-tagged
		vector<Jet> goodBJets4;
		vector<Jet> goodBJets8;

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

		// apply b-tag cut to jets4
		for (int i = 0; i < goodJets4.size(); i++)
		{
			if (goodJets4[i].BTag == 1)
			{
				goodBJets4.push_back(goodJets4[i]);
			}
		}

		// apply the pt/eta cut to jets8
		for (int i = 0; i < branchJet8->GetEntries(); i++)
		{
			Jet * jet8 = (Jet*) branchJet8->At(i);

			if (jet8->PT > jet8MinPt && fabs(jet8->Eta) < jet8MaxEta)
			{
				goodJets8.push_back(*jet8);
			}
		}

		// apply b-tag cut to jets8
		for (int i = 0; i < goodJets8.size(); i++)
		{
			if (goodJets8[i].BTag == 1)
			{
				goodBJets8.push_back(goodJets8[i]);
			}
		}

		// if there is at least 1 positive and 1 negative lepton (which passed the pt/eta cut)
		if (havePositiveLepton && haveNegativeLepton)
		{
			// electron channel
			// ee cut
			if (goodElectrons.size() >= 2)// && goodMuons.size() == 0)
			{

				ee_events += 1;

				double mee = (goodElectrons[0].P4() + goodElectrons[1].P4()).M();
				double ptee = (goodElectrons[0].P4() + goodElectrons[1].P4()).Pt();
				double dree = goodElectrons[0].P4().DeltaR(goodElectrons[1].P4());

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
						double ptbb = (goodBJets4[0].P4() + goodBJets4[1].P4()).Pt();
						double drbb = goodBJets4[0].P4().DeltaR(goodBJets4[1].P4());

						double mA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).M();
						double ptA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Pt();
						double rA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Rapidity();

						// meembb cut
						if (mee < llMaxMass && mee > llMinMass && mbb < bbMaxMass && mbb > bbMinMass)
						{
							meembb_events += 1;

							h_res_dielectron_mass->Fill(mee);
							h_res_dielectron_pt->Fill(ptee);
							h_res_electron_deltar->Fill(dree);

							h_ec_res_dijet4_mass->Fill(mbb);
							h_ec_res_dijet4_pt->Fill(ptbb);
							h_ec_res_jet4_deltar->Fill(drbb);
							h_ec_res_dijet4_pt_vs_deltar->Fill(ptbb, drbb);

							h_ec_res_A_mass->Fill(mA);
							h_ec_res_A_pt->Fill(ptA);
							h_ec_res_A_rapidity->Fill(rA);
						}
					}
				}
				// eeJ cut
				if (goodBJets4.size() < 2 && goodJets8.size() >= 1)
				{
					eeJ_events += 1;

					// meeJ cut
					if (mee < llMaxMass && mee > llMinMass)
					{
						meeJ_events += 1;
					}

					// eeB cut
					if (goodBJets8.size() >= 1)
					{
						eeB_events += 1;

						// meeB cut
						if (mee < llMaxMass && mee > llMinMass)
						{
							meeB_events += 1;
						}

						double mB = goodBJets8[0].P4().M();
						double ptB = goodBJets8[0].P4().Pt();

						double mA = (goodBJets8[0].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).M();
						double ptA = (goodBJets8[0].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Pt();
						double rA = (goodBJets8[0].P4() + goodElectrons[0].P4() + goodElectrons[1].P4()).Rapidity();

						// meemB cut
						if (mee < llMaxMass && mee > llMinMass && mB < BMaxMass && mB > BMinMass)
						{
							meemB_events += 1;

							h_mer_dielectron_mass->Fill(mee);
							h_mer_dielectron_pt->Fill(ptee);
							h_mer_electron_deltar->Fill(dree);

							h_ec_mer_Jet8_mass->Fill(mB);
							h_ec_mer_Jet8_pt->Fill(ptB);

							h_ec_mer_A_mass->Fill(mA);
							h_ec_mer_A_pt->Fill(ptA);
							h_ec_mer_A_rapidity->Fill(rA);
						}
					}
				}
			}

			// muon channel
			// mumu cut
			if (goodMuons.size() >= 2)// && goodElectrons.size() == 0)
			{

				mumu_events += 1;

				double mmumu = (goodMuons[0].P4() + goodMuons[1].P4()).M();
				double ptmumu = (goodMuons[0].P4() + goodMuons[1].P4()).Pt();
				double drmumu = goodMuons[0].P4().DeltaR(goodMuons[1].P4());

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
						double ptbb = (goodBJets4[0].P4() + goodBJets4[1].P4()).Pt();
						double drbb = goodBJets4[0].P4().DeltaR(goodBJets4[1].P4());

						double mA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).M();
						double ptA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Pt();
						double rA = (goodBJets4[0].P4() + goodBJets4[1].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Rapidity();

						// mmumumbb cut
						if (mmumu < llMaxMass && mmumu > llMinMass && mbb < bbMaxMass && mbb > bbMinMass)
						{
							mmumumbb_events += 1;

							h_res_dimuon_mass->Fill(mmumu);
							h_res_dimuon_pt->Fill(ptmumu);
							h_res_muon_deltar->Fill(drmumu);

							h_mc_res_dijet4_mass->Fill(mbb);
							h_mc_res_dijet4_pt->Fill(ptbb);
							h_mc_res_jet4_deltar->Fill(drbb);
							h_mc_res_dijet4_pt_vs_deltar->Fill(ptbb, drbb);

							h_mc_res_A_mass->Fill(mA);
							h_mc_res_A_pt->Fill(ptA);
							h_mc_res_A_rapidity->Fill(rA);
						}
					}
				}
				// mumuJ cut
				if (goodBJets4.size() < 2 && goodJets8.size() >= 1)
				{
					mumuJ_events += 1;

					// mmumuJ cut
					if (mmumu < llMaxMass && mmumu > llMinMass)
					{
						mmumuJ_events += 1;
					}

					// mumuB cut
					if (goodBJets8.size() >= 1)
					{
						mumuB_events += 1;

						// meeB cut
						if (mmumu < llMaxMass && mmumu > llMinMass)
						{
							mmumuB_events += 1;
						}

						double mB = goodBJets8[0].P4().M();
						double ptB = goodBJets8[0].P4().Pt();

						double mA = (goodBJets8[0].P4() + goodMuons[0].P4() + goodMuons[1].P4()).M();
						double ptA = (goodBJets8[0].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Pt();
						double rA = (goodBJets8[0].P4() + goodMuons[0].P4() + goodMuons[1].P4()).Rapidity();

						// mmumumB cut
						if (mmumu < llMaxMass && mmumu > llMinMass && mB < BMaxMass && mB > BMinMass)
						{
							mmumumB_events += 1;

							h_mer_dimuon_mass->Fill(mmumu);
							h_mer_dimuon_pt->Fill(ptmumu);
							h_mer_muon_deltar->Fill(drmumu);

							h_mc_mer_Jet8_mass->Fill(mB);
							h_mc_mer_Jet8_pt->Fill(ptB);

							h_mc_mer_A_mass->Fill(mA);
							h_mc_mer_A_pt->Fill(ptA);
							h_mc_mer_A_rapidity->Fill(rA);
						}
					}
				}
			}
		}
	}

	// make the directories for the plots
	mkdir("cut_plots", 0777);
	mkdir("cut_plots/electron_resolved_channel", 0777);
	mkdir("cut_plots/muon_resolved_channel", 0777);
	mkdir("cut_plots/electron_merged_channel", 0777);
	mkdir("cut_plots/muon_merged_channel", 0777);

	// output the histograms
	TCanvas * c1 = new TCanvas("c1", "c1", 640, 480);

	// electron resolved channel
	// mass
	h_res_dielectron_mass->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/res_dielectron_mass.eps");
	h_ec_res_dijet4_mass->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/ec_res_dijet4_mass.eps");
	h_ec_res_A_mass->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/ec_res_A_mass.eps");

	// pt
	h_res_dielectron_pt->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/res_dielectron_pt.eps");
	h_ec_res_dijet4_pt->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/ec_res_dijet4_pt.eps");
	h_ec_res_A_pt->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/ec_res_A_pt.eps");

	// rapidity
	h_ec_res_A_rapidity->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/ec_res_A_rapidity.eps");

	// deltar
	h_res_electron_deltar->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/res_electron_deltar.eps");
	h_ec_res_jet4_deltar->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/ec_res_jet4_deltar.eps");
	h_ec_res_dijet4_pt_vs_deltar->Draw();
	c1->SaveAs("cut_plots/electron_resolved_channel/ec_res_dijet4_pt_vs_deltar.eps");

	// electron merged channel
	// mass
	h_mer_dielectron_mass->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/mer_dielectron_mass.eps");
	h_ec_mer_Jet8_mass->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/ec_mer_Jet8_mass.eps");
	h_ec_mer_A_mass->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/ec_mer_A_mass.eps");

	// pt
	h_mer_dielectron_pt->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/mer_dielectron_pt.eps");
	h_ec_mer_Jet8_pt->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/ec_mer_Jet8_pt.eps");
	h_ec_mer_A_pt->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/ec_mer_A_pt.eps");

	// rapidity
	h_ec_mer_A_rapidity->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/ec_mer_A_rapidity.eps");

	// deltar
	h_mer_electron_deltar->Draw();
	c1->SaveAs("cut_plots/electron_merged_channel/mer_electron_deltar.eps");

	// muon resolved channel
	// mass
	h_res_dimuon_mass->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/res_dimuon_mass.eps");
	h_mc_res_dijet4_mass->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/mc_res_dijet4_mass.eps");
	h_mc_res_A_mass->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/mc_res_A_mass.eps");

	// pt
	h_res_dimuon_pt->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/res_dimuon_pt.eps");
	h_mc_res_dijet4_pt->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/mc_res_dijet4_pt.eps");
	h_mc_res_A_pt->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/mc_res_A_pt.eps");

	// rapidity
	h_mc_res_A_rapidity->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/mc_res_A_rapidity.eps");

	// deltar
	h_res_muon_deltar->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/res_muon_deltar.eps");
	h_mc_res_jet4_deltar->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/mc_res_jet4_deltar.eps");
	h_mc_res_dijet4_pt_vs_deltar->Draw();
	c1->SaveAs("cut_plots/muon_resolved_channel/mc_res_dijet4_pt_vs_deltar.eps");

	// muon merged channel
	// mass
	h_mer_dimuon_mass->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mer_dimuon_mass.eps");
	h_mc_mer_Jet8_mass->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mc_mer_Jet8_mass.eps");
	h_mc_mer_A_mass->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mc_mer_A_mass.eps");

	// pt
	h_mer_dimuon_pt->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mer_dimuon_pt.eps");
	h_mc_mer_Jet8_pt->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mc_mer_Jet8_pt.eps");
	h_mc_mer_A_pt->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mc_mer_A_pt.eps");

	// rapidity
	h_mc_mer_A_rapidity->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mc_mer_A_rapidity.eps");

	// deltar
	h_mer_muon_deltar->Draw();
	c1->SaveAs("cut_plots/muon_merged_channel/mer_muon_deltar.eps");

	// output event efficiency information
	int ll_events = ee_events + mumu_events;
	int mll_events = mee_events + mmumu_events;
	int lljj_events = eejj_events + mumujj_events;
	int mlljj_events = meejj_events + mmumujj_events;
	int llbb_events = eebb_events + mumubb_events;
	int mllbb_events = meebb_events + mmumubb_events;
	int mllmbb_events = meembb_events + mmumumbb_events;
	int llJ_events = eeJ_events + mumuJ_events;
	int mllJ_events = meeJ_events + mmumuJ_events;
	int llB_events = eeB_events + mumuB_events;
	int mllB_events = meeB_events + mmumuB_events;
	int mllmB_events = meemB_events + mmumumB_events;

	cout << total_events << "\t";

	cout << ee_events << "\t";
	cout << mee_events << "\t";
	cout << eejj_events << "\t";
	cout << meejj_events << "\t";
	cout << eebb_events << "\t";
	cout << meebb_events << "\t";
	cout << meembb_events << "\t";
	cout << eeJ_events << "\t";
	cout << meeJ_events << "\t";
	cout << eeB_events << "\t";
	cout << meeB_events << "\t";
	cout << meemB_events << "\t";

	cout << mumu_events << "\t";
	cout << mmumu_events << "\t";
	cout << mumujj_events << "\t";
	cout << mmumujj_events << "\t";
	cout << mumubb_events << "\t";
	cout << mmumubb_events << "\t";
	cout << mmumumbb_events << "\t";
	cout << mumuJ_events << "\t";
	cout << mmumuJ_events << "\t";
	cout << mumuB_events << "\t";
	cout << mmumuB_events << "\t";
	cout << mmumumB_events << "\t";

	cout << ll_events << "\t";
	cout << mll_events << "\t";
	cout << lljj_events << "\t";
	cout << mlljj_events << "\t";
	cout << llbb_events << "\t";
	cout << mllbb_events << "\t";
	cout << mllmbb_events << "\t";
	cout << llJ_events << "\t";
	cout << mllJ_events << "\t";
	cout << llB_events << "\t";
	cout << mllB_events << "\t";
	cout << mllmB_events << "\t";
	cout << endl;

	return 0;
}
