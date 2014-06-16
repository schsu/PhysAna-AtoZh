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

	// get pointer to particle, electron, muon, and jet branches
	TClonesArray * branchParticle = tr->UseBranch("Particle");
	TClonesArray * branchElectron = tr->UseBranch("Electron");
	TClonesArray * branchMuon = tr->UseBranch("Muon");
	TClonesArray * branchJet = tr->UseBranch("Jet");

	// counter for electrons, muons, and jets in particle, electron, muon, and jet branches
	int64_t electrons = 0;
	int64_t muons = 0;
	int64_t jets = 0;
	int64_t truth_electrons = 0;
	int64_t truth_muons = 0;

	// create the histograms
	TH1F * h_electron_pt = new TH1F("electron pt", "electron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_electron_e = new TH1F("electron E", "electron E; E (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_electron_eta = new TH1F("electron eta", "electron eta; eta (100 bins); count", 100, -3.4, 3.4);
	TH1F * h_electron_phi = new TH1F("electron phi", "electron phi; phi (100 bins); count", 100, -3.4, 3.4);

	TH1F * h_muon_pt = new TH1F("muon pt", "muon pt; pt (GeV 100 bins); count", 100, 0, 300);
	TH1F * h_muon_e = new TH1F("muon E", "muon E; E (GeV 100 bins); count", 100, 0, 300);
	TH1F * h_muon_eta = new TH1F("muon eta", "muon eta; eta (100 bins); count", 100, -3.4, 3.4);
	TH1F * h_muon_phi = new TH1F("muon phi", "muon phi; phi (100 bins); count", 100, -3.4, 3.4);

	TH1F * h_jet_pt = new TH1F("jet pt", "jet pt; pt (GeV 100 bins); count", 100, 0, 300);
	TH1F * h_jet_e = new TH1F("jet E", "jet E; E (GeV 100 bins); count", 100, 0, 300);
	TH1F * h_jet_eta = new TH1F("jet eta", "jet eta; eta (100 bins); count", 100, -3.4, 3.4);
	TH1F * h_jet_phi = new TH1F("jet phi", "jet phi; phi (100 bins); count", 100, -3.4, 3.4);

	TH1F * h_truth_electron_pt = new TH1F("truth electron pt", "truth electron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	TH1F * h_truth_muon_pt = new TH1F("truth muon pt", "truth muon pt; pt (GeV 100 bins); count", 100, 0, 300);

	TH1F * h_electron_multiplicity = new TH1F("electron multiplicity", "electron multiplicity; electrons per event; count", 4, 0, 4);
	TH1F * h_muon_multiplicity = new TH1F("muon multiplicity", "muon multiplicity; muons per event; count", 4, 0, 4);

	// loop over each event
	for(Int_t entry = 0; entry < numberofEntries; entry++)
	{
		// load selected branches
		tr->ReadEntry(entry);

		// add the number of electrons, muons, and jets in the event
		electrons += branchElectron->GetEntries();
		muons += branchMuon->GetEntries();
		jets += branchJet->GetEntries();

		// add the number of truth electrons and muons in the event
		for (int i = 0; i < branchParticle->GetEntries(); i++)
		{
			GenParticle * particle = (GenParticle*) branchParticle->At(i);
			if (particle->Status == 1)
			{
				if (particle->PID == 11 || particle->PID == -11)
				{
					truth_electrons += 1;
					h_truth_electron_pt->Fill(particle->PT);
				}
				if (particle->PID == 13 || particle->PID == -13)
				{
					truth_muons += 1;
					h_truth_muon_pt->Fill(particle->PT);
				}
			}
		}

		// fill the histograms
		h_electron_multiplicity->Fill(branchElectron->GetEntries());
		h_muon_multiplicity->Fill(branchMuon->GetEntries());

		for (int i = 0; i < branchElectron->GetEntries(); i++)
		{
			Electron * electron = (Electron*) branchElectron->At(i);
			h_electron_pt->Fill(electron->PT);
			h_electron_e->Fill(electron->P4().E());
			h_electron_eta->Fill(electron->Eta);
			h_electron_phi->Fill(electron->Phi);
		}
		for (int i = 0; i < branchMuon->GetEntries(); i++)
		{
			Muon * muon = (Muon*) branchMuon->At(i);
			h_muon_pt->Fill(muon->PT);
			h_muon_e->Fill(muon->P4().E());
			h_muon_eta->Fill(muon->Eta);
			h_muon_phi->Fill(muon->Phi);
		}
		for (int i = 0; i < branchJet->GetEntries(); i++)
		{
			Jet * jet = (Jet*) branchJet->At(i);
			h_jet_pt->Fill(jet->PT);
			h_jet_e->Fill(jet->P4().E());
			h_jet_eta->Fill(jet->Eta);
			h_jet_phi->Fill(jet->Phi);
		}
	}

	// output the number of electrons, muons, and jets
	cout.width(20);
	cout << "particle";
	cout.width(20);
	cout << "detected";
	cout.width(20);
	cout << "truth";
	cout.width(20);
	cout << "efficiency" << endl;

	cout.width(20);
	cout << "electrons";
	cout.width(20);
	cout << electrons;
	cout.width(20);
	cout << truth_electrons;
	cout.width(20);
	cout << (float) electrons / (float) truth_electrons << endl;

	cout.width(20);
	cout << "muons";
	cout.width(20);
	cout << muons;
	cout.width(20);
	cout << truth_muons;
	cout.width(20);
	cout << (float) muons / (float) truth_muons << endl;

	cout.width(20);
	cout << "jets";
	cout.width(20);
	cout << jets << endl;

	// make the directory for the plots
	mkdir("inclusive_plots", 0777);

	// output the histograms
	TCanvas * c1 = new TCanvas("c1", "c1", 640, 480);

	// electron, muon, and jet branch kinematics
	h_electron_pt->Draw();
	c1->SaveAs("inclusive_plots/electron_pt.eps");
	h_electron_e->Draw();
	c1->SaveAs("inclusive_plots/electron_E.eps");
	h_electron_eta->Draw();
	c1->SaveAs("inclusive_plots/electron_eta.eps");
	h_electron_phi->Draw();
	c1->SaveAs("inclusive_plots/electron_phi.eps");

	h_muon_pt->Draw();
	c1->SaveAs("inclusive_plots/muon_pt.eps");
	h_muon_e->Draw();
	c1->SaveAs("inclusive_plots/muon_E.eps");
	h_muon_eta->Draw();
	c1->SaveAs("inclusive_plots/muon_eta.eps");
	h_muon_phi->Draw();
	c1->SaveAs("inclusive_plots/muon_phi.eps");

	h_jet_pt->Draw();
	c1->SaveAs("inclusive_plots/jet_pt.eps");
	h_jet_e->Draw();
	c1->SaveAs("inclusive_plots/jet_E.eps");
	h_jet_eta->Draw();
	c1->SaveAs("inclusive_plots/jet_eta.eps");
	h_jet_phi->Draw();
	c1->SaveAs("inclusive_plots/jet_phi.eps");

	// particle branch pt plots
	h_truth_electron_pt->Draw();
	c1->SaveAs("inclusive_plots/truth_electron_pt.eps");
	h_truth_muon_pt->Draw();
	c1->SaveAs("inclusive_plots/truth_muon_pt.eps");

	// pt divide plots
	h_electron_pt->Sumw2();
	h_truth_electron_pt->Sumw2();
	TH1F * h_electron_pt_divide = (TH1F*) h_electron_pt->Clone("electron pt divide");
	h_electron_pt_divide->SetTitle("electron pt divide; pt (GeV, 100 bins); count");
	h_electron_pt_divide->Divide(h_truth_electron_pt);
	h_electron_pt_divide->Draw();
	c1->SaveAs("inclusive_plots/electron_pt_divide.eps");

	h_muon_pt->Sumw2();
	h_truth_muon_pt->Sumw2();
	TH1F * h_muon_pt_divide = (TH1F*) h_muon_pt->Clone("muon pt divide");
	h_muon_pt_divide->SetTitle("muon pt divide; pt (GeV, 100 bins); count");
	h_muon_pt_divide->Divide(h_truth_muon_pt);
	h_muon_pt_divide->Draw();
	c1->SaveAs("inclusive_plots/muon_pt_divide.eps");

	// multiplicity plots
	h_electron_multiplicity->Draw();
	c1->SaveAs("inclusive_plots/electron_multiplicity.eps");
	h_muon_multiplicity->Draw();
	c1->SaveAs("inclusive_plots/muon_multiplicity.eps");

	return 0;
}
