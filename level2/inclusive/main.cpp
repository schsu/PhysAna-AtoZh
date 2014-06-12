#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>

#include <errno.h>
#include <unistd.h>

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

	// get pointer to particle, electron, and muon branches
	TClonesArray * branchParticle = tr->UseBranch("Particle");
	TClonesArray * branchElectron = tr->UseBranch("Electron");
	TClonesArray * branchMuon = tr->UseBranch("Muon");

	// counter for electrons and muons in particle, electron, and muon branches
	int64_t electrons = 0;
	int64_t muons = 0;
	int64_t truth_electrons = 0;
	int64_t truth_muons = 0;

	// create the histograms
	TH1F * h_electron_pt = new TH1F("electron pt", "electron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	h_electron_pt->Sumw2();
	TH1F * h_truth_electron_pt = new TH1F("truth electron pt", "truth electron pt; pt (GeV, 100 bins); count", 100, 0, 300);
	h_truth_electron_pt->Sumw2();
	TH1F * h_electron_eta = new TH1F("electron eta", "electron eta; eta (100 bins); count", 100, -3.4, 3.4);
	h_electron_eta->Sumw2();
	TH1F * h_electron_phi = new TH1F("electron phi", "electron phi; phi (100 bins); count", 100, -3.4, 3.4);
	h_electron_phi->Sumw2();
	TH1F * h_electron_multiplicity = new TH1F("electron multiplicity", "electron multiplicity; electrons per event; count", 4, 0, 4);
	h_electron_multiplicity->Sumw2();

	TH1F * h_muon_pt = new TH1F("muon pt", "muon pt; pt (GeV 100 bins); count", 100, 0, 300);
	h_muon_pt->Sumw2();
	TH1F * h_truth_muon_pt = new TH1F("truth muon pt", "truth muon pt; pt (GeV 100 bins); count", 100, 0, 300);
	h_truth_muon_pt->Sumw2();
	TH1F * h_muon_eta = new TH1F("muon eta", "muon eta; eta (100 bins); count", 100, -3.4, 3.4);
	h_muon_eta->Sumw2();
	TH1F * h_muon_phi = new TH1F("muon phi", "muon phi; phi (100 bins); count", 100, -3.4, 3.4);
	h_muon_phi->Sumw2();
	TH1F * h_muon_multiplicity = new TH1F("muon multiplicity", "muon multiplicity; muons per event; count", 4, 0, 4);
	h_muon_multiplicity->Sumw2();

	// loop over each event
	for(Int_t entry = 0; entry < numberofEntries; entry++)
	{
		// load selected branches
		tr->ReadEntry(entry);

		// add the number of electrons and muons in the event
		electrons += branchElectron->GetEntries();
		muons += branchMuon->GetEntries();

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
			h_electron_eta->Fill(electron->Eta);
			h_electron_phi->Fill(electron->Phi);
		}
		for (int i = 0; i < branchMuon->GetEntries(); i++)
		{
			Muon * muon = (Muon*) branchMuon->At(i);
			h_muon_pt->Fill(muon->PT);
			h_muon_eta->Fill(muon->Eta);
			h_muon_phi->Fill(muon->Phi);
		}
	}

	// output the number of electrons and muons
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

	// output the histograms
	TCanvas * c1 = new TCanvas("c1", "c1", 640, 480);

	h_electron_pt->Draw();
	c1->SaveAs("electron_pt.eps");
	h_truth_electron_pt->Draw();
	c1->SaveAs("truth_electron_pt.eps");
	h_electron_eta->Draw();
	c1->SaveAs("electron_eta.eps");
	h_electron_phi->Draw();
	c1->SaveAs("electron_phi.eps");

	h_muon_pt->Draw();
	c1->SaveAs("muon_pt.eps");
	h_truth_muon_pt->Draw();
	c1->SaveAs("truth_muon_pt.eps");
	h_muon_eta->Draw();
	c1->SaveAs("muon_eta.eps");
	h_muon_phi->Draw();
	c1->SaveAs("muon_phi.eps");

	TH1F * h_electron_pt_divide = (TH1F*) h_electron_pt->Clone("electron pt divide");
	h_electron_pt_divide->SetTitle("electron pt divide; pt (GeV, 100 bins); count");
	h_electron_pt_divide->Divide(h_truth_electron_pt);
	h_electron_pt_divide->Draw();
	c1->SaveAs("electron_pt_divide.eps");

	TH1F * h_muon_pt_divide = (TH1F*) h_muon_pt->Clone("muon pt divide");
	h_muon_pt_divide->SetTitle("muon pt divide; pt (GeV, 100 bins); count");
	h_muon_pt_divide->Divide(h_truth_muon_pt);
	h_muon_pt_divide->Draw();
	c1->SaveAs("muon_pt_divide.eps");

	h_electron_multiplicity->Draw();
	c1->SaveAs("electron_multiplicity.eps");
	h_muon_multiplicity->Draw();
	c1->SaveAs("muon_multiplicity.eps");

	return 0;
}
