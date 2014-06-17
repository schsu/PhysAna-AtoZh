// C headers
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdint.h>

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
			fprintf(stderr, "%s does not exist\n", inputFile);
		}
		else if (errno == EACCES)
		{
			fprintf(stderr, "%s is not accessible\n", inputFile);
		}
		return 1;
	}

	// create object of class TChain
	TChain * tc = new TChain("Delphes");
	tc->Add(inputFile);

	// create object of class ExRootTreeReader
	ExRootTreeReader * tr = new ExRootTreeReader(tc);
	int64_t numberofEntries = tr->GetEntries();

	// get pointer to particle branch
	TClonesArray * branchParticle = tr->UseBranch("Particle");

	// create the histogram
	TH1F * h_truth_electron_pt = new TH1F("truth electron pt", "truth electron pt; pt (GeV, 100 bins); count", 100, 0, 300);

	// loop over each event
	for (int64_t entry = 0; entry < numberofEntries; entry++)
	{
		// load selected branches
		tr->ReadEntry(entry);

		// fill the histogram
		for (int i = 0; i < branchParticle->GetEntries(); i++)
		{
			GenParticle * particle = (GenParticle*) branchParticle->At(i);
			if (particle->Status == 1)
			{
				if (particle->PID == 11 || particle->PID == -11)
				{
					h_truth_electron_pt->Fill(particle->PT);
				}
			}
		}
	}

	// output the histogram
	TCanvas * c1 = new TCanvas("c1", "c1", 640, 480);

	h_truth_electron_pt->Draw();
	c1->SaveAs("truth_electron_pt.eps");

	return 0;
}
