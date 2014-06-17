// C++ headers
#include <iostream>

// C headers
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>

// ROOT headers
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TObject.h"

// Delphes headers (requires that Delphes' library path be added to your LD_LIBRARY_PATH)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

using namespace std;

int main (int argc, char * argv[])
{
	// command line input
	if (argc != 2 && argc != 4)
	{
		fprintf(stderr, "ERROR - wrong number of arguments\n");
		fprintf(stderr, "usage: %s inputFile [event] [numberofParticles]\n", argv[0]);
		return 1;
	}

	char inputFile[1024];
	sprintf(inputFile, argv[1]);

	// the event and number of particles to read
	int64_t event;
	int64_t numberofParticles;

	if (argc == 4)
	{
		event = atoi(argv[2]);
		numberofParticles = atoi(argv[3]);
	}

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
	int64_t numberofEvents = tr->GetEntries();

	// get pointer to particle branch
	TClonesArray * branchParticle = tr->UseBranch("Particle");

	if (argc != 4)
	{
		cout << "There are " <<  numberofEvents << " Events." << endl;
		cout << "Enter the number of an Event you would like to read: ";
		cin >> event;
		tr->ReadEntry(event);

		cout << "There are " <<  branchParticle->GetEntries() << " particles in this Event." << endl;
		cout << "Enter the number of particles you would like information on: ";
		cin >> numberofParticles;
	}
	else
	{
		if (event > numberofEvents)
		{
			fprintf(stderr, "The user has specified an Event not contained in the root file, aborting.\n");
		}
		else
		{
			tr->ReadEntry(event);
			if (numberofParticles > branchParticle->GetEntries())
			{
				fprintf(stderr, "There are %d particles in Event %d, the user has requested information on %d particles, only the information for %d particles will be output.\n", branchParticle->GetEntries(), event, numberofParticles, branchParticle->GetEntries());
				numberofParticles = branchParticle->GetEntries();
			}
		}
	}

	cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "   #  PID  Status  M1  M2  D1  D2  Charge    Mass        E        Px        Py        Pz        PT       Eta       Phi     Rapidity " << endl;
	cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;

	for (int64_t i = 0; i < numberofParticles; i++)
	{
		GenParticle * particle = (GenParticle*) branchParticle->At(i);
		cout.width(4);
		cout.fill(' ');

		cout << i;
		cout.width(5);
		cout << particle->PID;

		cout.width(8);
		cout << particle->Status;

		cout.width(4);
		cout << particle->M1;
		cout.width(4);
		cout << particle->M2;
		cout.width(4);
		cout << particle->D1;
		cout.width(4);
		cout << particle->D2;

		cout.width(8);
		cout << particle->Charge;

		cout.precision(2);
		cout << scientific;

		cout.width(10);
		cout << particle->Mass;
		cout.width(10);
		cout << particle->E;
		cout.width(10);
		cout << particle->Px;
		cout.width(10);
		cout << particle->Py;
		cout.width(10);
		cout << particle->Pz;
		cout.width(10);
		cout << particle->PT;
		cout.width(10);
		cout << particle->Eta;
		cout.width(10);
		cout << particle->Phi;
		cout.width(10);
		cout << particle->Rapidity;

		cout << endl;
	}

	cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;

	return 0;
}
