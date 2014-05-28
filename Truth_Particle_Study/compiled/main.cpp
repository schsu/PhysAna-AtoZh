#include <iostream>
#include <errno.h>
#include <unistd.h>

#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TObject.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

using namespace std;

int main() {
  char inputFile[1024];
  cout << "Enter the full path to a root file: ";
  cin >> inputFile;

  if (access(inputFile, F_OK) != 0) {
    if (errno == ENOENT) {
      cout << inputFile << " does not exist" << endl;
    }
    else if (errno == EACCES) {
      cout << inputFile << " is not accessible" << endl;
    }
    return 1;
  }

  TChain* tc = new TChain("Delphes");
  tc->Add(inputFile);
  ExRootTreeReader* tr = new ExRootTreeReader(tc);
  Long64_t numberofEntries = tr->GetEntries();
  Int_t entry;

  TClonesArray* branchParticle = tr->UseBranch("Particle");
  Int_t numberofParticle(0) ; //branchParticle->GetEntries();

  cout << "There are " <<  numberofEntries << " Events." << endl;
  cout << "Enter the number of an Event you would like to read: ";
  cin >> entry;
  tr->ReadEntry(entry);

  cout << "There are " <<  branchParticle->GetEntries() << " particles." << endl;
  cout << "Enter the number of particles you would like information on: ";
  cin >> numberofParticle;

  //print out truth

  cout << "Particle: " << branchParticle->GetEntries() << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "   #  PID  Status  M1  M2  D1  D2  Charge    Mass        E        Px        Py        Pz        PT       Eta       Phi     Rapidity " << endl;
  cout << "------------------------------------------------------------------------------------------------------------------------------------" << endl;

  for (Int_t i; i < numberofParticle; i++) {
    GenParticle* particle = (GenParticle*) branchParticle->At(i);
    cout.width(4);
    cout.fill(' ');

    cout << i ;
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
    cout << std::scientific;

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
