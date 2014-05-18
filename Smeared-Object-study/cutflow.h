#ifndef _CUTFLOW_H
#define _CUTFLOW_H

#include <vector>
#include "TClonesArray.h"
#include "HelperClasses.h"

void SelectGoodJets(std::vector<BaseParticle*>& goodJets, TClonesArray* branchJet, const std::vector<BaseParticle*>& goodElectrons);
void cutflow(char* inputFile);



#endif
