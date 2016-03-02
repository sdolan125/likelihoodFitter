//////////////////////////////////////////////////////////
//
//  Nucleon FSI modeling parameters
//
//
//
//  Created: Oct 2013
//  Modified:
//
//////////////////////////////////////////////////////////

#ifndef __NuclFSIParameters_hh__
#define __NuclFSIParameters_hh__

#include "AnaFitParameters.hh"
#include "XsecParameters.hh"
#include <TFile.h>
#include <TGraph.h>

struct FSIBin
{
  SampleTypes topology;
  ReactionTypes reaction;
  std::vector<TGraph*> respfuncs;
};

class NuclFSIParameters : public AnaFitParameters
{
public:
  NuclFSIParameters(const char *name = "par_NuclFSI");
  ~NuclFSIParameters();
  
  void StoreResponseFunctions(std::vector<TFile*> respfuncs);
  void InitEventMap(std::vector<AnaSample*> &sample); 
  void EventWeights(std::vector<AnaSample*> &sample, 
		    std::vector<double> &params);
  void ReWeight(AnaEvent *event, int nsample, int nevent,
		std::vector<double> &params);
private:
  int GetBinIndex(SampleTypes sampletype, ReactionTypes reactype);
  std::vector<FSIBin> m_bins;
};

#endif
