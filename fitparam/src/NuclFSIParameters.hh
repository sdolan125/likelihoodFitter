//////////////////////////////////////////////////////////
//
//  NuclFSI modeling parameters
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

struct NuclFSIBin
{
  double recoPlow, recoPhigh;
  double truePlow, truePhigh;
  double recoCTHlow, recoCTHhigh;
  double trueCTHlow, trueCTHhigh;
  SampleTypes topology;
  ReactionTypes reaction;
  std::vector<TGraph*> respfuncs;
};

class NuclFSIParameters : public AnaFitParameters
{
public:
  NuclFSIParameters(const char *name = "par_NuclFSI");
  ~NuclFSIParameters();
  
  void StoreResponseFunctions(std::vector<TFile*> respfuncs,
			      std::vector<std::pair <double,double> > v_pedges, 
			      std::vector<std::pair <double,double> > v_cthedges);
  void InitEventMap(std::vector<AnaSample*> &sample); 
  void EventWeights(std::vector<AnaSample*> &sample, 
		    std::vector<double> &params);
  void ReWeight(AnaEvent *event, int nsample, int nevent,
		std::vector<double> &params);
private:
  int GetBinIndex(SampleTypes sampletype, ReactionTypes reactype, 
		  double recoP, double trueP, double recoCTH, double trueCTH);
  std::vector<NuclFSIBin> m_bins;
};

#endif
