//////////////////////////////////////////////////////////
//
//  Xsec modeling parameters
//
//
//
//  Created: Oct 2013
//  Modified:
//
//////////////////////////////////////////////////////////

#ifndef __XsecParameters_hh__
#define __XsecParameters_hh__

#include "AnaFitParameters.hh"
#include <TFile.h>
#include <TGraph.h>

// Sample types
// 0 --> single track
// 1 --> mu+p in TPC
// 2 --> mu+p in FGD
// 3 --> p in TPC, mu in FGD
// 4 --> control sample CC1pi
// 5 --> control sample DIS
enum SampleTypes { mutrack = 0, 
		   mupTPC = 1,
		   mupFGD = 2,
		   muFGDpTPC = 3,
		   crCC1pi = 4,
		   crDIS = 5 };
		   
// Reaction and backgrounds
// the indices should match with ntuples
// 0 - CCQE
// 1 - CC1pi
// 2 - CCDIS
// 3 - NC & Anti-nu background
// 4 - Ouf-of-FGD background
enum ReactionTypes { ReCCQE = 0, 
		     ReCC1pi = 1, 
		     ReCCDIS = 2,
		     ReNCAntiNu = 3, 
		     OutFGD = 4 };


struct XsecBin
{
  double recoPlow, recoPhigh;
  double truePlow, truePhigh;
  double recoCTHlow, recoCTHhigh;
  double trueCTHlow, trueCTHhigh;
  SampleTypes topology;
  ReactionTypes reaction;
  std::vector<TGraph*> respfuncs;
};

class XsecParameters : public AnaFitParameters
{
public:
  XsecParameters(const char *name = "par_xsec");
  ~XsecParameters();
  
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
  std::vector<XsecBin> m_bins;
};

#endif
