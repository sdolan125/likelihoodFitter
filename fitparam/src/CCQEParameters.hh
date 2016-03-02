//////////////////////////////////////////////////////////
//
//  CCQE cross-section parameters
//
//
//
//  Created: Thu Jun 13 11:46:03 CEST 2013
//  Modified:
//
//////////////////////////////////////////////////////////

#ifndef __CCQEParameters_hh__
#define __CCQEParameters_hh__

#include "AnaFitParameters.hh"

struct CCQEBin
{
  double plow, phigh;
  double cthlow, cthhigh;
};

class CCQEParameters : public AnaFitParameters
{
public:
  CCQEParameters(const char *fname,
		 const char *name = "par_ccqe");
  ~CCQEParameters();
  
  void InitEventMap(std::vector<AnaSample*> &sample);
  void EventWeights(std::vector<AnaSample*> &sample, 
		    std::vector<double> &params);
  void ReWeight(AnaEvent *event, int nsample, int nevent,
		std::vector<double> &params);
  void SetBinning(const char *fname);

private:
  //binnig function
  int GetBinIndex(double p, double cth); 

  int ccqe_recode;              //reaction code
  std::vector<CCQEBin> m_bins; //binning
};

#endif
