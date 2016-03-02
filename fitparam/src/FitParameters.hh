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

#ifndef __FitParameters_hh__
#define __FitParameters_hh__

#include "AnaFitParameters.hh"

struct CCQEBin
{
  double D1low, D1high;
  double D2low, D2high;
};

class FitParameters : public AnaFitParameters
{
public:
  FitParameters(const char *fname,
		 const char *name = "par_ccqe");
  ~FitParameters();
  
  void InitEventMap(std::vector<AnaSample*> &sample);
  void EventWeights(std::vector<AnaSample*> &sample, 
		    std::vector<double> &params);
  void ReWeight(AnaEvent *event, int nsample, int nevent,
		std::vector<double> &params);
  void SetBinning(const char *fname);

private:
  //binnig function
  int GetBinIndex(double D1, double D2); 

  int ccqe_recode;              //reaction code
  std::vector<CCQEBin> m_bins; //binning
};

#endif
