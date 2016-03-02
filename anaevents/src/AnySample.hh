//////////////////////////////////////////////////////////
//
//  A class for event samples for for any analysis
//
//
//
//  Created: Nov 17 CEST 2015   
//
//////////////////////////////////////////////////////////
#ifndef __AnySample_hh__
#define __AnySample_hh__

#include <string>
#include <vector>

#include <TH2D.h>
#include <TDirectory.h>
#include <TRandom3.h>
#include <TTree.h>

#include "AnaEvent.hh"
#include "AnaSample.hh"

///////////////////////////////////////
// Class definition
///////////////////////////////////////
class AnySample : public AnaSample
{
public:
  AnySample(int sample_id, std::string name, 
       std::vector<std::pair <double,double> > v_d1edges, 
       std::vector<std::pair <double,double> > v_d2edges, TTree *data);
  ~AnySample();
  
  //binning for various histograms
  void SetPmuBinning(int nbins, double *bins);
  void SetCThmuBinning(int nbins, double *bins);
  void SetEnuBinning(int nbins, double *bins);
  void MakeHistos(); //must be called after binning is changed
  
  //histogram for event distributions
  void SetData(TObject *hdata);
  void FillEventHisto(int datatype);
  double CalcChi2();
  
  TH1D* GetPredHisto(){ return m_hpred; }
  TH1D* GetDataHisto(){ return m_hdata; }

  void GetSampleBreakdown(TDirectory *dirout, std::string tag);
  void Write(TDirectory *dirout, const char *bsname, int fititer);
  int GetSampleId(){ return sample_id; }

private:
  int sample_id;
  TH1D *m_hmc;
  TH1D *m_hpred; //n(pRec_mu, thetaRec_mu)
  TH1D *m_hdata;
  TTree * m_data_tree;
  int nbins_pmu, nbins_cthmu, nbins_enu, nAnybins;
  double *bins_pmu, *bins_cthmu, *bins_enu, *bins_Any; 
  std::vector<std::pair<double, double> > m_pedges;
  std::vector<std::pair<double, double> > m_cthedges;
};

#endif
