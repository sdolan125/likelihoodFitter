#ifndef ccqeTreeMC_h
#define ccqeTreeMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <vector>

#include "AnaEvent.hh"
#include "AnaSample.hh"

class ccqeTreeMC 
{
public :
  TChain         *fChain; //!pointer to the analyzed TTree or TChain
 
  // Declaration of leaf types
  Int_t           reaction;
  Double_t        trueMom;
  Double_t        trueCostheta;
  Int_t           qesampleFinal;
  Double_t        MainMomGlb;
  Double_t        MainCosTheta;
  Double_t        MainRecEneGlb;
  Double_t        TrueEnergy;
  Double_t        weight;

  // List of branches
  TBranch        *b_reaction;   //!
  TBranch        *b_trueMom;   //!
  TBranch        *b_trueCostheta;   //!
  TBranch        *b_qesampleFinal;   //!
  TBranch        *b_MainMomGlb;   //!
  TBranch        *b_MainCosTheta;   //!
  TBranch        *b_MainRecEneGlb;   //!
  TBranch        *b_TrueEnergy;   //!
  TBranch        *b_weight;   //!

  ccqeTreeMC(const char *fname);
  virtual ~ccqeTreeMC();
  virtual Int_t GetEntry(Long64_t entry);
  virtual void  Init();
  virtual void  GetEvents(std::vector<AnaSample*> ana_samples);
};

#endif

#ifdef ccqeTreeMC_cxx
ccqeTreeMC::ccqeTreeMC(const char *fname)
{
  const char *trname = "selectedEvents";
  //init TChain object
  fChain = new TChain(trname); 
  //add files to TChain
  fChain->Add(fname);
  //init tree leaf pntrs
  Init();
}

ccqeTreeMC::~ccqeTreeMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ccqeTreeMC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void ccqeTreeMC::Init()
{
   // Set branch addresses and branch pointers
  fChain->SetBranchAddress("reaction", &reaction, &b_reaction);
  fChain->SetBranchAddress("Pmutrue", &trueMom, &b_trueMom);
  fChain->SetBranchAddress("CTHmutrue", &trueCostheta, &b_trueCostheta);
  fChain->SetBranchAddress("topology", &qesampleFinal, &b_qesampleFinal);
  fChain->SetBranchAddress("Pmureco", &MainMomGlb, &b_MainMomGlb);
  fChain->SetBranchAddress("CTHmureco", &MainCosTheta, &b_MainCosTheta);
  fChain->SetBranchAddress("Enureco", &MainRecEneGlb, &b_MainRecEneGlb);
  fChain->SetBranchAddress("Enutrue", &TrueEnergy, &b_TrueEnergy);
  fChain->SetBranchAddress("weight", &weight, &b_weight);
}

#endif
