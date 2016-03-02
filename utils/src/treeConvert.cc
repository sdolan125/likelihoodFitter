/******************************************************

Code to convert a HL2 tree into the format required 
for the fitting code. C

Can't simply read HL2 tree directly since we don't
know what variables will be the tree

Author: Stephen Dolan
Date Created: November 2015

******************************************************/


#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TMath.h>


using namespace std;

int treeConvert(TString inFileName, TString inTreeName, TString outFileName, TString D1NameRec, TString D1NameTrue,
                TString D2NameRec="", TString D2NameTrue="")
{
  // You need to provide the number of branches in your HL2 tree
  // And the accum_level you want to cut each one at to get your selected events
  // i.e choosing n in accum_level[0][branch_i]>n
  const int nbranches = 10;
  //const int accumToCut[nbranches] = {5,6,7,7,6,4,3,5,7,6};
  const int accumToCut[nbranches] =   {6,7,8,8,7,5,4,6,8,7};

  TFile *infile = new TFile(inFileName);
  TTree *intree = (TTree*)infile->Get(inTreeName);

  TFile *outfile = new TFile(outFileName,"recreate");
  TTree *outtree = new TTree("selectedEvents", "selectedEvents");

  // Declaration of leaf types
  Int_t           accum_level[1500][50];
  Int_t           reaction;
  Int_t           cutBranch=-999;
  Int_t           mectopology;
  Float_t        D1true;
  Float_t        D2true;
  Float_t        D1Reco;
  Float_t        D2Reco;
  Float_t        pMomRec;
  Float_t        pThetaRec;
  Float_t        pMomTrue;
  Float_t        pThetaTrue;
  Float_t        muMomRec;
  Float_t        muThetaRec;
  Float_t        muMomTrue;
  Float_t        muThetaTrue;
  Float_t        muCosThetaRec;
  Float_t        muCosThetaTrue;
  Float_t        pCosThetaRec;
  Float_t        pCosThetaTrue;
  Float_t        RecoNuEnergy=0;
  Float_t        TrueNuEnergy=0;
  Float_t        weight;

  intree->SetBranchAddress("accum_level", &accum_level);
  intree->SetBranchAddress("reaction", &reaction);
  intree->SetBranchAddress("mectopology", &mectopology);
  intree->SetBranchAddress(D1NameTrue, &D1true);
  intree->SetBranchAddress(D2NameTrue, &D2true);
  intree->SetBranchAddress(D1NameRec, &D1Reco);
  intree->SetBranchAddress(D2NameRec, &D2Reco);
  intree->SetBranchAddress("selp_mom", &pMomRec);
  intree->SetBranchAddress("selp_theta" ,&pThetaRec);
  intree->SetBranchAddress("selp_truemom" ,&pMomTrue);
  intree->SetBranchAddress("selp_trueztheta" ,&pThetaTrue);
  intree->SetBranchAddress("selmu_mom", &muMomRec);
  intree->SetBranchAddress("selmu_theta", &muThetaRec);
  intree->SetBranchAddress("selmu_truemom", &muMomTrue);
  //intree->SetBranchAddress("selmu_trueztheta", &muThetaTrue);
  intree->SetBranchAddress("truemu_costheta", &muCosThetaTrue);
  //intree->SetBranchAddress("nu_trueE", &RecoNuEnergy);
  intree->SetBranchAddress("nu_trueE", &TrueNuEnergy);
  intree->SetBranchAddress("weight", &weight);

  outtree->Branch("reaction", &reaction, "reaction/I");
  outtree->Branch("cutBranch", &cutBranch, "cutBranch/I");
  outtree->Branch("mectopology", &mectopology, "mectopology/I");
  outtree->Branch("D1True", &D1true, ("D1True/F"));
  outtree->Branch("D1Rec", &D1Reco, ("D1Rec/F"));
  outtree->Branch("D2True", &D2true, ("D2True/F"));
  outtree->Branch("D2Rec", &D2Reco, ("D2Rec/F"));
  outtree->Branch("muMomRec", &muMomRec, ("muMomRec/F"));
  outtree->Branch("muMomTrue", &muMomTrue, ("muMomTrue/F"));
  outtree->Branch("muCosThetaRec", &muCosThetaRec, ("muCosThetaRec/F"));
  outtree->Branch("muCosThetaTrue", &muCosThetaTrue, ("muCosThetaTrue/F"));
  outtree->Branch("pMomRec", &pMomRec, ("pMomRec/F"));
  outtree->Branch("pMomTrue", &pMomTrue, ("pMomTrue/F"));
  outtree->Branch("pCosThetaRec", &pCosThetaRec, ("pCosThetaRec/F"));
  outtree->Branch("pCosThetaTrue", &pCosThetaTrue, ("pCosThetaTrue/F"));
  outtree->Branch("Enureco", &RecoNuEnergy, "Enureco/F");
  outtree->Branch("Enutrue", &TrueNuEnergy, "Enutrue/F");
  outtree->Branch("weight", &weight, "weight/F");

  Long64_t nentries = intree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int passCount=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = intree->GetEntry(jentry); nbytes += nb;
    passCount=0;
    RecoNuEnergy=TrueNuEnergy;
    pCosThetaRec   = TMath::Cos(pThetaRec);
    pCosThetaTrue  = TMath::Cos(pThetaTrue);
    muCosThetaRec  = TMath::Cos(muThetaRec);
    //muCosThetaTrue = TMath::Cos(muThetaTrue);
    int branches_passed[10]={0};
    for(int i=0; i<nbranches; i++){
      if(accum_level[0][i]>accumToCut[i]){
        cutBranch=i; passCount++; 
        branches_passed[i]++;
        outtree->Fill();
      }
    }
    if(passCount>1){
      printf("***Warning: More than one cut branch passed***\n");
      for(int j=0;j<10;j++){
        if(branches_passed[j]==1) printf("branch %d passed ...",j);
      }
      printf("\n");
    }
  }

  outtree->Print();
  outfile->Write();

  delete infile;
  delete outfile;
  return 0;
}
