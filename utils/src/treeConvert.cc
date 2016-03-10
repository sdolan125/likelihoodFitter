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
    const int nbranches = 9;
    //const int accumToCut[nbranches] = {5,6,7,7,6,4,3,5,7,6};
    const int accumToCut[nbranches] =   {7,8,7,6,7,7,8,6,7};
    
    TFile *infile = new TFile(inFileName);
    TTree *intree = (TTree*)infile->Get(inTreeName);
    
    TFile *outfile = new TFile(outFileName,"recreate");
    TTree *outtree = new TTree("selectedEvents", "selectedEvents");
    
    // Declaration of leaf types
    Int_t          accum_level[1500][50];
    Int_t          reaction;
    Int_t          cutBranch=-999;
    Int_t          mectopology;
    Float_t        D1true;//truth var 1 you are fitting for 1D diff xsec
    Float_t        D2true;//truth var 2 you are fitting for 1D diff xsec
    Float_t        D1Reco;//recon var 1 you are fitting for 12 diff xsec
    Float_t        D2Reco;//recon var 2 you are fitting for 12 diff xsec
    Float_t        pMomRec;
    Float_t        pThetaRec;
    Float_t        pMomTrue;
    Float_t        pThetaTrue;
    Float_t        muMomRec;
    Float_t        muThetaRec;
    Float_t        muMomTrue;
    //Float_t        muThetaTrue;
    //Float_t        muCosThetaRec;
    Float_t        muCosThetaTrue;
    //Float_t        pCosThetaRec;
    //Float_t        pCosThetaTrue;
    Float_t        RecoNuEnergy=0;
    Float_t        TrueNuEnergy=0;
    Float_t        weight;
    
    //**********  Coplowe vars **********//
    Int_t           nu_truereac;
    Int_t           target;
    Float_t         selpi_mom;
    Float_t         selpi_truemom;
    Float_t         selpi_costheta;
    Float_t         selpi_truecostheta;
    //***********************************//
    
    printf("Setting Branch Inofrmation\n");
    intree->SetBranchAddress("accum_level", &accum_level);//in coplowe's HL2
    intree->SetBranchAddress("reaction", &reaction);//in coplowe's HL2
    intree->SetBranchAddress("mectopology", &mectopology);//in coplowe's HL2
    intree->SetBranchAddress(D1NameTrue, &D1true);
    intree->SetBranchAddress(D2NameTrue, &D2true);
    intree->SetBranchAddress(D1NameRec, &D1Reco);
    intree->SetBranchAddress(D2NameRec, &D2Reco);
    intree->SetBranchAddress("selp_mom", &pMomRec);//in coplowe's HL2
    intree->SetBranchAddress("selp_costheta" ,&pThetaRec);//in coplowe's HL2
    intree->SetBranchAddress("selp_truemom" ,&pMomTrue);//in coplowe's HL2
    intree->SetBranchAddress("selp_truecostheta" ,&pThetaTrue);//in coplowe's HL2
    intree->SetBranchAddress("selmu_mom", &muMomRec);//in coplowe's HL2
    intree->SetBranchAddress("selmu_costheta", &muThetaRec);//in coplowe's HL2
    intree->SetBranchAddress("selmu_truemom", &muMomTrue);//in coplowe's HL2
    //intree->SetBranchAddress("selmu_trueztheta", &muThetaTrue);
    intree->SetBranchAddress("truemu_costheta", &muCosThetaTrue);//in coplowe's HL2
    //intree->SetBranchAddress("nu_trueE", &RecoNuEnergy);
    intree->SetBranchAddress("nu_trueE", &TrueNuEnergy);//in coplowe's HL2
    intree->SetBranchAddress("weight", &weight);//in coplowe's HL2
    
    //**********  Coplowe vars **********//
    intree->SetBranchAddress("nu_truereac", &nu_truereac);
    intree->SetBranchAddress("target", &target);
    intree->SetBranchAddress("selpi_mom", &selpi_mom);
    intree->SetBranchAddress("selpi_truemom", &selpi_truemom);
    intree->SetBranchAddress("selpi_costheta", &selpi_costheta);
    intree->SetBranchAddress("selpi_truecostheta", &selpi_truecostheta);
    //***********************************//
    
    outtree->Branch("reaction", &reaction, "reaction/I");
    outtree->Branch("cutBranch", &cutBranch, "cutBranch/I");
    outtree->Branch("mectopology", &mectopology, "mectopology/I");
    outtree->Branch("D1True", &D1true, ("D1True/F"));
    outtree->Branch("D1Rec", &D1Reco, ("D1Rec/F"));
    outtree->Branch("D2True", &D2true, ("D2True/F"));
    outtree->Branch("D2Rec", &D2Reco, ("D2Rec/F"));
    outtree->Branch("muMomRec", &muMomRec, ("muMomRec/F"));
    outtree->Branch("muMomTrue", &muMomTrue, ("muMomTrue/F"));
    outtree->Branch("muCosThetaRec", &muThetaRec, ("muCosThetaRec/F"));
    outtree->Branch("muCosThetaTrue", &muCosThetaTrue, ("muCosThetaTrue/F"));
    outtree->Branch("pMomRec", &pMomRec, ("pMomRec/F"));
    outtree->Branch("pMomTrue", &pMomTrue, ("pMomTrue/F"));
    outtree->Branch("pCosThetaRec", &pThetaRec, ("pCosThetaRec/F"));
    outtree->Branch("pCosThetaTrue", &pThetaTrue, ("pCosThetaTrue/F"));
    outtree->Branch("Enureco", &RecoNuEnergy, "Enureco/F");
    outtree->Branch("Enutrue", &TrueNuEnergy, "Enutrue/F");
    outtree->Branch("weight", &weight, "weight/F");
    
    //**********  Coplowe vars **********//
    outtree->Branch("nu_truereac", &nu_truereac, "nu_truereac/I");
    outtree->Branch("target", &target, "target/I");
    outtree->Branch("selpi_mom", &selpi_mom, "selpi_mom/F");
    outtree->Branch("selpi_truemom", &selpi_truemom, "selpi_truemom/F");
    outtree->Branch("selpi_costheta", &selpi_costheta, "selpi_costheta/F");
    outtree->Branch("selpi_truecostheta", &selpi_truecostheta, "selpi_truecostheta/F");
    //***********************************//
    
    Long64_t nentries = intree->GetEntriesFast();
    printf("Coverting tree with %d entries.\n",nentries);
    Long64_t nbytes = 0, nb = 0;
    int passCount=0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nb = intree->GetEntry(jentry); nbytes += nb;
        passCount=0;
        RecoNuEnergy=TrueNuEnergy;
        
        if(jentry%5000==0) printf("%d/%d complete\n",jentry, nentries);
        
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
