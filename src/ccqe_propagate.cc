#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom3.h>

#include "FluxParameters.hh"
#include "FluxParameters_norm.hh"
#include "DetParameters.hh"
#include "FSIParameters.hh"
#include "NuclFSIParameters.hh"
#include "XsecParameters.hh"
#include "CCQEParameters.hh"
#include "XsecFitter.hh"
#include "ccqeTreeMC.hh"
#include "CCQESample.hh"


using namespace std;

int main(int argc, char *argv[])
{
  string fsel     = "../inputs/CCQE_neut.root";
  string fdetcov = "../inputs/cov_matrix_old/cov_matrix_det_newbin.root";  //USING NEW BINNING DET MATRIX
  string fccqebin = "../inputs/ccqebins.txt";
  string fnameout = "ccqe_fiterrors.root";
  double potD     = 60;   //in units of 10^19
  double potMC    = 389.5;//204.95; //in units 10^19 
  int ntoys        = 100;
  int seed        = 1019;

  //get command line options
  char cc;
  while((cc = getopt(argc, argv, "i:e:b:o:n:N:s:h:f:")) != -1)
    {
      switch(cc)
        {
        case 'i': //selected events
	  fsel = optarg;
          break; 
        case 'd': //det covariance matrix
          fdetcov = optarg;
          break;
	case 'b': //binning for xsec weights
	  fccqebin = optarg;
	  break;
	case 'o': //output file
	  fnameout = optarg;
	  break;
	case 'n': //data POT
	  potD = atof(optarg);
	  break;
	case 'N': //MC POT
	  potMC = atof(optarg);
	  break;
	case 't': //number of toys
	  ntoys = atof(optarg);
	  break;
	case 's': //random seed
	  seed = atoi(optarg);
	  break;
	case 'h': //help 
          std::cout << "USAGE: " 
                    << argv[0] << " OPTIONS:" << std::endl
                    << "-i : \tset event selection" << std::endl
                    << "-e : \tset prefit covariances" << std::endl
		    << "-d : \tset detector covariances" << std::endl
		    << "-b : \tset ccqe xsec bin definitions" << std::endl
  		    << "-o : \tset name of output file" << std::endl
		    << "-n : \tset POT for data in units 10**19" << std::endl
		    << "-N : \tset POT for MC in units 10**19" << std::endl
		    << "-s : \tset random seed" << std::endl;
          return 0;
          break;
	default:
	  return 1;
	}
    }

 
  /********************************COV.MATRIX FOR DETECTOR SYSTEMATICS*********************/
  //setup cth,p bins, starting param values and covm for det syst --------------
  TFile *findetcov = TFile::Open(fdetcov.c_str()); //contains det. systematics info
  TVectorD* det_weights = (TVectorD*)findetcov->Get("det_weights");
 
  TMatrixDSym *cov_det_in   = (TMatrixDSym*)findetcov->Get("cov_det");
  TMatrixDSym cov_det(cov_det_in->GetNrows());
  for(size_t m=0; m<cov_det_in->GetNrows(); m++){
    for(size_t k=0; k<cov_det_in->GetNrows(); k++){
      cov_det(m, k) = (*cov_det_in)(m,k);
    }
  }

  findetcov->Close();

  /***********************************************************************************/

  /********************************XSEC RESPONSE FUNCTIONS****************************/
  //Signal modeling has ~0 impact on results -> commented out
  vector<TFile*> responsefunctions;
  //TFile* MAQErespfunc = new TFile("../inputs/responsefunc/MAQEshape_all_variation_shapeNorm.root");
  //responsefunctions.push_back(MAQErespfunc);
  //TFile* PFrespfunc = new TFile("../inputs/responsefunc/PF_all_variation.root");
  //responsefunctions.push_back(PFrespfunc);
  //TFile* SFrespfunc = new TFile("../inputs/responsefunc/SF_all_variation.root");
  //responsefunctions.push_back(SFrespfunc);
  TFile* MAResrespfunc = new TFile("../inputs/responsefunc/MARes_all_variation.root");
  responsefunctions.push_back(MAResrespfunc);
  TFile* CCOthrespfunc = new TFile("../inputs/responsefunc/CCOther_all_variation.root");
  responsefunctions.push_back(CCOthrespfunc);
  TFile* CC1piE0respfunc = new TFile("../inputs/responsefunc/CC1piE0_all_variation.root");
  responsefunctions.push_back(CC1piE0respfunc);
  TFile* CC1piE1respfunc = new TFile("../inputs/responsefunc/CC1piE1_all_variation.root");
  responsefunctions.push_back(CC1piE1respfunc);
  TFile* NC1pi0E0respfunc = new TFile("../inputs/responsefunc/NC1pi0E0_all_variation.root");
  responsefunctions.push_back(NC1pi0E0respfunc);

  TMatrixDSym cov_xsec(5);
  cov_xsec(0,0) = 0.03344; //MARes
  cov_xsec(0,1) = 0;
  cov_xsec(0,2) = -0.01549;//MARes-CC1piE0
  cov_xsec(0,3) = 0;
  cov_xsec(0,4) = -0.01841;//MARes-NC1pi0
  cov_xsec(1,1) = 0.16; //CCother
  cov_xsec(1,2) = 0;
  cov_xsec(1,3) = 0;
  cov_xsec(1,4) = 0;
  cov_xsec(2,2) = 0.1003; //CC1piE0
  cov_xsec(2,3) = 0;
  cov_xsec(2,4) = 0.07703; //CC1piE0-NC1pi0
  cov_xsec(3,3) = 0.16; //CC1piE1
  cov_xsec(3,4) = 0;
  cov_xsec(4,4) = 0.1075; //NC1pi0

  /*****************************************************************************************/


  /************************FSI RESPONSE FUNCTIONS*****************************************/
  vector<TFile*> responsefunctions_FSI;
  TFile* FrInelLowrespfunc = new TFile("../inputs/responsefunc/FrInelLow_all_variation.root");
  responsefunctions_FSI.push_back(FrInelLowrespfunc);
  TFile* FrInelHighrespfunc = new TFile("../inputs/responsefunc/FrInelHigh_all_variation.root");
  responsefunctions_FSI.push_back(FrInelHighrespfunc);
  TFile* FrPiProdrespfunc = new TFile("../inputs/responsefunc/FrPiProd_all_variation.root");
  responsefunctions_FSI.push_back(FrPiProdrespfunc);
  TFile* FrAbsrespfunc = new TFile("../inputs/responsefunc/FrAbs_all_variation.root");
  responsefunctions_FSI.push_back(FrAbsrespfunc);
  TFile* FrCExLowrespfunc = new TFile("../inputs/responsefunc/FrCExLow_all_variation.root");
  responsefunctions_FSI.push_back(FrCExLowrespfunc);
  TFile* FrCExHighrespfunc = new TFile("../inputs/responsefunc/FrCExHigh_all_variation.root");
  responsefunctions_FSI.push_back(FrCExHighrespfunc);

  //FSI parameters are correlated (new from iRODS file)
  TMatrixDSym cov_fsi(6);
  cov_fsi(0,0)= 0.17;
  cov_fsi(0,1)= -0.002778;
  cov_fsi(0,2)= 0;
  cov_fsi(0,3)= 0.02273;
  cov_fsi(0,4)= 0.005;
  cov_fsi(0,5)= 0;
  cov_fsi(1,1)= 0.1142;
  cov_fsi(1,2)= -0.1667;
  cov_fsi(1,3)= -0.001263;
  cov_fsi(1,4)= -0.002083;
  cov_fsi(1,5)=  -0.09259;
  cov_fsi(2,2)= 0.25;
  cov_fsi(2,3)= -5.204e-18;
  cov_fsi(2,4)= 0;
  cov_fsi(2,5)= 0.1389;
  cov_fsi(3,3)= 0.1694;
  cov_fsi(3,4)= -0.002273;
  cov_fsi(3,5)= -3.469e-18;
  cov_fsi(4,4)= 0.3213;
  cov_fsi(4,5)= 1.735e-18;
  cov_fsi(5,5)= 0.07716;
  // TMatrixDSym cov_fsi(6);
  // cov_fsi(0,0)= 0.412*0.412;
  // cov_fsi(0,1)= -(0.053*0.053);
  // cov_fsi(0,2)= 0;
  // cov_fsi(0,3)= 0.151*0.151;
  // cov_fsi(0,4)= 0.071*0.071;
  // cov_fsi(0,5)= 0;
  // cov_fsi(1,1)= 0.338*0.338;
  // cov_fsi(1,2)= -(0.408*0.408);
  // cov_fsi(1,3)= -(0.036*0.036);
  // cov_fsi(1,4)= -(0.046*0.046);
  // cov_fsi(1,5)= -(0.304*0.304);
  // cov_fsi(2,2)= 0.500*0.500;
  // cov_fsi(2,3)= 0;
  // cov_fsi(2,4)= 0;
  // cov_fsi(2,5)= 0.373*0.373;
  // cov_fsi(3,3)= 0.412*0.412;
  // cov_fsi(3,4)= -(0.048*0.048);
  // cov_fsi(3,5)= 0;
  // cov_fsi(4,4)= 0.567*0.567;
  // cov_fsi(4,5)= 0;
  // cov_fsi(5,5)= 0.278*0.278;
  
  // //Need to add 0.001 to the diagonal, otherwise not positive definite
   cov_fsi(0,0)= cov_fsi(0,0)+0.001;
   cov_fsi(1,1)= cov_fsi(1,1)+0.001;
   cov_fsi(2,2)= cov_fsi(2,2)+0.001;
   cov_fsi(3,3)= cov_fsi(3,3)+0.001;
   cov_fsi(4,4)= cov_fsi(4,4)+0.001;
   cov_fsi(5,5)= cov_fsi(5,5)+0.001;
  
  /*********************************************************************************/
  
  /************************Nucleon FSI RESPONSE FUNCTIONS*****************************/
  vector<TFile*> responsefunctions_NuclFSI;
  TFile* NMFPrespfunc = new TFile("../inputs/responsefunc/N_MFP_all_variation.root");
  responsefunctions_NuclFSI.push_back(NMFPrespfunc);
  TFile* NFrElasrespfunc = new TFile("../inputs/responsefunc/N_FrElas_all_variation.root");
  responsefunctions_NuclFSI.push_back(NFrElasrespfunc);
  TFile* NFrAbsrespfunc = new TFile("../inputs/responsefunc/N_FrAbs_all_variation.root");
  responsefunctions_NuclFSI.push_back(NFrAbsrespfunc);

  //NuclFSI parameters are not correlated for now 
  TMatrixDSym cov_nuclfsi(3);
  cov_nuclfsi(0,0)= 0.20*0.20;
  cov_nuclfsi(0,1)= 0;
  cov_nuclfsi(0,2)= 0;
  cov_nuclfsi(1,1)= 0.3*0.3;
  cov_nuclfsi(1,2)= 0;
  cov_nuclfsi(2,2)= 0.20*0.20;
 
  /*********************************************************************************/
  

   TFile *fdata = new TFile(TString("../inputs/CCQE_neut.root")); 
   TTree *tdata = (TTree*)(fdata->Get("selectedEvents"));

   std::vector<std::pair<double, double> > v_pedges;
   std::vector<std::pair<double, double> > v_cthedges;
   ifstream fin(fccqebin.c_str());
   assert(fin.is_open());
   string line;
   while (getline(fin, line))
     {
       stringstream ss(line);
       double p1, p2, cth1, cth2;
       if(!(ss>>cth1>>cth2>>p1>>p2))
	 {
	   cerr<<"Bad line format: "<<endl
	       <<"     "<<line<<endl;
	   continue;
	 }
       v_pedges.push_back(make_pair(p1,p2));
       v_cthedges.push_back(make_pair(cth1,cth2));
     }
   fin.close();

  TFile *fout = TFile::Open(fnameout.c_str(), "RECREATE");
  cout<<"file open"<<endl;
  vector<AnaSample*> samples;
  CCQESample sam1(0, "MuTrack",v_pedges, v_cthedges,tdata);
  sam1.SetNorm(potD/potMC);
  samples.push_back(&sam1);

  CCQESample sam2(1, "MuPtpc",v_pedges, v_cthedges,tdata);
  sam2.SetNorm(potD/potMC);
  samples.push_back(&sam2);
  
  CCQESample sam3(2, "MuPfgd",v_pedges, v_cthedges,tdata);
  sam3.SetNorm(potD/potMC);
  samples.push_back(&sam3);
  
  CCQESample sam4(3, "MufgdPtpc",v_pedges, v_cthedges,tdata);
  sam4.SetNorm(potD/potMC);
  samples.push_back(&sam4);

  CCQESample sam5(4, "CC1pi",v_pedges, v_cthedges,tdata);
  sam5.SetNorm(potD/potMC);
  samples.push_back(&sam5);
  
  CCQESample sam6(5, "DIS",v_pedges, v_cthedges,tdata);
  sam6.SetNorm(potD/potMC);
  samples.push_back(&sam6);

  //--
  //read MC events
  ccqeTreeMC selTree(fsel.c_str());
  selTree.GetEvents(samples);
  //get brakdown by reaction
  for(size_t s=0;s<samples.size();s++){
    ((CCQESample*)(samples[s]))->GetSampleBreakdown(fout,"nominal");
  }

  
  //CCQE parameters
  CCQEParameters ccqepara(fccqebin.c_str());
  ccqepara.InitEventMap(samples);

  //Det parameters
  DetParameters detpara("../inputs/ccqebins_det.txt",det_weights);
  //DetParameters detpara("../inputs/ccqebins.txt",det_weights);
  detpara.InitEventMap(samples);

  //Xsec parameters
  XsecParameters xsecpara;
  xsecpara.StoreResponseFunctions(responsefunctions, v_pedges, v_cthedges);
  xsecpara.InitEventMap(samples);

  //FSI parameters
  FSIParameters fsipara;
  fsipara.StoreResponseFunctions(responsefunctions_FSI, v_pedges, v_cthedges);
  fsipara.InitEventMap(samples);

  //Nucleon FSI parameters
  NuclFSIParameters nuclfsipara;
  nuclfsipara.StoreResponseFunctions(responsefunctions_NuclFSI, v_pedges, v_cthedges);
  nuclfsipara.InitEventMap(samples);

  vector<double> allparams;
  TFile *finput= new TFile("SaveCovMatrix/ccqe_results_13.root","READ");
  TMatrixDSym *covariance   = (TMatrixDSym*)finput->Get("res_cov_matrix");
  TVectorD *priorVec = ((TVectorD*)finput->Get("res_vector"));
  if(covariance->GetNrows() != priorVec->GetNrows() ){
    cout<<"Error dimension vector vs matrix"<<endl;
    abort();
  }
  for(uint i=0; i<covariance->GetNrows(); i++){
    allparams.push_back((*priorVec)[i]);
    cout<<i<<" "<<((*priorVec)[i])<<endl;
  }

  ThrowParms *throwParms = new ThrowParms((*priorVec),(*covariance));
  TRandom3 *rand = new TRandom3(seed);
  //set gRandom to our rand
  gRandom = rand;

  //to be redone for each toy!////////////////////////////////////////
  for(int t=0; t<ntoys; t++){
    vector<double> allparams_throw;
    throwParms->ThrowSet(allparams_throw);
    vector<double>::const_iterator first = allparams_throw.begin() + 0;
    vector<double>::const_iterator last =  allparams_throw.begin() + ccqepara.Npar;
    vector <double> ccqeparams_throw (first,last);
    first = allparams_throw.begin() + ccqepara.Npar;
    last =  allparams_throw.begin() + ccqepara.Npar + detpara.Npar;
    vector <double> detparams_throw (first,last);
    first = allparams_throw.begin() + ccqepara.Npar + detpara.Npar;
    last =  allparams_throw.begin() + ccqepara.Npar + detpara.Npar + xsecpara.Npar;
    vector <double> xsecparams_throw (first,last);
    first = allparams_throw.begin() + ccqepara.Npar + detpara.Npar + xsecpara.Npar;
    last =  allparams_throw.begin() + ccqepara.Npar + detpara.Npar + xsecpara.Npar + fsipara.Npar;
    vector <double> fsiparams_throw (first,last);
    first = allparams_throw.begin() + ccqepara.Npar + detpara.Npar + xsecpara.Npar + fsipara.Npar;
    last =  allparams_throw.begin() + ccqepara.Npar + detpara.Npar + xsecpara.Npar + fsipara.Npar + nuclfsipara.Npar;
    vector <double> nuclfsiparams_throw (first,last);
 
    cout<<"CHECK "<<endl;
    cout<<"ccqe size "<<ccqeparams_throw.size()<<endl;
    for(int j=0; j<ccqeparams_throw.size(); j++){
      cout<<j<<" "<<ccqeparams_throw[j]<<endl;
    }
    cout<<"det size "<<detparams_throw.size()<<endl;
    for(int j=0; j<detparams_throw.size(); j++){
      cout<<detparams_throw[j]<<endl;
    }
    cout<<"xsec size "<<xsecparams_throw.size()<<endl;
    for(int j=0; j<xsecparams_throw.size(); j++){
      cout<<xsecparams_throw[j]<<endl;
    }
    cout<<"fsi size "<<fsiparams_throw.size()<<endl;
    for(int j=0; j<fsiparams_throw.size(); j++){
      cout<<fsiparams_throw[j]<<endl;
    }
    cout<<"nuclfsi size "<<nuclfsiparams_throw.size()<<endl;
    for(int j=0; j<nuclfsiparams_throw.size(); j++){
      cout<<nuclfsiparams_throw[j]<<endl;
    }

    for(size_t s=0;s<samples.size();s++)
      {
	//loop over events
	for(int i=0;i<samples[s]->GetN();i++)
	  {
	    AnaEvent* ev = samples[s]->GetEvent(i);
	    ev->SetEvWght(ev->GetEvWghtMC()); 
	    //do weights for each AnaFitParameters obj
	    ccqepara.ReWeight(ev, s, i, ccqeparams_throw);
	    detpara.ReWeight(ev, s, i, detparams_throw);
	    xsecpara.ReWeight(ev, s, i, xsecparams_throw);
	    fsipara.ReWeight(ev, s, i, fsiparams_throw);
	    nuclfsipara.ReWeight(ev, s, i, nuclfsiparams_throw);
	  }
      }

    for(size_t s=0;s<samples.size();s++){
      ((CCQESample*)(samples[s]))->GetSampleBreakdown(fout,Form("toy_%d",t));
    }
    
    fout->cd();
    TVectorD throwVec(ccqeparams_throw.size());
    for(uint i=0; i<ccqeparams_throw.size(); i++)
      throwVec(i) = ccqeparams_throw[i];
    throwVec.Write(Form("ccqeparam_toy_%d",t));

    throwVec.ResizeTo(detparams_throw.size());
    for(uint i=0; i<detparams_throw.size(); i++)
      throwVec(i) = detparams_throw[i];
    throwVec.Write(Form("detparam_toy_%d",t));

    throwVec.ResizeTo(xsecparams_throw.size());
    for(uint i=0; i<xsecparams_throw.size(); i++)
      throwVec(i) = xsecparams_throw[i];
    throwVec.Write(Form("xsecparam_toy_%d",t));

    throwVec.ResizeTo(fsiparams_throw.size());
    for(uint i=0; i<fsiparams_throw.size(); i++)
      throwVec(i) = fsiparams_throw[i];
    throwVec.Write(Form("fsiparam_toy_%d",t));

    throwVec.ResizeTo(nuclfsiparams_throw.size());
    for(uint i=0; i<nuclfsiparams_throw.size(); i++)
      throwVec(i) = nuclfsiparams_throw[i];
    throwVec.Write(Form("nuclfsiparam_toy_%d",t));

    }
   ///////////////////////////////////////////////

  return 0;
}
