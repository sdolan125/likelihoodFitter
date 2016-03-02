#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <assert.h>

#include <TCanvas.h>
#include <TH1F.h>

#include "FluxParameters.hh"
#include "FluxParameters_norm.hh"
#include "DetParameters.hh"
#include "FSIParameters.hh"
#include "NuclFSIParameters.hh"
#include "XsecParameters.hh"
//#include "CCQEParameters.hh"
#include "FitParameters.hh"
#include "XsecFitter.hh"
#include "anyTreeMC.hh"
#include "AnySample.hh"


using namespace std;

int main(int argc, char *argv[])
{
  //string fsel     = "../inputs/CCQE_neut_h2.root"; 
  //string fsel     = "../inputs/ConvertedCC0PiV2.root";
  //string fsel     = "../inputs/CC0Pi_6B_dpTvsMuMom_NeutAir.root";
  string fsel     = "../inputs/CC0Pi_6B_dpTvsMuMom_kineVars_NeutAir.root";
  //string fsel     = "../inputs/CCQE_neut_h2_sigOppCrazy.root";
  //string fsel     = "../inputs/CCQE_neut_h2_sigCrazy.root";
  //string fsel     = "../inputs/CCQE_neut_h2_pionversionBiased.root";
  //string fsel     = "../inputs/CCQE_nuwro_h2.root";
  //string fakeData     = "../inputs/CCQE_neut_h2_sigOppCrazy.root";
  //string fakeData     = "../inputs/CCQE_nuwro_h2.root";
  //string fakeData     = "../inputs/CCQE_genie_h2.root"; 
  //string fakeData     = "../inputs/CCQE_data_h2.root";
  //string fakeData     = "../inputs/ConvertedCC0PiV2.root";
  //string fakeData     = "../inputs/CC0Pi_6B_dpTvsMuMom_NeutWater.root";
  //string fakeData      = "../inputs/CC0Pi_6B_dpTvsMuMom_Genie.root";
  string fakeData      = "../inputs/CC0Pi_6B_dpTvsMuMom_kineVars_Genie.root";
  //string fakeData      = "../inputs/CC0Pi_6B_dpTvsMuMom_kineVars_NeutAir.root";
  //string fakeData      = "../inputs/CC0Pi_6B_dpTvsMuMom_kineVars_NeutWater.root";
  //string fakeData     = "../inputs/CCQE_neut_h2_pionversionBiased.root";
  //string fakeData       = "../inputs/CCQE_neut_h2_bkgBias.root";
  //string fakeData       = "../inputs/CCQE_neut_h2.root";
  //string ffluxcov     = "../inputs/flux_covariance_banff_run1to4_v1.root";  
  string ffluxcov     = "../inputs/flux_covariance_banff_13av1.1.root";  
  string fdetcov_fine = "../inputs/cov_matrix_det_finebins.root";  //USING NEW BINNING DET MATRIX 
  string fdetcov = "../inputs/cov_matrix_det_averagebins.root";  //USING NEW BINNING DET MATRIX 
  string fcovFSI     = "../inputs/fromBANFF/fsi_cov.root";   //not used anymore
  //string fccqebin = "../inputs/ccqebins_Dec2014.txt";
  //string fccqebin = "../inputs/dptbins_V1.txt";
  string fccqebin = "../inputs/dptbins_coarseV6_8binRec.txt";
  //string fccqebin = "../inputs/OneBin.txt";
  //string fccqebin = "../inputs/ccqebins_Andy.txt";
  string fnameout = "ccqe_fitresults.root";
  //double potD     = 60;   //in units of 10^19
  //double potMC    = 389.5;//204.95; //in units 10^19 
  double potD     = 57.34;   //in units of 10^19
  double potMC    = 57.34; //in units 10^19 
  int seed        = 1019;
  double regparam = 0.0;

  //get command line options
  char cc;
  while((cc = getopt(argc, argv, "i:e:b:o:n:N:s:h:f:r:")) != -1)
    {
      switch(cc)
        {
        case 'i': //selected events
	  fsel = optarg;
          break; 
        case 'e': //flux covariance matrix
          ffluxcov = optarg;
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
	case 's': //random seed
	  seed = atoi(optarg);
	  break;
	case 'f': //fake data file
	  fakeData = optarg;
	  break;
  case 'r': //regularisation param
    regparam = atof(optarg);
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
        << "-s : \tset random seed" << std::endl
        << "-r : \tset regularisation parameter" << std::endl;
          return 0;
          break;
	default:
	  return 1;
	}
    }

  //get info on the priors 
  TFile *finfluxcov = TFile::Open(ffluxcov.c_str()); //contains flux and det. systematics info
  TFile *findetcov = TFile::Open(fdetcov.c_str()); //contains flux and det. systematics info
  TFile *findetcov_fine = TFile::Open(fdetcov_fine.c_str()); //contains flux and det. systematics info


  /*************************************** FLUX *****************************************/
  
  //setup enu bins and covm for flux 
  //TAxis *nd_numu_bins = (TAxis*)finfluxcov->Get("nd_numu_bins");
  TAxis *nd_numu_bins = (TAxis*)finfluxcov->Get("nd5_numode_numu_bins");
  //TMatrixDSym *covm   = (TMatrixDSym*)fincov->Get("prefit_cov");
  TMatrixDSym *cov_flux_in   = (TMatrixDSym*)finfluxcov->Get("total_flux_cov");
  TMatrixDSym cov_flux(nd_numu_bins->GetNbins());
  vector<double> enubins;
  enubins.push_back(nd_numu_bins->GetBinLowEdge(1));
   for(int i=0;i<nd_numu_bins->GetNbins();i++)
    {
      enubins.push_back(nd_numu_bins->GetBinUpEdge(i+1));
      for(int j=0;j<nd_numu_bins->GetNbins();j++)
	{
	  cov_flux(i, j) = (*cov_flux_in)(i,j);
	}
    }
  
    finfluxcov->Close();
  /*****************************************************************************************/

  /********************************COV.MATRIX FOR DETECTOR SYSTEMATICS*********************/
  //setup D1,D2 bins, starting param values and covm for det syst --------------
  TVectorD* det_weights = (TVectorD*)findetcov->Get("det_weights");
 
  TMatrixDSym *cov_det_in   = (TMatrixDSym*)findetcov->Get("cov_det");
  TMatrixDSym cov_det(cov_det_in->GetNrows());
  for(size_t m=0; m<cov_det_in->GetNrows(); m++){
    for(size_t k=0; k<cov_det_in->GetNrows(); k++){
      cov_det(m, k) = (*cov_det_in)(m,k);
    }
  }

  findetcov->Close();

  TVectorD* det_weights_fine = (TVectorD*)findetcov_fine->Get("det_weights");
  TMatrixDSym *cov_det_in_fine   = (TMatrixDSym*)findetcov_fine->Get("cov_det");
  TMatrixDSym cov_det_fine(cov_det_in_fine->GetNrows());
  for(size_t m=0; m<cov_det_in_fine->GetNrows(); m++){
    for(size_t k=0; k<cov_det_in_fine->GetNrows(); k++){
      cov_det_fine(m, k) = (*cov_det_in_fine)(m,k);
    }
  }
  
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
  /*  TFile* MAResrespfunc = new TFile("../inputs/responsefunc/MARes_all_variation.root");
  responsefunctions.push_back(MAResrespfunc);
  TFile* CCOthrespfunc = new TFile("../inputs/responsefunc/CCOther_all_variation.root");
  responsefunctions.push_back(CCOthrespfunc);
  TFile* PilessDcyrespfunc = new TFile("../inputs/responsefunc/PilessDcy_all_variation.root");
  responsefunctions.push_back(PilessDcyrespfunc);
  TFile* CC1piE0respfunc = new TFile("../inputs/responsefunc/CC1piE0_all_variation.root");
  responsefunctions.push_back(CC1piE0respfunc);
  TFile* CC1piE1respfunc = new TFile("../inputs/responsefunc/CC1piE1_all_variation.root");
  responsefunctions.push_back(CC1piE1respfunc);
  TFile* CCCohE0respfunc = new TFile("../inputs/responsefunc/CCCohE0_all_variation.root");
  responsefunctions.push_back(CCCohE0respfunc);
  TFile* NCOthrespfunc = new TFile("../inputs/responsefunc/NCOtherE0_all_variation.root");
  responsefunctions.push_back(NCOthrespfunc);
  TFile* NC1pi0E0respfunc = new TFile("../inputs/responsefunc/NC1pi0E0_all_variation.root");
  responsefunctions.push_back(NC1pi0E0respfunc);
  TFile* NC1piE0respfunc = new TFile("../inputs/responsefunc/NC1piE0_all_variation.root");
  responsefunctions.push_back(NC1piE0respfunc);
  */
  //TMatrixDSym cov_xsec(5);
  //cov_xsec(0,0) = 0.03344; //MARes
  //cov_xsec(0,1) = 0;
  //cov_xsec(0,2) = -0.01549;//MARes-CC1piE0
  //cov_xsec(0,3) = 0;
  //cov_xsec(0,4) = -0.01841;//MARes-NC1pi0
  //cov_xsec(1,1) = 0.16; //CCother
  //cov_xsec(1,2) = 0;
  //cov_xsec(1,3) = 0;
  //cov_xsec(1,4) = 0;
  //cov_xsec(2,2) = 0.1003; //CC1piE0
  //cov_xsec(2,3) = 0;
  //cov_xsec(2,4) = 0.07703; //CC1piE0-NC1pi0
  //cov_xsec(3,3) = 0.16; //CC1piE1
  //cov_xsec(3,4) = 0;
  //cov_xsec(4,4) = 0.1075; //NC1pi0
 TMatrixDSym cov_xsec(9);
  cov_xsec(0,0) = 0.03344; //MARes
  cov_xsec(0,1) = 0;
  cov_xsec(0,2) = 0;
  cov_xsec(0,3) = -0.01549;//MARes-CC1piE0
  cov_xsec(0,4) = 0;
  cov_xsec(0,5) = 0;
  cov_xsec(0,6) = 0;
  cov_xsec(0,7) = -0.01841;//MARes-NC1pi0
  cov_xsec(0,8) = 0;

  cov_xsec(1,1) = 0.16; //CCother
  cov_xsec(1,2) = 0;
  cov_xsec(1,3) = 0;
  cov_xsec(1,4) = 0;
  cov_xsec(1,5) = 0;
  cov_xsec(1,6) = 0;
  cov_xsec(1,7) = 0;
  cov_xsec(1,8) = 0;

  cov_xsec(2,2) = 0.04; //PilessDcy
  cov_xsec(2,3) = 0;
  cov_xsec(2,4) = 0;
  cov_xsec(2,5) = 0;
  cov_xsec(2,6) = 0;
  cov_xsec(2,7) = 0;
  cov_xsec(2,8) = 0;

  cov_xsec(3,3) = 0.1003; //CC1piE0
  cov_xsec(3,4) = 0;
  cov_xsec(3,5) = 0;
  cov_xsec(3,6) = 0;
  cov_xsec(3,7) = 0.07703; //CC1piE0-NC1pi0
  cov_xsec(3,8) = 0;

  cov_xsec(4,4) = 0.16; //CC1piE1
  cov_xsec(4,5) = 0;
  cov_xsec(4,6) = 0;
  cov_xsec(4,7) = 0;
  cov_xsec(4,8) = 0;

  cov_xsec(5,5) = 1; //CCCoh
  cov_xsec(5,6) = 0;
  cov_xsec(5,7) = 0;
  cov_xsec(5,8) = 0;

  cov_xsec(6,6) = 0.09; //NCOther
  cov_xsec(6,7) = 0;
  cov_xsec(6,8) = 0;
 
  cov_xsec(7,7) = 0.1075; //NC1pi0
  cov_xsec(7,8) = 0;

  cov_xsec(8,8) = 0.09; //NC1pi (not present in 2013 BANFF matrix, added)
          
  /*****************************************************************************************/


  /************************FSI RESPONSE FUNCTIONS*****************************************/
  vector<TFile*> responsefunctions_FSI;
  /*TFile* FrInelLowrespfunc = new TFile("../inputs/responsefunc/FrInelLow_all_variation.root");
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
  */
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
  /*TFile* NMFPrespfunc = new TFile("../inputs/responsefunc/N_MFP_all_variation.root");
  responsefunctions_NuclFSI.push_back(NMFPrespfunc);
  TFile* NFrElasrespfunc = new TFile("../inputs/responsefunc/N_FrElas_all_variation.root");
  responsefunctions_NuclFSI.push_back(NFrElasrespfunc);
  TFile* NFrAbsrespfunc = new TFile("../inputs/responsefunc/N_FrAbs_all_variation.root");
  responsefunctions_NuclFSI.push_back(NFrAbsrespfunc);*/

  //NuclFSI parameters are not correlated for now 
  TMatrixDSym cov_nuclfsi(3);
  cov_nuclfsi(0,0)= 0.20*0.20;
  cov_nuclfsi(0,1)= 0;
  cov_nuclfsi(0,2)= 0;
  cov_nuclfsi(1,1)= 0.3*0.3;
  cov_nuclfsi(1,2)= 0;
  cov_nuclfsi(2,2)= 0.20*0.20;
  
  /*********************************************************************************/
  

   TFile *fdata = new TFile(TString(fakeData)); 
   TTree *tdata = (TTree*)(fdata->Get("selectedEvents"));

   std::vector<std::pair<double, double> > v_D1edges;
   std::vector<std::pair<double, double> > v_D2edges;
   ifstream fin(fccqebin.c_str());
   assert(fin.is_open());
   string line;
   while (getline(fin, line))
     {
       stringstream ss(line);
       double D1_1, D1_2, D2_1, D2_2;
       if(!(ss>>D2_1>>D2_2>>D1_1>>D1_2))
	 {
	   cerr<<"Bad line format: "<<endl
	       <<"     "<<line<<endl;
	   continue;
	 }
       v_D1edges.push_back(make_pair(D1_1,D1_2));
       v_D2edges.push_back(make_pair(D2_1,D2_2));
     }
   fin.close();

  TFile *fout = TFile::Open(fnameout.c_str(), "RECREATE");
  cout<<"file open"<<endl;
  /*Add this for data too (true variables should be just empty....)
    and not by reaction!!!!!!!!!!!!!!!!!*/
  //Fake samples from external source (eg from Genie)
  /*vector<AnaSample*> fake_samples;
  //double potMC_genie=384.762;
  //double potMC_genie=389.5; //neut
  //double potMC_genie=380.0; //nuwro
  double potMC_genie=57.34; //neut
  CCQESample fake_sam1(0, "genie_MuTrack",v_D1edges, v_D2edges,tdata);
  fake_sam1.SetNorm(potD/potMC_genie);
  fake_samples.push_back(&fake_sam1);

  CCQESample fake_sam2(1, "genie_MuPtpc",v_D1edges, v_D2edges,tdata);
  fake_sam2.SetNorm(potD/potMC_genie);
  fake_samples.push_back(&fake_sam2);
  
  CCQESample fake_sam3(2, "genie_MuPfgd",v_D1edges, v_D2edges,tdata);
  fake_sam3.SetNorm(potD/potMC_genie);
  fake_samples.push_back(&fake_sam3);
  
  CCQESample fake_sam4(3, "genie_MufgdPtpc",v_D1edges, v_D2edges,tdata);
  fake_sam4.SetNorm(potD/potMC_genie);
  fake_samples.push_back(&fake_sam4);
  
  CCQESample fake_sam5(4, "genie_CC1pi",v_D1edges, v_D2edges,tdata);
  fake_sam5.SetNorm(potD/potMC_genie);
  fake_samples.push_back(&fake_sam5);
  
  CCQESample fake_sam6(5, "genie_DIS",v_D1edges, v_D2edges,tdata);
  fake_sam6.SetNorm(potD/potMC_genie);
  fake_samples.push_back(&fake_sam6);

  //read MC fake events
  ccqeTreeMC selTreeFake(fakeData.c_str());
  selTreeFake.GetEvents(fake_samples);
  for(size_t s=0;s<fake_samples.size();s++)
    {
      ((CCQESample*)fake_samples[s])->GetSampleBreakdown(fout,"fake");
    }
  
  cout<<"fake data done"<<endl;
  */  
  
  vector<AnaSample*> samples;
  /*
  AnySample sam1(0, "MuTPC",v_D1edges, v_D2edges,tdata);
  sam1.SetNorm(potD/potMC);
  samples.push_back(&sam1);
  */
  AnySample sam2(1, "MuTPCpTPC",v_D1edges, v_D2edges,tdata);
  sam2.SetNorm(potD/potMC);
  samples.push_back(&sam2);
  
  AnySample sam3(2, "MuTPCpFGD",v_D1edges, v_D2edges,tdata);
  sam3.SetNorm(potD/potMC);
  samples.push_back(&sam3);
  
  AnySample sam4(3, "MuFGDPTPC",v_D1edges, v_D2edges,tdata);
  sam4.SetNorm(potD/potMC);
  samples.push_back(&sam4);

  /*
  AnySample sam5(4, "MuFGD",v_D1edges, v_D2edges,tdata);
  sam5.SetNorm(potD/potMC);
  samples.push_back(&sam5);
  
  AnySample sam6(5, "CC1pi",v_D1edges, v_D2edges,tdata);
  sam6.SetNorm(potD/potMC);
  samples.push_back(&sam6);
    
  AnySample sam7(6, "DIS",v_D1edges, v_D2edges,tdata);
  sam7.SetNorm(potD/potMC);
  samples.push_back(&sam7);
  */
  
  AnySample sam8(7, "muTPC_Np",v_D1edges, v_D2edges,tdata);
  sam8.SetNorm(potD/potMC);
  samples.push_back(&sam8);

  AnySample sam9(8, "muFGDpTPC_Np",v_D1edges, v_D2edges,tdata);
  sam9.SetNorm(potD/potMC);
  samples.push_back(&sam9);

  AnySample sam10(9, "muFGD_Np",v_D1edges, v_D2edges,tdata);
  sam10.SetNorm(potD/potMC);
  samples.push_back(&sam10);
  
  //--
  //read MC events
  anyTreeMC selTree(fsel.c_str());
  cout << "Reading and collecting events" << endl;
  selTree.GetEvents(samples);
  //get brakdown by reaction
  cout << "Getting sample breakdown by reaction" << endl;
  for(size_t s=0;s<samples.size();s++){
    ((AnySample*)(samples[s]))->GetSampleBreakdown(fout,"nominal");
  }

  cout<<"nominal data done"<<endl;

  //define fit param classes
  vector<AnaFitParameters*> fitpara;
  //fitpara.SetFluxHisto(h_flux);
  
  //CCQE parameters
  FitParameters ccqepara(fccqebin.c_str());
  ccqepara.InitEventMap(samples);
  fitpara.push_back(&ccqepara);
  
  cout<<"ccqe parameters done"<<endl;

  //Flux parameters
  /*FluxParameters fluxpara(enubins);
  fluxpara.SetCovarianceMatrix(&cov_flux);
  fluxpara.InitEventMap(samples);
  fitpara.push_back(&fluxpara);
  */
  /*
    //Flux shape
  FluxParameters fluxpara(enubins,"par_flux_shape");
  fluxpara.SetCovarianceMatrix(&cov_flux);
  TFile *fflux= new TFile("../inputs/c_flux.root","READ");
  TH1F *flux=(TH1F*)(((TCanvas*)(fflux->Get("c_flux")))
		     ->GetPrimitive("flux"));
  fluxpara.SetFluxHisto(flux);
  fluxpara.InitEventMap(samples);
  fitpara.push_back(&fluxpara);
  */
  //Flux parameters normalization only
  /*vector<double> enubins_norm;
  enubins_norm.push_back(0);
  enubins_norm.push_back(9999999);
  FluxParameters fluxpara_norm(enubins_norm,"par_flux_norm");
  TMatrixDSym cov_flux_norm(1);
  (cov_flux_norm)(0,0) = (0.11*0.11);
  fluxpara_norm.SetCovarianceMatrix(&cov_flux_norm);
  fluxpara_norm.InitEventMap(samples);
  fitpara.push_back(&fluxpara_norm);
  */
  //Det parameters
  /*DetParameters detpara("../inputs/ccqebins_det_Dec2014.txt",det_weights,"par_detAve");
  detpara.SetCovarianceMatrix(&cov_det);
  detpara.InitEventMap(samples);
  fitpara.push_back(&detpara);*/
  /*DetParameters detpara_fine("../inputs/ccqebins_Dec2014.txt",det_weights_fine,"par_detFine");
  detpara_fine.SetCovarianceMatrix(&cov_det_fine);
  detpara_fine.InitEventMap(samples);
  fitpara.push_back(&detpara_fine);
  */
  //Xsec parameters
  /*XsecParameters xsecpara;
  xsecpara.SetCovarianceMatrix(&cov_xsec);
  xsecpara.StoreResponseFunctions(responsefunctions, v_D1edges, v_D2edges);
  xsecpara.InitEventMap(samples);
  fitpara.push_back(&xsecpara);
  */
  //FSI parameters
  /*FSIParameters fsipara;
  fsipara.SetCovarianceMatrix(&cov_fsi);
  fsipara.StoreResponseFunctions(responsefunctions_FSI, v_D1edges, v_D2edges);
  fsipara.InitEventMap(samples);
  fitpara.push_back(&fsipara);
  */
  //Nucleon FSI parameters
  /*NuclFSIParameters nuclfsipara;
  nuclfsipara.SetCovarianceMatrix(&cov_nuclfsi);
  nuclfsipara.StoreResponseFunctions(responsefunctions_NuclFSI, v_D1edges, v_D2edges);
  nuclfsipara.InitEventMap(samples);
  fitpara.push_back(&nuclfsipara);
  */
  //Instantiate fitter obj
  XsecFitter xsecfit(seed);
  //init w/ para vector
  std::cout << "initialising fitter with regularisaion param " << regparam << std::endl;
  xsecfit.InitFitter(fitpara, regparam);
   
  //fix parameters
  //xsecfit.FixParameter("flux_norm",1);
  //xsecfit.FixParameter("flux_shape",1);
  //xsecfit.FixParameter("flux",1);
  //xsecfit.FixParameter("detFine",1);
  /*xsecfit.FixParameter("MAres",1.41);
  xsecfit.FixParameter("CCoth",0.0);
  xsecfit.FixParameter("Piless",0.0);
  xsecfit.FixParameter("CC1piE0",1.1);
  xsecfit.FixParameter("CC1piE1",1.0);
  xsecfit.FixParameter("CCCoh",1.0);
  xsecfit.FixParameter("NCoth",1.0);
  xsecfit.FixParameter("NC1pi0E0",0.96);
  xsecfit.FixParameter("NC1piE0",1.0);
  xsecfit.FixParameter("PionFSI",1);
  */
  //set frequency to save output
  xsecfit.SetSaveMode(fout, 1);

  //do fit: 1 = generate toy dataset from nuisances (WITH stat fluct)
  //        2 = fake data from MC or real data (+ stat fluct)
  //        3 = no nuisance sampling only stat fluctuation 
  xsecfit.Fit(samples, 2);
  
  fout->Close();

  return 0;
}
