//////////////////////////////////////////////////////////
//
//  A class for event samples for for Any analysis
//
//
//
//  Created: Nov 17 2015   
//
//////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>

#include "AnySample.hh"

using namespace std;

// ctor
AnySample::AnySample(int sample_id, string name, 
		       std::vector<std::pair <double,double> > v_d1edges, 
		       std::vector<std::pair <double,double> > v_d2edges, 
		       TTree* data)
{  
  m_sampleid = sample_id; //unique id
  m_name     = name;      //some comprehensible name
  //cout<<"NEW SAMPLE with name "<<name<<" sample id "<<sample_id<<endl;
  m_data_tree = data;
  m_pedges = v_d1edges;
  m_cthedges = v_d2edges;

  for(int i=0;i<v_d1edges.size(); i++){
    cout<<v_d2edges[i].first<<"  "<<v_d2edges[i].second<<"  "<<v_d1edges[i].first<<"  "<<v_d1edges[i].second<<endl;
  }

  //Default binning choices
  nbins_enu = 28;
  bins_enu = new double[nbins_enu + 1];
  for(int i=0;i<=nbins_enu;i++)
    {
      if(i<10) bins_enu[i] = i*0.1;
      else if(i>=10 && i<18) bins_enu[i] = 1.0 + (i-10)*0.2;
      else if(i==18) bins_enu[i] = 2.7;
      else if(i>=19 && i<23) bins_enu[i] = 3.0 + (i-19)*0.5;
      else if(i>=23) bins_enu[i] = 5.0 + (i-23)*1.0;
    }
  /*
  nbins_pmu = 15;
  bins_pmu = new double[nbins_pmu + 1];
  for(int i=0;i<=nbins_pmu;i++)
    {
      if(i==0) bins_pmu[i] = 0.0;
      else if(i>=1 && i<=9) bins_pmu[i] = 0.2 + (i-1)*0.1;
      else if(i==10) bins_pmu[i] = 1.25;
      else if(i==11) bins_pmu[i] = 1.5;
      else if(i==12) bins_pmu[i] = 2.0;
      else if(i==13) bins_pmu[i] = 3.0;
      else if(i==14) bins_pmu[i] = 5.0;
      else if(i==15) bins_pmu[i] = 30.0;
    }
  
  nbins_cthmu = 11;
  bins_cthmu = new double[nbins_cthmu + 1];
  bins_cthmu[0]  = -1.0;
  bins_cthmu[1]  =  0.6;
  bins_cthmu[2]  =  0.7;
  bins_cthmu[3]  =  0.8;
  bins_cthmu[4]  =  0.85;
  bins_cthmu[5]  =  0.9;
  bins_cthmu[6]  =  0.92;
  bins_cthmu[7]  =  0.94;
  bins_cthmu[8]  =  0.96;
  bins_cthmu[9]  =  0.98;
  bins_cthmu[10] =  0.99;
  bins_cthmu[11] =  1.00;
  */

  // New way of choosing default binning based on input binning, 
  // relies on later bins always having larger values than earlier bins!

  int nbins_temp = v_d1edges.size();
  vector<double> bins_d1_vector;
  for(int i=0;i<nbins_temp;i++){
    if (bins_d1_vector.size()==0) bins_d1_vector.push_back(v_d1edges[i].first);
    else{
      if (bins_d1_vector.back()!=v_d1edges[i].first) bins_d1_vector.push_back(v_d1edges[i].first);
    }
  }
  bins_d1_vector.push_back(v_d1edges[nbins_temp-1].second);

  nbins_pmu = bins_d1_vector.size()-1;
  cout << "There are " << nbins_pmu << " d1 bins" << endl;

  bins_pmu = new double[nbins_pmu + 1];

  for(int i=0;i<=nbins_pmu;i++){
    bins_pmu[i]  =  bins_d1_vector[i];
    cout << "bins_pmu " << i << " is " << bins_d1_vector[i] << endl;
  }


  nbins_temp = v_d1edges.size();
  vector<double> bins_d2_vector;
  for(int i=0;i<nbins_temp;i++){
    if (bins_d2_vector.size()==0) bins_d2_vector.push_back(v_d2edges[i].first);
    else{
      if (bins_d2_vector.back()!=v_d2edges[i].first) bins_d2_vector.push_back(v_d2edges[i].first);
    }
  }
  bins_d2_vector.push_back(v_d2edges[nbins_temp-1].second);

  nbins_cthmu = bins_d2_vector.size()-1;
  cout << "There are " << nbins_cthmu << " d2 bins" << endl;

  bins_cthmu = new double[nbins_cthmu + 1];

  for(int i=0;i<=nbins_cthmu;i++){
    bins_cthmu[i]  =  bins_d2_vector[i];
    cout << "bins_cthmu " << i << " is " << bins_d2_vector[i] << endl;
  }


  nAnybins=m_pedges.size();
  bins_Any = new double[nAnybins+1];
  for (int i=0; i<=nAnybins; i++){
    bins_Any[i]=i;
  }
  cout<<"Any bins defined"<<endl;
  //event distribution histo
  m_hpred = NULL;
  m_hmc   = NULL;
  MakeHistos(); //with default binning
  
  cout<<"MakeHistos called"<<endl;
  //data (or toy) histo
  m_hdata = NULL;

  m_norm  = 1.0;
}

// dtor
AnySample::~AnySample()
{
  m_hpred->Delete();
  m_hmc->Delete();
  if(m_hdata != NULL) m_hdata->Delete();
  delete [] bins_pmu;
  delete [] bins_cthmu;
  delete [] bins_enu;
}

// MakeEventHisto
void AnySample::MakeHistos()
{
  if(m_hpred != NULL) m_hpred->Delete();
  m_hpred = new TH1D(Form("%s_pred_recpmucthmu", m_name.c_str()),
		     Form("%s_pred_recpmucthmu", m_name.c_str()),
		     nAnybins, bins_Any);
  /*new TH2D(Form("%s_pred_recpmucthmu", m_name.c_str()),
		       Form("%s_pred_recpmucthmu", m_name.c_str()),
		       nbins_pmu, bins_pmu, nbins_cthmu, bins_cthmu);
  m_hpred->GetXaxis()->SetTitle("p_{#mu}^{Rec} (GeV/c)");
  m_hpred->GetYaxis()->SetTitle("cos#theta_{#mu}^{Rec}");*/
  m_hpred->SetDirectory(0);

  if(m_hmc != NULL) m_hmc->Delete();
  m_hmc = new TH1D(Form("%s_mc_recpmucthmu", m_name.c_str()),
		   Form("%s_mc_recpmucthmu", m_name.c_str()),
		   nAnybins, bins_Any);
  /*new TH2D(Form("%s_mc_recpmucthmu", m_name.c_str()),
		       Form("%s_mc_recpmucthmu", m_name.c_str()),
		       nbins_pmu, bins_pmu, nbins_cthmu, bins_cthmu);
  m_hmc->GetXaxis()->SetTitle("p_{#mu}^{Rec} (GeV/c)");
  m_hmc->GetYaxis()->SetTitle("cos#theta_{#mu}^{Rec}");*/
  m_hmc->SetDirectory(0);
  cout<<nAnybins<<" bins inside MakeHistos"<<endl;
}

void AnySample::SetData(TObject *hdata)
{
  //clone the data histogram internally
  if(m_hdata != NULL) m_hdata->Delete();
  //m_hdata = (TH2D*)hdata->Clone(Form("%s_data", m_name.c_str()));
  m_hdata = (TH1D*)hdata->Clone(Form("%s_data", m_name.c_str()));
  m_hdata->SetDirectory(0);
}

// SetPmuBinning
void AnySample::SetPmuBinning(int nbins, double *bins)
{
  nbins_pmu = nbins;
  delete [] bins_pmu;
  bins_pmu = new double[nbins_pmu + 1];
  for(int i=0;i<=nbins_pmu;i++)
    bins_pmu[i] = bins[i];
}

// SetCThmuBinning
void AnySample::SetCThmuBinning(int nbins, double *bins)
{
  nbins_cthmu = nbins;
  delete [] bins_cthmu;
  bins_cthmu = new double[nbins_cthmu + 1];
  for(int i=0;i<=nbins_cthmu;i++)
    bins_cthmu[i] = bins[i];
}

// SetEnuBinning
void AnySample::SetEnuBinning(int nbins, double *bins)
{
  nbins_enu = nbins;
  delete [] bins_enu;
  bins_enu = new double[nbins_enu + 1];
  for(int i=0;i<=nbins_enu;i++)
    bins_enu[i] = bins[i];
}

// FillEventHist
void AnySample::FillEventHisto(int datatype)
{
  if(m_hpred) m_hpred->Reset();
  if(m_hmc) m_hmc->Reset();
  for(size_t i=0;i<m_events.size();i++)
    {
      double pmu_rec   = m_events[i].GetRecPtrk();
      double cthmu_rec = m_events[i].GetRecCThtrk();
      double wght      = m_events[i].GetEvWght();
      for(int j=0; j<nAnybins; j++){
	if((pmu_rec > m_pedges[j].first) && (pmu_rec  < m_pedges[j].second)  &&
	  (cthmu_rec  > m_cthedges[j].first) && (cthmu_rec  < m_cthedges[j].second)){
	  m_hpred->Fill(j+0.5,wght);
	  m_hmc->Fill(j+0.5,wght);
	  break;
	}
      }
    }
  m_hpred->Scale(m_norm);
  m_hmc->Scale(m_norm);
  
  //data without stat variation: useful when nuisance parameters 
  //varied in the toys
  if(datatype==1) 
    {
      SetData(m_hpred);
      m_hdata->Reset();
      for(int j=1;j<=m_hpred->GetNbinsX();j++)
	{
	  double val = m_hpred->GetBinContent(j);
	  //cout<<"bin "<<j<<" entry "<<val<<endl;
	  if(val == 0.0) {
	    cout<<"AnySample:"<<m_sampleid<<" bin "<<j<<" with 0 entries may cause proble on chi2 computations"<<endl;
	    continue;
	  }
	  int binc = gRandom->Poisson(val);
	  //m_hdata->SetBinContent(j,binc); //with statistical fluctuations
	  m_hdata->SetBinContent(j,val);  //without statistical fluctuations
	}
    }

  //data from external (fake) dataset 
  else if(datatype==2) {
    SetData(m_hpred);
    m_hdata->Reset();
    double potD = 57.34;   //in units of 10^19
    //double potMC_genie=384.762;
    //double potMC_genie=389.5; //neut
    //double potMC_genie=380.0; //nuwro
    double potMC_genie = 57.34; 
    /*
    Float_t pmu_rec_tree,cthmu_rec_tree,wght; 
    Float_t recDphiT,recDpT,recDalphaT; 
    Int_t topology;
    int** accum_level;
    m_data_tree->SetBranchAddress("weight",&wght); 
    m_data_tree->SetBranchAddress("selmu_mom",&pmu_rec_tree);
    m_data_tree->SetBranchAddress("selmu_costheta",&cthmu_rec_tree);
    m_data_tree->SetBranchAddress("recDphiT",&recDphiT);
    m_data_tree->SetBranchAddress("recDpT",&recDpT);
    m_data_tree->SetBranchAddress("recDalphaT",&recDalphaT);
    m_data_tree->SetBranchAddress("topology",&topology);  
    m_data_tree->SetBranchAddress("accum_level",&accum_level);  
    */

    Float_t pmu_rec_tree,cthmu_rec_tree,wght; 
    Int_t topology;
    m_data_tree->SetBranchAddress("weight",&wght); 
    m_data_tree->SetBranchAddress("D1Rec",&pmu_rec_tree);
    m_data_tree->SetBranchAddress("D2Rec",&cthmu_rec_tree);
    m_data_tree->SetBranchAddress("cutBranch",&topology);        
   
    for(size_t i=0;i<m_data_tree->GetEntries();i++){
      m_data_tree->GetEntry(i);
      if(topology != m_sampleid) continue;
      for(int j=0; j<nAnybins; j++){
	if((pmu_rec_tree/1000 > m_pedges[j].first) && (pmu_rec_tree/1000  < m_pedges[j].second)  &&
	   (cthmu_rec_tree  > m_cthedges[j].first) && (cthmu_rec_tree  < m_cthedges[j].second)){
	  m_hdata->Fill(j+0.5,wght);
	  //cout<<"filling "<<pmu_rec_tree/1000<<" "<<wght_syst<<endl;
	  break;
	}
      }  
    }
    //m_hdata->Scale(potD/potMC_genie);
    //add MC or data (!!!!) statistical variations also to genie dataset to evaluate genie MC stat uncert
    /*for(int j=1;j<=m_hdata->GetNbinsX();j++)
	{
	  double val = m_hdata->GetBinContent(j);
	  //cout<<"bin "<<j<<" entry "<<val<<endl;
	  if(val == 0.0) {
	    cout<<"AnySample:"<<m_sampleid<<" bin "<<j<<" with 0 entries may cause proble on chi2 computations"<<endl;
	    continue;
	  }
	  int binc = gRandom->Poisson(val);
	  m_hdata->SetBinContent(j,binc);  //with statistical fluctuations
	}
    //m_hdata->Scale(potD/potMC_genie);
    */
  }
  //data with statistical variation 
  //(used when no nuisance sampling but nuisances are fitted)
   else if(datatype==3) 
    {
      SetData(m_hpred);
      m_hdata->Reset();
      for(int j=1;j<=m_hpred->GetNbinsX();j++)
	{
	  double val = m_hpred->GetBinContent(j);
	  //cout<<"bin "<<j<<" entry "<<val<<endl;
	  if(val == 0.0) {
	    cout<<"AnySample:"<<m_sampleid<<" bin "<<j<<" with 0 entries may cause proble on chi2 computations"<<endl;
	    continue;
	  }
	  int binc = gRandom->Poisson(val);
	  m_hdata->SetBinContent(j,binc); //with statistical fluctuations
	  //m_hdata->SetBinContent(j,val);  //without statistical fluctuations
	}
    }

}

double AnySample::CalcChi2()
{
  if(m_hdata == NULL) 
    {
      cerr<<"ERROR: need to define data histogram"<<endl;
      return 0.0;
    }
  
  int nx = m_hpred->GetNbinsX();
  //int ny = m_hpred->GetNbinsY();
  
  if(nx != m_hdata->GetNbinsX())// || ny != m_hdata->GetNbinsY())
    {
      cerr<<"ERROR: binning mismatch between data and mc"<<endl;
      return 0.0;
    }

  double chi2 = 0.0;
  for(int j=1;j<=nx;j++)
    {
      double obs = m_hdata->GetBinContent(j);
      double exp = m_hpred->GetBinContent(j);
      if(exp>0.0){  //added when external fake datasets (you cannot reweight when simply 0)
                    // this didn't happen when all from same MC since if exp=0 then obs =0
	chi2 += 2*(exp - obs);
	if(obs>0.0)
	  chi2 += 2*obs*TMath::Log(obs/exp);
	//else
	//cout<<"bin "<<j<<" has observation "<<obs<<" PROBLEM???"<<endl;
      }
    }
  
  if(chi2 != chi2)
    {
      cerr<<"ERROR: stat chi2 is nan"<<endl;
      chi2 = 0.0;
    }
 
  return chi2;
}

// GetSampleBreakdown 
void AnySample::GetSampleBreakdown(TDirectory *dirout, string tag)
{
  int nreac = 7;
  const char *names[] = {"cc0pi0p", "cc0pi1p", "cc0pinp", "cc1pi+", "ccother", 
			 "backg", "OOFV"};
  TH1D *henu_rec[nreac];
  TH1D *hpmu_rec[nreac];
  TH1D *hcth_rec[nreac];
  TH2D *hpmucth_rec[nreac];
  TH1D *hpmu_true[nreac];
  TH1D *hcth_true[nreac];
  TH1D *hAnybin_true[nreac];
  TH1D *hAnybin_rec[nreac];
  int compos[nreac];
  
  cout<<"AnySample::GetSampleBreakdown - Inializing histos of reactions" << endl;

  for(int i=0;i<nreac;i++)
    {
      compos[i] = 0;
      henu_rec[i] = new TH1D(Form("%s_RecEnu_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     Form("%s_RecEnu_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     nbins_enu, bins_enu);
      henu_rec[i]->SetDirectory(0);
      henu_rec[i]->GetXaxis()->SetTitle("Recon E_{#nu} (GeV)");

      hpmu_rec[i] = new TH1D(Form("%s_RecPmu_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     Form("%s_RecPmu_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     nbins_pmu, bins_pmu);
      hpmu_rec[i]->SetDirectory(0);
      hpmu_rec[i]->GetXaxis()->SetTitle("Recon p_{#mu} (GeV/c)");

      hcth_rec[i] = new TH1D(Form("%s_RecCTh_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     Form("%s_RecCTh_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     nbins_cthmu, bins_cthmu);
      hcth_rec[i]->SetDirectory(0);
      hcth_rec[i]->GetXaxis()->SetTitle("Recon cos#theta_{#mu}");
      
      hpmucth_rec[i] = new TH2D(Form("%s_RecPmuCTh_%s_%s",
				     m_name.c_str(),names[i],tag.c_str()),
				Form("%s_RecPmuCTh_%s_%s",
				     m_name.c_str(),names[i],tag.c_str()),
				nbins_pmu, bins_pmu, 
				nbins_cthmu, bins_cthmu);
      hpmucth_rec[i]->SetDirectory(0);
      hpmucth_rec[i]->GetXaxis()->SetTitle("Recon p_{#mu} (GeV/c)");
      hpmucth_rec[i]->GetXaxis()->SetTitle("Recon cos#theta_{#mu}");

      hpmu_true[i] = new TH1D(Form("%s_TruePmu_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     Form("%s_TruePmu_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     nbins_pmu, bins_pmu);
      hpmu_true[i]->SetDirectory(0);
      hpmu_true[i]->GetXaxis()->SetTitle("True p_{#mu} (GeV/c)");

      hcth_true[i] = new TH1D(Form("%s_TrueCTh_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     Form("%s_TrueCTh_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     nbins_cthmu, bins_cthmu);
      hcth_true[i]->SetDirectory(0);
      hcth_true[i]->GetXaxis()->SetTitle("True cos#theta_{#mu}");

      hAnybin_true[i] = new TH1D(Form("%s_Anybins_true_%s_%s",
				       m_name.c_str(),names[i],tag.c_str()),
				  Form("%s_Anybins_true_%s_%s",
				       m_name.c_str(),names[i],tag.c_str()),
			     nAnybins, bins_Any);
      hAnybin_true[i]->SetDirectory(0);
      hAnybin_true[i]->GetXaxis()->SetTitle("Any bins");

      hAnybin_rec[i] = new TH1D(Form("%s_Anybins_rec_%s_%s",
				      m_name.c_str(),names[i],tag.c_str()),
			     Form("%s_Anybins_rec_%s_%s",
				  m_name.c_str(),names[i],tag.c_str()),
			     nAnybins, bins_Any);
      hAnybin_rec[i]->SetDirectory(0);
      hAnybin_rec[i]->GetXaxis()->SetTitle("Any bins");
    }

  //loop over the events
  // ////////////////////////////////////////
  // double Enutrue,Pmutrue,CTHmutrue;
  // double Enureco,Pmureco,CTHmureco;
  // double weight;
  // int reaction;
  // TTree *tree = new TTree(Form("finalMC_%s",m_name.c_str()),Form("finalMC_%s",m_name.c_str()));
  // tree->Branch("reaction",&reaction,"reaction/I");
  // tree->Branch("Pmureco",&Pmureco,"Pmureco/D");
  // tree->Branch("CTHmureco",&CTHmureco,"CTHmureco/D");
  // tree->Branch("Pmutrue",&Pmutrue,"Pmutrue/D");
  // tree->Branch("CTHmutrue",&CTHmutrue,"CTHmutrue/D");
  // tree->Branch("weight",&weight,"weight/D");
  // ////////////////////////////////////////


  cout<<"AnySample::GetSampleBreakdown - Collecting events" << endl;
  int Ntot = GetN();
  for(size_t i=0;i<m_events.size();i++)
    {
      //cout<<"AnySample::GetSampleBreakdown - In event loop iteration " << i << " out of " << m_events.size() << endl;
      double enu_rec, pmu_rec, cth_rec, pmu_true, cth_true, wght;
      enu_rec = m_events[i].GetRecEnu();
      pmu_rec = m_events[i].GetRecPtrk();
      cth_rec = m_events[i].GetRecCThtrk();
      pmu_true = m_events[i].GetTruePtrk();
      cth_true = m_events[i].GetTrueCThtrk();
      wght    = m_events[i].GetEvWght();
      int rtype = m_events[i].GetReaction();
      //cout<< "AnySample::GetSampleBreakdown - rtype is: " << rtype << endl;
      if((rtype==7)) rtype=6; //BKG is 5 then OOFV is 7, 6 is skipped causing array to overrun

      /*if((rtype==0)||(rtype==1)||(rtype==2)) rtype=0; //cc0pi
      else if(rtype==3) rtype=1; //cc1pi
      else if(rtype==4) rtype=2; //ccother
      else if(rtype==5) rtype=3; //BKG
      else if(rtype>5) rtype=4; //OOFV/noTruth*/


      //cout<< "AnySample::GetSampleBreakdown - Event breakdown:" << endl;
      //m_events[i].Print();

      // /////////////////
      // reaction=rtype;
      // Pmutrue=pmu_true ;
      // CTHmutrue= cth_true;
      // weight=wght;
      // tree->Fill();
      // ////////////////////////////////////////
      
      //cout<< "AnySample::GetSampleBreakdown - Filling histos" << endl;

      //if(rtype<0 && rtype>6) continue;
      compos[rtype]++;
      henu_rec[rtype]->Fill(enu_rec, wght);
      hpmu_rec[rtype]->Fill(pmu_rec, wght);
      hcth_rec[rtype]->Fill(cth_rec, wght);
      hpmu_true[rtype]->Fill(pmu_true, wght);
      hcth_true[rtype]->Fill(cth_true, wght);
      hpmucth_rec[rtype]->Fill(pmu_rec, cth_rec, wght);


      //cout<< "AnySample::GetSampleBreakdown - Filling histos with analysis binning" << endl;

      for(int j=0; j<nAnybins; j++){
	      if((pmu_true > m_pedges[j].first) && (pmu_true < m_pedges[j].second)  &&
	         (cth_true > m_cthedges[j].first) && (cth_true < m_cthedges[j].second)){
	        hAnybin_true[rtype]->Fill(j+0.5,wght);
	        break;
	      }
      }
      for(int j=0; j<nAnybins; j++){
      	if((pmu_rec > m_pedges[j].first) && (pmu_rec < m_pedges[j].second)  &&
	         (cth_rec > m_cthedges[j].first) && (cth_rec < m_cthedges[j].second)){
	           hAnybin_rec[rtype]->Fill(j+0.5,wght);
	           break;
	      }
      } 
    }
  
  cout<<"AnySample::GetSampleBreakdown - Wrapping up" << endl;

  dirout->cd();
  //tree->Write();
  for(int i=0;i<nreac;i++)
    {
      henu_rec[i]->Scale(m_norm);
      hpmu_rec[i]->Scale(m_norm);
      hcth_rec[i]->Scale(m_norm);
      hpmu_true[i]->Scale(m_norm);
      hcth_true[i]->Scale(m_norm);
      hpmucth_rec[i]->Scale(m_norm);
      hAnybin_true[i]->Scale(m_norm);
      hAnybin_rec[i]->Scale(m_norm);
 
      henu_rec[i]->Write();
      hpmu_rec[i]->Write();
      hcth_rec[i]->Write();
      hpmu_true[i]->Write();
      hcth_true[i]->Write();
      hpmucth_rec[i]->Write();
      hAnybin_true[i]->Write();
      hAnybin_rec[i]->Write();
 
      henu_rec[i]->Delete();
      hpmu_true[i]->Delete();
      hcth_true[i]->Delete();
      hpmu_rec[i]->Delete();
      hcth_rec[i]->Delete();
      hpmucth_rec[i]->Delete();
      hAnybin_true[i]->Delete();
      hAnybin_rec[i]->Delete();
    }
  cout<<"============> Sample "<<m_name<<" <============"<<endl;
  for(int j=0;j<nreac;j++)
    cout<<setw(10)<<names[j]<<setw(5)<<j<<setw(10)<<compos[j]
	<<setw(10)<<(float)(compos[j])/Ntot*100.0<<"%"<<endl;
}

// Write
void AnySample::Write(TDirectory *dirout, const char *bsname, int fititer)
{
  dirout->cd();
  m_hpred->Write(Form("%s_pred", bsname));
  if(fititer==0){
    m_hmc->Write(Form("%s_mc", bsname));
    if(m_hdata != NULL) m_hdata->Write(Form("%s_data", bsname));
  }
}
