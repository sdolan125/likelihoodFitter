#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>

#include "FSIParameters.hh"

using namespace std;

//ctor
FSIParameters::FSIParameters(const char *name)
{
  m_name = name;
 
  //parameter number, names & prior values
  //(Npar, values and ordering should be in agreement with input TFiles)
  //Npar=5;
  Npar=6;
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"InelLow"));
  pars_prior.push_back(1.0);
  pars_step.push_back(0.412);
  pars_limlow.push_back(-0.236);
  pars_limhigh.push_back(2.236);
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"InelHigh"));
  pars_prior.push_back(1.0);
  pars_step.push_back(0.338);
  pars_limlow.push_back(-0.014);
  pars_limhigh.push_back(2.014);
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"PiAbs"));
  pars_prior.push_back(1.0);
  pars_step.push_back(0.412);
  pars_limlow.push_back(-0.236);
  pars_limhigh.push_back(2.236);
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"PiProd"));
  pars_prior.push_back(1.0);
  pars_step.push_back(0.5);
  pars_limlow.push_back(-0.5);
  pars_limhigh.push_back(2.5);
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"CExLow"));
  pars_prior.push_back(1.0);
  pars_step.push_back(0.567);
  pars_limlow.push_back(-0.701);
  pars_limhigh.push_back(2.701);
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"CExHigh"));
  pars_prior.push_back(1.0);
  pars_step.push_back(0.278);
  pars_limlow.push_back(0);
  pars_limhigh.push_back(2);

 }

//dtor
FSIParameters::~FSIParameters()
{;}

// store response functions in vector of FSI "bins" (Ereco, Etrue, reac, topo)
void FSIParameters::StoreResponseFunctions(vector<TFile*> respfuncs, std::vector<std::pair <double,double> > v_pedges, 
					    std::vector<std::pair <double,double> > v_cthedges)
{
  for ( int stInt = mutrack; stInt != crDIS+1; stInt++ ){
    SampleTypes sampletype = static_cast <SampleTypes> (stInt);
    for ( int rtInt = ReCCQE; rtInt != OutFGD+1; rtInt++){
      ReactionTypes reactype = static_cast<ReactionTypes>(rtInt);
      cout<<"reading response functions for topology "<<stInt<<"  reaction "<<rtInt<<endl;
      int nccqebins=v_pedges.size();
      for(int br=0;br<nccqebins;br++){//reco kinematics bin
	//cout<<"reading rewighting function for reco bin "<<br<<endl;
	for(int bt=0;bt<nccqebins;bt++){//true kinematics bin
	  //cout<<"reading rewighting function for true bin "<<bt<<endl;
	  FSIBin bin;
	  bin.recoPlow = v_pedges[br].first;
	  bin.recoPhigh = v_pedges[br].second;
	  bin.truePlow = v_pedges[bt].first; //same binning for reco and true
	  bin.truePhigh = v_pedges[bt].second;
	  bin.recoCTHlow = v_cthedges[br].first;
	  bin.recoCTHhigh = v_cthedges[br].second;
	  bin.trueCTHlow = v_cthedges[bt].first; //same binning for reco and true
	  bin.trueCTHhigh = v_cthedges[bt].second;
	  bin.topology = sampletype; 
	  bin.reaction = reactype; 
	  if(fabs(br-bt)<21) {  //save memory if reco bin and true bin very far away
	    for(uint i=0; i<Npar; i++){
	      char name[200];
	      sprintf(name,"topology_%d/RecBin_%d_trueBin_%d_topology_%d_reac_%d",stInt,br,bt,stInt,rtInt);
	      TGraph* g=(TGraph*)respfuncs[i]->Get(name);
	      g->SetName(name);
	      bin.respfuncs.push_back(g);
	    }
	  }
	  m_bins.push_back(bin);
	}
      }
    }
  }

  /*for(size_t j=0; j<m_bins.size();j++){
    cout<<j<<" topology: "<<m_bins[j].topology<<"  reaction: "<<m_bins[j].reaction
	<<"  recoP: "<<m_bins[j].recoPlow<<"-"<<m_bins[j].recoPhigh
	<<"  trueP: "<<m_bins[j].truePlow<<"-"<<m_bins[j].truePhigh
	<<"  recoCTH: "<<m_bins[j].recoCTHlow<<"-"<<m_bins[j].recoCTHhigh
	<<"  trueCTH: "<<m_bins[j].trueCTHlow<<"-"<<m_bins[j].trueCTHhigh<<endl;
    if(m_bins[j].respfuncs.size()>0)
	cout<<" response function name "<<m_bins[j].respfuncs[0]->GetName()<<endl;
      else
	cout<<" no response function"<<endl;
	}*/
  
}

// --
int FSIParameters::GetBinIndex(SampleTypes sampletype, ReactionTypes reactype, 
				double Pmureco, double Pmutrue, double CTHmureco, double CTHmutrue)
{
  int binn = BADBIN;
  for(size_t i=0;i<m_bins.size();i++)
    {
      if(m_bins[i].topology == sampletype && m_bins[i].reaction == reactype &&
	 (Pmureco > m_bins[i].recoPlow) && (Pmureco  < m_bins[i].recoPhigh)  &&
	 (CTHmureco  > m_bins[i].recoCTHlow) && (CTHmureco  < m_bins[i].recoCTHhigh) &&
	 (Pmutrue > m_bins[i].truePlow) && (Pmutrue  < m_bins[i].truePhigh)  &&
	 (CTHmutrue  > m_bins[i].trueCTHlow) && (CTHmutrue  < m_bins[i].trueCTHhigh)){
	binn = (int)i;
	break;
	}
    }
  /*cout<<"topology "<<sampletype<<"  reaction "<<reactype<<endl;
  cout<<"recoP "<<Pmureco<<"  trueP "<<Pmutrue<<"    recoCTH "<<CTHmureco<<"  trueCTH "<<CTHmutrue<<endl;
  cout<<"BIN "<<binn<<endl<<endl;*/
  return binn;
}

// initEventMap
void FSIParameters::InitEventMap(std::vector<AnaSample*> &sample) 
{
  if(m_bins.empty()) 
    {
      cout<<"Need to build map of response functions for "<<m_name<<" ... exiting ..."<<endl;
      exit(-1);
    }
  m_evmap.clear();
  
  //loop over events to build index map
  for(size_t s=0;s<sample.size();s++)
    {
      vector<int> row;
      for(int i=0;i<sample[s]->GetN();i++)
	{
	  AnaEvent *ev = sample[s]->GetEvent(i);
	  //skip reactions not prepared in response function
	  /*int code = PASSEVENT;
	  if(ev->GetReaction() == AntiNu || 
	     ev->GetReaction() == OutFGD) 
	    {
	      row.push_back(code);
	      continue;
	      }*/	  
	  //get event info
	  int binn = GetBinIndex(static_cast<SampleTypes>(ev->GetSampleType()),
				 static_cast<ReactionTypes>(ev->GetReaction()),
				 ev->GetRecPtrk(),ev->GetTruePtrk(),ev->GetRecCThtrk(),ev->GetTrueCThtrk());
	  if(binn == BADBIN) 
	    {
	      cout<<"WARNING: "<<m_name<<" event "
		  <<" fall outside bin ranges"<<endl;
	      cout<<"        This event will be ignored in analysis."
	      <<endl;
	      ev->Print();
	    }
	  row.push_back(binn); 
	}//event loop
      m_evmap.push_back(row);
    }//sample loop
}


// EventWeghts
void FSIParameters::EventWeights(std::vector<AnaSample*> &sample, 
				  std::vector<double> &params)
{
  
  for(size_t s=0;s<sample.size();s++)
    {
      for(int i=0;i<sample[s]->GetN();i++)
	{
	  AnaEvent *ev = sample[s]->GetEvent(i);
	  ReWeight(ev, s, i, params);
	}
    }
}

// ReWeight
void FSIParameters::ReWeight(AnaEvent *event, int nsample, int nevent,
			      std::vector<double> &params)
{
  if(m_evmap.empty()) //need to build an event map first
    {
      cout<<"Need to build event map index for "<<m_name<<endl;
      return;
    }

  int binn = m_evmap[nsample][nevent];

  if(binn == PASSEVENT) return;
  if(binn == BADBIN)
    event->AddEvWght(0.0); //skip!!!!
  else
    {
      vector <TGraph*> respfuncs = m_bins[binn].respfuncs;
      double weight=1;
      if(respfuncs.size()>0){ //needed because there are missing reponse functions when reco very different from true (to save memory)
	for(uint i=0; i<Npar; i++){
	  weight = weight*(respfuncs[i]->Eval(params[i]));
	  //if(weight!=1)
	    //cout<<"reweighting using weight "<<weight<<"  from bin "<<binn<<endl;
	}
      }
       event->AddEvWght(weight);
    }
}


