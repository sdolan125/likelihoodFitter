#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>

#include "NuclFSIParameters.hh"

using namespace std;

//ctor
NuclFSIParameters::NuclFSIParameters(const char *name)
{
  m_name = name;
 
  Npar=3;
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"MFP"));
  pars_prior.push_back(1);
  pars_step.push_back(0.2);
  pars_limlow.push_back(0);
  pars_limhigh.push_back(2);
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"FrElas"));
  pars_prior.push_back(1);
  pars_step.push_back(0.3);
  pars_limlow.push_back(0);
  pars_limhigh.push_back(2);
  pars_name.push_back(Form("%s%d%s", m_name.c_str(), 1,"FrAbs"));
  pars_prior.push_back(1);
  pars_step.push_back(0.2);
  pars_limlow.push_back(0);
  pars_limhigh.push_back(2);
}

//dtor
NuclFSIParameters::~NuclFSIParameters()
{;}

// store response functions in vector of FSI "bins" (Ereco, Etrue, reac, topo)
void NuclFSIParameters::StoreResponseFunctions(vector<TFile*> respfuncs)
{
  for ( int stInt = mutrack; stInt != crDIS+1; stInt++ ){
    SampleTypes sampletype = static_cast <SampleTypes> (stInt);
    for ( int rtInt = ReCCQE; rtInt != AntiNu; rtInt++){
      ReactionTypes reactype = static_cast<ReactionTypes>(rtInt);
      cout<<"reading response functions for topology "<<stInt<<"  reaction "<<rtInt<<endl;
      bin.topology = sampletype; 
      bin.reaction = reactype; 
      for(uint i=0; i<Npar; i++){
	char name[200];
	sprintf(name,"topology_%d/topology_%d_reac_%d",stInt,stInt,rtInt);
	TGraph* g=(TGraph*)respfuncs[i]->Get(name);
	g->SetName(name);
	bin.respfuncs.push_back(g);
      }
      m_bins.push_back(bin);
    }
  }

  /*for(size_t j=0; j<m_bins.size();j++){
    cout<<j<<" topology: "<<m_bins[j].topology<<"  reaction: "<<m_bins[j].reaction
    if(m_bins[j].respfuncs.size()>0)
	cout<<" response function name "<<m_bins[j].respfuncs[0]->GetName()<<endl;
      else
	cout<<" no response function"<<endl;
	}*/
}

// --
int NuclFSIParameters::GetBinIndex(SampleTypes sampletype, ReactionTypes reactype, 
				double Pmureco, double Pmutrue, double CTHmureco, double CTHmutrue)
{
  int binn = BADBIN;
  for(size_t i=0;i<m_bins.size();i++)
    {
      if(m_bins[i].topology == sampletype && m_bins[i].reaction == reactype ){
	binn = (int)i;
	break;
      }
    }
  /*cout<<"topology "<<sampletype<<"  reaction "<<reactype<<endl;
  cout<<"BIN "<<binn<<endl<<endl;*/
  return binn;
}

// initEventMap
void NuclFSIParameters::InitEventMap(std::vector<AnaSample*> &sample) 
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
	  int code = PASSEVENT;
	  if(ev->GetReaction() == AntiNu || 
	     ev->GetReaction() == OutFGD) 
	    {
	      row.push_back(code);
	      continue;
	    }	  
	  //get event info
	  int binn = GetBinIndex(static_cast<SampleTypes>(ev->GetSampleType()),
				 static_cast<ReactionTypes>(ev->GetReaction()));
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
void NuclFSIParameters::EventWeights(std::vector<AnaSample*> &sample, 
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
void NuclFSIParameters::ReWeight(AnaEvent *event, int nsample, int nevent,
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


