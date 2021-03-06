#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>

#include "CCQEParameters.hh"

using namespace std;

CCQEParameters::CCQEParameters(const char *fname, const char *name)
{
  m_name    = name;
  hasCovMat = false;

  ccqe_recode = 0; //reaction code for CCQE
  
  //get the binning from a file
  SetBinning(fname);

  //parameter names & prior values
  for(size_t i=0;i<Npar;i++)
    {
      pars_name.push_back(Form("%s%d", m_name.c_str(), (int)i));
      pars_prior.push_back(1.0); //all weights are 1.0 a priori 
      pars_step.push_back(0.1);
      pars_limlow.push_back(0.0);
      pars_limhigh.push_back(10.0);
    }
}

CCQEParameters::~CCQEParameters()
{;}

void CCQEParameters::SetBinning(const char *fname)
{
  ifstream fin(fname);
  assert(fin.is_open());
  string line;
  CCQEBin bin;
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
      bin.plow    = p1;
      bin.phigh   = p2;
      bin.cthlow  = cth1;
      bin.cthhigh = cth2;
      m_bins.push_back(bin);
    }
  fin.close();
  Npar = m_bins.size();
  
  cout<<endl<<"CCQE binning:"<<endl;
  for(size_t i = 0;i<m_bins.size();i++)
    {
      cout<<setw(3)<<i
	  <<setw(5)<<m_bins[i].cthlow
	  <<setw(5)<<m_bins[i].cthhigh
	  <<setw(5)<<m_bins[i].plow
	  <<setw(5)<<m_bins[i].phigh<<endl;
    }
  cout<<endl;
}

int CCQEParameters::GetBinIndex(double p, double cth)
{
  int binn = BADBIN;
  for(size_t i=0;i<m_bins.size();i++)
    {
      if(p>=m_bins[i].plow && p<m_bins[i].phigh &&
	 cth>=m_bins[i].cthlow && cth<m_bins[i].cthhigh)
	binn = (int)i;
    }

  return binn;
}
      

// initEventMap
void CCQEParameters::InitEventMap(std::vector<AnaSample*> &sample) 
{
  cout << "CCQEParameters:InitEventMap: SHOULD NOT BE HERE WHY AM I HERE????" << endl;
  getchar();
  m_evmap.clear();
  
  //loop over events to build index map
  for(size_t s=0;s<sample.size();s++)
    {
      vector<int> row;
      for(int i=0;i<sample[s]->GetN();i++)
	{
	  AnaEvent *ev = sample[s]->GetEvent(i);
	  //
	  int code = PASSEVENT; // -1 by default
   //if(ev->GetReaction() != ccqe_recode) //pass if not CCQE
   //if(( (ev->GetReaction()!=1)&&(ev->GetReaction()!=2) ) ||
   //   ( (ev->GetpMomTrue()<460)||(ev->GetmuMomTrue()<460)||(ev->GetmuCosThetaTrue()<0.4)||(ev->GetmuCosThetaTrue()<0.4) )) //pass if not CC0Pi+np (n>=1) in mectopology cat
	  if(true)   
      {
	      row.push_back(code);
        cout << "CCQEParameters:InitEventMap: event not CCQE, therefore do not reweight" << endl;
	      continue;
	    }
	  //get event true p and cth
	  double p   = ev->GetTruePtrk();
	  double cth = ev->GetTrueCThtrk();
	  int binn   = GetBinIndex(p, cth);
	  if(binn == BADBIN)
	    {
	      cout<<"WARNING: "<<m_name<<" p = "<<p<<" cth = "
		  <<cth<<" fall outside bin ranges"<<endl;
	      cout<<"        This event will be ignored in analysis."
		  <<endl;
	    }
    //cout << "CCQEParameters:InitEventMap: binn is " << binn << endl;  
	  row.push_back(binn);
	}
      m_evmap.push_back(row);
    }
}

// EventWeghts
void CCQEParameters::EventWeights(std::vector<AnaSample*> &sample, 
				  std::vector<double> &params)
{
  if(m_evmap.empty()) //build an event map
    {
      cout<<"Need to build event map index for "<<m_name<<endl;
      InitEventMap(sample);
    }
  
  for(size_t s=0;s<sample.size();s++)
    {
      for(int i=0;i<sample[s]->GetN();i++)
	{
	  AnaEvent *ev = sample[s]->GetEvent(i);
	  ReWeight(ev, s, i, params);
	}
    }
}


void CCQEParameters::ReWeight(AnaEvent *event, int nsample, int nevent,
			      std::vector<double> &params)
{
  if(m_evmap.empty()) //need to build an event map first
    {
      cout<<"Need to build event map index for "<<m_name<<endl;
      return;
    }

  int binn = m_evmap[nsample][nevent];
  
  //skip event if not CCQE
  if(binn == PASSEVENT) return;

  /*
  double p   = event->GetTruePtrk();
  double cth = event->GetTrueCThtrk();
  cout<<"Event: "<<event->GetEvId()<<endl;
  cout<<event->GetReaction()<<endl;
  cout<<p<<" "<<cth<<" "<<" "<<binn<<endl;
  */
  
  if(binn == BADBIN)
    event->AddEvWght(0.0); //skip!!!!
  else
    {
      if(binn>(int)params.size())
	{
	  cerr<<"ERROR: number of bins "<<m_name
	      <<" does not match num of param"<<endl;
	  event->AddEvWght(0.0);
	}
      event->AddEvWght(params[binn]);
      //cout << "ReWeight param " << binn << endl;
      //cout << "Weight is " << params[binn] << endl;
    }
}
