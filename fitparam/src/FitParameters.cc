#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>

#include "FitParameters.hh"

using namespace std;

FitParameters::FitParameters(const char *fname, const char *name)
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

FitParameters::~FitParameters()
{;}

void FitParameters::SetBinning(const char *fname)
{
  ifstream fin(fname);
  assert(fin.is_open());
  string line;
  CCQEBin bin;
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
      bin.D1low    = D1_1;
      bin.D1high   = D1_2;
      bin.D2low  = D2_1;
      bin.D2high = D2_2;
      m_bins.push_back(bin);
    }
  fin.close();
  Npar = m_bins.size();
  
  cout<<endl<<"CCQE binning:"<<endl;
  for(size_t i = 0;i<m_bins.size();i++)
    {
      cout<<setw(3)<<i
	  <<setw(5)<<m_bins[i].D2low
	  <<setw(5)<<m_bins[i].D2high
	  <<setw(5)<<m_bins[i].D1low
	  <<setw(5)<<m_bins[i].D1high<<endl;
    }
  cout<<endl;
}

int FitParameters::GetBinIndex(double D1, double D2)
{
  int binn = BADBIN;
  for(size_t i=0;i<m_bins.size();i++)
    {
      if(D1>=m_bins[i].D1low && D1<m_bins[i].D1high &&
	 D2>=m_bins[i].D2low && D2<m_bins[i].D2high)
	binn = (int)i;
    }

  return binn;
}
      

// initEventMap
void FitParameters::InitEventMap(std::vector<AnaSample*> &sample) 
{
  m_evmap.clear();
  
  //loop over events to build index map
  for(size_t s=0;s<sample.size();s++)
    {
      vector<int> row;
      for(int i=0;i<sample[s]->GetN();i++)
	    {
	      AnaEvent *ev = sample[s]->GetEvent(i);
        //DEBUG TIME:
        //cout << "*** out of if ***" << endl;
        //if( (s==0 && i==0) && ((ev->GetReaction()==1) || (ev->GetReaction()==2) ) && ((ev->GetpMomTrue()<460)||(ev->GetmuMomTrue()<460)||(ev->GetpCosThetaTrue()<0.4)||(ev->GetmuCosThetaTrue()<0.4))) ev->Print();
	      //
	      int code = PASSEVENT; // -1 by default
        //if(ev->GetReaction() != ccqe_recode) //pass if not CCQE
        //if((ev->GetReaction()!=1)&&(ev->GetReaction()!=2)) //pass if not CC0Pi in mectopology cat
        if(( (ev->GetReaction()!=1)&&(ev->GetReaction()!=2) )  ||
           ( (ev->GetpMomTrue()<450)||(ev->GetmuMomTrue()<200)||(ev->GetpCosThetaTrue()<-1.0)||(ev->GetmuCosThetaTrue()<0.0) )) //pass if not CC0Pi+np (n>=1) in mectopology cat
	      {
	        row.push_back(code);

          //DEBUG TIME:
          //cout << "*** In if ***" << endl;
          //if( ((ev->GetReaction()==1) || (ev->GetReaction()==2) ) && ((ev->GetpMomTrue()<460)||(ev->GetmuMomTrue()<460)||(ev->GetpCosThetaTrue()<0.4)||(ev->GetmuCosThetaTrue()<0.4))) ev->Print();
          //
	        continue;
	      }
	      //get event true D1 and D2
	      double D1   = ev->GetTruePtrk();
	      double D2 = ev->GetTrueCThtrk();
	      int binn   = GetBinIndex(D1, D2);
	      if(binn == BADBIN)
	      {
	        cout<<"WARNING: "<<m_name<<" D1 = "<<D1<<" D2 = "
		          <<D2<<" fall outside bin ranges"<<endl;
	        cout<<"        This event will be ignored in analysis."
		          <<endl;
	      }
        //cout << "FitParameters:InitEventMap: binn is " << binn << endl;  
	      row.push_back(binn);
	    }
      m_evmap.push_back(row);
    }
}

// EventWeghts
void FitParameters::EventWeights(std::vector<AnaSample*> &sample, 
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


void FitParameters::ReWeight(AnaEvent *event, int nsample, int nevent,
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
  double D1   = event->GetTrueD1trk();
  double D2 = event->GetTrueD2trk();
  cout<<"Event: "<<event->GetEvId()<<endl;
  cout<<event->GetReaction()<<endl;
  cout<<D1<<" "<<D2<<" "<<" "<<binn<<endl;
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
