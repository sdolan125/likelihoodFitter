#define ccqeTreeMC_cxx

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "ccqeTreeMC.hh"

using namespace std;

void ccqeTreeMC::GetEvents(std::vector<AnaSample*> ana_samples)
{
  if (fChain == 0) return;
  if (ana_samples.empty()) return;
  
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  //nentries = 500;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if(jentry % (int)1e+5 == 0)
	cout<<"Processing "<<jentry<<" out of "<<nentries<<endl;
      nb = fChain->GetEntry(jentry); nbytes += nb;
      //create and fill event structure
      AnaEvent ev(jentry);
      ev.SetSampleType(qesampleFinal);
      int reac=reaction; //0 CC0pi, 1 CC1pi, 2 CCOther, 3 backg(NC+antinu), 7 OOFV
      if(reac==7)
	reac=4;
      ev.SetReaction(reac);
      ev.SetTrueEnu(TrueEnergy/1000.0);   //MeV --> GeV
      ev.SetRecEnu(MainRecEneGlb/1000.0); //MeV --> GeV
      ev.SetTruePtrk(trueMom/1000.0);     //MeV/c --> GeV/c
      ev.SetRecPtrk(MainMomGlb/1000.0);   //MeV/c --> GeV/c
      ev.SetTrueCThtrk(trueCostheta);
      ev.SetRecCThtrk(MainCosTheta);
      ev.SetEvWght(weight);
      ev.SetEvWghtMC(weight);
      //ev.Print();
      
      //Loop over AnaSample objs and fill appropriate ones
      for(size_t i=0;i<ana_samples.size();i++)
	{
	  if(ana_samples[i]->GetSampleType() == qesampleFinal)
	    {
	      ana_samples[i]->AddEvent(ev);
	      break;
	    }
	}
    } // jentry loop

  //Print some stats
  for(size_t i=0;i<ana_samples.size();i++)
    ana_samples[i]->PrintStats();
}
