#define anyTreeMC_cxx

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "anyTreeMC.hh"

using namespace std;

void anyTreeMC::GetEvents(std::vector<AnaSample*> ana_samples)
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
      int evtTopo=evtTopology; //0 CC0pi, 1 CC1pi, 2 CCOther, 3 backg(NC+antinu), 7 OOFV
      //if(reac==7) reac=4;
      //cout << "Evt Topology is " << evtTopo << endl;
      if(evtTopo>7) cout << "*** Warning: evtTopology>7, evt Topology is " << evtTopo << endl;
      ev.SetReaction(evtTopology);
      ev.SetTrueEnu(TrueEnergy/1000.0);   //MeV --> GeV
      ev.SetRecEnu(MainRecEneGlb/1000.0); //MeV --> GeV
      ev.SetTruePtrk(trueD1/1000.0);     //MeV/c --> GeV/c
      ev.SetRecPtrk(MainD1Glb/1000.0);   //MeV/c --> GeV/c
      ev.SetTrueCThtrk(trueD2);
      ev.SetRecCThtrk(MainD2);
      ev.SetEvWght(weight);
      ev.SetEvWghtMC(weight);

      ev.SetpMomRec(pMomRec);
      ev.SetpMomTrue(pMomTrue);
      ev.SetmuMomRec(muMomRec);
      ev.SetmuMomTrue(muMomTrue);
      ev.SetmuCosThetaRec(muCosThetaRec);
      ev.SetmuCosThetaTrue(muCosThetaTrue);
      ev.SetpCosThetaRec(pCosThetaRec);
      ev.SetpCosThetaTrue(pCosThetaTrue);
      //ev.Print();

      //DEBUG TIME
      //if(ev.GetEvId()==28) ev.Print();




      
      //Loop over AnaSample objs and fill appropriate ones
      for(size_t i=0;i<ana_samples.size();i++)
	    {
	      if(ana_samples[i]->GetSampleType() == qesampleFinal)
	      {
          //DEBUG TIME
          //if( ((ev.GetReaction()==1) || (ev.GetReaction()==2) ) && ((ev.GetpMomTrue()<460)||(ev.GetmuMomTrue()<460)||(ev.GetpCosThetaTrue()<0.4)||(ev.GetmuCosThetaTrue()<0.4))) ev.Print();
          //
	        ana_samples[i]->AddEvent(ev);
	        break;
	      }  
	    }  
  } // jentry loop

  //Print some stats
  for(size_t i=0;i<ana_samples.size();i++)
    ana_samples[i]->PrintStats();
}
