//////////////////////////////////////////////////////////
//
//  A container with subset of necessary event 
//  information for CCQE analysis
//
//
//  Created: Thu Jun  6 11:28:13 CEST 2013   
//  Modified:
//
//////////////////////////////////////////////////////////
#ifndef __AnaEvent_hh__
#define __AnaEvent_hh__

#include <iostream>

#include <TMath.h>

class AnaEvent
{
 public:
  AnaEvent(Long64_t evid) //unique event id
  {
    m_evid       = evid;
    m_evtype     = 0; 
    m_reaction   = -1;
    m_sample     = -1;
    m_trueEnu    = -999.0;
    m_recEnu     = -999.0;
    m_truePtrk   = -999.0;
    m_trueThtrk  = -999.0;
    m_trueCThtrk = -999.0;
    m_recPtrk    = -999.0;
    m_recThtrk   = -999.0;
    m_recCThtrk  = -999.0;
    m_wght       = 1.0;
    m_wghtMC       = 1.0;

    // New kinematic variables always included for phase space cuts
    m_pMomRec = -999.0;
    m_pMomTrue = -999.0;
    m_muMomRec = -999.0;
    m_muMomTrue = -999.0;
    m_muCosThetaRec = -999.0;
    m_muCosThetaTrue = -999.0;
    m_pCosThetaRec = -999.0;
    m_pCosThetaTrue = -999.0;

  }
  ~AnaEvent(){;}

  //Set/Get methods
  void SetEvType(int val){ m_evtype = val; }
  int GetEvType(){ return m_evtype; }
  
  void SetReaction(int val){ m_reaction = val; }
  int GetReaction(){ return m_reaction; }
  
 void SetSampleType(int val){ m_sample = val; }
  int GetSampleType(){ return m_sample; }
  
  Long64_t GetEvId(){ return m_evid; }

  void SetTrueEnu(double val) {m_trueEnu = val;}
  double GetTrueEnu(){ return m_trueEnu; }

  void SetRecEnu(double val){ m_recEnu = val; }
  double GetRecEnu(){ return m_recEnu; }

  void SetTruePtrk(double val){ m_truePtrk = val; }
  double GetTruePtrk(){ return m_truePtrk; }

  void SetRecPtrk(double val){ m_recPtrk = val; }
  double GetRecPtrk(){ return m_recPtrk; }

  void SetTrueThtrk(double val)
  { 
    m_trueThtrk  = val; 
    m_trueCThtrk = TMath::Cos(m_trueThtrk);
  }
  void SetTrueCThtrk(double val)
  { 
    m_trueCThtrk = val;
    m_trueThtrk  = TMath::ACos(m_trueCThtrk);
  }
  double GetTrueThtrk(){ return m_trueThtrk; }
  double GetTrueCThtrk(){ return m_trueCThtrk; }
  
  void SetRecThtrk(double val)
  { 
    m_recThtrk  = val; 
    m_recCThtrk = TMath::Cos(m_recThtrk);
  }
  void SetRecCThtrk(double val)
  { 
    m_recCThtrk = val;
    m_recThtrk  = TMath::ACos(m_recCThtrk);
  }
  double GetRecThtrk(){ return m_recThtrk; }
  double GetRecCThtrk(){ return m_recCThtrk; }
  
  void SetEvWght(double val){ m_wght  = val; }
  void SetEvWghtMC(double val){ m_wghtMC  = val; }
  void AddEvWght(double val){ m_wght *= val; }
  double GetEvWght(){ return m_wght; }
  double GetEvWghtMC(){ return m_wghtMC; }

  // New kinematic variables always included for phase space cuts

  void SetpMomRec(double val){ m_pMomRec = val; }
  void SetpMomTrue(double val){ m_pMomTrue = val; }
  void SetmuMomRec(double val){ m_muMomRec = val; }
  void SetmuMomTrue(double val){ m_muMomTrue = val; }
  void SetmuCosThetaRec(double val){ m_muCosThetaRec = val; }
  void SetmuCosThetaTrue(double val){ m_muCosThetaTrue = val; }
  void SetpCosThetaRec(double val){ m_pCosThetaRec = val; }
  void SetpCosThetaTrue(double val){ m_pCosThetaTrue = val; }

  double GetpMomRec(){ return m_pMomRec; }
  double GetpMomTrue(){ return m_pMomTrue; }
  double GetmuMomRec(){ return m_muMomRec; }
  double GetmuMomTrue(){ return m_muMomTrue; }
  double GetmuCosThetaRec(){ return m_muCosThetaRec; }
  double GetmuCosThetaTrue(){ return m_muCosThetaTrue; }
  double GetpCosThetaRec(){ return m_pCosThetaRec; }
  double GetpCosThetaTrue(){ return m_pCosThetaTrue; }


  void Print()
  {
    std::cout<<"Event ID "<<m_evid<<std::endl
	     <<"  Reaction        "<<GetReaction()<<std::endl
	     <<"  Sample          "<<GetSampleType()<<std::endl
	     <<"  True energy     "<<GetTrueEnu()<<std::endl
	     <<"  Recon energy    "<<GetRecEnu()<<std::endl
	     <<"  True track mom  "<<GetTruePtrk()<<std::endl
	     <<"  Recon track mom "<<GetRecPtrk()<<std::endl
       <<"  True track cth  "<<GetTrueCThtrk()<<std::endl
	     <<"  Recon track cth "<<GetRecCThtrk()<<std::endl
	     <<"  True track th   "<<GetTrueThtrk()<<std::endl
	     <<"  Recon track th  "<<GetRecThtrk()<<std::endl
	     <<"  Event weight    "<<GetEvWght()<<std::endl
	     <<"  Event weight MC  "<<GetEvWghtMC()<<std::endl
       <<"  Recon proton momentum  " <<GetpMomRec() <<std::endl
       <<"  True proton momentum   " <<GetpMomTrue() <<std::endl
       <<"  Recon muon momentum    " <<GetmuMomRec() <<std::endl
       <<"  True muon momentum     " <<GetmuMomTrue() <<std::endl
       <<"  Recon muon cos theta   " <<GetmuCosThetaRec() <<std::endl
       <<"  True muon cos theta    " <<GetmuCosThetaTrue() <<std::endl
       <<"  Recon Proton Cos Theta " <<GetpCosThetaRec() <<std::endl
       <<"  True Proton Cos Theta  " <<GetpCosThetaTrue() <<std::endl;
  }
 
private:
  Long64_t m_evid;     //unique event id
  int m_evtype;        //0 - MC, 1 - Data event
  int m_reaction;      //reaction type
  int m_sample;        //sample type (aka topology)
  double m_trueEnu;    //true nu energy
  double m_recEnu;     //recon nu energy
  double m_truePtrk;   //true track momentum
  double m_trueThtrk;  //true track angle
  double m_trueCThtrk; //true cos track angle
  double m_recPtrk;    //recon track momentum
  double m_recThtrk;   //recon track angle
  double m_recCThtrk;  //recon cos track angle
  double m_wght;       //event weight
  double m_wghtMC;       //event weight from original MC

  // New kinematic variables always included for phase space cuts
  double        m_pMomRec;
  double        m_pMomTrue;
  double        m_muMomRec;
  double        m_muMomTrue;
  double        m_muCosThetaRec;
  double        m_muCosThetaTrue;
  double        m_pCosThetaRec;
  double        m_pCosThetaTrue;

};

#endif
