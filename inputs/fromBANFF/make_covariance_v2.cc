{

   //File with 2000 detector systematic throws
   TFile *fi = new TFile("throws_det_mcstat_ccfix.root");

   //File with no throws (only 1 toy)
   TFile *fin = new TFile("throws_nom_det_ccfix.root");

   TH2D *mean[3];
   TH2D *nom[3];
   mean[0] = (TH2D*)fi->Get("cc0pi_hist_0");
   mean[1] = (TH2D*)fi->Get("cc1pi_hist_0");
   mean[2] = (TH2D*)fi->Get("ccNpi_hist_0");

   for(int i=1; i<2000; i++){
     TH2D *tmp = (TH2D*)fi->Get(Form("cc0pi_hist_%d",i));
     mean[0]->Add(tmp);
     tmp = (TH2D*)fi->Get(Form("cc1pi_hist_%d",i));
     mean[1]->Add(tmp);
     tmp = (TH2D*)fi->Get(Form("ccNpi_hist_%d",i));
     mean[2]->Add(tmp);
   }

   mean[0]->Scale(1./2000.);
   mean[1]->Scale(1./2000.);
   mean[2]->Scale(1./2000.);

   nom[0] = (TH2D*)fin->Get("cc0pi_hist_0");
   nom[1] = (TH2D*)fin->Get("cc1pi_hist_0");
   nom[2] = (TH2D*)fin->Get("ccNpi_hist_0");

   /*mean[0] = (TH2D*)fin->Get("cc0pi_hist_0");
   mean[1] = (TH2D*)fin->Get("cc1pi_hist_0");
   mean[2] = (TH2D*)fin->Get("ccNpi_hist_0");*/
  
  
   TMatrixDSym *cov = new TMatrixDSym(mean[0]->GetNbinsX()*mean[0]->GetNbinsY()+
                                      mean[1]->GetNbinsX()*mean[1]->GetNbinsY()+
                                      mean[2]->GetNbinsX()*mean[2]->GetNbinsY() );

   TH1D *hist_weights = new TH1D("hist_weights","",50,0.5,1.5);

   TVectorD *weights = new TVectorD(cov->GetNrows());
   std::cout << weights->GetNrows() << std::endl;

   char names[3][50] = {"cc0pi_hist_","cc1pi_hist_","ccNpi_hist_"};

   for(int t=0; t<2000; t++){
     int iter1=0;
     for(int s1=0; s1<3; s1++){
      TH2D *tmp1 = (TH2D*)fi->Get(Form("%s%d",names[s1],t));
      for(int p1=1; p1<=mean[s1]->GetNbinsX(); p1++)
        for(int c1=1; c1<=mean[s1]->GetNbinsY(); c1++){
          double delta1 = (1.0-tmp1->GetBinContent(p1,c1)/mean[s1]->GetBinContent(p1,c1));
          //std::cout << "Set Weight " << s1 << " " <<  p1 << " " << c1  << " " << iter1 << std::endl;
          //std::cout << mean[s1]->GetBinContent(p1,c1) << std::endl;
          //std::cout << nom[s1]->GetBinContent(p1,c1) << std::endl;
          if(t==0){
           (*weights)(iter1) = mean[s1]->GetBinContent(p1,c1)/nom[s1]->GetBinContent(p1,c1);
           hist_weights->Fill((*weights)(iter1));
          }
          int iter2=0;
          for(int s2=0; s2<3; s2++){
            TH2D *tmp2 = (TH2D*)fi->Get(Form("%s%d",names[s2],t));
            for(int p2=1; p2<=mean[s2]->GetNbinsX(); p2++)
              for(int c2=1; c2<=mean[s2]->GetNbinsY(); c2++){
                double delta2 = (1.0-tmp2->GetBinContent(p2,c2)/mean[s2]->GetBinContent(p2,c2));
                (*cov)(iter1,iter2)+=delta1*delta2/2000.;
                iter2++;
              }
           }    
           iter1++;
         }
     }
   }         

   cov->Draw("colz");

   TFile *fout = new TFile("nd280_det_covariance.root","RECREATE");
   cov->Write("nddet_cov");
   mean[0]->GetXaxis()->Write("p_axis");
   mean[0]->GetYaxis()->Write("th_axis");
   weights->Write("det_weights");
   hist_weights->Write();
   fout->Close();

}







