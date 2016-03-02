
{
  //Configure the parameters
  //                     MAQE   MARES     MPISHP SF     Eb     pF     PDD    QEN1   QEN2   QEN3   RESN1    RESN2  COHN   NCN    NC1Pi0N   FSIINELLO FSIINELHI FSIPIPROD FSIPIABS FSICEXLO FSICEXHI
  //Nominal value
  double nom_val[21] =  {0.0,   0.0,      0.0,   0.0,   0.0,   0.0,   0.0,   0.163267, 0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   1.0,   1.0,   1.15402, 1.0,   1.0,   1.0,   0.962761};
  //NEUT default
  double neut_val[21] = {0.,       0.,       0.,       0.,      0.,      0., 0.0,   0.0,      0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   1.0,   1.0,   1.0,     1.0,   1.0,   1.0,   1.0};
  //Fit parameter?
  //bool do_param[15] =   {true,  true,     false, false, false, false, false, true,  true,  true,  true,    true,  false, false, true};
  bool do_param[21] =   {true,     true,     true,     true,   true,    true, true,  true,     true,  true,  true,  true,  true,  true,  true,  true,  true,    true,  true,  true,  true};
  //bool do_param[21] =   {true,  true,     true,  true,  true,  true,  true,  true,  true,  true,  true,    true,  true,  true,  true,    false,     false,      false,     false,   false,    false};
  //Throw parameter?
  //bool do_throw[15] =   {true,  true,     false, false, false, false, false, true,  true,  true,  true,    true,  false, false, true};
  bool do_throw[21] =   {true,  true,     true,  true,  true,  true,  true,  true,  true,  true,  true,    true,  true,  true,  true,    true,     true,      true,     true,    true,    true};
  //bool do_throw[21] =   {true,  true,     true,  true,  true,  true,  true,  true,  true,  true,  true,    true,  true,  true,  true,    false,     false,      false,     false,   false,    false};
  //Prior constraint?
  bool do_prior[21] =   {true,  true,     true,  true,  true,  true,  true,  true,  true,  true,  true,    true,  true,  true,  true,    true,     true,      true,     true,    true,    true};
  //Parameter error
  //double err_val[15] =  {0.3719,0.18286,  0.4,   1.0,   0.36,  0.14,  0.2,   0.11,  0.3,   0.3,   0.31675, 0.4,   1.0,   0.3,   0.327907};
  double err_val[21] =  {0.2,       0.2,      0.2,      0.2,    0.2, 0.2, 0.3719,   0.18286, 0.4,    1.0,   0.36,  0.14,  0.2,   0.11,  0.3,   0.3,   0.31675, 0.4,   1.0,   0.3,   0.327907};
  //Lower limit
  double low_li[21] =   {-0.8,   -0.8,     -0.8,    -0.8,    -0.8,  -0.8, -1.0,  -1.0,     -999., -999., -1.0,  -1.0,  -0.2,   0.0,   0.0,   0.0,   0.0,     0.0,   0.0,   0.0,   0.0};
  double up_li[23] =   {999.,   999.,     999.,    999.,    999.,  999., 999.,  999.,     999.,  999., 999.,  999.,  999.,   999.,   999.,   999.,   999.,     999.,   999.,   999.,   999.};
  //Parameter ID
  int parm_id[21] =     {7,       8,         9,        10,       11,      12,  0,     1,        2,     3,     4,     5,     6,     -1,    -1,    -1,    -1,      -1,    -1,    -1,    -1};
  //Interaction mode
  int int_mode[21] =    {1,       1,         1,          1,       1,       1,  0,     1,        4,     0,     0,     0,     1,      0,     0,     0,     1,      1,     2,     6,     5};
  //Use spline (1) or linear (0)
  int do_spline[21] =   {1,       1,         1,         1,        1,       1,   1,     1,        1,     0,     0,     1,     0,      -1,   -1,    -1,    -1,      -1,    -1,    -1,    -1};
  //Index for loading MB results
  int nind[21] =        { -1,       -1,        -1,       -1,       -1,       -1,   -1,    0,        -1,    -1,    -1,    -1,    -1,     -1,   -1,    -1,     1,      -1,    -1,    -1,    2};
  //Parameter name
  char name[21][21] =   {"FSIINELLO","FSIINELHI","FSIPIPROD","FSIPIABS","FSICEXLO","FSICEXHI","MAQE","MARES","MPISHP","SF","Eb","pF",  "PDD",  "NQE1","NQE2","NQE3","NRES1","NRES2","NCCCOH","NNCOTHER","NNC1PI0"}; 

  //MB fit results
  TFile *fi_niwg = new TFile("niwg_1pi_output_disappearance.root");
  TVectorD *xsec_1pi_fit_params = (TVectorD*)(fi_niwg->Get("xsec_1pi_fit_params"));
  TMatrixDSym *xsec_1pi_total_cov = (TMatrixDSym*)(fi_niwg->Get("xsec_1pi_total_cov"));

  //Use the MB fit results
  TMatrixDSym corr(21);
  for(int i=0; i<21; i++) 
    for(int j=0; j<21; j++){
      if(i==j) corr(i,j)=1.0;
      else if( nind[i]>=0 && nind[j]>=0 ){
       err_val[i] = sqrt((*xsec_1pi_total_cov)(nind[i],nind[i]));
       err_val[j] = sqrt((*xsec_1pi_total_cov)(nind[j],nind[j]));
       corr(i,j) = (*xsec_1pi_total_cov)(nind[i],nind[j])/err_val[i]/err_val[j];
       nom_val[i] = (*xsec_1pi_fit_params)(nind[i]);
       if(nind[i]==0) nom_val[i] = nom_val[i]-1.0;
      } 
      else corr(i,j) = 0.0;
    } 
       


  TMatrixDSym xsec_cov(21);

  for(int i=0; i<xsec_cov.GetNrows(); i++)
    for(int j=0; j<xsec_cov.GetNcols(); j++)
       xsec_cov(i,j) = 0.;


  TVectorD nominal(21);
  TVectorD neut(21);
  TVectorD low_limit(21);
  TVectorD up_limit(21);
  TVectorD param_id(21);
  TVectorD do_param_throw(21);
  TVectorD do_param_fit(21);
  TVectorD do_param_prior(21);
  TVectorD mode(21);
  TVectorD spline(21);
  TList name_strings;
  TObjString *name_string[21];

  int piter = 0;
  for(int i=0; i<21; i++){
    xsec_cov(piter,piter) = pow(err_val[i],2);
    nominal(piter) = nom_val[i];
    neut(piter) = neut_val[i];
    low_limit(piter) = low_li[i];
    up_limit(piter) = up_li[i];
    param_id(piter) = parm_id[i];
    mode(piter) = int_mode[i];
    do_param_throw(piter) = (do_throw[i] ? 1 : -1);
    do_param_fit(piter) = (do_param[i] ? 1 : -1);
    do_param_prior(piter) = (do_prior[i] ? 1 : -1);
    spline(piter) = do_spline[i];
    name_string[piter] = new TObjString(name[i]);
    int jiter=0;
    for(int j=0; j<21; j++) {
        xsec_cov(piter,jiter) = corr(i,j)*err_val[i]*err_val[j];
        jiter++;
    }   
    piter++;
  }  

  TFile *ffsi = new TFile("fsi_cov.root");
  TMatrixDSym *fsi_cov= (TMatrixDSym*)ffsi->Get("fsi_cov");
  for(int i=0; i<6; i++)
    for(int j=0; j<6; j++)
      xsec_cov(i,j) = (*fsi_cov)(i,j) + (i==j ? pow(0.001,2) : 0);
  ffsi->Close();
   

  Int_t nbins = 3;
  Double_t bins[4] = {0.0,1.5,3.5,30.0};

  Int_t nbinsR = 2;
  Double_t binsR[3] = {0.0,2.5,30.0};

  Int_t nbins1 = 1;
  Double_t bins1[2] = {0.0,30.0};

  TAxis ccqe_ebins(nbins,bins);
  TAxis ccres_ebins(nbinsR,binsR);
  TAxis cccoh_ebins(nbins1,bins1);
  TAxis ncdis_ebins(nbins1,bins1);
  TAxis ncres_ebins(nbins1,bins1);

  TFile *fout = new TFile("xsec_parameters_fsifirst.root","RECREATE");
  xsec_cov.Write("xsec_cov");
  nominal.Write("xsec_param_prior");
  neut.Write("xsec_param_nom");
  low_limit.Write("xsec_param_lb");
  up_limit.Write("xsec_param_ub");
  param_id.Write("xsec_param_id");
  mode.Write("xsec_int_mode");
  spline.Write("xsec_param_spline");
  do_param_throw.Write("xsec_param_throw");
  do_param_fit.Write("xsec_param_fit");
  do_param_prior.Write("xsec_param_constrain");
  ccqe_ebins.Write("ccqe_ebins");
  ccres_ebins.Write("cc1pi_ebins");
  ncres_ebins.Write("nc1pi0_ebins");
  cccoh_ebins.Write("cccoh_ebins");
  ncdis_ebins.Write("ncother_ebins");
  for(int i=0; i<piter; i++) (name_string[i])->Write(Form("xsec_param_name_%d",i));
  fout->Close();

}

