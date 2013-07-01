////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

// Correction model on impuls error with experimental data (only proton !!!)
// and create new root file (to same as original only correct proton impuls)
// All is LAB system

void correct_on_landau(UShort_t mod=1)
{
  // mod=1 for Moscow(default), mod==2 for Musulmanbekov
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // change only this
  // Correction factor from experiment (from Landau)
  Float_t mpv=0.0167;
  Float_t sigma=0.002527;
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  UShort_t const MAX_TRACK = 50;
  UInt_t nfrag=0;
  Float_t px[MAX_TRACK];
  Float_t py[MAX_TRACK]; 
  Float_t pz[MAX_TRACK];
  Float_t mass[MAX_TRACK];
  Float_t charge[MAX_TRACK];
  //  Float_t source[MAX_TRACK];
  
  Text_t tmp[100], tmp_c[100];
  if (mod==1) // Moscow model
    {
      cout << "Calculate from Moscow model !!!" << endl;
      sprintf(tmp,"%s","moscow_model.root");
      sprintf(tmp_c,"%s","moscow_model_c.root");
    }
  else if (mod==2) // Musulmanbekov model
    {
      cout << "Calculate from Musulmanbekov model !!!" << endl;
      sprintf(tmp,"%s","musul_model.root");
      sprintf(tmp_c,"%s","musul_model_c.root");
    }
  else
    cout << "Sorry, I dont have this model !!!" << endl;

  TFile *fin = new TFile(tmp,"READ");
  TTree *oxy_model = (TTree*)fin->Get("oxy_model");
  oxy_model->SetBranchAddress("nfrag",&nfrag);
  oxy_model->SetBranchAddress("px",px);
  oxy_model->SetBranchAddress("py",py);
  oxy_model->SetBranchAddress("pz",pz); 
  oxy_model->SetBranchAddress("mass",mass);
  oxy_model->SetBranchAddress("charge",charge);
  //  oxy_model->SetBranchAddress("source",source);

  TFile *newfile = new TFile(tmp_c,"recreate");
  TTree *oxy_model_correct = oxy_model->CloneTree(0);

  Float_t tot_p=0; // total impuls
  Float_t new_delp_vs_p=0;
  Float_t temp_gaus=0;
  Float_t new_tot_p=0;
  Float_t const MIN_P=1.75-0.0;
  Float_t const MAX_P=4.75+0.0;

  
  UInt_t nevents=0;
  nevents=oxy_model->GetEntries();
  
  for (UInt_t i=0; i<nevents; i++)
    {
      oxy_model->GetEntry(i);
      for (UShort_t j=0; j<nfrag; j++)
	{
	  // correct only proton in interval (MIN_P, MAX_P)
	  if (charge[j]==1 && mass[j]>0.93 && mass[j]<0.95) // proton
	    {
	      tot_p = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
	      if ( tot_p > MIN_P && tot_p < MAX_P ) // only our interval
		{
		  new_delp_vs_p = gRandom->Landau(mpv,sigma);
		  // rozmyvat aj na gauss alebo len +- 1 !!!
		  temp_gaus=gRandom->Gaus(0,1);
		  if (temp_gaus>0) temp_gaus=1.0;
		  else temp_gaus=-1.0;
		  new_tot_p = tot_p + (temp_gaus*tot_p*new_delp_vs_p);
		  px[j] = (px[j]/tot_p) * new_tot_p;
		  py[j] = (py[j]/tot_p) * new_tot_p;
		  pz[j] = (pz[j]/tot_p) * new_tot_p;

		} // if impuls interval
	    } // if proton
	} // every fragment
      oxy_model_correct->Fill();
    } // every event

  oxy_model_correct->Write();
  newfile->Close();
  fin->Close();
}
