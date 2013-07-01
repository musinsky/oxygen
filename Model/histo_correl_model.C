////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

void histo_correl_model(const Text_t *input_f, TH1F *true_h, TH1F *back_h)
{
  TFile *fin = new TFile(input_f);
  cout << "Open -> " << input_f << endl;
  TTree *oxy_model = (TTree*)fin->Get("oxy_model");
  
  UShort_t const MAX_TRACK = 50;
  UInt_t nfrag=0;
  Float_t px[MAX_TRACK];
  Float_t py[MAX_TRACK]; 
  Float_t pz[MAX_TRACK];
  Float_t mass[MAX_TRACK];
  Float_t charge[MAX_TRACK];
  //  Float_t source[MAX_TRACK];
  oxy_model->SetBranchAddress("nfrag",&nfrag);
  oxy_model->SetBranchAddress("px",px);
  oxy_model->SetBranchAddress("py",py);
  oxy_model->SetBranchAddress("pz",pz);
  oxy_model->SetBranchAddress("mass",mass);
  oxy_model->SetBranchAddress("charge",charge);
  //  oxy_model->SetBranchAddress("source",source);
  
  // tree->GetEntry(random) bez virtual size velmi pomale
  oxy_model->SetMaxVirtualSize(64000000); // 64 MB

  if (nevents==0)
    nevents = oxy_model->GetEntries();
  const UShort_t MAX_NUM_P=MAX_TRACK;
  UShort_t count=0;
  UShort_t index[MAX_NUM_P]={0};
  Float_t tot_p=0;
  Int_t multi[MAX_NUM_P]={0};
  
  TStopwatch timer;
  timer.Start();
  
  // mixing true pairs fragments
  for (UInt_t i=0; i<nevents; i++)
    {
      oxy_model->GetEntry(i);
      count=0; // reset
      for (UShort_t j=0; j<nfrag; j++) // every fragment
	{
 	  if (charge[j]==1 && mass[j]>0.93 && mass[j]<0.95) // proton
	    {
	      tot_p = sqrt(px[j]**2+py[j]**2+pz[j]**2);
	      if ( (tot_p > MIN_IMP) && (tot_p < MAX_IMP) )
		index[count++]=j;
	    }
	} // every fragment
      
      if ( (count > 1) && (count <= MAX_NUM_FRAG) ) // how many fragment
	{
	  // mixing pairs of fragments
	  for (UShort_t k=0; k < count-1; k++)
	    for (UShort_t l=k+1; l <= count-1; l++)
	      true_h->Fill( mix(px[index[k]],py[index[k]],pz[index[k]],px[index[l]],py[index[l]],pz[index[l]]) );
	} // end how many fragment
      multi[count]++;
    } // every event

  UInt_t nevents_back=0;
  UInt_t count_back_events=0;
  UInt_t n_b_event=0;
  Float_t tot_p_back[MAX_NUM_P]={0};
  UShort_t topology_e=0;
  Int_t sort_index[MAX_NUM_P]={0}; // musi byt Int_t (TMath::Sort)
  Float_t px_back[MAX_NUM_P]={0}; // for background pairs
  Float_t py_back[MAX_NUM_P]={0};
  Float_t pz_back[MAX_NUM_P]={0};
  Int_t verify_e[MAX_NUM_P]={0};
  UShort_t to_same_event=0;
  
  // mixing background pairs fragments
  for (UShort_t m=2; m <= MAX_NUM_FRAG; m++)
    {
      nevents_back = multi[m];
      if (nevents_back==0) continue;
      cout << "Multiplicity = " << m << " -> " << nevents_back << " ";
      cout << " * " << NUM_BACK_EVENTS << " = ";
      nevents_back = nevents_back*NUM_BACK_EVENTS;
      cout << nevents_back;
      if (m >= (nevents_back/NUM_BACK_EVENTS) )
	{
	  cout << " !!! LOOP !!! h_true filled with this multi" << endl;
	  continue;
	}
      cout << endl;
      
      count_back_events=0;
      while (count_back_events < nevents_back)
	{
	  topology_e = 0;
	  while (topology_e < m)
	    {
	      n_b_event = (gRandom->Rndm())*nevents;
	      oxy_model->GetEntry( n_b_event );
	      count=0;
	      for (UShort_t j=0; j<nfrag; j++) // every fragment
		{
		  if (charge[j]==1 && mass[j]>0.93 && mass[j]<0.95) // proton
		    {
		      tot_p = sqrt(px[j]**2+py[j]**2+pz[j]**2);
		      if ( (tot_p > MIN_IMP) && (tot_p < MAX_IMP) )
			{
			  tot_p_back[count] = tot_p;
			  index[count++]=j;
			}
		    }
		} // every fragment
	      
	      // je ta topologia ktoru teraz pocitame ?
	      if (count == m)
		{
		  // potrebny len i-ty fragment (zoradeny)
		  TMath::Sort(count,tot_p_back,sort_index);
		  px_back[topology_e] = px[index[sort_index[topology_e]]];
		  py_back[topology_e] = py[index[sort_index[topology_e]]];
		  pz_back[topology_e] = pz[index[sort_index[topology_e]]];
		  
		  // potrebny len i-ty fragment (!!! nezoradeny !!!)
		  //px_back[topology_e] = px[topology_e];
		  //py_back[topology_e] = py[topology_e];
		  //pz_back[topology_e] = pz[topology_e];
		  
		  verify_e[topology_e]=n_b_event;
		  topology_e++;
		}
	    } // pokial pocet eventov nie je taky isty ako je multiplicita
	  
	  // preverenie aby kazdy fragment bol z ineho eventu
	  to_same_event=0;
	  for (UShort_t kk=0; kk < m-1; kk++)
	    for (UShort_t ll=kk+1; ll <= m-1; ll++)
	      if ( verify_e[kk] == verify_e[ll] ) to_same_event++;
	  if (to_same_event == 0)
	    {
	      
	      // mame tolko fragmentov (eventov) kolko je multiplicita
	      // t.j. akoby mame jeden event s danym poctom fragmentov
	      for (UShort_t k=0; k < m-1; k++)
		for (UShort_t l=k+1; l <= m-1; l++)
		  back_h->Fill( mix(px_back[k],py_back[k],pz_back[k],px_back[l],py_back[l],pz_back[l]) );
	      
	      count_back_events++;
	      
	    }
	  
	} // pokial pocet eventov v danej multiplicite nie je ten isty ako true
      
    } // every multiplicity
  
  cout << "Calculating time is = " << timer.RealTime() << " sec" << endl;
  cout << endl << "Total events = " << i << endl;
}
////////////////////////////////////////
Float_t mix(Float_t px1,Float_t py1,Float_t pz1,Float_t px2,Float_t py2,Float_t pz2)
{
  // mixing pairs of fragment
  Float_t temp_p = (px1-px2)**2+(py1-py2)**2+(pz1-pz2)**2;
  Float_t temp_e1= sqrt(px1**2+py1**2+pz1**2+MASS**2);
  Float_t temp_e2= sqrt(px2**2+py2**2+pz2**2+MASS**2);
  Float_t temp_e= (temp_e1-temp_e2)**2;
  return ( sqrt(temp_p-temp_e) / 2 );
}
