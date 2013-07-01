////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

void histo_correl_DST_old()
{
  const UShort_t TOPOLOGY=1; // !!! (to same topology)
  TFile *fin = new TFile("oxygen_DST.root");
  TTree *oxy_DST = (TTree*)fin->Get("oxy_DST");
  
  UShort_t const MAX_TRACK = 50;
  UInt_t nchar=0;
  UShort_t ident[MAX_TRACK];
  Float_t px[MAX_TRACK];
  Float_t py[MAX_TRACK]; 
  Float_t pz[MAX_TRACK];
  oxy_DST->SetBranchAddress("nchar",&nchar);
  oxy_DST->SetBranchAddress("ident",ident);
  oxy_DST->SetBranchAddress("px",px);
  oxy_DST->SetBranchAddress("py",py);
  oxy_DST->SetBranchAddress("pz",pz); 
  oxy_DST->SetBranchStatus("*",0);
  oxy_DST->SetBranchStatus("nchar",1);
  oxy_DST->SetBranchStatus("ident",1);
  oxy_DST->SetBranchStatus("px",1);
  oxy_DST->SetBranchStatus("py",1);
  oxy_DST->SetBranchStatus("pz",1);
  // tree->GetEntry(random) bez virtual size velmi pomale
  oxy_DST->SetMaxVirtualSize(64000000); // 64 MB
  
  Float_t tot_p=0;
  const UShort_t MAX_NUM_P=MAX_TRACK; // max. num. protons
  Float_t px_true[MAX_NUM_P]={0}; // for true pairs
  Float_t py_true[MAX_NUM_P]={0};
  Float_t pz_true[MAX_NUM_P]={0};
  Float_t px_back[MAX_NUM_P]={0}; // for background pairs
  Float_t py_back[MAX_NUM_P]={0};
  Float_t pz_back[MAX_NUM_P]={0};
  
  UInt_t back_event=0;
  UInt_t nevents = oxy_DST->GetEntries();

  TStopwatch timer;
  timer.Start();
  
  for (UInt_t i=0; i<nevents; i++)
    {
      if (i%1000==0) cout << "Calculating = " << i << " events" << endl;
      oxy_DST->GetEntry(i);
      UShort_t count=0; // reset

      for (UShort_t j=0; j<nchar; j++) // every fragment
	{
	  if (ident[j] == FRAGMENT)
	    {
	      px_true[count]=px[j];
	      py_true[count]=py[j];
	      pz_true[count]=pz[j];
	      count++;
	    }
	  
	  if (ident[j]==FRAGMENT_IMP) // nejednoznacne jedno-naboje fragmenty
	    {
	      tot_p = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
	      if ( (tot_p > MIN_IMP) && (tot_p < MAX_IMP) ) // proton
		{
		  px_true[count]=px[j];
		  py_true[count]=py[j];
		  pz_true[count]=pz[j];
		  count++;
		}
	    }
	} // every fragment
      
      if ( count > 1 ) // two or more protons
	{
	  // mixing true pairs protons
	  for (UShort_t k=0; k < count-1; k++)
	    for (UShort_t l=k+1; l <= count-1; l++)
	      h_true->Fill( mix(px_true[k],py_true[k],pz_true[k],px_true[l],py_true[l],pz_true[l]) );
	} // more then one proton
      
      // finding background events
      for (UShort_t nb=0; nb < NUM_BACK_EVENTS; nb++)
	{
	  back_event=(gRandom->Rndm())*nevents;
	  // not to same event
	  if (back_event == i) back_event=(gRandom->Rndm())*nevents;
	  oxy_DST->GetEntry(back_event);
	  UShort_t count_b=0; // reset

	  for (UShort_t j=0; j<nchar; j++) // every fragment
	    {
	      if (ident[j]==FRAGMENT)
		{
		  px_back[count_b]=px[j];
		  py_back[count_b]=py[j];
		  pz_back[count_b]=pz[j];
		  count_b++;
		}
	      
	      if (ident[j]==FRAGMENT_IMP) // nejednoznacne jedno-naboje fragmenty
		{
		  tot_p = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
		  if ( (tot_p > MIN_IMP) && (tot_p < MAX_IMP) ) // proton
		    {
		      px_back[count_b]=px[j];
		      py_back[count_b]=py[j];
		      pz_back[count_b]=pz[j];
		      count_b++;
		    }
		}
	    } // every fragment
	  
	  if (TOPOLOGY) // to same topology
	    if ( count_b != count ) continue;
	  if ( (count_b>1) && (count>1) ) // two or more protons
	    {
	      // mixing background pairs protons
	      for (UShort_t k=0; k < count; k++) // every true proton
		for (UShort_t l=0; l < count_b; l++) // every back proton
		  h_back->Fill( mix(px_true[k],py_true[k],pz_true[k],px_back[l],py_back[l],pz_back[l]) );
	    } // more then one proton (back)
	} // end num back protons

    } // end evey event
  cout << "Calculating time is = " << timer.RealTime() << " sec" << endl;
  cout << endl << "Total events = " << i << endl;
}
////////////////////////////////////////
Float_t mix(Float_t px1,Float_t py1,Float_t pz1,Float_t px2,Float_t py2,Float_t pz2)
{
  // mixing pairs of protons
  Float_t temp_p = (px1-px2)**2+(py1-py2)**2+(pz1-pz2)**2;
  Float_t temp_e1= sqrt(px1*px1+py1*py1+pz1*pz1+MASS*MASS);
  Float_t temp_e2= sqrt(px2*px2+py2*py2+pz2*pz2+MASS*MASS);
  Float_t temp_e= (temp_e1-temp_e2)**2;
  return ( sqrt(temp_p-temp_e) / 2 );
}
