////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

// Make model_root file from "moscow" model (in bin format)

void moscow_model()
{
  UShort_t const MAX_TRACK = 50;
  UInt_t nfrag=0;
  Float_t px[MAX_TRACK];
  Float_t py[MAX_TRACK];
  Float_t pz[MAX_TRACK];
  Float_t mass[MAX_TRACK];
  Float_t charge[MAX_TRACK];
  
  TFile fout("moscow_model.root","RECREATE");
  TTree *oxy_model = new TTree("oxy_model","Oxygen model Tree");
  oxy_model->Branch("nfrag",&nfrag,"nfrag/I");
  oxy_model->Branch("px",px,"px[nfrag]/F");
  oxy_model->Branch("py",py,"py[nfrag]/F");
  oxy_model->Branch("pz",pz,"pz[nfrag]/F");
  oxy_model->Branch("mass",mass,"mass[nfrag]/F");
  oxy_model->Branch("charge",charge,"charge[nfrag]/F");
  
  // for reading integer from binary
  union mix_integer {
    Int_t value;
    char buffer[4]; // 4 bytes
  } my_integer;
  // for reading float from binary
  union mix_float {
    Float_t value;
    char buffer[4]; // 4 bytes
  } my_float;

  ifstream fin("moscow_model.bin"); // input binary DST file
  
  UInt_t nevent=0;
  Float_t p_tot, py_temp;
  Float_t cos_theta, sin_phi, cos_phi;
  
  do
    {
      if (!fin.good()) break;
      // first word, number of fragments
      fin.read (my_integer.buffer,4);
      nfrag = (my_integer.value )/6;
      // every fragment
      for (UShort_t j=0; j<nfrag; j++)
	{
	  fin.read (my_float.buffer,4); // 1 - impuls p
	  p_tot = my_float.value;
	  fin.read (my_float.buffer,4); // 2 - cos(theta)
	  cos_theta = my_float.value;
	  fin.read (my_float.buffer,4); // 3 - sin(phi)
	  sin_phi = my_float.value;
	  fin.read (my_float.buffer,4); // 4 - cos(phi)
	  cos_phi = my_float.value;
	  fin.read (my_float.buffer,4); // 5 - charge
	  charge[j] = my_float.value;
	  fin.read (my_float.buffer,4); // 6 - mass
	  mass[j] = my_float.value;
	  // transform to px, py, pz
	  px[j] = p_tot * cos_theta;
	  py_temp = p_tot * sqrt(1-cos_theta**2);
	  py[j] = py_temp * cos_phi;
	  pz[j] = py_temp * sin_phi;
	}
      nevent++;
      oxy_model->Fill();
    }
  while (1);
  oxy_model->Write();
  fout.Close();
  fin.close();
  cout << "Number of events = " << nevent << endl;
}
