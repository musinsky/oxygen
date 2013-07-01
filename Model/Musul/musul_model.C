////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

// Make oxygen ROOT files from Musulmanbekov model

void musul_model()
{
  UInt_t nevent=0;
  UShort_t const MAX_TRACK = 50;
  UInt_t nfrag=0; // !!! pokial je to index of array in TTree musi byt integer !!!
  Float_t px[MAX_TRACK];
  Float_t py[MAX_TRACK];
  Float_t pz[MAX_TRACK];
  Float_t mass[MAX_TRACK];
  Float_t charge[MAX_TRACK];
  Float_t source[MAX_TRACK];

  ifstream fin("musul_model.asc");
  if (!fin)
    {
      cout << "Sorry can not open your DST-model file !!!" << endl;
      break;
    }

  TFile fout("musul_model.root","RECREATE");
  TTree *oxy_model = new TTree("oxy_model","Oxygen model Tree"); // pozor s * je OK !!!

  oxy_model->Branch("nfrag",&nfrag,"nfrag/I"); //  tento index musi byt integer
  oxy_model->Branch("px",px,"px[nfrag]/F");
  oxy_model->Branch("py",py,"py[nfrag]/F");
  oxy_model->Branch("pz",pz,"pz[nfrag]/F");
  oxy_model->Branch("mass",mass,"mass[nfrag]/F");
  oxy_model->Branch("charge",charge,"charge[nfrag]/F");
  oxy_model->Branch("source",source,"source[nfrag]/F"); 

  // Read header
  Float_t temp[9]=0;
  for (UShort_t tem=0; tem<9; tem++)
    fin >> temp[tem];
  cout << "Your model: " << endl;
  cout << "----------------------------------------------------------" << endl;
  if (temp[1]==0) cout << "FERMI-GAS MODEL " << endl;
  else if (temp[1]==1) cout << "FCC LATTICE MODEL " << endl;
  else cout << " !!! " << endl;
  if (temp[8]>=0) cout << "INPUT VALUE OF IMPACT PARAMETER [FM]" << endl;
  else cout << "RANDOM VALUE OF IMPACT PARAMETER THAT IS GIVEN BY PROGRAM" << endl;
  cout << "----------------------------------------------------------" << endl;
  cout << "Number of events from model file is = " << temp[7] << endl;
  cout << "Beam kinetic energy per nucleon is = " << temp[6] << " GeV/c" << endl;
  cout << "Beam     Z = " << temp[4] << "   A = " << temp[2] << endl;
  cout << "Target   Z = " << temp[5] << "   A = " << temp[3] << endl;

  //  transformation from AntiLAB to LAB
  //  collison axis is Z, beam is proton
  //  Pbeam[] = (0,0,3.257,3.39)
  //  Pbeam[4] = 0.940+Ekin    ; Ekin=2.45 GeV/c
  //  Pbeam[3] = sqrt(Pbeam[4]*Pbeam[4] - 0.94*0.94)
  Float_t Ekin=temp[6];
  Float_t Pbeam4=0.940+Ekin;
  Float_t Pbeam3=sqrt(Pbeam4*Pbeam4-0.940*0.940);
  Float_t beta=Pbeam3/Pbeam4;
  Float_t gamma=1/sqrt(1-beta*beta);
  Float_t anti_energy=0;
  
  // Event by event
  UShort_t ncas;
  do
    {
      fin >> ncas;
      if (!fin.good()) break;
      fin >> temp[1]; fin >> temp[2]; fin >> temp[3]; // pimpac, mvs1, mvs2
      nfrag=ncas/6;
      for (UShort_t i=0; i< nfrag; i++)
	{
	  fin >> px[i];
	  fin >> py[i];
	  fin >> pz[i];
	  fin >> mass[i];
	  fin >> charge[i];
	  fin >> source[i];
	  // ALABtoLAB
	  anti_energy = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]+mass[i]*mass[i]);
	  pz[i] = gamma*(pz[i]-beta*anti_energy);
	}
      oxy_model->Fill();
      nevent++;
    }
  while (1);
  oxy_model->Write();
  cout << "Find " << nevent << " events" << endl;
  if (temp[7]!=nevent) cout << "Num events is not to same !!!" << endl;
  fin.close();
  fout.Close();
}
