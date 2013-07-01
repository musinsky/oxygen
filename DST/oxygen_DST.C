////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

// Make oxygen ROOT files with input values

void oxygen_DST()
{
  UInt_t nevent=0;
  UShort_t nwords=0;
  UShort_t const SIZE = 500;
  Float_t events[SIZE];
  UShort_t const MAX_TRACK = 50;
  Float_t ptot,vx0,vy0,vz0,delptot;
  UShort_t flag,nsec,nfrag1m,np;
  UShort_t nfrag2,nfrag3,nfrag4,nfrag5,nfrag6,nfrag7,nfrag8;
  UShort_t npip,npim,nfrag1;
  UInt_t nchar=0; // !!! pokial je to index of array in TTree musi byt integer !!!
  UShort_t ident[MAX_TRACK];
  Float_t len[MAX_TRACK];
  Float_t px[MAX_TRACK];
  Float_t py[MAX_TRACK]; 
  Float_t pz[MAX_TRACK];
  Float_t delp[MAX_TRACK];

  ifstream fin("oxygen_DST.asc");
  if (!fin)
    {
      cout << "Sorry can not open your DST file !!!" << endl;
      break;
    }

  TFile fout("oxygen_DST.root","RECREATE");
  TTree *oxy_DST = new TTree("oxy_DST","Oxygen Tree"); // pozor takto s * je OK !!!

  oxy_DST->Branch("ptot",&ptot,"ptot/F");
  oxy_DST->Branch("delptot",&delptot,"delptot/F");
  oxy_DST->Branch("flag",&flag,"flag/s");
  oxy_DST->Branch("nsec",&nsec,"nsec/s");
  oxy_DST->Branch("vx0",&vx0,"vx0/F");
  oxy_DST->Branch("vy0",&vy0,"vy0/F");
  oxy_DST->Branch("vz0",&vz0,"vz0/F");
  oxy_DST->Branch("nfrag1m",&nfrag1m,"nfrag1m/s");
  oxy_DST->Branch("np",&np,"np/s");
  oxy_DST->Branch("nfrag2",&nfrag2,"nfrag2/s");
  oxy_DST->Branch("nfrag3",&nfrag3,"nfrag3/s");
  oxy_DST->Branch("nfrag4",&nfrag4,"nfrag4/s");
  oxy_DST->Branch("nfrag5",&nfrag5,"nfrag5/s");
  oxy_DST->Branch("nfrag6",&nfrag6,"nfrag6/s");
  oxy_DST->Branch("nfrag7",&nfrag7,"nfrag7/s");
  oxy_DST->Branch("nfrag8",&nfrag8,"nfrag8/s");
  oxy_DST->Branch("npip",&npip,"npip/s");
  oxy_DST->Branch("npim",&npim,"npim/s");
  oxy_DST->Branch("nfrag1",&nfrag1,"nfrag1/s");
  oxy_DST->Branch("nchar",&nchar,"nchar/I"); //  tento index musi byt integer
  oxy_DST->Branch("ident",ident,"ident[nchar]/s");
  oxy_DST->Branch("len",len,"len[nchar]/F");
  oxy_DST->Branch("px",px,"px[nchar]/F");
  oxy_DST->Branch("py",py,"py[nchar]/F");
  oxy_DST->Branch("pz",pz,"pz[nchar]/F");
  oxy_DST->Branch("delp",delp,"delp[nchar]/F");
  
  do
    {
      nwords=get_event(fin,events);
      if (nwords==0) break;
      ptot=events[17];
      delptot=events[18];
      flag=events[0];
      nsec=events[1];
      vx0=events[2];
      vy0=events[3];
      vz0=events[4];
      nfrag1m=events[5];
      np=events[6];
      nfrag2=events[7];
      nfrag3=events[8];
      nfrag4=events[9];
      nfrag5=events[10];
      nfrag6=events[11];
      nfrag7=events[12];
      nfrag8=events[13];
      npip=events[15];
      npim=events[14];
      nfrag1=events[16];
      nchar=events[19];
      for (UShort_t j=0; j < nchar; j++)
	{
	  ident[j]=events[20+(6*j)];
	  len[j]=events[21+(6*j)];
	  px[j]=events[22+(6*j)];
	  py[j]=events[23+(6*j)];
	  pz[j]=events[24+(6*j)];
	  delp[j]=events[25+(6*j)];
	}
      oxy_DST->Fill();
      nevent++;
    }
  while (1);
  oxy_DST->Write();
  cout << "Find " << nevent << " events" << endl;
  fin.close();
  fout.Close();
}

////////////////////////////////////////
UInt_t get_event(ifstream &fin, Float_t *values)
{
  UShort_t nrec=0;
  fin >> nrec;
  for (UShort_t i=0;i < nrec;i++)
    {
      fin >> values[i];
    }
  return nrec;
}
