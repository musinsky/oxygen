// Jan Musinsky
// 13/09/2005 (01/07/2013)

void dst2root()
{
  // for reading integer from binary DST
  union dst_integer {
    Int_t value;
    char buffer[4]; // 4 bytes
  } my_integer;

  // for reading float from binary DST
  union dst_float {
    Float_t value;
    char buffer[4]; // 4 bytes
  } my_float;

  // input
  const UShort_t nfiles = 5;
  TString fnames[nfiles] = {"DUB57.DST", "DUB61OK.DST",
                            "KIS57.DST", "KIS61.DST", "KIS65.DST"};
  ifstream fin;

  // output
  TFile fout("oxygen_DST.root", "RECREATE");
  TTree *oxy_DST = new TTree("oxy_DST", "Oxygen DST");
  const UShort_t MAX_TRACK = 50;
  Int_t   ntracks;
  Int_t   flag[MAX_TRACK];
  Float_t length[MAX_TRACK];
  Float_t deltap[MAX_TRACK];
  Float_t px[MAX_TRACK];
  Float_t py[MAX_TRACK];
  Float_t pz[MAX_TRACK];
  oxy_DST->Branch("ntracks",&ntracks,"ntracks/I");
  oxy_DST->Branch("flag",flag,"flag[ntracks]/I");
  oxy_DST->Branch("length",length,"length[ntracks]/F");
  oxy_DST->Branch("deltap",deltap,"deltap[ntracks]/F");
  oxy_DST->Branch("px",px,"px[ntracks]/F");
  oxy_DST->Branch("py",py,"py[ntracks]/F");
  oxy_DST->Branch("pz",pz,"pz[ntracks]/F");

  Int_t nw;
  Float_t b[5+5*MAX_TRACK]; // header info + tracks info
  Int_t total_event = 0;
  Int_t file_event  = 0;
  Float_t pmod, dip, azim;

  for (UShort_t f = 0; f < nfiles; f++) { // loop over all files
    if (gSystem->AccessPathName(fnames[f].Data(), kReadPermission)) {
      Printf("File %s can not be opened", fnames[f].Data());
      continue;
    }

    fin.open(fnames[f].Data());
    file_event = 0;

    do { // loop over all events in file
      fin.read(my_integer.buffer,4);
      if (!fin.good()) break;

      nw = my_integer.value;
      ntracks = (nw-5)/5;
      if (nw%5 != 0 || ntracks >= MAX_TRACK) {
        Printf("Problem with reading DST");
        return;
      }

      for (Int_t i = 0; i < nw; i++) {
        fin.read(my_float.buffer,4);
        b[i] = my_float.value;
      }
      file_event++;

      for (Int_t j = 0; j < ntracks; j++) {
        flag[j]   = b[9+5*j]/1000;
        length[j] = b[9+5*j]-1000*flag[j];
        deltap[j] = b[8+5*j];
        pmod      = b[5+5*j];
        dip       = b[6+5*j];
        azim      = b[7+5*j];
        px[j]     = pmod*TMath::Cos(dip)*TMath::Cos(azim);
        py[j]     = pmod*TMath::Cos(dip)*TMath::Sin(azim);
        pz[j]     = pmod*TMath::Sin(dip);
        if (j == 0) { // primary track (beam)
          px[j] = -px[j];
          py[j] = -py[j];
          pz[j] = -pz[j];
        }
      }

      oxy_DST->Fill();
    } while (1);

    fin.close();
    fin.clear();

    total_event += file_event;
    Printf("in file %s number of events = %d", fnames[f].Data(), file_event);
  }

  oxy_DST->Write();
  fout.Close();
  Printf("total number of events = %d", total_event);
}
