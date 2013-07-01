// Jan Musinsky
// 12/09/2005 (01/07/2013)

void dst2convert(Bool_t ascii = kTRUE)
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

  Printf("creating %s files", (ascii) ? "ascii" : "binary");
  const UShort_t nfiles = 5;
  TString fnames[nfiles] = {"DUB57.DST", "DUB61OK.DST",
                            "KIS57.DST", "KIS61.DST", "KIS65.DST"};
  ifstream fin;  // input binary DST file
  ofstream fout; // output file (ascii or binary)

  Int_t total_event = 0;
  Int_t file_event  = 0;
  Int_t nw;  // nw = 5 + 5*ntracks;
  Float_t b;

  for (UShort_t f = 0; f < nfiles; f++) { // loop over all files
    if (gSystem->AccessPathName(fnames[f].Data(), kReadPermission)) {
      Printf("file %s can not be opened", fnames[f].Data());
      continue;
    }

    fin.open(fnames[f].Data());
    fnames[f].ToLower();
    if (ascii) {
      fnames[f].ReplaceAll("ok", "");
      fnames[f].Replace(6,3,"asc");
      fout.open(fnames[f].Data(), ios::out );
    }
    else {
      fnames[f].ReplaceAll("ok", "");
      fnames[f].Replace(6,3,"bin");
      fout.open(fnames[f].Data(), ios::out | ios::binary );
    }
    file_event = 0;

    do { // loop over all events in file
      fin.read(my_integer.buffer,4);
      if (!fin.good()) break;

      nw = my_integer.value;
      if (nw%5 != 0) Printf("Size of record is corrupted");
      if (ascii) fout << nw << " ";
      else       fout.write((char*)&nw,sizeof(nw));

      for (Int_t i = 0; i < nw; i++) {
        fin.read (my_float.buffer,4);
        b = my_float.value;
        if (ascii) fout << b << " ";
        else       fout.write((char*)&b,sizeof(b));
      }

      if (ascii) fout << endl;
      file_event++;
    } while (1);

    fin.close();
    fin.clear();
    fout.close();
    fout.clear();

    total_event += file_event;
    Printf("in file %s number of events = %d", fnames[f].Data(), file_event);
  }

  Printf("total number of events = %d", total_event);
}
