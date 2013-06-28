// -*- mode: c++ -*-
// @(#) 12 Sep 2005

void convertDST(Bool_t ascii = kTRUE)
{
  // for reading integer from binary
  union dst_integer {
    Int_t value;
    char buffer[4]; // 4 bytes
  } my_integer;
  
  // for reading float from binary
  union dst_float {
    Float_t value;
    char buffer[4]; // 4 bytes
  } my_float;
  
  printf("create %s files\n", (ascii) ? "ascii" : "binary");
  const UShort_t nfiles = 5;
  TString fnames[nfiles] = {"DUB57.DST", "DUB61.DST",
			    "KIS57.DST", "KIS61.DST", "KIS65.DST"};
  ifstream fin;  // input binary DST file
  ofstream fout; // output file (binary or ascii)
  Int_t total_event = 0;
  Int_t file_event  = 0;
  Int_t nw; // nw = 5 + 5*ntracks;
  Float_t b;

  for (UShort_t f = 0; f < nfiles; f++) { // loop over all files
    if (gSystem->AccessPathName(fnames[f].Data(),kReadPermission)) {
      printf("File %s can not be opened\n",fnames[f].Data());
      continue;
    }
    fin.open(fnames[f].Data());
    fnames[f].ToLower();
    if (ascii) {
      fnames[f].Replace(6,3,"asc");
      fout.open(fnames[f].Data(), ios::out );
    }
    else {
      fnames[f].Replace(6,3,"bin");
      fout.open(fnames[f].Data(), ios::out | ios::binary );
    }
    
    file_event = 0;
    do { // loop over all events in file
      fin.read(my_integer.buffer,4);
      if (!fin.good()) break;
      nw = my_integer.value;
      if (nw%5 != 0) printf("Size of record is corrupted\n");
      if (ascii) fout << nw << " ";
      else fout.write((char*)&nw,sizeof(nw));
      
      for (Int_t i = 0; i < nw; i++) {
	fin.read (my_float.buffer,4);
	b = my_float.value;
	if (ascii) fout << b << " ";
	else fout.write((char*)&b,sizeof(b));
      }
      if (ascii) fout << endl;
      file_event++;
    }
    while (1);
    fin.close(); fin.clear();
    fout.close(); fout.clear();
    total_event += file_event;
    printf("in file %s number of events = %d\n",fnames[f].Data(),file_event);
  }
  printf("total number of events = %d\n",total_event);
}
