// 05-apr-2004

void DST2ASCII(Bool_t bin=1)
{
  // for reading integer from binary
  union dst_integer {
    Int_t value;
    char buffer[4]; // 4 bytes
  } my_integer;
  
  // for reading float from binary
  union dst_float {
    Float_t value;
    char buffer[8]; // 4 bytes
  } my_float;
  
  ifstream fin("DUB57.DST"); // input binary DST file
  
  ofstream fout("dub57.bin",ios::out | ios::binary);
  
  //  FILE *fout;
  //  fout = fopen("dub57.bin","w");
  
  UInt_t nevent=0;
  int nw;  // see notes about OXYGEN DST
  float b; // see notes about OXYGEN DST
  
  do {
    if (nevent == 1 ) break; // only first 10 events

    // dlina rekorda
    fin.read (my_integer.buffer,4); 
    
    if (!fin.good()) break;
    cout << "event = " << nevent << endl;
    
    nw = my_integer.value;
    
    cout << "nw  = " << nw << endl << endl;
    //    fout.write((char*)&nw,sizeof(nw));
    //fout << "AAAA" << endl;
    //cout << sizeof(nw) << endl;
    fout.write((char*)&nw,sizeof(nw));
    for (Int_t i=0; i<nw; i++) {
      fin.read (my_float.buffer,4);
      b = my_float.value;
      cout << "b(1) = " << b << endl;
      fout.write((char*)&b,sizeof(b));
      //cout << sizeof(b) << endl;
      //fout.write((char*)&b,4);
    }
    
    
    
    nevent++;
    cout << "--------------------------------------" << endl;
  }
  while (1);
  fin.close();
  fout.close();
  cout << "Number of total events = " << nevent << endl;
}
