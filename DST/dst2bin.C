// Jan Musinsky
// 05/04/2004 (28/06/2013)

void dst2bin(UInt_t max = 99999)
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

  ifstream fin("DUB57.DST"); // input binary DST file
  ofstream fout("dub57.bin", ios::out | ios::binary); // output binary file

  UInt_t nevent = 0;
  Int_t nw;  // see notes about OXYGEN DST
  Float_t b; // see notes about OXYGEN DST

  do {
    // dlina rekorda
    fin.read (my_integer.buffer,4);

    if (!fin.good()) break;
    //    cout << "event = " << nevent << endl;

    nw = my_integer.value;
    //    cout << "nw  = " << nw << endl << endl;
    fout.write((char*)&nw,sizeof(nw));
    // fout.write((char*)&nw,4);

    for (Int_t i = 0; i < nw; i++) {
      fin.read (my_float.buffer,4);
      b = my_float.value;
      //      cout << Form("b(%2d) = ", i) << b << endl;
      fout.write((char*)&b,sizeof(b));
      // fout.write((char*)&b,4);
    }

    nevent++;
    //    cout << "##############################" << endl;
  } while (nevent < max);

  fin.close();
  fout.close();
  cout << "Number of total events = " << nevent << endl;
}
