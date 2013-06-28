// 02-apr-2004

void DST2ASCII_full()
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
  
  ifstream fin("DUB57.DST"); // input binary DST file
  
  UInt_t nevent=0;
  Int_t nw;  // see notes about OXYGEN DST
  Float_t b; // see notes about OXYGEN DST
  
  do {
    //    if (nevent == 10 ) break; // only first 10 events

    // dlina rekorda
    fin.read (my_integer.buffer,4); 

    if (!fin.good()) break;
    cout << "event = " << nevent << endl;

    nw = my_integer.value;
    cout << "nw  = " << nw << endl << endl;
    // b(1)   ( !!! je to float a nie integer !!!)
    fin.read (my_float.buffer,4);
    b = my_float.value;
    cout << "b(1) = " << b << endl;
    // b(2)   ( !!! je to float a nie integer !!!)
    fin.read (my_float.buffer,4);
    b = my_float.value;
    cout << "b(2) = " << b << endl;
    // b(3), b(4), b(5)
    fin.read (my_float.buffer,4);
    b = my_float.value;
    cout << "b(3) = " << b << endl;
    fin.read (my_float.buffer,4);
    b = my_float.value;
    cout << "b(4) = " << b << endl;
    fin.read (my_float.buffer,4);
    b = my_float.value;
    cout << "b(5) = " << b << endl;
      
    // dal'she s b(6) do b(nw) po 5 slov na kazhdogo treka
    for (Int_t j=0; j<(nw-5)/5; j++) {
      cout << endl;
      // b(6+5*j)
      fin.read (my_float.buffer,4);
      b = my_float.value;
      cout << Form("b( 6+5*%u) = ",j)  << b << endl;
      // b(7+5*j)
      fin.read (my_float.buffer,4);
      b = my_float.value;
      cout << Form("b( 7+5*%u) = ",j)  << b << endl;
      // b(8+5*j)
      fin.read (my_float.buffer,4);
      b = my_float.value;
      cout << Form("b( 8+5*%u) = ",j)  << b << endl;
      // b(9+5*j)
      fin.read (my_float.buffer,4);
      b = my_float.value;
      cout << Form("b( 9+5*%u) = ",j)  << b << endl;
      // b(10+6*j)
      fin.read (my_float.buffer,4);
      b = my_float.value;
      cout << Form("b(10+5*%u) = ",j)  << b << endl;
    }
    nevent++;
    cout << "--------------------------------------" << endl;
  }
  while (1);
  fin.close();
  cout << "Number of total events = " << nevent << endl;
}
