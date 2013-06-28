// 06-apr-2004

void DST2FOR()
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
  
  const UShort_t NF=5;
  Text_t *input_names[NF]={"DUB57.DST","DUB61OK.DST",
			   "KIS57.DST","KIS61.DST","KIS65.DST"};
  Text_t *output_names[NF]={"dub57.bin","dub61.bin",
			    "kis57.bin","kis61.bin","kis65.bin"};
  
  ifstream fin;  // input binary DST files
  ofstream fout; // output binary bin file
  
  UInt_t total_event=0;
  UInt_t current_event=0;
  Int_t nw;
  
  for (UShort_t f=0; f<NF; f++) {         // loop over all files
    fin.open(input_names[f]);
    fout.open(output_names[f], ios::out | ios::binary );
    current_event=0;
    do {
      //    if (nevent == 10 ) break;
      // size of record
      fin.read(my_integer.buffer,4); 
      if (!fin.good()) break;
      nw = my_integer.value;
      my_integer.value = (nw+1)*4;        // byte counter for FORTRAN binary
      fout.write(my_integer.buffer,4);    // goes first
      my_integer.value = nw ;
      fout.write(my_integer.buffer,4);
      for (UShort_t i=0; i < nw ; i++) {
	fin.read(my_float.buffer,4);
	fout.write(my_float.buffer,4);
      }
      my_integer.value = (nw+1)*4;        // byte counter for FORTRAN binary
      fout.write(my_integer.buffer,4);    // goes last
      current_event++;
    }
    while (1);
    fin.close();
    fin.clear();
    fout.close();
    fout.clear();
    total_event+=current_event;
    cout << "In file " << input_names[f] << " number of events = "
	 << current_event << endl;
  }
  cout << "Number of total events = " << total_event << endl;
}
