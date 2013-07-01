// root macro
{
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // change only this
  Float_t ZR_min = 1.0;
  Float_t ZR_max = 3.0;
  Float_t ZR_step= 0.1;
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Text_t tmp[100];
  for (Float_t r = ZR_min; r <= ZR_max; r += ZR_step)
    {
      sprintf(tmp,"echo '%0.2f' | ./lednic.exe",r);
      cout << tmp << endl;
      gSystem->Exec(tmp);
      sprintf(tmp,"mv out.dat led_ZR=%0.2f",r);
      cout << tmp << endl;
      gSystem->Exec(tmp);
    }
}
