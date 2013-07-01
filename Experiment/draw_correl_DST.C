////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

void draw_correl_DST()
{
  // pocitanie normalizacneho faktor
  UInt_t sum_true=0;
  UInt_t sum_back=0;
  Int_t begin_norm = h_true->FindBin(norm); // od ktoreho binu normujeme
  for (UShort_t b=begin_norm; b <= h_true->GetNbinsX(); b++)
    {
      sum_true += h_true->GetBinContent(b);
      sum_back += h_back->GetBinContent(b);
    }
  Float_t ratio = (float(sum_back)) / (float(sum_true)); // normovaci faktor
  
  // podiel dvoch funkcii
  Float_t fraction=0;
  Float_t error=0;
  Float_t x=0;
  UShort_t count=0;
  os_x = h_true->GetXaxis();
  for (UShort_t b=1; b<= h_true->GetNbinsX(); b++) // zacinat od 1 !!!
    {
      fraction = (h_true->GetBinContent(b)*ratio)/h_back->GetBinContent(b);
      error = fraction / sqrt( h_true->GetBinContent(b) );
      x = os_x->GetBinCenter(b);
      correl_graph->SetPoint(count,x,fraction);
      correl_graph->SetPointError(count,0,error);
      count++;
    }
  correl_graph->SetTitle("p+p correlation function, ^{16}O+p interactions");
  correl_graph->SetMarkerStyle(8);
  correl_graph->SetMinimum(0);
  correl_graph->GetXaxis()->SetTitle("k [GeV/c]");
  correl_graph->GetYaxis()->SetTitle("R");
  correl_graph->GetYaxis()->CenterTitle();
}
////////////////////////////////////////
void draw_chi2(TGraph *chi2_gra,TGraph *led_gra,Float_t *RR,Int_t correct=0)
{
  // calculate chi2 square
  Double_t temp_x=0, temp_y=0, temp_error=0; // musi byt Double_t
  Float_t chi2_square=0;
  Text_t tmp[250];
  ifstream fin;
  UShort_t count_chi2=0;
  UShort_t ndf=0;
  Float_t x=0, y=0;
  
  if (correct)
    {
      const UShort_t max_c=bining;
      Float_t x_c[max_c]={0};
      Float_t y_c[max_c]={0};
      ifstream fin_c;
      fin_c.open(file_c);
      UShort_t count_c=0;
      // nacitanie korekcnych faktorov
      while (fin_c.good())
	fin_c >> x_c[count_c] >> y_c[count_c++];
      fin_c.close();
      fin_c.clear();
    }
  
  for (Float_t r0=first_r0; r0<=last_r0; r0+=0.1)
    {
      sprintf(tmp,"%s/led_ZR=%1.2f",gPath,r0);
      fin.open(tmp);
      chi2_square=0;
      ndf=0;
      while(fin.good())
	{
	  fin >> x >> y;

	  // oprava na correct factor
	  if (correct)
	    for (UShort_t jj=0; jj < correl_graph->GetN(); jj++)
	      if ( int(floor(x*1000+0.5)) == int(floor(x_c[jj]*1000+0.5)) )
		y = y * y_c[jj];
	  
	  for (UShort_t j=0; j < correl_graph->GetN(); j++)
	    {
	      correl_graph->GetPoint(j,temp_x,temp_y);
	      temp_error = correl_graph->GetErrorY(j);
	      if (j==0) continue; // !!! without first point
	      if (temp_x >= MAX_k) continue; // pocitanie chi2 len do MAX_k
	      
	      // len pre tie iste x v Lednic aj v Graph (zaokruhlenie)
	      if ( int(floor(x*1000+0.5)) == int(floor(temp_x*1000+0.5)) )
		
		{
		  chi2_square += ( (temp_y-y)/temp_error )**2;
		  ndf++;
		}
	    }
	}      
      
      fin.close();
      fin.clear();
      chi2_gra->SetPoint(count_chi2++,r0,chi2_square/float(ndf));
    }
  sprintf(tmp,"#chi^{ 2} distribution, ndf = %i",ndf);
  chi2_gra->SetTitle(tmp);
  chi2_gra->SetMarkerStyle(8);
  chi2_gra->SetMinimum(0);
  chi2_gra->GetXaxis()->SetTitle("r_{ 0} [fm]");
  chi2_gra->GetYaxis()->SetTitle("#chi^{ 2} / ndf");
  chi2_gra->GetYaxis()->CenterTitle();
  // fiting and find R0
  chi2_gra->Fit("pol6","Q0");
  pol6->SetLineStyle(2);
  pol6->SetLineWidth(1);
  RR[0] = pol6->GetMinimumX();
  Float_t y_R0 = pol6->GetMinimum()+1.0; // +(1.0/ndf)
  RR[1] = pol6->GetX(y_R0,RR[0],RR[0]+10); // after minimum
  RR[2] = pol6->GetX(y_R0,0,RR[0]); // before minimum
  
  // drawing Lednic correlation function
  sprintf(tmp,"%s/led_ZR=%1.2f",gPath,which_lednic);
  fin.open(tmp);
  UShort_t count=0;
  while(fin.good())
    {
      fin >> x >> y;

      // oprava na correct factor
      if (correct)
	{
	  for (UShort_t jj=0; jj < correl_graph->GetN(); jj++)
	    if ( int(floor(x*1000+0.5)) == int(floor(x_c[jj]*1000+0.5)) )
	      {
		y = y * y_c[jj];
		led_gra->SetPoint(count++,x,y);
	      }
	}
      else led_gra->SetPoint(count++,x,y);  
    }

  led_gra->SetMarkerStyle(4);
  led_gra->SetMinimum(0);
  led_gra->SetMinimum(0);

}
