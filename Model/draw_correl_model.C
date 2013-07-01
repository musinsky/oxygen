////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

void draw_correl_model(const TH1F *his_true, const TH1F *his_back, TGraphErrors *graph)
{
  // pocitanie normalizacneho faktor
  UInt_t sum_true=0;
  UInt_t sum_back=0;
  Int_t begin_norm = his_true->FindBin(norm); // od ktoreho binu normujeme
  for (UShort_t b=begin_norm; b <= his_true->GetNbinsX(); b++)
    {
      sum_true += his_true->GetBinContent(b);
      sum_back += his_back->GetBinContent(b);
    }
  Float_t ratio = (float(sum_back)) / (float(sum_true)); // normovaci faktor

  // podiel dvoch funkcii
  Float_t fraction=0;
  Float_t error=0;
  Float_t x=0;
  UShort_t count=0;
  os_x = his_true->GetXaxis();
  for (UShort_t b=1; b<= his_true->GetNbinsX(); b++) // zacinat od 1 !!!
    {
      fraction = (his_true->GetBinContent(b)*ratio)/his_back->GetBinContent(b);
      error = fraction / sqrt( his_true->GetBinContent(b) );
      x = os_x->GetBinCenter(b);
      graph->SetPoint(count,x,fraction);
      graph->SetPointError(count,0,error);
      count++;
    }
  graph->SetTitle("p+p correlation function MODEL, ^{16}O+p interactions");
  graph->SetMarkerStyle(8);
  graph->SetMinimum(0);
  graph->GetXaxis()->SetTitle("k [GeV/c]");
  graph->GetYaxis()->SetTitle("R");
  graph->GetYaxis()->CenterTitle();
}
////////////////////////////////////////
void find_ratio(const TGraphErrors *g, const TGraphErrors *g_c)
{
  ofstream info_out(out_info);
  Double_t x,y,x_c,y_c,rat; // must be Double_t !
  UShort_t n=g->GetN();
  for (UShort_t j=0; j<n; j++)
    {
      g  ->GetPoint(j,x  ,y);
      g_c->GetPoint(j,x_c,y_c);
      if ( int(floor(x*1000+0.5)) != int(floor(x_c*1000+0.5)) )
	cout << " !!! other bin !!!" << endl;
      rat = y_c / y;
      info_out << x << " " << rat << endl;
    }
}
