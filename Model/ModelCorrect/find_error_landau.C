////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///     last changes: 26 oct 2003    ///
///                                  ///
////////////////////////////////////////

// Landau fit (is very good)

{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // change only this
  Float_t min_p=1.75; // min impuls
  Float_t max_p=4.75; // max impuls
  Float_t max_bin=0.05;
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TFile *fin = TFile::Open("oxygen_DST.root");
  h1 = new TH1F("h1","",max_bin*10000,0,max_bin);
  
  Text_t temp_p[250]="sqrt(px*px+py*py+pz*pz)>%f && sqrt(px*px+py*py+pz*pz)<%f";
  Text_t cut_p[250];
  sprintf(cut_p,temp_p,min_p,max_p);
  TCut cut1=cut_p;
  TCut cut2="ident==11 || ident==1";
  
  oxy_DST->Draw("delp/sqrt(px*px+py*py+pz*pz)>>h1",cut1&&cut2);
  h1->Fit("landau","","",0,max_bin);

}
