////////////////////////////////////////
///                                  ///
///           Jan Musinsky           ///
///                                  ///
///        mucha@sunhe.jinr.ru       ///
///                                  ///
///      last changes: 26 oct 2003   ///
///                                  ///
////////////////////////////////////////

{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gROOT->LoadMacro("histo_correl_DST.C");
  gROOT->LoadMacro("draw_correl_DST.C");
  
  // mass of correlation particle (pre protony MASS = 0.938)
  const Float_t MASS=0.938272;
  // kolko max. fragmentov v jednom evente (max = MAX_NUM_P = 50 )
  const UShort_t MAX_NUM_FRAG=8;
  // kolko krat viac pozadie ako true pairs
  const UInt_t NUM_BACK_EVENTS=1;
  // aky fragment (jednoznacne protony=1)
  const UShort_t FRAGMENT=1;
  // + aky fragment v akom impulsnom intervale (nejednoznacne p=11)
  const UShort_t FRAGMENT_IMP=11;
  // jeho impulsny interval ( pre FRAGMENT_IMP )
  const Float_t MIN_IMP=1.75, MAX_IMP=4.75;
  
  UShort_t bining=30;
  Float_t max_bin=0.3;
  h_true = new TH1F("h_true","mix true protons",bining,0,max_bin);
  h_back = new TH1F("h_back","mix background protons",bining,0,max_bin);

  //////////////////////////////////////
  histo_correl_DST(); // spocitnaie korelacnych histogramov
  //////////////////////////////////////

  // od ktoreho k* zacinat normovat
  Float_t norm=0.25;
  correl_graph = new TGraphErrors;

  //////////////////////////////////////  
  draw_correl_DST(); // draw korelacnej funkcie
  //////////////////////////////////////

  // cesta pre subory s korelacnymi funkciami lednickeho
  Text_t gPath[]="/home/mucha/OXYGEN/LednicCorrel";
  // pre ktore r0 zobarzit lednickeho funkciu + korigovanu na chyby
  Float_t which_lednic=1.7; // 0 for do not draw
  // pre ktore vsetky lednickeho r0 pocitat chi2
  Float_t first_r0=1.0;
  Float_t last_r0=3.0;
  // pocitanie chi2 len do urciteho k
  Float_t MAX_k = 0.06;
  chi2_graph = new TGraph();
  lednic_graph = new TGraph();
  // Vysledne polomer a jeho chyba R0[0]=r0, R0[1]=r0+, R0[2]=r0-
  Float_t R0[3]={0};

  //////////////////////////////////////
  draw_chi2(chi2_graph,lednic_graph,R0,0); // spocitanie chi2
  //////////////////////////////////////
  cout << endl;
  cout <<"R0= "<<R0[0]<<" + "<<R0[1]-R0[0]<<" - "<<R0[0]-R0[2]<<endl;
  
  // ceste pre subor s korekciami lednickeho funkcii na exp. chybu
  Text_t file_c[]="/home/mucha/OXYGEN/Model/correct_on_moscow.asc";
  chi2_graph_c = new TGraph();
  lednic_graph_c = new TGraph();
  // Vysledne polomer a jeho chyba R0[0]=r0, R0[1]=r0+, R0[2]=r0-
  Float_t R0_c[3]={0};
  
  //////////////////////////////////////
  draw_chi2(chi2_graph_c,lednic_graph_c,R0_c,1); // spocitanie chi2 (corrected)
  //////////////////////////////////////
  cout << endl;
  cout << "After correct on experimental errors" << endl;
  cout <<"R0= "<<R0_c[0]<<" + "<<R0_c[1]-R0_c[0]<<" - "<<R0_c[0]-R0_c[2]<<endl;
 
  c1 = new TCanvas();
  c1->Divide(2);
  c1_1->cd();
  correl_graph->Draw("AP");
  c1_2->cd();
  chi2_graph->Draw("AP");
}
