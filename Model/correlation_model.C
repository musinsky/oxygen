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
  gROOT->LoadMacro("histo_correl_model.C");
  gROOT->LoadMacro("draw_correl_model.C");
  
  // !!! korekcia sa robi iba pre protony v intervale MIN_IMP, MAX_IMP !!!
  // t.j. pre jednoznacne protony z experimentu nie korekcia na chybu
  
  // Moscow model
  Text_t in_file[]  ="moscow_model.root";
  Text_t in_file_c[]="moscow_model_c.root";
  Text_t out_info[]="correct_on_moscow.asc";
  
  //  // Musulmanbekov model
  //  Text_t in_file[]  ="musul_model.root";
  //  Text_t in_file_c[]="musul_model_c.root";
  //  Text_t out_info[]="correct_on_musul.asc";
  
  // How many events for calculat correlation, 0 for all events
  UInt_t nevents=0;
  
  // mass of correlation particle (pre protony MASS = 0.938)
  const Float_t MASS=0.938272;
  // kolko max. fragmentov v jednom evente (max = MAX_NUM_P = 50 )
  const UShort_t MAX_NUM_FRAG=8;
  // kolko krat viac pozadie ako true pairs
  const UInt_t NUM_BACK_EVENTS=1;
  // jeho impulsny interval ( pre FRAGMENT_IMP )
  const Float_t MIN_IMP=1.75, MAX_IMP=4.75;

  UShort_t bining=30;
  Float_t max_bin=0.3;
  // od ktoreho k* zacinat normovat
  Float_t norm=0.25;
  
  // for model
  h_true   = new TH1F("h_true","true",bining,0,max_bin);
  h_back   = new TH1F("h_back","background",bining,0,max_bin);
  histo_correl_model(in_file,h_true, h_back); // model
  
  correl_graph = new TGraphErrors;
  draw_correl_model(h_true, h_back, correl_graph);
  
  // for correct model
  h_true_c = new TH1F("h_true_c","true_correct",bining,0,max_bin);
  h_back_c = new TH1F("h_back_c","background_correct",bining,0,max_bin);
  histo_correl_model(in_file_c, h_true_c, h_back_c); // correct model

  correl_graph_c = new TGraphErrors;
  draw_correl_model(h_true_c, h_back_c, correl_graph_c);

  // find ratio between model and model correct on errors
  find_ratio(correl_graph,correl_graph_c);
}
