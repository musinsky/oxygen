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
  gROOT->SetStyle("Plain");
  // default lednic program use ZT = 0.01 a V = 0.001
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // change only this
  Float_t ZR_min=0.5;
  Float_t ZR_max=2.5;
  Float_t ZR_step=0.1;

  Float_t ZT_min=0.5;
  Float_t ZT_max=6.0;
  // graf ZR od R pre rozne V ( pritom hodnota ZT = ZT_const )
  Float_t ZT_const=1.0;
  Float_t ZT_step=1.0;
  
  Float_t V_min=0.2;
  Float_t V_max=0.8;
  // graf ZR od R pre rozne ZT ( pritom hodnota V = V_const )
  Float_t V_const=0.2;
  Float_t V_step=0.1;
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Text_t tmp[200];
  ifstream fin;
  const UShort_t MAX=100;
  Float_t x_ZR[MAX]={0};
  Float_t x_ZT[MAX]={0};
  Float_t x_V[MAX] ={0};
  Float_t ZR_ZT_V[MAX][MAX][MAX]={0};
  
  // calculate
  UShort_t count_ZR=0;
  for (Float_t zr = ZR_min; zr <= ZR_max; zr += ZR_step)
    {
      UShort_t count_ZT=0;
      x_ZR[count_ZR]=zr;
      for (Float_t zt = ZT_min; zt <= ZT_max; zt += ZT_step)
	{
	  UShort_t count_V=0;
	  x_ZT[count_ZT]=zt;
	  for (Float_t v = V_min; v <= V_max; v += V_step)
	    {
	      x_V[count_V]=v;
	      
	      sprintf(tmp,"echo '%0.5f,%0.5f,%0.5f'|./lednic_input.exe>out",zr,zt,v);
	      gSystem->Exec(tmp);
	      fin.open("out");
	      fin >> ZR_ZT_V[count_ZR][count_ZT][count_V];
	      fin.close();
	      fin.clear();
	      gSystem->Exec("rm out");
	      
	      count_V++;
	    } // v
	  count_ZT++;
	} // zt
      count_ZR++;
    } // zr
  
  // drawing
  // ZR od R pre rozne ZT, ZT =1.0, 2.0, 3.0, ...
  // to vsetko pri V = V_const, hladanie indexa pre V_const
  UShort_t V_index=0;
  Float_t round=0;
  if ( (V_const<=V_max) && (V_const>=V_min) ) // ak je V v pocitanom intervale
    round = ( V_const - V_min ) / V_step;
  // inak index = 0; => pre V = V_min
  V_index = short( floor( (round*100+0.5)/100 ) );

  TGraph *g_ZR_cZT[count_ZT];
  leg_cZT = new TLegend(0.7,0.6,0.89,0.89);
  c1 = new TCanvas();
  c1->cd();
  for (UShort_t i = 0; i < count_ZT; i++)
    {
      g_ZR_cZT[i] = new TGraph;
      for (UShort_t j = 0; j < count_ZR; j++)
	g_ZR_cZT[i]->SetPoint(j,x_ZR[j],ZR_ZT_V[j][i][V_index]);
      g_ZR_cZT[i]->SetMarkerStyle(8);
      g_ZR_cZT[i]->SetMarkerColor(i+1);
      sprintf(tmp,"v_{const} = %0.2f",x_V[V_index]);
      g_ZR_cZT[i]->SetTitle(tmp);
      g_ZR_cZT[i]->GetXaxis()->SetTitle("r_{ 0} [fm]");
      g_ZR_cZT[i]->GetYaxis()->SetTitle("R( k=0.02 [GeV/c] )");
      g_ZR_cZT[i]->GetYaxis()->CenterTitle();
      if (i==0)
	g_ZR_cZT[i]->Draw("AP");
      g_ZR_cZT[i]->Draw("P");
      sprintf(tmp,"t_{ 0} = %0.2f",x_ZT[i]);
      leg_cZT->AddEntry(g_ZR_cZT[i],tmp,"P");
    }
  leg_cZT->Draw();

  // ZR od R pre rzone V, V =0.2, 0.2, 0.4, ...
  // to vsetko pri ZT = ZT_const, hladanie indexa pre ZT_const
  UShort_t ZT_index=0;
  round=0;
  if ( (ZT_const<=ZT_max) && (ZT_const>=ZT_min) ) // ak je ZT v pocitanom intervale
    round = ( ZT_const - ZT_min ) / ZT_step;
  // inak index = 0; => pre ZT_min
  ZT_index = short( floor( (round*100+0.5)/100 ) );
  
  TGraph *g_ZR_cV[count_V];
  leg_cV = new TLegend(0.7,0.6,0.89,0.89);
  c2 = new TCanvas();
  c2->cd();
  for (UShort_t i = 0; i < count_V; i++)
    {
      g_ZR_cV[i] = new TGraph;
      for (UShort_t j = 0; j < count_ZR; j++)
	g_ZR_cV[i]->SetPoint(j,x_ZR[j],ZR_ZT_V[j][ZT_index][i]);
      g_ZR_cV[i]->SetMarkerStyle(8);
      g_ZR_cV[i]->SetMarkerColor(i+1);
      sprintf(tmp,"t_{ 0,const} = %0.2f [fm]",x_ZT[ZT_index]);
      g_ZR_cV[i]->SetTitle(tmp);
      g_ZR_cV[i]->GetXaxis()->SetTitle("r_{ 0} [fm]");
      g_ZR_cV[i]->GetYaxis()->SetTitle("R( k=0.02 [GeV/c] )");
      g_ZR_cV[i]->GetYaxis()->CenterTitle();
      if (i==0)
	g_ZR_cV[i]->Draw("AP");
      g_ZR_cV[i]->Draw("P");
      sprintf(tmp,"v = %0.2f",x_V[i]);
      leg_cV->AddEntry(g_ZR_cV[i],tmp,"P");
    }
  leg_cV->Draw();
}
