// Author: Jan Musinsky
// 15/10/2010

Double_t CalcLednicR(Double_t zr, Double_t zt, Double_t v)
{
  TString cmd = Form("echo '%0.5f,%0.5f,%0.5f'", zr, zt, v);
  cmd+= " | ./lednic_input.exe";
  Printf("%s", cmd.Data());
  cmd = gSystem->GetFromPipe(cmd.Data());
  if (!cmd.IsFloat()) {
    Printf("Problem with program 'lednic_input.exe'");
    return 0;
  }
  return cmd.Atof();
}

void auto_lednic_input_new()
{
  gROOT->SetStyle("Plain");
  Double_t ZR_min = 0.5;
  Double_t ZR_max = 3.0;
  const Int_t n   = 100;

  // lednic_input default values (ZT = 0.01, V = 0.001)
  Double_t ZT_min = 0.001;
  Double_t ZT_max = 2.0;
  Double_t V_min  = 0.0001;
  Double_t V_max  = 0.3;

  Double_t x[n], y_min[n], y_max[n], step = (ZR_max - ZR_min)/n;
  for (Int_t i = 0; i < n; i++) {
    x[i]     = ZR_min + i*step;
    y_min[i] = CalcLednicR(x[i], ZT_min, V_min);
    y_max[i] = CalcLednicR(x[i], ZT_max, V_max);
  }

  TGraph *g_min = new TGraph(n, x, y_min);
  TGraph *g_max = new TGraph(n, x, y_max);
  TGraph *gg    = new TGraph(2*n);
  for (Int_t i = 0; i < n; i++) {
    gg->SetPoint(i, x[i], y_max[i]);
    gg->SetPoint(n+i, x[n-i-1], y_min[n-i-1]);
  }

  gg->SetTitle(";r_{0} [fm];R (k = 0.02 GeV/c)");
  gg->GetYaxis()->CenterTitle(); gg->SetMinimum(0); gg->SetFillStyle(3003);
  gg->Draw("AF");
  g_min->Draw("C");
  g_max->SetLineStyle(2);
  g_max->Draw("C");
  TLegend *leg = new TLegend(0.50, 0.70, 0.85, 0.85);
  leg->AddEntry(g_min, Form("t_{0} = %.3f fm, v = %.4fc", ZT_min, V_min), "L");
  leg->AddEntry(g_max, Form("t_{0} = %.3f fm, v = %.4fc", ZT_max, V_max), "L");
  leg->SetFillColor(0);
  leg->Draw();
  gPad->Print("pp_correl_r0.pdf");
}
