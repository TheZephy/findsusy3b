#include <iostream>
#include <sstream>

using namespace std;

#include "plot.h"
#include "fakerate.h"
#include "TColor.h"
#include "TLatex.h"
#include "TFile.h"
#include "TF1.h"

void setNiceColorPalette () {
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;


    Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 }; // 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 1.00, 1.00, 0.50, 0.00 }; // 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 1.00, 1.00, 0.55, 0.00, 0.00 }; // 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.00, 0.00, 0.00, 0.00 }; // 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void default_selection() {
  selection("fullrun75");
}

void vtx_n() {
  default_selection();
  MakeRatioCanvas();
  plot("bSkim_vtx_n");
  logy();
  zoom(0,47);
  min(10);
  max(1e8);
  drawoverflow();
  legend(3e-6, 1, -1, 2);
  lumi();
  ratio();
}

void vtx_n_nopu() {
  selection("fullrun76_nopu");
  MakeRatioCanvas();
  plot("bSkim_vtx_n");
  logy();
  zoom(0,47);
  min(10);
  max(1e8);
  drawoverflow();
  legend(3e-6, 1, -1, 2);
  lumi();
  ratio();
}

void jer_smear() {
  selection("fullrun75_jercalc");
  MakeCanvas(1,1);
  gMoveOverflow = kFALSE;
  plot("JER_deltae");
  legend(1e-2);
  lumi();
  gMoveOverflow = kTRUE;
}

void jer_smear_2d() {
  selection("fullrun75_jercalc");
  // default_selection();
  MakeCanvas(1,1);
  cd(1);
  plot3("JER_deltaepteta");

  TH3D * mc_sum = 0;
  for (Int_t i = gMaxSignal; i < gMaxProcess-1; i++) {
    Int_t process = gOrder[gPadNr][i];

    if (i == gMaxSignal) {
      mc_sum = gHisto3[0][process];
      continue;
    }
    mc_sum->Add(gHisto3[0][process], 1.);
  }
  
  TH2D * pt_projection = (TH2D *) mc_sum->Project3D("yxNUFNOF");
  TH1D * pt_slice = 0;
  gStyle->SetOptStat(1111);
  for (int i = 1; i < pt_projection->GetYaxis()->GetNbins()-1; i=i+10) {
    cout << i << endl;
    // if (i == 1) {
      pt_slice = (TH1D *) pt_projection->ProjectionX("_py", i, i+10, "e");

      new TCanvas(Form("c%d", i), Form("c%d", i), 800, 600);
      pt_slice->Draw();

      TF1 * func = new TF1("fit", "gaus", -20, 20);      
      pt_slice->Fit(func, "R");
      gPad->Print(Form("proj/c%d.pdf", i));
    // }
  }  

  TH2D * eta_projection = (TH2D *) mc_sum->Project3D("zxNUFNOF");
  TH1D * eta_slice = 0;
  for (int i = 1; i < eta_projection->GetYaxis()->GetNbins()-1; i=i+2) {
    cout << i << endl;
    // if (i == 1) {
      eta_slice = (TH1D *) eta_projection->ProjectionX("_pz", i, i+1, "e");

      new TCanvas(Form("d%d", i), Form("d%d", i), 800, 600);
      eta_slice->Draw();

      TF1 * func = new TF1("fit", "gaus", -20, 20);      
      eta_slice->Fit(func, "R");
      gPad->Print(Form("proj/d%d.pdf", i));
    // }
  }  

  
  // setopt(gStyle);
  // setNiceColorPalette();
 
  // gPad->SetRightMargin(0.155);
  // setopt(pt_projection);
  // pt_projection->GetXaxis()->SetTitle("p_{T, gen} - p_{T}");
  // pt_projection->GetYaxis()->SetTitle("jet transverse momentum [GeV]");
  // pt_projection->GetZaxis()->SetTitle("Number of events");
  // pt_projection->Draw("COLZ");
  // pt_projection->SetFillColor(kOrange+1);
}


void pfmet() {
  default_selection();
  MakeRatioCanvas();
  plot("pfmet");
  logy();
  rebin(4);
  zoom(0,330);
  min(0.8);
  legend(1e-6);
  lumi();
  ratio();
}

void pfmet_nojer() {
  selection("fullrun75_jercalc");
  MakeRatioCanvas();
  plot("pfmet");
  logy();
  rebin(4);
  zoom(0,330);
  min(0.8);
  legend(1e-6);
  lumi();
  ratio();
}

void apfmet() {
  default_selection();
  MakeRatioCanvas();
  plot("pfmet");
  logy();
  rebin(4);
  zoom(0,330);
  min(0.8);
  arrow(50., 0, "l");
  legend(1e-6);
  lumi();
  ratio();
}

void jet_n() {
  selection("fullrun76");
  MakeCanvas(1,1);
  plot("check_njets");
  logy();
  min(5e1);
  max(1e8);
  legend(1e-6, 1, -1, 2);
  lumi();
  arrow(1.5, 0, "r");
}

void jet_pt() {
  default_selection();
  MakeCanvas(1,1);
  plot("Jet_pt0");
  rebin(8);
  logy();
  zoom(20, 600);
  min(8e-2);
  max(1e7);
  legend(1e-7, 1, -1, 2);
  lumi();
  arrow(30., 0, "r");
}

void rel_iso() {
  default_selection();
  MakeCanvas(1,1);
  plot("nMuon_relIso");
  zoom(0,1.5);
  logy();
  min(0.9);
  max(1e11);
  legend(1e-6, 1, -1, 2);
  lumi();
  arrow(0.12, 0, "l");
}

void dz_mumu() {
  default_selection();
  MakeCanvas(1,1);
  plot("mumudz");
  logy();
  rebin(50);
  zoom(0,0.2);
  max(1e8);
  legend(1e-6, 1, -1, 2);
  lumi();
  arrow(0.08, 0, "l", 0.6);
}

void m_mumu_zpeak() {
  default_selection();
  MakeRatioCanvas();
  cd(1); plot("m_mumu_zpeak");
  logy();
  zoom(20,200);
  min(8);
  legend(1e-3);
  lumi();
  ratio();
}

void m_mumu_zpeak_noscale() {
  setup("plot_noscale.cfg");
  default_selection();
  MakeRatioCanvas();
  cd(1); plot("m_mumu_zpeak");
  logy();
  zoom(20,200);
  min(8);
  legend(1e-3);
  lumi();
  ratio();
  // setup("plot.cfg");
}

void m_mumu_zpeak_zoom() {
  default_selection();
  MakeCanvas(1,1);
  gMoveOverflow=kFALSE;
  cd(1); plot("m_mumu_zpeak");
  logy();
  zoom(0, 20);
  min(0.8e-1);
  max(1e8);
  legend(1e-3, 1, -1, 2);
  lumi();
  arrow(15., 0 , "r");
  gMoveOverflow=kTRUE;
}

void btag_discriminator() {
  default_selection();
  MakeCanvas(1,1);
  cd(1); plot("1st_btag");
  logy();
  min(1e-1);
  max(1e7);
  legend(1e-6, 1, -1, 2);
  lumi();
  arrow(0.679, 0, "r");
}

// calculate_b_efficiency() 

void crbt() {
  default_selection();
  MakeRatioCanvas();
  cd(1); plot("CR2_m_mumu");
  rebin(2);
  logy();
  zoom(20, 300);
  min(0.8);
  max(1e5);
  legend(1e-3, 1, -1, 2);
  lumi();
  ratio();
}

void crbv() {
  default_selection();
  MakeRatioCanvas();
  cd(1); plot("CR3_m_mumu");
  rebin(2);
  logy();
  zoom(20, 300);
  min(0.8);
  max(1e5);
  legend(1e-6, 1, -1, 2);
  lumi();
  ratio();
}

void crbvc_m_mumu() {
  default_selection();
  MakeCanvas(1,1);
  plot("CR6_m_mumu");
  logy();
  rebin(20);
  max(1e4);
  legend(1e-6, 1, -1, 2);
  lumi();
}

void crbvc_m_smuon() {
  setup("plot_def_nofakes.cfg");
  default_selection();
  MakeCanvas(1,1);
  plot("CR6_m_smuon");
  logy();
  rebin(20);
  max(1e4);
  legend(1e-6, 1, -1, 2);
  lumi();
}

void crbvoc_m_smuon() {
  setup("plot_def_nofakes.cfg");
  default_selection();
  MakeCanvas(1,1);
  plot("CR5_m_smuon");
  logy();
  rebin(20);
  max(1e4);
  legend(1e-6, 1, -1, 2);
  lumi();
}

void crbvc_m_gaugino() {
  default_selection();
  MakeCanvas(1,1);
  plot("CR6_m_gaugino");
  logy();
  rebin(20);
  max(1e4);
  legend(1e-6, 1, -1, 2);
  lumi();
}

void ntl_plots(int i) {
  default_selection();
  if (i == 1) {
    MakeCanvas(1,1);
    plot("nTL_nloose");
    logy();
    min(0.1);
    max(1e9);
    arrow(0.5, 0, "r");
    arrow(1.5, 0, "l");
    legend(1e-6);
    lumi();
  } else if (i == 2) {
    MakeCanvas(1,1);
    plot("nTL_njets");
    logy();
    min(2);
    max(1e7);
    arrow(1.5, 0, "r");
    legend(1e-6);
    lumi();
  } else if (i == 3) {
    MakeCanvas(1,1);
    plot("nTL_mupt");
    rebin(2);
    logy();
    zoom(0, 100);
    min(11);
    max(1e6);
    arrow(20., 0, "r");
    arrow(10., 0, "r");
    legend(1e-4);
    lumi();
  } else if (i == 4) {
    MakeCanvas(1,1);
    plot("nTL_jetpt");
    logy();
    zoom(20, 600);
    min(5e-1);
    max(1e6);
    arrow(50., 0, "r");
    legend(1e-4);
    lumi();
  } else if (i == 5) {
    MakeCanvas(1,1);
    plot("nTL_zmass");
    rebin(2);
    logy();
    zoom(0, 150);
    min(11);
    max(5e5);
    arrow(10., 0, "r");
    arrow(80., 0, "l");
    legend(1e-4);
    lumi();
  } else if (i == 6) {
    MakeCanvas(1,1);
    plot("nTL_met");
    logy();
    zoom(0, 300);
    min(1);
    max(5e5);
    arrow(50., 0, "l");
    legend(1.5e-5);
    lumi();
  } else if (i == 7) {
    MakeCanvas(1,1);
    plot("nTL_mt");
    rebin(2);
    zoom(0, 150);
    logy();
    min(1);
    max(1e6);
    arrow(40., 0, "l", 0.55);
    legend(1e-4);
    lumi();
  } else if (i == 8) {
    MakeCanvas(1,1);
    plot("nTL_jetdphi");
    logy();
    rebin(5);
    min(11);
    max(5e5);
    arrow(1., 0, "r", 0.4);
    legend(1e-6, 1, -1, 2);
    lumi();
  }
}

void cr6_singlefake() {
  MakeCanvas(1,1);
  selection("fullrun75_singlefake");
  get_fakes_1d("CR6_m_smuon");
  rebin(20);
  logy();
  max(1e3);
  legend(1e-6);
  lumi();
}

void cr6_doublefake() {
  MakeCanvas(1,1);
  selection("fullrun75_doublefake");
  get_fakes_1d("CR6_m_smuon");
  rebin(20);
  logy();
  max(1e3);
  legend(1e-6);
  lumi();
}

void cr6_m_smuon() {
  MakeCanvas(1,1);
  setup("plot_fakes.cfg");
  default_selection();
  plot("CR6_m_smuon");
  rebin(20);
  logy();
  min(8e-3);
  max(1e3);
  legend(1e-7, 1, -1, 2);
  lumi();
}

void cr6_m_gaugino() {
  MakeCanvas(1,1);
  setup("plot_fakes.cfg");
  default_selection();
  plot("CR6_m_gaugino");
  rebin(20);
  logy();
  min(8e-3);
  max(1e3);
  legend(1e-7, 1, -1, 2);
  lumi();
}

void m_smuon () {
  setup("plot_fakes.cfg");
  MakeCanvas(1,1);
  default_selection();
  logy();
  plot("m_smuon");
  rebin(80);
  max(1e3);
  legend(1e-7, 1, -1 ,2);
  lumi();
}

void m_gaugino () {
  setup("plot_fakes.cfg");
  MakeCanvas(1,1);
  default_selection();
  logy();
  plot("m_gaugino");
  rebin(20);
  max(1e3);
  min(8e-3);
  legend(1e-7, 1, -1 ,2);
  lumi();
}

void trigger_efficiency (const char * data_file = "jet2012all_Mu17TkMu8.root", const char * mc_file = "ttbar_Mu17TkMu8.root") {
  TH1D * data = (TH1D *) get_object(Form("../config/%s", data_file), "h1_trig_eff_s");
  TH1D * mc = (TH1D *) get_object(Form("../config/%s", mc_file), "h1_trig_eff_s");

  MakeCanvas(1,1);
  setopt(data);
  setopt(mc);

  data->SetLineColor(4);
  data->SetMarkerColor(4);
  data->SetMarkerStyle(20);
  data->GetXaxis()->SetTitle("transverse momentum of muons [GeV]");
  data->GetXaxis()->CenterTitle(true);
  data->GetXaxis()->SetRangeUser(0, 1);
  data->GetYaxis()->SetTitle("Trigger Efficiency");

  mc->SetLineColor(2);
  mc->SetMarkerColor(2);
  mc->SetMarkerStyle(23);

  data->Draw("P");
  mc->Draw("Psame");

  TLegend * t = new TLegend(0.72, 0.74, 0.92, 0.92);
  setopt(t);
  t->AddEntry(data, "Jet Data", "lp");
  t->AddEntry(mc, "t#bar{t} Monte Carlo", "lp");
  t->Draw();
  lumi();
  drawperiod();
}


void m_smu_chi() {
  setup("plot_fakes.cfg");
  default_selection();
  plot2("m_smu_chi");

  TH2D * mc_sum = 0;
  TH2D * data = gHisto2[gMaxProcess-1];
  
  for (Int_t i = gMaxSignal; i < gMaxProcess-1; i++) {
    Int_t process = gOrder[gPadNr][i];

    if (i == gMaxSignal) {
      mc_sum = gHisto2[process];
      continue;
    }
    mc_sum->Add(gHisto2[process]);
  }
  
  MakeCanvas(1,1);
  setopt(gStyle);
  setNiceColorPalette();
  gPad->SetRightMargin(0.155);
  setopt(mc_sum);
  mc_sum->GetXaxis()->SetTitle("m_{#tilde{#mu}} (#mu_{1}^{#pm}, #mu_{2}^{#pm}, jets) [GeV]");
  mc_sum->GetYaxis()->SetTitle("m_{#tilde{#chi}} (#mu_{2}^{#pm}, jets) [GeV]");
  mc_sum->GetZaxis()->SetTitle("Number of events");
  mc_sum->Draw("COLZ");
  mc_sum->SetFillColor(kOrange+1);

  data->SetMarkerSize(1.3);
  data->SetMarkerStyle(34);
  data->SetMarkerColor(kWhite);
  data->DrawCopy("same");
  data->SetMarkerStyle(28);
  data->SetMarkerColor(kBlack);
  data->Draw("same");

  mc_sum->GetXaxis()->SetRangeUser(0, 1400);
  mc_sum->GetYaxis()->SetRangeUser(0, 1000);
  TLine * y300 = new TLine(mc_sum->GetXaxis()->GetBinLowEdge(mc_sum->GetXaxis()->GetFirst()), 300.,
			   mc_sum->GetXaxis()->GetBinUpEdge(mc_sum->GetXaxis()->GetLast()), 300.);
  TLine * x300 = new TLine(300., mc_sum->GetYaxis()->GetBinLowEdge(mc_sum->GetYaxis()->GetFirst()),
			   300., 700.);
  TLine * y700 = new TLine(300., 700.,
			   mc_sum->GetXaxis()->GetBinUpEdge(mc_sum->GetXaxis()->GetLast()), 700.);
  TLine * x700 = new TLine(700., mc_sum->GetYaxis()->GetBinLowEdge(mc_sum->GetYaxis()->GetFirst()),
			   700., mc_sum->GetYaxis()->GetBinUpEdge(mc_sum->GetYaxis()->GetLast()));
  
  y300->Draw();
  x300->Draw();
  y700->Draw();
  x700->Draw();

  TLatex * t2 = new TLatex(0.190, 0.323, "SR 1");
  t2->SetTextSize(0.06);
  t2->SetNDC();
  t2->Draw();

  TLatex * t3 = new TLatex(0.367, 0.154, "SR 2");
  t3->SetTextSize(0.06);
  t3->SetNDC();
  t3->Draw();

  TLatex * t4 = new TLatex(0.367, 0.625, "SR 3");
  t4->SetTextSize(0.06);
  t4->SetNDC();
  t4->Draw();

  TLatex * t5 = new TLatex(0.668, 0.220, "SR 4");
  t5->SetTextSize(0.06);
  t5->SetNDC();
  t5->Draw();

  TLatex * t6 = new TLatex(0.668, 0.488, "SR 5");
  t6->SetTextSize(0.06);
  t6->SetNDC();
  t6->Draw();

  TLatex * t7 = new TLatex(0.531, 0.867, "SR 6");
  t7->SetTextSize(0.06);
  t7->SetNDC();
  t7->Draw();

  TLegend * t = new TLegend(0.168, 0.671, 0.368, 0.851);
  setopt(t);
  t->AddEntry(data, "Data", "p");
  t->AddEntry(mc_sum, "Background", "f");
  t->Draw();


  lumi();
  drawperiod();
}

void signal_smu_chi(int i = 0) {
  setup("plot_fakes.cfg");
  default_selection();
  plot2("m_smu_chi");

  TH2D * signal = 0;
  
  for (Int_t j = 0; j < gMaxSignal; j++) {
    Int_t process = j; //gOrder[gPadNr][j];

    INFO("name " << gProcess[process].fname);
    if (i == 0 && !strcmp(gProcess[process].fname, "signal_m0_1000_m12_1200")) {
      signal = gHisto2[process];
      cout << "hi" << endl;
    }
    if (i == 1 && !strcmp(gProcess[process].fname, "signal_m0_1000_m12_200"))
      signal = gHisto2[process];
    
  }

  if (signal == 0)
    return;
  
  MakeCanvas(1,1);
  setopt(gStyle);
  setNiceColorPalette();
  gPad->SetRightMargin(0.155);
  setopt(signal);
  signal->GetXaxis()->SetTitle("m_{#tilde{#mu}} (#mu_{1}^{#pm}, #mu_{2}^{#pm}, jets) [GeV]");
  signal->GetYaxis()->SetTitle("m_{#tilde{#chi}} (#mu_{2}^{#pm}, jets) [GeV]");
  signal->GetZaxis()->SetTitle("Number of events");
  signal->Draw("COLZ");
  signal->SetFillColor(kOrange+1);

  signal->GetXaxis()->SetRangeUser(0, 1400);
  signal->GetYaxis()->SetRangeUser(0, 1000);
  TLine * y300 = new TLine(signal->GetXaxis()->GetBinLowEdge(signal->GetXaxis()->GetFirst()), 300.,
			   signal->GetXaxis()->GetBinUpEdge(signal->GetXaxis()->GetLast()), 300.);
  TLine * x300 = new TLine(300., signal->GetYaxis()->GetBinLowEdge(signal->GetYaxis()->GetFirst()),
			   300., 700.);
  TLine * y700 = new TLine(300., 700.,
			   signal->GetXaxis()->GetBinUpEdge(signal->GetXaxis()->GetLast()), 700.);
  TLine * x700 = new TLine(700., signal->GetYaxis()->GetBinLowEdge(signal->GetYaxis()->GetFirst()),
			   700., signal->GetYaxis()->GetBinUpEdge(signal->GetYaxis()->GetLast()));
  
  y300->Draw();
  x300->Draw();
  y700->Draw();
  x700->Draw();

  TLatex * t2 = new TLatex(0.190, 0.323, "SR 1");
  t2->SetTextSize(0.06);
  t2->SetNDC();
  t2->Draw();

  TLatex * t3 = new TLatex(0.367, 0.154, "SR 2");
  t3->SetTextSize(0.06);
  t3->SetNDC();
  t3->Draw();

  TLatex * t4 = new TLatex(0.367, 0.625, "SR 3");
  t4->SetTextSize(0.06);
  t4->SetNDC();
  t4->Draw();

  TLatex * t5 = new TLatex(0.668, 0.220, "SR 4");
  t5->SetTextSize(0.06);
  t5->SetNDC();
  t5->Draw();

  TLatex * t6 = new TLatex(0.668, 0.488, "SR 5");
  t6->SetTextSize(0.06);
  t6->SetNDC();
  t6->Draw();

  TLatex * t7 = new TLatex(0.531, 0.867, "SR 6");
  t7->SetTextSize(0.06);
  t7->SetNDC();
  t7->Draw();

  TLegend * t = new TLegend(0.182, 0.625, 0.382, 0.805);
  setopt(t);
  t->SetTextSize(0.0338947);
  if (i == 0)
    t->AddEntry(signal, "m_{0} = 1000, m_{1/2} = 1200", "f");
  if (i == 1)
    t->AddEntry(signal, "m_{0} = 1000, m_{1/2} = 200", "f");
  t->Draw();


  lumi();
  drawperiod();
}

void matched_triggers() {
  default_selection();
  MakeCanvas(1,1);
  setup("plot_data.cfg");
  plot("matched");
  logy();
  max(1e8);
  min(0.5);
  legend(1e-4);
  lumi();
}

void closure_test() {
  setup("plot_ttjets.cfg");
  TH1D * data_fakes = closure_estimate_1d("fullrun75", "m_smuon");  
  TH1D * qcd_fakes = closure_estimate_1d("fullrun75_qcdclosure", "m_smuon");

  default_selection();
  MakeCanvas(1,1);
  plot("m_smuon");
  rebin(80);
  logy();
  max(1e2);

  data_fakes->Rebin(80);
  data_fakes->SetLineColor(kRed);
  data_fakes->SetLineWidth(2);
  data_fakes->Draw("HISTsame");
  qcd_fakes->Rebin(80);
  qcd_fakes->SetLineColor(kBlue);
  qcd_fakes->SetLineWidth(2);
  qcd_fakes->Draw("HISTsame");

  TLegend * t = new TLegend(0.61, 0.62, 0.88, 0.90);
  setopt(t);
  t->AddEntry(gStack[0][0], "t#bar{t} (MC)", "f");
  t->AddEntry(data_fakes, "t#bar{t} (Data T/L)", "l");
  t->AddEntry(qcd_fakes, "t#bar{t} (QCD T/L)", "l");
  t->Draw();  
  
  lumi();
  drawperiod();
}

void closure_proof() {
  setup("plot_ttjets.cfg");
  TH1D * data_fakes = closure_estimate_1d("fullrun75", "CR4_m_smuon");  
  TH1D * qcd_fakes = closure_estimate_1d("fullrun75_qcdclosure", "CR4_m_smuon");

  default_selection();
  MakeCanvas(1,1);
  plot("CR4_m_smuon");
  rebin(20);
  logy();
  max(1e3);

  data_fakes->Rebin(20);
  data_fakes->SetLineColor(kRed);
  data_fakes->SetLineWidth(2);
  data_fakes->Draw("HISTsame");
  qcd_fakes->Rebin(20);
  qcd_fakes->SetLineColor(kBlue);
  qcd_fakes->SetLineWidth(2);
  qcd_fakes->Draw("HISTsame");
  
  TLegend * t = new TLegend(0.61, 0.62, 0.88, 0.90);
  setopt(t);
  t->AddEntry(gStack[0][0], "t#bar{t} (MC)", "f");
  t->AddEntry(data_fakes, "t#bar{t} (Data T/L)", "l");
  t->AddEntry(qcd_fakes, "t#bar{t} (QCD T/L)", "l");
  t->Draw();

  lumi();
  drawperiod();
}

void signal_plots(int i) {
  setup("plot_sig.cfg");
  default_selection();
  MakeCanvas(1,1);
  if (i == 0) {
    plot("bSkim_muo_n");
    logy();
    min(1e-2);
    max(1e5);
    legend(1e-5);
    lumi();
  } else if (i == 1) {
    plot("Sig_ptMu1");
    rebin(5);
    logy();
    min(5e-2);
    max(3e3);
    legend(1e-5);
    lumi();
  } else if (i == 2) {
    plot("Sig_ptMu2");
    rebin(5);
    zoom(0, 500);
    logy();
    min(1e-2);
    max(3e3);
    legend(1e-5);
    lumi();
  } else if (i == 3) {
    plot("Sig_SleptonCharge");
    logy();
    // min(1e-2);
    // max(3e3);
    legend(1e-5);
    lumi();
  }
}

void lambda_ratio() {
  TFile * f = new TFile("../limits/limit_histos.root", "READ");
  TH2D * h2_LambdaRatio = (TH2D*) f->Get("h2_LambdaRatio");

  MakeCanvas(1,1);
  gPad->SetRightMargin(0.10);
  
  setopt(h2_LambdaRatio);
  h2_LambdaRatio->GetYaxis()->SetTitle("m_{1/2} [GeV]");
  h2_LambdaRatio->GetXaxis()->SetTitle("m_{0} [GeV]");
  h2_LambdaRatio->GetZaxis()->SetRangeUser(0.55, 1.45);
  h2_LambdaRatio->Draw("COLZ");
  
  drawperiod();
  lumi();
}

void lambda_1d_limit(int j) {
  TFile * f = new TFile("../limits/limit_histos.root", "READ");
  MakeCanvas(1,1);
  
  TGraph * graphs[6];
  for (int n = 0; n < 6; n++) {
    if (j == 0) {
      graphs[n] = (TGraph *) f->Get(Form("m12_200_%d", n));
    } else if (j == 1) {
      graphs[n] = (TGraph *) f->Get(Form("m12_500_%d", n));
    } else if (j == 2) {
      graphs[n] = (TGraph *) f->Get(Form("m12_800_%d", n));
    } else if (j == 3) {
      graphs[n] = (TGraph *) f->Get(Form("m12_1200_%d", n));
    } else if (j == 4) {
      graphs[n] = (TGraph *) f->Get(Form("m12_1600_%d", n));
    } else if (j == 5) {
      graphs[n] = (TGraph *) f->Get(Form("m12_2000_%d", n));
    } 

  }

  const int n_pnts = graphs[0]->GetN();
  TGraph * updown = new TGraph(2*n_pnts);
  TGraph * updown2 = new TGraph(2*n_pnts);

  double first_pnt = 0;
  for (int i = 0; i < n_pnts; i++) {
    double x[4], y[4];
    graphs[2]->GetPoint(i, x[0], y[0]);
    cout << "i " << i << " " << x[0] << " " << y[0] << endl;
    graphs[3]->GetPoint(n_pnts-i-1, x[1], y[1]);
    graphs[4]->GetPoint(i, x[2], y[2]);
    graphs[5]->GetPoint(n_pnts-i-1, x[3], y[3]);

    if (x[0] == 0)
      first_pnt = i;

    updown->SetPoint(i, x[0], y[0]);
    updown->SetPoint(n_pnts+i, x[1], y[1]);
    updown2->SetPoint(i, x[2], y[2]);
    updown2->SetPoint(n_pnts+i, x[3], y[3]);
  }
  setopt(updown2);
  
  updown2->SetFillColor(kYellow);
  updown2->Draw("AF");
  updown2->GetXaxis()->SetTitle("m_{0} [GeV]");
  updown2->GetYaxis()->SetTitle("95% CL upper limit on #lambda'_{211}");
  updown2->GetYaxis()->SetTitleOffset(1.55);
  cout << first_pnt << endl;
  updown2->GetXaxis()->SetRangeUser(200 + first_pnt * 100, 1500);
  
  updown->SetFillColor(kGreen);
  updown->Draw("f");
  
  graphs[1]->SetLineWidth(2);
  graphs[1]->SetLineStyle(2);
  graphs[1]->Draw("l");

  graphs[0]->SetLineWidth(2);
  graphs[0]->Draw("l");

  TLegend * t = new TLegend(0.192, 0.63, 0.463, 0.91);
  setopt(t);
  t->AddEntry(graphs[0], "Observed", "l");
  t->AddEntry(graphs[1], "Expected", "l");
  t->AddEntry(updown, "#pm 1#sigma", "f");
  t->AddEntry(updown2, "#pm 2#sigma", "f");
  t->Draw();

  drawperiod();
  lumi();
 
  gCanvas->Update(); 
}

void muon_n_1(int i) {
  default_selection();
  MakeCanvas(1,1);
 
  if (i == 0) {
    plot("nMuon_pt");
    rebin(4);
    logy();
    zoom(0, 300);
    min(2e-1);
    max(3e7);
    legend(1e-7, 1, -1, 2);
    lumi();
    max(1e9);
    arrow(20., 0, "r");
  } else if (i == 1) {
    plot("nMuon_eta");
    logy();
    zoom(0, 3);
    min(1);
    max(1e10);
    legend(1e-7, 1, -1 ,2);
    lumi();
    arrow(2.1, 0, "l");
  } else if (i == 2) {
    plot("nMuon_d0Tk");
    rebin(10);
    logy();
    zoom(0, 0.4);
    min(1e-1);
    max(1e9);
    legend(1e-7, 1, -1 ,2);
    lumi();
    arrow(0.2, 0, "l");
  } else if (i == 3) {
    plot("nMuon_TrkChiNormCm");
    logy();
    zoom(0, 20);
    min(3);
    max(1e9);
    legend(1e-7, 1, -1, 2);
    lumi();
    arrow(10., 0, "l", 0.6);
  } else if (i == 4) {
    plot("nMuon_dzTk");
    logy();
    zoom(0, 3);
    max(1e10);
    min(10);
    min(1);
    legend(1e-7, 1, -1, 2);
    arrow(0.5, 0, "l");
    lumi();
  } else if (i == 5) {
    plot("nMuon_TrackerLayersMeasCm");
    logy();
    zoom(1, 14);
    max(1e10);
    min(1);
    legend(1e-7, 1, -1, 2);
    arrow(4.5, 0, "r");
    lumi();
  } else if (i == 6) {
    plot("nMuon_StationsMatched");
    logy();
    zoom(0, 5);
    min(1);
    max(1e12);
    arrow(2, 0, "r");
    legend(1e-7, 1, -1, 2);
    lumi();
  } else if (i == 7) {
    plot("nMuon_ValidPixelHitsCm");
    logy();
    min(2e2);
    max(1e11);
    arrow(0.5, 0, "r");
    legend(1e-7, 1, -1, 2);
    lumi();
  }
}
