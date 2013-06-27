#include "plot.h"
#include "TFile.h"
#include "string.h"

void calculate_b_efficiency()
{
  TH2D * h2_btag_num_b = 0;
  plot2("btag_num_b");
  for (int i = 0; i < gMaxProcess; i++) {
    if (strstr(gProcess[i].fname, "ttjets") != 0) {
      h2_btag_num_b = new TH2D(*gHisto2[i]);
    }
  }
  TH2D * h2_btag_denom_b = 0;
  plot2("btag_denom_b");
  for (int i = 0; i < gMaxProcess; i++) {
    if (strstr(gProcess[i].fname, "ttjets") != 0) {
      h2_btag_denom_b = new TH2D(*gHisto2[i]);
    }
  }

  TH2D * h2_b_eff = new TH2D(*h2_btag_num_b);
  h2_b_eff->Divide(h2_btag_num_b, h2_btag_denom_b, 1., 1., "B");
  new TCanvas;
  //h2_b_eff->GetXaxis()->SetRangeUser(0,500);
  //h2_b_eff->GetYaxis()->SetRangeUser(-2.4,2.4);
  h2_b_eff->Draw();

  TFile * f = new TFile("btag_eff_map.root", "UPDATE");
  h2_b_eff->Smooth();
  h2_b_eff->Write();
  f->Close();
  INFO("Written efficiency map for b-quarks to btag_eff_map.root");

  
  TH1D * h1_btag_num_b = h2_btag_num_b->ProjectionX();
  TH1D * h1_btag_denom_b = h2_btag_denom_b->ProjectionX();

  TH1D * h1_b_eff = new TH1D(*h1_btag_num_b);
  h1_b_eff->Divide(h1_btag_num_b, h1_btag_denom_b, 1., 1., "B");
  
  new TCanvas;
  h1_b_eff->Draw();
}



void calculate_c_efficiency()
{
  TH2D * h2_btag_num_c = 0;
  plot2("btag_num_c");
  for (int i = 0; i < gMaxProcess; i++) {
    if (strstr(gProcess[i].fname, "ttjets") != 0) {
      h2_btag_num_c = new TH2D(*gHisto2[i]);
    }
  }
  TH2D * h2_btag_denom_c = 0;
  plot2("btag_denom_c");
  for (int i = 0; i < gMaxProcess; i++) {
    if (strstr(gProcess[i].fname, "ttjets") != 0) {
      h2_btag_denom_c = new TH2D(*gHisto2[i]);
    }
  }

  TH2D * h2_c_eff = new TH2D(*h2_btag_num_c);
  h2_c_eff->Divide(h2_btag_num_c, h2_btag_denom_c, 1., 1., "B");
  new TCanvas;
  //h2_c_eff->GetXaxis()->SetRangeUser(0,500);
  //h2_c_eff->GetYaxis()->SetRangeUser(-2.4,2.4);
  h2_c_eff->DrawCopy();

  TH1D * h1_btag_num_c = h2_btag_num_c->ProjectionX();
  TH1D * h1_btag_denom_c = h2_btag_denom_c->ProjectionX();

  TH1D * h1_c_eff = new TH1D(*h1_btag_num_c);
  h1_c_eff->Divide(h1_btag_num_c, h1_btag_denom_c, 1., 1., "B");
  new TCanvas;
  h1_c_eff->Draw();

  TFile * f = new TFile("btag_eff_map.root", "UPDATE");
  h2_c_eff->Smooth();
  h2_c_eff->Write();
  f->Close();
  INFO("Written efficiency map for c-quarks to btag_eff_map.root");
}


void calculate_b_mistag_efficiency()
{
  TH2D * h2_btag_num_l = 0;
  plot2("btag_num_l");
  for (int i = 0; i < gMaxProcess; i++) {
    if (strstr(gProcess[i].fname, "ttjets") != 0) {
      h2_btag_num_l = new TH2D(*gHisto2[i]);
    }
  }
  TH2D * h2_btag_denom_l = 0;
  plot2("btag_denom_l");
  for (int i = 0; i < gMaxProcess; i++) {
    if (strstr(gProcess[i].fname, "ttjets") != 0) {
      h2_btag_denom_l = new TH2D(*gHisto2[i]);
    }
  }

  TH2D * h2_l_eff = new TH2D(*h2_btag_num_l);
  h2_l_eff->Divide(h2_btag_num_l, h2_btag_denom_l, 1., 1., "B");
  new TCanvas;
  //h2_l_eff->GetXaxis()->SetRangeUser(0,500);
  //h2_l_eff->GetYaxis()->SetRangeUser(-2.4,2.4);
  h2_l_eff->Draw();

  TH1D * h1_btag_num_l = h2_btag_num_l->ProjectionX();
  TH1D * h1_btag_denom_l = h2_btag_denom_l->ProjectionX();

  TH1D * h1_l_eff = new TH1D(*h1_btag_num_l);
  h1_l_eff->Divide(h1_btag_num_l, h1_btag_denom_l, 1., 1., "B");
  new TCanvas;
  h1_l_eff->Draw();

  TFile * f = new TFile("btag_eff_map.root", "UPDATE");
  h2_l_eff->Smooth();
  h2_l_eff->Write();
  f->Close();
  INFO("Written efficiency map for c-quarks to btag_eff_map.root");
}