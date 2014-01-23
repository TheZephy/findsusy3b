#include "TMath.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "plot.h"
#include "rpv.h"
#include "fakerate.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

struct proc {
  double N;          // number of events
  double staterr;    // MC statistical error
  double systfactor; // systematical uncertainty factor (percent/100) e.g. cross-section
  const char * name; // process name
};

const int nMax = 19;
const int nbins = 6;
const int firststage = 7; // first stage that contains the first bin

proc procs[nbins][nMax];

proc sumprocs[nbins][nMax];

int readprocs(const char * fname)
{
  FILE * infile = fopen(fname, "r");
  if (infile == 0) {
    ERROR("could not open file " << fname);
    return 0;
  }
  char buffer[256];
  int stage = 0;
  int n = 0;
  int i = 0;
  while (fgets(buffer, 256, infile)) {
    if (!strncmp(buffer, " Stage ", 7)) {
      // found a new mass bin, extract stage
      if (sscanf(buffer, " Stage [%d]", & stage) != 1) {
	ERROR("could not get stage: " << buffer);
      }
      INFO("Stage: " << stage);
      // reset counter for processes in this mass bin
      n = 0;
      i = stage-firststage;
      if (i < 0) {
	ERROR("first stage in file must be at least " << firststage);
	return 0;
      }
      else if (i >= nbins) {
	ERROR("last stage must be less than " << firststage + nbins);
	return 0;
      }
    }
    else if (!strncmp(buffer, " +", 2) || !strncmp(buffer, " *", 2)) {
      buffer[1] = '+'; // fix for sscanf below, expects +
      // found a MC or data
      if (stage == 0) {
	ERROR("stage is null!");
	return 0;
      }
      char name[256];
      if (sscanf(buffer, " + %lf +/- %lf %s", 
		 & procs[i][n].N, 
		 & procs[i][n].staterr,
		 name) != 3) {
	ERROR("could not scan line: " << buffer);
	return 0;
      }
      procs[i][n].name = strdup_new(name);
      // assign cross-section error
      if (!strcmp(name, "WW")) {
	procs[i][n].systfactor = 1.5/43.;
      }
      else if (!strcmp(name, "WZ_Q")) {
	procs[i][n].systfactor = 0.7/18.2;
      }
      else if (!strcmp(name, "WZ_Nu")) {
	procs[i][n].systfactor = 0.7/18.2;
      }
      else if (!strcmp(name, "ZZ_nu")) {
	procs[i][n].systfactor = 0.15/5.9;
      }
      else if (!strcmp(name, "ZZ_Q")) {
	procs[i][n].systfactor = 0.15/5.9;
      }
      else if (!strcmp(name, "ZZ_L")) {
	procs[i][n].systfactor = 0.15/5.9;
      }
      else if (!strcmp(name, "fakes")) {
	procs[i][n].systfactor = 0.3;
      }
      else if (!strcmp(name, "data")) {
	procs[i][n].systfactor = 0;
	procs[i][n].staterr = 0; // fix statistical error
      }
      else {
	procs[i][n].systfactor = 0.5;
      }
      INFO(n << " -> " << procs[i][n].name << ": " << procs[i][n].N << " +/- " 
	   << procs[i][n].staterr << " (stat), syst: " << procs[i][n].systfactor);
      n++;
      if (n > nMax) {
	ERROR("Maximum number of processes reached, increase nMax");
      }
    }
    else {
      WARNING("unknown line: " << buffer);
    }
  }
  fclose(infile);
  return n;
}


double fakes[nbins][2];

void fill_fakes(const char * sel = "default29", const char * hname = "btag_jjmm_m")
{
  // six bins with value and error in each bin
  get2dstatistics(fake_estimate_2d(sel, hname), fakes);
}

void summarize_procs(const char * sel = "default29", const char * hname = "btag_jjmm_m")
{
  fill_fakes(sel, hname);
  int maxProc = readprocs("DoublePrompt_MC.txt");
  ofstream out("sum_MC.txt");
  for (int i = 0; i < nbins; i++) {
    out << " Stage [" << i+firststage << "]" << endl;
    double sumN = 0;
    double sumErr2 = 0;
    const char * name = 0;
    for (int j = 0; j < maxProc; j++) {
      name = 0;
      // give a summary at these names
      if (!strcmp(procs[i][j].name, "TTWplus")) {
	name = "VVV";
      }
      else if (!strcmp(procs[i][j].name, "DoublePartonWW")) {
	name = "tt+V";
      }
      else if (!strcmp(procs[i][j].name, "WW")) {
	name = "rare";
      }
      else if (!strcmp(procs[i][j].name, "data")) {
	name = "VV";
      }
      else
	name = 0;
      if (name != 0) {
	out << " + " << sumN << " +/- " << TMath::Sqrt(sumErr2) << " " << name << endl;
	sumN = 0;
	sumErr2 = 0;
      }
      sumN += procs[i][j].N;
      sumErr2 += TMath::Power(procs[i][j].staterr, 2);
    }
    // output fake line
    out << " + " << fakes[i][0] << " +/- " << fakes[i][1] << " fakes" << endl;
    // output data line
    out << " * " << sumN << " +/- " << TMath::Sqrt(sumErr2) << " data" << endl;
    out << "------------------------------------------------------" << endl;
  }
  out.close();
}

void background_systematics()
{
  const int precision = 3;
  cout.precision(precision);
  int maxProc = readprocs("sum_MC.txt");
  // global systematic error, Table 5 of AN
  const double systerrglobal = TMath::Sqrt(1.*1.+2.*2+1.*1.+1.*1.+1.*1.+6.*6.)/100.;
  cout << "global systematic error in %: " << systerrglobal << endl;
  
  double globN = 0;
  double globErrorUncorrelated2 = 0;
  double globErrorCorrelated = 0;
  // walk through table bins and columns and summarize
  for (int bin = 0; bin < nbins; bin++) {
    INFO("bin: " << bin);
    double binN = 0;
    double binErrorUncorrelated2 = 0; 
    double binErrorCorrelated2 = 0; // correlated with other bins
    for (int n = 0; n < maxProc; n++) {
      double locN = 0;
      double locErrorCorrelated2 = 0;
      double locErrorUncorrelated2 = 0;
      if (strcmp(procs[bin][n].name, "data")) {
	locN = procs[bin][n].N;
	// MC cross-section
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*procs[bin][n].systfactor, 2.);
	// MC statistics
	locErrorUncorrelated2 += TMath::Power(procs[bin][n].staterr, 2.);
	// global systematic error, see above
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*systerrglobal, 2.);
      }
      INFO(setw(7) << procs[bin][n].name 
	   << ": " << setw(precision+3) << procs[bin][n].N
	   << " +/- " << setw(precision+3) << TMath::Sqrt(locErrorCorrelated2+locErrorUncorrelated2));

      binN += locN;
      binErrorCorrelated2 += locErrorCorrelated2;
      binErrorUncorrelated2 += locErrorUncorrelated2;
    }
    INFO(setw(7) << "sum" 
	 << ": " << setw(precision+3) << binN 
	 << " +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorUncorrelated2+binErrorCorrelated2)
	 << " [ +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorUncorrelated2) << " (uncorr) "
	 << " +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorCorrelated2) << " (corr) ]");
    globN += binN;
    globErrorUncorrelated2 += binErrorUncorrelated2;
    globErrorCorrelated += TMath::Sqrt(binErrorUncorrelated2);
  }
  INFO(setw(7) << "sum" 
       << ": " << setw(precision+3) << globN 
       << " +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2+globErrorCorrelated*globErrorCorrelated)
       << " [ +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2) << " (uncorr) "
       << " +/- " << setw(precision+3)
       << globErrorCorrelated << " (corr) ]");

  // compute sum over rows with correlated systematics
  globN = 0;
  globErrorUncorrelated2 = 0;
  globErrorCorrelated = 0;
  for (int n = 0; n < maxProc; n++) {
    double procN = 0;
    double procErrorUncorrelated2 = 0; 
    double procErrorCorrelated = 0; // correlated with other bins
    for (int bin = 0; bin < nbins; bin++) {
      double locN = 0;
      double locErrorCorrelated2 = 0;
      double locErrorUncorrelated2 = 0;
      locN += procs[bin][n].N;						       
      if (strcmp(procs[bin][n].name, "data")) {
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*procs[bin][n].systfactor, 2);  
	locErrorUncorrelated2 += TMath::Power(procs[bin][n].staterr, 2);                 
	locErrorCorrelated2 += TMath::Power(procs[bin][n].N*systerrglobal, 2.);
      }
      procN += locN;
      procErrorCorrelated += TMath::Sqrt(locErrorCorrelated2);
      procErrorUncorrelated2 += locErrorUncorrelated2;
    }
    INFO(setw(7) << procs[0][n].name 
	 << ": " << setw(precision+3) << procN 
	 << " +/- " << setw(precision+3)
	 << TMath::Sqrt(procErrorUncorrelated2+procErrorCorrelated*procErrorCorrelated)
	 << " [ +/- " << setw(precision+3)
	 << TMath::Sqrt(procErrorUncorrelated2) << " (uncorr) "
	 << " +/- " << setw(precision+3)
	 << procErrorCorrelated << " (corr) ]");
    if (strcmp(procs[0][n].name, "data")) {
      globN += procN;
      globErrorCorrelated = TMath::Sqrt(procErrorCorrelated*procErrorCorrelated+procErrorUncorrelated2);
      globErrorUncorrelated2 += procErrorUncorrelated2;
    }
  }
  INFO(setw(7) << "sum" 
       << ": " << setw(precision+3) << globN 
       << " +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2+globErrorCorrelated*globErrorCorrelated)
       << " [ +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2) << " (uncorr) "
       << " +/- " << setw(precision+3)
       << globErrorCorrelated << " (corr) ]");
}

void drawPDFUncertainties () {
  MakeCanvas(1,1);
  cd(1);
  stage("0");
  plot("m_smuon");
  rebin(80);

  TH1D * h_U = new TH1D(*gStack[0][gOrder[0][0]]);
  h_U->SetName("h_U");
  h_U->SetTitle("h_U");
  h_U->Reset();

  TH1D * h_L = new TH1D(*gStack[0][gOrder[0][0]]);
  h_L->SetName("h_L");
  h_L->SetTitle("h_L");
  h_L->Reset();

  TH1D * h_M = new TH1D(*gStack[0][gOrder[0][0]]);
  h_M->SetName("h_M");
  h_M->SetTitle("h_M");
  h_M->Reset();

  TH1D * h_EP = new TH1D(*gStack[0][gOrder[0][0]]);
  h_EP->SetName("h_EP");
  h_EP->SetTitle("h_EP");
  h_EP->Reset();

  TH1D * h_EM = new TH1D(*gStack[0][gOrder[0][0]]);
  h_EM->SetName("h_EM");
  h_EM->SetTitle("h_EM");
  h_EM->Reset();
  
  TH1D * h_pdf[6] = { 0 };
  const char * pdf[] = {
    "m_smuon_sys_pdf_CT10_up",
    "m_smuon_sys_pdf_CT10_down",
    "m_smuon_sys_pdf_MSTW_up",
    "m_smuon_sys_pdf_MSTW_down",
    "m_smuon_sys_pdf_NNPDF_up",
    "m_smuon_sys_pdf_NNPDF_down"
  };

  TCanvas * c2 = new TCanvas("c2", "PDF plots", 849, 600);
  TCanvas * c3 = new TCanvas("c3", "PDF ratio plots", 849, 600);
  // TCanvas * c4 = new TCanvas("c4", "PDF ratio plots (nice)", 849, 600);

  for (int i = 0; i < 3; i++) {
    gCanvas->cd(gPadNr);
    plot(pdf[2*i]);
    rebin(80);
    h_pdf[2*i] = new TH1D(*gStack[0][gOrder[0][0]]);
    c2->cd();
    h_pdf[2*i]->SetLineColor(i+2);
    h_pdf[2*i]->SetLineWidth(2);
    h_pdf[2*i]->SetFillStyle(0);
    INFO("title: " << h_pdf[2*i]->GetTitle());
    h_pdf[2*i]->DrawCopy(Form("hist%s", i == 0 ? "" : "same"));

    gCanvas->cd(gPadNr);
    plot(pdf[2*i+1]);
    rebin(80);
    h_pdf[2*i+1] = new TH1D(*gStack[0][gOrder[0][0]]);
    c2->cd();
    h_pdf[2*i+1]->SetLineColor(i+2);
    h_pdf[2*i+1]->SetLineWidth(2);
    h_pdf[2*i+1]->SetFillStyle(0);
    h_pdf[2*i+1]->DrawCopy(Form("hist%s", i == 0 ? "" : "same"));
  }

  for (int bin = 0; bin < h_pdf[0]->GetNbinsX()+2; bin++) {
    Double_t u = TMath::Max(h_pdf[0]->GetBinContent(bin), TMath::Max(h_pdf[2]->GetBinContent(bin), h_pdf[4]->GetBinContent(bin)));
    Double_t l = TMath::Min(h_pdf[1]->GetBinContent(bin), TMath::Min(h_pdf[3]->GetBinContent(bin), h_pdf[5]->GetBinContent(bin)));
    Double_t m = (u+l)/2.;

    h_U->SetBinContent(bin, u);
    h_L->SetBinContent(bin, l);
    h_M->SetBinContent(bin, m);

    if (m != 0) {
      h_EP->SetBinContent(bin, (u-m)/m);
      h_EM->SetBinContent(bin, (l-m)/m);

      for (int j = 0; j < 6; j++) {
	h_pdf[j]->SetBinContent(bin, (h_pdf[j]->GetBinContent(bin)-m)/m);
      }
    }
    else {
      h_EP->SetBinContent(bin, 0.);
      h_EM->SetBinContent(bin, 0.);

      for (int j = 0; j < 6; j++) {
	h_pdf[j]->SetBinContent(bin, 0.);
      }
    }
  }
  
  c2->cd(); // weighted distributions
  h_U->SetLineColor(kOrange);
  h_U->SetLineStyle(kDotted);
  h_U->SetFillStyle(0);
  h_U->Draw("same");
  h_L->SetLineColor(kOrange);
  h_L->SetLineStyle(kDotted);
  h_L->SetFillStyle(0);
  h_L->Draw("same");

  c3->cd(); // percentual deviations
  // coloured areas
  const int maxXBins = h_U->GetNbinsX()+2;
  Double_t x[3][maxXBins];
  Double_t y[3][maxXBins];
  Double_t ex[3][maxXBins];
  Double_t ey[3][maxXBins];

  for (int bin = 0; bin < h_pdf[0]->GetNbinsX()+2; bin++) {
    for (int i = 0; i < 3; i++) {
      x[i][bin] = h_pdf[2*i]->GetBinCenter(bin);
      ex[i][bin] = h_pdf[2*i]->GetBinWidth(bin)/2;
      y[i][bin] = ( h_pdf[2*i]->GetBinContent(bin) + h_pdf[2*i+1]->GetBinContent(bin) )/2;
      ey[i][bin] = (h_pdf[2*i]->GetBinContent(bin) - y[i][bin]);
    }
  }

  // c4->cd();
  TGraphErrors * h_fillpdf[3];
  TMultiGraph * mg = new TMultiGraph();
  for (int i = 0; i < 3; i++) { 
    h_fillpdf[i] = new TGraphErrors(maxXBins, x[i], y[i], ex[i], ey[i]);
    h_fillpdf[i]->SetFillStyle(3004);
    h_fillpdf[i]->SetFillColor(i+2);
    // h_fillpdf[i]->SetRangeUser(0, maxXBins);
    mg->Add(h_fillpdf[i]); 
  }  


  h_EP->SetLineColor(kBlack);
  h_EP->SetLineStyle(kDashed);
  h_EP->SetFillStyle(0);
  h_EP->GetYaxis()->SetTitle("#frac{Variation - Best fit}{Best fit}");
  h_EP->SetMaximum(0.2);
  h_EP->SetMinimum(-0.15);

  h_EP->Draw("hist");
  mg->Draw("2");


  // outlines
  for (int i = 0; i < 6; i++) { 
    h_pdf[i]->Draw("histsame");
    // Form("hist%s", i == 0 ? "" : "same")
  }

  // redraw over the lines
  // h_EP->Draw("same");

  h_EM->SetLineColor(kBlack);
  h_EM->SetLineStyle(kDashed);
  h_EM->SetFillStyle(0);
  h_EM->Draw("same");



  TLegend * t = new TLegend(0.72, 0.74, 0.92, 0.92);
  setopt(t);
  t->AddEntry(h_fillpdf[0], "CT10", "f");
  t->AddEntry(h_fillpdf[1], "MSTW2008", "f");
  t->AddEntry(h_fillpdf[2], "NNPDF 2.3", "f");
  t->Draw();
  
  c3->Modified();
  c3->Update();
}
