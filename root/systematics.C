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

const int nbins = 6; // number of bins
vector<vector<proc> > procs;


double fakes[nbins][2];
void fill_fakes(const char * sel = "fullrun76", const char * hname = "jjmm_m")
{
  // six bins with value and error in each bin
  get2dstatistics(fake_estimate_2d(sel, hname), fakes);
}

int readprocs(const char * fname)
{
  FILE * infile = fopen(fname, "r");
  if (infile == 0) {
    ERROR("could not open file " << fname);
    return 0;
  }
  double fakes_syst_factor = TMath::Sqrt(45.*45. + 10.*10. + 8.*8.)/100;
  cout << "systematic factor for fakes: " << fakes_syst_factor << endl;

  char buffer[256];

  while (fgets(buffer, 256, infile)) {
    char name[256];
    vector<proc> current_procs;
    for (int i = 0; i < nbins; i++) {
      proc temp;
      current_procs.push_back(temp);
    }
      

    if (sscanf(buffer, "%s %lf +/- %lf | %lf +/- %lf | %lf +/- %lf | %lf +/- %lf | %lf +/- %lf | %lf +/- %lf |",
	       name,
	       & current_procs[0].N,
	       & current_procs[0].staterr,
	       & current_procs[1].N,
	       & current_procs[1].staterr,
	       & current_procs[2].N,
	       & current_procs[2].staterr,
	       & current_procs[3].N,
	       & current_procs[3].staterr,
	       & current_procs[4].N,
	       & current_procs[4].staterr,
	       & current_procs[5].N,
	       & current_procs[5].staterr) != 13) {
      ERROR("could not scan line: " << buffer);
    }
    for (int i = 0; i < nbins; i++) {
      current_procs[i].name = strdup_new(name);   
      // cout << setw(12) << current_procs[i].N << " +/- " << setw(12) << current_procs[i].staterr << "  ";
    }
    // cout << endl;


    // assign cross section errors
    if (!strcmp(name, "fakes")) {
      fill_fakes("fullrun75", "jjmm_m");
      for (int i = 0; i < nbins; i++) {
	current_procs[i].staterr = fakes[i][1]; // the fakes error 
	current_procs[i].systfactor = fakes_syst_factor;
      }
    }
    else if (!strcmp(name, "ttwjets")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.067/0.232;
      }
    }
    else if (!strcmp(name, "ttzjets")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.024/0.2057;
      }
    }
    else if (!strcmp(name, "wwjetsto2l2nu")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.396/5.885;
      }
    }
    else if (!strcmp(name, "wzjetsto2l2q")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 1.68/33.60;
      }
    }
    else if (!strcmp(name, "wzjetsto2qlnu")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.4552/7.649;
      }
    }
    else if (!strcmp(name, "wzjetsto3lnu")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.0661/1.1050;
      }
    }
    else if (!strcmp(name, "zzjetsto2l2nu")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.0194/0.3579;
      }
    }
    else if (!strcmp(name, "zzjetsto2l2q")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.0650/1.5213;
      }
    }
    else if (!strcmp(name, "zzjetsto4l")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.00944/0.18076;
      }
    }
    else if (!strcmp(name, "www")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.047;
      }
    }
    else if (!strcmp(name, "wwz")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.056;
      }
    }
    else if (!strcmp(name, "wzz")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.060;
      }
    }
    else if (!strcmp(name, "zzz")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.027;
      }
    }
    else if (!strcmp(name, "ttwjets") ||
	     !strcmp(name, "dyll50")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.05; // higher order Xs without error
      }
    }
    else if (!strcmp(name, "data_doublemu")) {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].staterr = 0.; // data has no stat error
	current_procs[i].systfactor = 0.; // data has no syst factor
      }
    }
    else {
      for (int i = 0; i < nbins; i++) {
	current_procs[i].systfactor = 0.5; // leading order Xs without error
      }
    }

    procs.push_back(current_procs);
    current_procs.clear();
  }

  fclose(infile);
  cout << "Read " << procs.size() << " processes" << endl;
  return procs.size();
}

void compute_xsec_uncertaintiy() {
  int maxProc = readprocs("final_binning_resolved.txt");

  double N_sum = 0; // Sum of all events
  double N_unc = 0; // Sum of events per process times xsec uncertainty

  for (int j = 0; j < maxProc; j++) {

    if (!strcmp(procs[j][0].name, "fakes") ||
	!strcmp(procs[j][0].name, "data_doublemu"))
	continue;
    
    double N_proc = 0;
    for (int i = 0; i < nbins; i++) {
      N_proc += procs[j][i].N;
    }
    N_unc += N_proc * procs[j][0].systfactor;
    N_sum += N_proc;
    INFO(procs[j][0].name << " " << N_proc *procs[j][0].systfactor );
  }

  INFO("Number of events: " << N_sum);
  INFO("Sum of all processes times xsec uncertainty: " << N_unc);
  INFO("Average xsec uncertainty: " << N_unc/N_sum);
}

void summarize_procs(const char * sel = "fullrun76", const char * hname = "jjmm_m")
{
  fill_fakes(sel, hname);
  int maxProc = readprocs("final_binning_resolved.txt");
  const double systerrglobal = TMath::Sqrt(0.6*0.6 + 0.1*0.1 + 0.6*0.6 + 3.6*3.6 + 0.2*0.2 + 0.2*0.2 + 1.0*1.0 + 1.5*1.5 + 2.6*2.6 + 5.0*5.0 + 6.0*6.0)/100.;
  double fakes_syst_factor = TMath::Sqrt(45.*45. + 10.*10. + 8.*8.)/100;
  

  ofstream out("sum_MC.txt");
  for (int i = 0; i < nbins; i++) {
    out << " Stage [" << i << "]" << endl;
    double sumN = 0;
    double sumErr2 = 0; // statistical errors
    double sumErr = 0; // systematic errors
    const char * name = 0;
    for (int j = 0; j < maxProc; j++) {
      name = 0;
      // give a summary at these names
      if (!strcmp(procs[j][i].name, "fakes")) 
	continue;

      if (!strcmp(procs[j][i].name, "ttwjets")) {
	name = "dyll";
      }
      else if (!strcmp(procs[j][i].name, "wwjetsto2l2nu")) {
	name = "tt+V";
      }
      else if (!strcmp(procs[j][i].name, "www")) {
	name = "VV";
      }
      else if (!strcmp(procs[j][i].name, "wminuswminus")) {
	name = "VVV";
      }
      else if (!strcmp(procs[j][i].name, "data_doublemu")) {
	name = "rare";
      }
      else
	name = 0;
      if (name != 0) {
	out << " + " << sumN << " +/- " << TMath::Sqrt(sumErr2 + sumErr*sumErr) << " " << name << endl;
	sumN = 0;
	sumErr2 = 0;
	sumErr = 0;
      }
      sumN += procs[j][i].N;
      sumErr2 += TMath::Power(procs[j][i].staterr, 2);
      sumErr += TMath::Sqrt(TMath::Power(procs[j][i].N*procs[j][i].systfactor, 2.) // xs uncertainty
					  + TMath::Power(procs[j][i].N*systerrglobal, 2.)); // global uncertainty
    }
    // output fake line
    out << " + " << fakes[i][0] << " +/- " << TMath::Sqrt(fakes[i][1] * fakes[i][1] + fakes[i][0] * fakes[i][0] * fakes_syst_factor * fakes_syst_factor) << " fakes" << endl;
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
  int maxProc = readprocs("final_binning_resolved.txt");
  // global systematic error, Table 5 of AN
  const double systerrglobal = TMath::Sqrt(0.6*0.6 + 0.1*0.1 + 0.6*0.6 + 3.6*3.6 + 0.2*0.2 + 0.2*0.2 + 1.0*1.0 + 1.5*1.5 + 2.6*2.6 + 5.0*5.0 + 6.0*6.0)/100.;
  cout << "global systematic error: " << systerrglobal << endl;
  
  double globN = 0;
  double globErrorCorrelated = 0; // correlated with other bins -> add linearly
  double globErrorUncorrelated2 = 0;
  // walk through table bins and columns and summarize
  for (int bin = 0; bin < nbins; bin++) {
    INFO("bin: " << bin);
    double binN = 0;
    double binErrorCorrelated = 0; // correlated with other bins -> add linearly
    double binErrorUncorrelated2 = 0; 
    for (int n = 0; n < maxProc; n++) {
      double locN = 0;
      double locErrorCorrelated = 0; // correlated with other bins -> add linearly
      double locErrorUncorrelated2 = 0;
      if (strcmp(procs[n][bin].name, "data_doublemu")) {
	locN = procs[n][bin].N;
	// MC cross-section
	locErrorCorrelated += TMath::Sqrt(TMath::Power(procs[n][bin].N*procs[n][bin].systfactor, 2.) // xs uncertainty
					  + TMath::Power(procs[n][bin].N*systerrglobal, 2.)); // global uncertainty
	// MC statistics
	locErrorUncorrelated2 += TMath::Power(procs[n][bin].staterr, 2.);
      }
      INFO(setw(16) << procs[n][bin].name 
	   << ": " << setw(precision+5) << procs[n][bin].N
	   << " +/- " << setw(precision+3) << TMath::Sqrt(locErrorCorrelated*locErrorCorrelated+locErrorUncorrelated2));

      binN += locN;
      binErrorCorrelated += locErrorCorrelated;
      binErrorUncorrelated2 += locErrorUncorrelated2;
    }
    INFO(setw(16) << "sum" 
	 << ": " << setw(precision+5) << binN 
	 << " +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorUncorrelated2+binErrorCorrelated*binErrorCorrelated)
	 << " [ +/- " << setw(precision+3) 
	 << TMath::Sqrt(binErrorUncorrelated2) << " (uncorr) "
	 << " +/- " << setw(precision+3) 
	 << binErrorCorrelated << " (corr) ]");
    globN += binN;
    globErrorCorrelated += binErrorCorrelated;
    globErrorUncorrelated2 += binErrorUncorrelated2;
  }
  INFO(setw(16) << "sum" 
       << ": " << setw(precision+5) << globN 
       << " +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2+globErrorCorrelated*globErrorCorrelated)
       << " [ +/- " << setw(precision+3) 
       << TMath::Sqrt(globErrorUncorrelated2) << " (uncorr) "
       << " +/- " << setw(precision+3)
       << globErrorCorrelated << " (corr) ]");

  // compute sum over rows with correlated systematics
  globN = 0;
  globErrorCorrelated = 0;
  globErrorUncorrelated2 = 0;
  for (int n = 0; n < maxProc; n++) {
    double procN = 0;
    double procErrorCorrelated = 0; // correlated with other bins
    double procErrorUncorrelated2 = 0;
    for (int bin = 0; bin < nbins; bin++) {
      double locN = 0;
      double locErrorCorrelated = 0;
      double locErrorUncorrelated2 = 0;
      locN += procs[n][bin].N;						       
      if (strcmp(procs[n][bin].name, "data_doublemu")) {
	locErrorCorrelated += TMath::Sqrt(TMath::Power(procs[n][bin].N*procs[n][bin].systfactor, 2)
					   + TMath::Power(procs[n][bin].N*systerrglobal, 2.));
	locErrorUncorrelated2 += TMath::Power(procs[n][bin].staterr, 2);
      }
      procN += locN;
      procErrorCorrelated += locErrorCorrelated;
      procErrorUncorrelated2 += locErrorUncorrelated2;
    }
    INFO(setw(16) << procs[n][0].name 
	 << ": " << setw(precision+5) << procN 
	 << " +/- " << setw(precision+3)
	 << TMath::Sqrt(procErrorUncorrelated2+procErrorCorrelated*procErrorCorrelated)
	 << " [ +/- " << setw(precision+3)
	 << TMath::Sqrt(procErrorUncorrelated2) << " (uncorr) "
	 << " +/- " << setw(precision+3)
	 << procErrorCorrelated << " (corr) ]");
    if (strcmp(procs[n][0].name, "data_doublemu")) {
      globN += procN;
      globErrorCorrelated += procErrorCorrelated;
      globErrorUncorrelated2 += procErrorUncorrelated2;
    }
  }
  INFO(setw(16) << "sum" 
       << ": " << setw(precision+5) << globN
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
