#include "BTagging.h"
#include "TFile.h"

#include "Utilities.h"

// Choose Track Counting High Efficiency with working point medium (CVSM)
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
const double btag_cut = 0.679;

BTagging::BTagging(const char * EffMapFile)
{
  // initialize private random number generator with given seed
  rand_ = new TRandom3(35535);

  // load efficiency maps
  TFile * f = new TFile(EffMapFile, "READ");
  btag_eff_map = (TH2D*) f->Get("h2_0_btag_num_b");
  ctag_eff_map = (TH2D*) f->Get("h2_0_btag_num_c");
  bmistag_eff_map = (TH2D*) f->Get("h2_0_btag_num_l");
}

BTagging::~BTagging()
{
  // delete private random number generator
  delete rand_;
}

bool BTagging::isBJet(double btag, int pdgIdPart, double pt, double eta, double phi)
{
  eta = TMath::Abs(eta);
  pdgIdPart = TMath::Abs(pdgIdPart);
  UInt_t seed = UInt_t(TMath::Abs(phi)/TMath::Pi()*100000) % 10000;
  rand_->SetSeed(seed);
  bool isBTagged = (btag > btag_cut);
  if (eta > 2.4) {
    WARNING("jet properties out of allowed range in BTagging::isBJet(): pt = " << pt 
	    << ", eta = " << eta);
    return isBTagged;
  }
  modifyBTagsWithSF(isBTagged,
		    pdgIdPart,
		    GetBTagScaleFactor(pt),
		    GetBTagScaleFactorError(pt),
		    GetBTagEfficiency(pt, eta),
		    GetCTagEfficiency(pt, eta),
		    GetBMisTagScaleFactor(pt, eta),
		    GetBMisTagEfficiency(pt, eta)
		    );
  return isBTagged;
}

double BTagging::GetBTagScaleFactor(double pt)
{
  return (0.938887+(0.00017124*pt))+(-2.76366e-07*(pt*pt));
}


double BTagging::GetBTagScaleFactorError(double pt)
{
  float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
  float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  double  SFb_error[] = {
    0.0415707,
    0.0204209,
    0.0223227,
    0.0206655,
    0.0199325,
    0.0174121,
    0.0202332,
    0.0182446,
    0.0159777,
    0.0218531,
    0.0204688,
    0.0265191,
    0.0313175,
    0.0415417,
    0.0740446,
    0.0596716
  };

  for (unsigned int i = 0; i < sizeof(ptmin)/sizeof(double); i++) {
    if (pt >= ptmin[i] && pt < ptmax[i]) {
      return SFb_error[i];
    }
  }
  // as stated on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2011_Data_and_MC
  if (pt < ptmin[0]) {
    // return constant error if pT < 30 GeV
    return 0.12;
  }
  else {
    // above 800 GeV, return twice the error from 800 GeV
    return 2*SFb_error[sizeof(ptmin)/sizeof(double)-1];
  }
}

double BTagging::GetBMisTagScaleFactor(double pt, double eta)
{
  // from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_EPS2013.C
  eta = TMath::Abs(eta);
  if (eta < 0.8) {
    if (systematic == "up") {
      return ((1.18638+(0.00314148*pt))+(-6.68993e-06*(pt*pt)))+(3.89288e-09*(pt*(pt*pt)));
    } else if (systematic == "down") {
      return ((0.964527+(0.00149055*pt))+(-2.78338e-06*(pt*pt)))+(1.51771e-09*(pt*(pt*pt)));
    }

    return ((1.07541+(0.00231827*pt))+(-4.74249e-06*(pt*pt)))+(2.70862e-09*(pt*(pt*pt)));
  }
  else if (eta < 1.6) {
    if (systematic == "up") {
      return ((1.16624+(0.00151884*pt))+(-3.59041e-06*(pt*pt)))+(2.38681e-09*(pt*(pt*pt)));
    } else if (systematic == "down") {
      return ((0.946051+(0.000759584*pt))+(-1.52491e-06*(pt*pt)))+(9.65822e-10*(pt*(pt*pt)));
    }

    return ((1.05613+(0.00114031*pt))+(-2.56066e-06*(pt*pt)))+(1.67792e-09*(pt*(pt*pt)));
  }
  else if (eta < 2.4) {
    if (systematic == "up") {
      return ((1.15575+(0.000693344*pt))+(-3.02661e-06*(pt*pt)))+(2.39752e-09*(pt*(pt*pt)));
    } else if (systematic == "down") {
      return ((0.956736+(0.000280197*pt))+(-1.42739e-06*(pt*pt)))+(1.0085e-09*(pt*(pt*pt)));
    }

    return ((1.05625+(0.000487231*pt))+(-2.22792e-06*(pt*pt)))+(1.70262e-09*(pt*(pt*pt)));
  }
  else {
    // same value as for eta < 2.4, but twice the uncertainty
    if (systematic == "up") {
      return ((1.15575+(0.000693344*pt))+(-3.02661e-06*(pt*pt)))+(2.39752e-09*(pt*(pt*pt)));
    } else if (systematic == "down") {
      return ((0.956736+(0.000280197*pt))+(-1.42739e-06*(pt*pt)))+(1.0085e-09*(pt*(pt*pt)));
    }

    return ((1.05625+(0.000487231*pt))+(-2.22792e-06*(pt*pt)))+(1.70262e-09*(pt*(pt*pt)));
  }
}

double BTagging::GetBMisTagEfficiency(double pt, double eta)
{
  if (pt > 800.) pt = 800.;
  if (eta > 2.4) eta = 2.4;
  return bmistag_eff_map->GetBinContent(bmistag_eff_map->GetXaxis()->FindBin(pt), bmistag_eff_map->GetYaxis()->FindBin(eta));
}

double BTagging::GetBTagEfficiency(double pt, double eta)
{
  if (pt > 800.) pt = 800.;
  if (eta > 2.4) eta = 2.4;
  return btag_eff_map->GetBinContent(btag_eff_map->GetXaxis()->FindBin(pt), btag_eff_map->GetYaxis()->FindBin(eta));
}

double BTagging::GetCTagEfficiency(double pt, double eta)
{
  if (pt > 800.) pt = 800.;
  if (eta > 2.4) eta = 2.4;
  return ctag_eff_map->GetBinContent(ctag_eff_map->GetXaxis()->FindBin(pt), ctag_eff_map->GetYaxis()->FindBin(eta));
}

void BTagging::modifyBTagsWithSF(bool & isBTagged, int pdgIdPart,
				 double Btag_SF, double Btag_SFerr, double Btag_eff, double Ctag_eff,
				 double Bmistag_SF, double Bmistag_eff)
{
  // assume no modification
  pdgIdPart = TMath::Abs(pdgIdPart);
  if (systematic == "up") {
    Btag_SF += Btag_SFerr;
  } else if (systematic == "down") {
    Btag_SF -= Btag_SFerr;
  }

  if(pdgIdPart == 5) {
    // b quarks
    isBTagged = applySF(isBTagged, Btag_SF, Btag_eff);
  }
  else if (pdgIdPart == 4)  {
    // c quarks
    isBTagged = applySF(isBTagged, Btag_SF, Ctag_eff);
  }
  else if((pdgIdPart >= 1 && pdgIdPart <= 3) || pdgIdPart == 21) {
    // light quarks, gluons
    isBTagged = applySF(isBTagged, Bmistag_SF, Bmistag_eff);
  }
  else if (pdgIdPart == 0) {
    // in data it is 0, do nothing
  }
  else {
    // treat everything else as unaffected by b-tagging
    DEBUG("unidentified particle " << pdgIdPart << " given to BTagging::modifyBTagsWithSF()");
  }
}

bool BTagging::applySF(bool isBTagged, double Btag_SF, double Btag_eff)
{
  bool newBTag = isBTagged;

  if (Btag_SF == 1)
    return newBTag; // no correction needed

  // throw dice
  double coin = rand_->Uniform(1.);

  if(Btag_SF > 1){  // use this if SF>1
    if( !isBTagged ) {
      // fraction of jets that need to be upgraded
      double mistagPercent = (1.0 - Btag_SF) / (1.0 - (1.0/Btag_eff) );

      // upgrade to tagged
      if( coin < mistagPercent ) {
	newBTag = true;
      }
    }
  }
  else {  // use this if SF<1
    // downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {
      newBTag = false;
    }
  }

  return newBTag;
}

void BTagging::SetSystematic(string sys) {
  systematic = sys;
}
