#include "BTagging.h"

#include "Utilities.h"

// Choose Track Counting High Efficiency with working point medium (CVSM)
const double btag_cut = 0.679;

BTagging::BTagging(int seed)
{
  // initialize private random number generator with given seed
  rand_ = new TRandom3(seed);
}

BTagging::~BTagging()
{
  // delete private random number generator
  delete rand_;
}

bool BTagging::isBJet(double btag, int pdgIdPart, double pt, double eta)
{
  eta = TMath::Abs(eta);
  pdgIdPart = TMath::Abs(pdgIdPart);
  bool isBTagged = (btag > btag_cut);
  if (eta > 2.4) {
    WARNING("jet properties out of allowed range in BTagging::isBJet(): pt = " << pt 
	    << ", eta = " << eta);
    return isBTagged;
  }
  modifyBTagsWithSF(isBTagged, 
		    pdgIdPart, 
		    GetBTagScaleFactor(btag),
		    GetBTagEfficiencyMC(btag),
		    GetCTagEfficiencyMC(btag),
		    GetBMisTagScaleFactor(pt, eta),
		    GetBMisTagEfficiency(pt, eta)
    );
  return isBTagged;
}

double BTagging::GetBTagScaleFactor(double pt)
{
  return 0.726981*((1.+(0.253238*pt))/(1.+(0.188389*pt)));
}

double BTagging::GetBTagScaleFactorError(double pt)
{
  float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
  float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  double SFb_error[] = {
    0.0554504,
    0.0209663,
    0.0207019,
    0.0230073,
    0.0208719,
    0.0200453,
    0.0264232,
    0.0240102,
    0.0229375,
    0.0184615,
    0.0216242,
    0.0248119,
    0.0465748,
    0.0474666,
    0.0718173,
    0.0717567
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

double BTagging::GetBMisTagScaleFactor(double pt, double eta, bool finebin)
{
  // from https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_Moriond2013.C
  eta = TMath::Abs(eta);
  if (finebin) {
    if (eta < 0.8) {
      return ((1.06238+(0.00198635*pt))+(-4.89082e-06*(pt*pt)))+(3.29312e-09*(pt*(pt*pt)));
    }
    else if (eta < 1.6) {
      return ((1.08048+(0.00110831*pt))+(-2.96189e-06*(pt*pt)))+(2.16266e-09*(pt*(pt*pt)));
    }
    else if (eta < 2.4) {
      return ((1.09145+(0.000687171*pt))+(-2.45054e-06*(pt*pt)))+(1.7844e-09*(pt*(pt*pt)));
    }
    else {
      // same value as for eta < 2.4, but twice the uncertainty
      return ((1.09145+(0.000687171*pt))+(-2.45054e-06*(pt*pt)))+(1.7844e-09*(pt*(pt*pt)));
    }
  }
  else {
    // has been disabled for 2012 CSVM
    return ((1.09145+(0.000687171*pt))+(-2.45054e-06*(pt*pt)))+(1.7844e-09*(pt*(pt*pt)));
  }
}
    
double BTagging::GetBMisTagEfficiency(double pt, double eta, bool finebin)
{
  eta = TMath::Abs(eta);
  if (finebin) {
    if (eta < 0.8) {
      return (0.00967751+(2.54564e-05*pt))+(-6.92256e-10*(pt*pt)) * (1.10422 + -0.000523856*pt + 1.14251e-06*pt*pt);
    }
    else if (eta < 1.6) {
      return (0.00974141+(5.09503e-05*pt))+(2.0641e-08*(pt*pt)) * (1.10422 + -0.000523856*pt + 1.14251e-06*pt*pt);
    }
    else if (eta < 2.4) {
      return (0.013595+(0.000104538*pt))+(-1.36087e-08*(pt*pt)) * (1.10422 + -0.000523856*pt + 1.14251e-06*pt*pt);
    }
    else {
      // same value as for eta < 2.4, but twice the uncertainty
      return (0.013595+(0.000104538*pt))+(-1.36087e-08*(pt*pt)) * (1.10422 + -0.000523856*pt + 1.14251e-06*pt*pt);
    }
  }
  else {
    // covering full range in eta
    return (0.0113428+(5.18983e-05*pt))+(-2.59881e-08*(pt*pt)) * (1.10422 + -0.000523856*pt + 1.14251e-06*pt*pt);
  }
}

double BTagging::GetBTagEfficiencyData(double btag)
{
  return -3.67153247396e-07*btag*btag*btag*btag +  -2.81599797034e-05*btag*btag*btag
    + 0.00293190163243*btag*btag +  -0.0849600849778*btag +  0.928524440715;
}

double BTagging::GetBTagEfficiencyDataError(double btag)
{
  return 3.03337430722e-06*btag*btag*btag*btag + -0.000171604835897*btag*btag*btag
    + 0.00474711667943*btag*btag + -0.0929933040514*btag + 0.978347619293
    - GetBTagEfficiencyDataError(btag);
}

double BTagging::GetBTagEfficiencyMC(double btag)
{
  return -1.73338329789*btag*btag*btag*btag +  1.26161794785*btag*btag*btag +  0.784721653518*btag*btag +  -1.03328577451*btag +  1.04305075822;
}

double BTagging::GetCTagEfficiencyMC(double btag)
{
  return -1.5734604211*btag*btag*btag*btag +  1.52798999269*btag*btag*btag +  0.866697059943*btag*btag +  -1.66657942274*btag +  0.780639301724;
}

void BTagging::modifyBTagsWithSF(bool & isBTagged, int pdgIdPart,
				 double Btag_SF, double Btag_eff, double Ctag_eff,
				 double Bmistag_SF, double Bmistag_eff)
{
  // assume no modification
  pdgIdPart = TMath::Abs(pdgIdPart);

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
      double mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );

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
