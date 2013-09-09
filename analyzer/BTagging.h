#ifndef _BTagging_h_
#define _BTagging_h_

// code partly borrowed from https://twiki.cern.ch/twiki/pub/CMS/BTagSFUtil/BTagSFUtil.h
// explanation on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG

#include "TRandom3.h"
#include "TMath.h"
#include "TH2.h"

using namespace std;

class BTagging 
{
 public: 
  BTagging(const char * EffMapFile);
  ~BTagging();

  // modify the b-tag of the jet according to its pt and eta
  bool isBJet(double btag, int pdgIdPart, double pt, double eta, double phi);
  void SetSystematic(string sys);
 protected:
  TH2D * btag_eff_map;
  TH2D * ctag_eff_map;
  TH2D * bmistag_eff_map;

  string systematic;

  
 private:
  double GetBTagScaleFactor(double pt);
  double GetBTagScaleFactorError(double pt);
  double GetBMisTagScaleFactor(double pt, double eta);
  double GetBMisTagEfficiency(double pt, double eta);
  double GetBTagEfficiency(double pt, double eta);
  double GetCTagEfficiency(double pt, double eta);

  void modifyBTagsWithSF(bool & isBTagged, int pdgIdPart,
			 double Btag_SF, double Btag_SFerr, double Btag_eff, double Ctag_eff,
			 double Bmistag_SF, double Bmistag_eff);
  bool applySF(bool isBTagged, double Btag_SF, double Btag_eff);


  TRandom3* rand_;
};

#endif
