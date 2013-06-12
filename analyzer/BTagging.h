#ifndef _BTagging_h_
#define _BTagging_h_

// code partly borrowed from https://twiki.cern.ch/twiki/pub/CMS/BTagSFUtil/BTagSFUtil.h
// explanation on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG

#include "TRandom3.h"
#include "TMath.h"

class BTagging 
{
public: 
  BTagging(int seed = 0);
  ~BTagging();

  // modify the b-tag of the jet according to its pt and eta
  bool isBJet(double btag, int pdgIdPart, double pt, double eta);

private:
  double GetBTagScaleFactor(double pt);
  double GetBTagScaleFactorError(double pt);
  double GetBMisTagScaleFactor(double pt, double eta, bool finebin = true);
  double GetBMisTagEfficiency(double pt, double eta, bool finebin = true);
  double GetBTagEfficiencyData(double btag);
  double GetBTagEfficiencyDataError(double btag);
  double GetBTagEfficiencyMC(double btag);
  double GetCTagEfficiencyMC(double btag);

  void modifyBTagsWithSF(bool & isBTagged, int pdgIdPart,
			 double Btag_SF, double Btag_eff, double Ctag_eff,
			 double Bmistag_SF, double Bmistag_eff);
  bool applySF(bool isBTagged, double Btag_SF, double Btag_eff);


  TRandom3* rand_;
};

#endif
