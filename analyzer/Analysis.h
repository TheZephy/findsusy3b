#ifndef __Analyis_h
#define __Analyis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TEnv.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "BTagging.h"
#include "RunLumiRanges.h"
#include "EventFilter.h"
#include "CfgParser.h"

#include "TreeContent.h"
#include "LumiReweightingStandAlone.h"

#ifndef __CINT__
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#else
class JetCorrectorParameters;
class JetCorrectionUncertainty;
#endif /* __CINT __ */

using namespace std;

struct ssys {
ssys(string * v, string s, string t) : var(v), sys(s), tag(t) {}
  string * var;
string sys;
string tag;
};
  
  
struct dsys {
dsys(double * v, double s, string t) : var(v), sys(s), tag(t) 
{}
  double * var;
double sys;
string tag;
};

struct bsys {
bsys(bool * v, bool s, string t) : var(v), sys(s), tag(t) {}
  bool * var;
bool sys;
string tag;
};


class Analysis : public TreeContent {
protected:
  TTree &   fInputTree;
  TTree &   fOutputTree;
  CfgParser  &   fCfgFile;

  //////////////////////////////////////////////////////////////////////
  // configuration options

  bool      fFill;      // fill output tree?
  string    fRunTag;     // name of the run
  string    fInputType; // type of input data
  string    fSample;    // the sample name
  string    fAnalysisType; // type of analysis to be performed
  Long64_t  fMaxEvents; // max events to process
  Long64_t  fMaxTreeSize; // maximum tree size (output)
  bool      fFindDuplicates; // find duplicate events?
  bool      fDumpAll;   // Dump all event information
  bool      fDumpBasic; // Dump basic information
  bool      fDumpTruth; // Dump MC truth information
  string    fDumpTrigger; // Dump Trigger information
  bool      fPileupReweighting; // use pileup reweighting?
  bool      fSkimActive; // Skimming Active?
  Int_t     fSkimMuons;  // how many muons to require
  Double_t  fSkimMuoptfirst; // cut on muon with highest pt
  Double_t  fSkimMuoptother; // cut on muon with second highest pt
  Double_t  fLooseMuonRelIso; // relative isolation cut for loose muon
  string    fFakeRateMethod; // fake rate method, "lastbin" or "zero"
  Double_t     fFakeRateDimensions; // number of dimensions for fake rate calculation
  vector<string> fTrigger; // trigger selection
  bool      fForceUnprescaledTriggers; // Force unprescaled triggers
  string    fLumiRanges; // lumi ranges (JSON file) for data processing
  string    fEventFilter; // zipped event filter file for data processing
  Double_t  fDeltaPhiMin; // max delta phi (leading mu, gaugino)
  bool fFirstWarnWeight;
  TString fCurrentTrigger;

  // stage histograms
  vector<const char *>	stageHistoNames;
  vector<const char *>	stageHistoXTitles;
  vector<Int_t>		stageHistoNBinsX;
  vector<Double_t>	stageHistoXLow;
  vector<Double_t>	stageHistoXUp;
  vector<Int_t>		stageHistoNBinsY;
  vector<Double_t>	stageHistoYLow;
  vector<Double_t>	stageHistoYUp;
  vector<Double_t *>	stageHistoFillValuesX;
  vector<Double_t *>	stageHistoFillValuesY;

  
  // cuts for T/L ratio
  Double_t fTL_met_max;
  Double_t fTL_ht_min;
  Double_t fTL_jetpt_min;
  Double_t fTL_nloose_min;
  Double_t fTL_nloose_max;
  Double_t fTL_mumudz_min;
  Double_t fTL_zmass_min;
  Double_t fTL_zmass_max;
  Double_t fTL_firstmupt_min;
  Double_t fTL_mupt_min;
  Double_t fTL_jetdphi_min;
  Double_t fTL_mumudphi_max;
  Double_t fTL_mt_max;
  Double_t fTL_njets_min;

  // values for b-tagging (BTAG)
  Double_t fBTAG_threshold;
  string fBTAG_systematic;

  // jet energy scale systematics (JES)
  JetCorrectionUncertainty * fJES_JetCorrectionUncertainty;
  string fJES_systematic;

  // values for smearing jet energies (JER)
  Bool_t   fJER_calculation;
  string   fJER_systematic;
  Double_t fJER_center;
  Double_t fJER_smear;
  Double_t fJER_met_old;

  vector<Double_t> fJER_jet_E;
  vector<Double_t> fJER_jet_Et;
  vector<Double_t> fJER_jet_p;
  vector<Double_t> fJER_jet_pt;
  vector<Double_t> fJER_jet_px;
  vector<Double_t> fJER_jet_py;
  vector<Double_t> fJER_jet_pz;

  Double_t fJER_met_et;
  Double_t fJER_met_sumet;
  Double_t fJER_met_ex;
  Double_t fJER_met_ey;
  Double_t fJER_met_phi;

  // muon energy resolution systematics (MER)
  Bool_t fMER_systematic;

  // muon energy resolution systematics (MER)
  string fMES_systematic;

  //////////////////////////////////////////////////////////////////////
  // analysis variables

  // signal study
  bool             fIsSignal;
  TLorentzVector * fSigMu0; // only in case of signal: generated muons
  TLorentzVector * fSigMu1;
  TLorentzVector * fSigJet0; // only in case of signal: generated jets
  TLorentzVector * fSigJet1; 

  // fake rate
  TH3D           * fFakeRateHisto3D;
  TH2D           * fFakeRateHisto2D;
  Double_t         fSingleFakeWeight;  // fake rate weight for one mu fluctuation 
  Double_t         fFakeRate[2];       // fake rate for the two muons

  // save particles for later analysis
  Int_t            fMuoId[2]; // ID's of selected muons
  TLorentzVector * fMuon0;
  TLorentzVector * fMuon1;
  Int_t            fJets;
  Int_t            fJetId[40]; // ID's of selected main jets
  TLorentzVector * fJet0;
  TLorentzVector * fJet1;
  Bool_t           fIsBTagged;
  TLorentzVector * fSmuon;
  TLorentzVector * fGaugino;

public:
  Analysis(TTree & inputTree, TTree & outputTree, CfgParser & cfgFile);
  virtual ~Analysis();
  
  void SetBranchAddresses();
  void CreateBranches();
  void CreateVariables();
  void StartLoop();
  void Loop();

protected:
  //////////////////////////////////////////////////////////////////////
  // config functions
  void GetValue(string &, cfgEntry &, string);
  void GetValue(bool &, cfgEntry &, string);
  void GetValue(double &, cfgEntry &, string);

  vector<ssys> sysStrings;
  vector<dsys> sysDoubles;

  vector<bsys> sysBools;
  

  //////////////////////////////////////////////////////////////////////
  // analysis functions
  void ReconstructFinalState(map<char, vector<int> > & particles, int vertex);
  void GetSimplifiedModelParticles(map<string, vector<int> > & particles, int vertex);
  int IsSimplifiedModel(map<string, vector<int> > & particles);
  void SignalStudy(int & charge);
  void TriggerMatchingComparison( int, const char *);
  void TightLooseRatioCalculation(const vector<int> & loose_muons,
				  const vector<int> & tight_muons,
				  const vector<int> & jets,
				  const double HT);
  double GetFakeRate(double muopt, double eta, double jetpt);
  double GetJERScale(double eta, string sys = "none");
  void PFJetSmearing();
  void PFJetSmearingCalculation();
  void JESandRecalculateMET(TString);
  void MuonEnergySmearing();
  void MuonEnergyScale(TString);
  void BTagEfficiencyMap();

  Bool_t TriggerMatched(Int_t muonIterator, vector <TString> triggerFilters);
  Bool_t dRMatched(Int_t muonIterator, Int_t triggerID);
  Bool_t MuonCuts(Int_t i);
  Bool_t JetCuts(Int_t i, vector <int> muons);
  Bool_t ElectronCuts(Int_t i);

  // helper functions
  void Analyze (Long64_t &, EventFilter &, lumi::RunLumiRanges &);
  void CreateHistograms();
  void CreateTaggedHistograms();
  void AddStageHisto(const char * name, const char * xtitle, Int_t nbinsx, Double_t xup, Double_t xlow, Double_t * fillvaluex);
  void AddStageHisto(const char * name, const char * xtitle, Int_t nbinsx, Double_t xup, Double_t xlow, Int_t nbinsy, Double_t yup, Double_t ylow, Double_t * fillvaluex, Double_t * fillvaluey);
  void RemoveStageHisto (string name);
  void FillStage(const char * stageTag);
  void ClearStages();
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup, 
		   Int_t nbinsy, Double_t ylow, Double_t yup);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, const Double_t * xbins, 
		   Int_t nbinsy, Double_t ylow, Double_t yup);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, const Double_t * xbins, 
		   Int_t nbinsy, const Double_t * ybins);
  void CreateHisto(const char * name, const char * title, 
		   Int_t nbinsx, Double_t xlow, Double_t xup, 
		   Int_t nbinsy, Double_t ylow, Double_t yup, 
		   Int_t nbinsz, Double_t zlow, Double_t zup);
  void Fill(const char * name, double value);
  void Fill(const char * name, const char * bin);
  void FillWithWeight(const char * name, double value, double weight);
  void Fill(const char * name, double x, double y);
  void Fill(const char * name, double x, double y, double z);

  void StoreSmearingObjects();
  void RestoreSmearingObjects();

  // taken over from ACSUSYAna
  void BasicDump(int i);
  void TriggerDump(TString sel);
  void MuonDump(bool full=1);
  void CaloJetDump();
  void PFJetDump();
  void TruthJetDump();
  void TruthDump();
  void VertexDump();
  void METDump();
  void SCDump();
  void EleDump(bool full=1);
  void PFEleDump(bool full=1);

  bool FindDuplicates(int run, int evt, double x1, double x2);
  Double_t DeltaPhi(double a, double b);
  
  typedef std::pair< pair<int,int> , pair<double,double> > Key;
  typedef std::set<Key> KeySet;
  typedef KeySet::const_iterator KeyIter;
  KeySet _keys;

  // this is the magic for pileup reweighting
  reweight::LumiReWeighting LumiWeights_;
  /* reweight::PoissonMeanShifter PShiftUp_; */
  /* reweight::PoissonMeanShifter PShiftDown_; */ 

  // b-tagging correction
  BTagging fBTagging;

  // histogram store
  map<string, TH1D *> histo;
  map<string, TH2D *> histo2;
  map<string, TH3D *> histo3;

  friend class truth_pt_comparator;
};

class truth_pt_comparator 
{
protected:
  const Analysis & fAnalysis;
public:
  truth_pt_comparator(const Analysis & analysis);
  bool operator()(int i, int j);
};

#endif
