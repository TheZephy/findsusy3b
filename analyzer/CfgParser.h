#ifndef _CfgParser_h_
#define _CfgParser_h_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
	
using namespace std;

struct cfgEntry {
  map<string, double> v;
  map<string, vector<double> > vv;
  map<string, string> s;
  map<string, vector<string> > ss;
  map<string, bool> b;
  map<string, vector<bool> > bb;
  map<string, cfgEntry> sub;
  struct cfgEntry * top;
};

class CfgParser {
 protected:
  bool IsEven(int);
  bool IsOdd(int);
  void ReadLines(std::ifstream &);
  void StoreLine(cfgEntry &,string, string);

  void FindCharNotInString(size_t &, char, string &, size_t startPos = 0);
  int RetrieveDepthLevel(string);
  string RetrieveEntryTitle(string, int &);
  
  
 public:
  CfgParser(const char *);
  ~CfgParser();
  void PrintCfgMap (map<string, cfgEntry> &, string);
  map<string, cfgEntry> entries;
  
};
#endif
