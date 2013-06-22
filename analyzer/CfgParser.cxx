#include "CfgParser.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

// trim functions
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

using namespace std;

CfgParser::CfgParser (const char * configPath) {
  std::ifstream ifs(configPath);
  if( !ifs.good() ) {
    throw "cannot open file";
  }
  
  ReadLines(ifs);
  //PrintCfgMap(entries, "All entries");
}

CfgParser::~CfgParser () {}

void CfgParser::ReadLines(ifstream & ifs) {
  // keep track of lines
  int lineCount = 0;
  try {
    // read lines one after the other
    string line;
    cfgEntry * currentEntry = 0;
    int currentDepthLevel = 1;


    while( std::getline(ifs, line) ) {
      lineCount++;

      // enable comments leading with a '#'
      size_t pos;
      FindCharNotInString(pos, '#', line);

      // trim off outcommented parts
      if ( pos != string::npos ) {      
	line = line.substr(0,pos);
      }

      // check for entry creation lines
      int depthLevel = RetrieveDepthLevel(line);
      // 0 indicates that this is not a new entry
      if (depthLevel != 0) {
	string title = RetrieveEntryTitle(line, depthLevel);
	// check if this is a subentry
	if ( depthLevel > currentDepthLevel ) {
	  if (currentEntry == 0 ) throw "Syntax Error - Wrong amount of '['?";

	  (*currentEntry).sub[title].top = currentEntry;
	  currentEntry = &((*currentEntry).sub[title]);
	  currentDepthLevel = depthLevel;

	} else {
	  // if it is not a subentry check which entry it is a sub to
	  if ( currentEntry != 0 ) {

	    // go to that level
	    for (int i = currentDepthLevel; i >= depthLevel; i--) {
	      currentEntry = (*currentEntry).top;
	    }
	  }
	  

	  if ( currentEntry != 0 ) {

	    // unless it is a completely new branch, enter it there	   
	    (*currentEntry).sub[title].top = currentEntry;
	    currentEntry = &((*currentEntry).sub[title]);
	    currentDepthLevel = depthLevel;

	  } else {
	    
	    // or if it is an entirely new branch of entries, create it
	    currentEntry =  &(entries[title]);
	    (*currentEntry).top = 0;
	    currentDepthLevel = depthLevel;
	  }
	}
      } else {

	 // if its not a new entry, its probably a value
	std::istringstream iss(line);
	string key;
	if( std::getline(iss, key, '=') ) {
	  if (!key.empty()) {
	    //cout << key << endl;
	    string value;
	    if( std::getline(iss, value) ) {
	      StoreLine((*currentEntry), key, value);
	    }
	  }
	}
      }
    }

  } catch (const char * s) {
    cout << lineCount << ": " << s << endl;
    exit(1);
  }
}


// trim from start
static inline string & ltrim(string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline string & rtrim(string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline string & trim(string &s) {
  return ltrim(rtrim(s));
}

void CfgParser::StoreLine(cfgEntry & currentEntry, string key, string value) {
  key = trim(key);
  value = trim(value);

  size_t pos, oldPos = -1;
  FindCharNotInString(pos, ',', value);
  
  // determine whether it is an array or just a single value
  if ( pos != string::npos ) {
    size_t valStart, valLength;
    valStart = oldPos+1; 
    valLength = pos-valStart;
    string trimmedVal = value.substr(valStart, valLength);
    trimmedVal = trim(trimmedVal);
    // array
    if ( trimmedVal[0] == '"' ) {
      // iterate over values
      while ( pos != string::npos ) {
	trimmedVal = value.substr(valStart, valLength);

	trimmedVal.erase(0, trimmedVal.find('"')+1);
	trimmedVal.erase(trimmedVal.rfind('"'), trimmedVal.length());
	currentEntry.ss[key].push_back(trimmedVal);

	oldPos = pos;
	FindCharNotInString(pos, ',', value, pos+1);
	valStart = oldPos+1; 
	valLength = pos-valStart;
	
      }

      // check for trailing value
      if ( value.length() > valStart ) {
	cout << value.length() << " " << valStart << " " << valLength << "   " << valStart+valLength+2 << endl;
	valStart = oldPos+1;
	valLength = value.length()-1;
	cout << value.length() << " " << valStart << " " << valLength << endl;
	trimmedVal = value.substr(valStart, valLength);
	trimmedVal.erase(0, trimmedVal.find('"')+1);
	trimmedVal.erase(trimmedVal.rfind('"'), trimmedVal.length());
	
	currentEntry.ss[key].push_back(trimmedVal);
      }
    } else if ( trimmedVal == "true" || trimmedVal == "false" ) { 
      // bool
      while ( pos != string::npos ) {
	valStart = oldPos+1; 
	valLength = pos-valStart;

	trimmedVal = value.substr(valStart, valLength);	
	trimmedVal = trim(trimmedVal);
	currentEntry.bb[key].push_back((trimmedVal == "true") ? true : false);

	oldPos = pos;
	FindCharNotInString(pos, ',', value, pos+1);	
      }

      // check for trailing value
      if ( value.length() > valStart ) {
	valStart = oldPos+1;
	valLength = value.length()-1;
	
	trimmedVal = value.substr(valStart, valLength);	
	trimmedVal = trim(trimmedVal);
	currentEntry.bb[key].push_back((trimmedVal == "true") ? true : false);
      }

    } else {
      
      // doubles
      // iterate over values
      while ( pos != string::npos ) {
	valStart = oldPos+1; 
	valLength = pos-valStart;

	std::stringstream ss;
	ss << value.substr(valStart, valLength);
	double doubleVal;
	ss >> doubleVal;
	currentEntry.vv[key].push_back(doubleVal);
	
	oldPos = pos;
	FindCharNotInString(pos, ',', value, pos+1);	
      }

      // check for trailing value
      if ( value.length() > valStart ) {
	valStart = oldPos+1;
	valLength = value.length()-1;

	std::stringstream ss;
	ss << value.substr(valStart, valLength);
	double doubleVal;
	ss >> doubleVal;
	currentEntry.vv[key].push_back(doubleVal);
      }
    }
    
  } else {
    // single value
    if ( value[0] == '"' ) {
      
      // string
      string trimmedVal = value;
      trimmedVal.erase(0, trimmedVal.find('"')+1);
      trimmedVal.erase(trimmedVal.rfind('"'), trimmedVal.length());
      currentEntry.s[key] = trimmedVal;
    } else if ( value == "true" ) { 
      // bool
      currentEntry.b[key] = true;      
    }
    else if ( value == "false" ) {
      currentEntry.b[key] = false;
    } else {
      // double
      std::stringstream ss;
      ss << value;
      double doubleVal;
      ss >> doubleVal;
      currentEntry.v[key] = doubleVal;
    }
  }
}

void CfgParser::FindCharNotInString(size_t & pos, char charToFind, string & line, size_t startPos) {
  
  // firstQuotePos = startPos
  size_t charPos = line.find(charToFind, startPos);
  if ( charPos == string::npos ) {
    pos = string::npos;
    return;
  }
  
  size_t firstQuotePos = line.find('"', startPos);
  if ( charPos < firstQuotePos ) {
    pos = charPos;
    return;
  }

  size_t secondQuotePos = line.find('"', firstQuotePos+1);
  if ( secondQuotePos == string::npos ) throw "Syntax error - Missing quotation mark";
  FindCharNotInString(pos, charToFind, line, secondQuotePos+1);
}


int CfgParser::RetrieveDepthLevel (string line) {
  int depthLevel = 0;
  line = trim(line);
  string::iterator it=line.begin();
  string::reverse_iterator rit=line.rbegin();
  while ( it!=line.end() && *it == '[' ) {
    if (*rit == ']') {
      depthLevel++;
    } else {
      throw "Syntax error - Missing parenthesis?";
    }
    ++it; ++rit;
  }
  return depthLevel;
}

string CfgParser::RetrieveEntryTitle (string line, int & depthLevel) {
  line = trim(line);
  line.erase(0, depthLevel);
  line.erase(line.length()-depthLevel, line.length());
  line = trim(line);
  return line;
}

void CfgParser::PrintCfgMap (map<string, cfgEntry> & mymap, string maptitle) {
  if ( !(mymap.empty()) ) {
    cout << "START of " << maptitle << endl;
    for (map<string, cfgEntry>::iterator it=mymap.begin(); it!=mymap.end(); ++it) {
      cout << it->first  << endl;
      
      PrintCfgMap( it->second.sub, it->first );
    }

    cout << "END of " << maptitle << endl;
  }
}




bool CfgParser::IsEven(int a) { return (a%2 == 0); }
bool CfgParser::IsOdd(int a) { return (a%2 == 1); }



