#include <fstream>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include "tools.h"

using namespace std;

void Tokenize(const string& str, vector<string>& tokens,const string& delimiters){
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos){
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
}

void Tokenize(const string& str, vector<int>& tokens,const string& delimiters){
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos){
      string str1 = str.substr(lastPos, pos - lastPos);
      tokens.push_back(atoi(str1.c_str()));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
}

void Tokenize(const string& str, vector<float>& tokens,const string& delimiters){
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos){
      string str1 = str.substr(lastPos, pos - lastPos);
      tokens.push_back(atof(str1.c_str()));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
}

int readSequences(vector<vector<int> > &sequences,ifstream &in){
  int lines = 0;
  string bufstr;
  getline(in,bufstr);
  while(!in.fail()){
    if((bufstr.length() == 0) || (bufstr[0] == '#')){
      getline(in,bufstr);
      continue;
    }
    vector<int> vec(0);
    Tokenize(bufstr,vec);
    sequences.push_back(vec);
    lines++;
    getline(in,bufstr);
  }
  return lines;
}

int readSequences(vector<vector<float> > &sequences,ifstream &in){
  int lines = 0;
  string bufstr;
  getline(in,bufstr);
  while(!in.fail()){
    if((bufstr.length() == 0) || (bufstr[0] == '#')){
      getline(in,bufstr);
      continue;
    }
    vector<float> vec(0);
    Tokenize(bufstr,vec);
    sequences.push_back(vec);
    lines++;
    getline(in,bufstr);
  }
  return lines;
}

void printSequence(vector<int> &sequence, ofstream &out){
  for(unsigned int j = 0; j < sequence.size()-1; j++) out<<sequence[j]<<" ";
  out<<sequence[sequence.size()-1];
}

void printSequence(vector<float> &sequence, ofstream &out){
  for(unsigned int j = 0; j < sequence.size()-1; j++)  out<<sequence[j]<<" ";
  out<<sequence[sequence.size()-1];
}

void printSequences(vector<vector<int> > &sequences, ofstream &out){
  for(unsigned int i = 0; i < sequences.size(); i++){printSequence(sequences[i],out);out<<endl;}
}

void printSequences(vector<vector<float> > &sequences, ofstream &out){
  for(unsigned int i = 0; i < sequences.size(); i++){ printSequence(sequences[i],out);out<<endl;}
}

void loadMap(map<pair<int,int>,float > &alphaMap, ifstream &in){
  string bufstr;
  getline(in,bufstr);
  while(!in.fail()){
    if((bufstr.length() == 0) || (bufstr[0] == '#')){
      getline(in,bufstr);
      continue;
    }
    //tokenize this string
    vector<float> tokens;
    Tokenize(bufstr, tokens,",");
    int a = (int) tokens[0];
    int b = (int) tokens[1];
    float val = tokens[2];
    alphaMap[pair<int,int>(a,b)] = val;
    getline(in,bufstr);
  }
}

void printMap(map<pair<int,int>,float > &alphaMap, ofstream &out){
  for(map<pair<int,int>,float >::iterator it=alphaMap.begin(); it != alphaMap.end(); it++)
    out<<(it->first).first<<","<<(it->first).second<<","<<it->second<<"\n";  
}

//Calculates the maximum value occuring in a vector of vector of integer sequences and stores in max_symbol and calcuates the frequency of the maximum occurring integer value and stores it in max_occ_any_symbol. 
//Note, if max_symbol and max_occ_any_symbol, already contain values, they are not changed if the new values are lesser.
void findMax(vector<vector<int> > &sequences, int *max_symbol, int* max_occ_any_symbol){
  map<int,int> alphaMap;
  for(unsigned int i = 0; i < sequences.size(); i++){
    for(unsigned int j = 0; j < sequences[i].size(); j++){
      if(alphaMap.find(sequences[i][j]) == alphaMap.end()) alphaMap[sequences[i][j]] = 1;
      else alphaMap[sequences[i][j]] += 1;
    }
  }
  for(map<int,int>::iterator it = alphaMap.begin(); it != alphaMap.end(); it++){
    if(it->first > *max_symbol) *max_symbol = it->first;
    if(it->second > *max_occ_any_symbol) *max_occ_any_symbol = it->second;
  }  
}

void getWindows(vector<vector<int> > &sequences,vector<vector<int> > &windows, int len){
  windows.resize(0);
  for(unsigned int i = 0; i < sequences.size(); i++){
    if(sequences[i].size() < (unsigned int) len){
      windows.push_back(sequences[i]);
    }else{
      for(unsigned int j = 0; j < sequences[i].size()-len+1;j++){
	vector<int>::iterator itbegin = sequences[i].begin() + j;
	vector<int>::iterator itend = sequences[i].begin()+ j + len;
	vector<int> tmp(itbegin,itend);
	windows.push_back(vector<int> (itbegin,itend));
      }
    }
  }
}

void getWindows(vector<vector<float> > &sequences,vector<vector<float> > &windows, int len){
  windows.resize(0);
  for(unsigned int i = 0; i < sequences.size(); i++){
    if(sequences[i].size() < (unsigned int) len){
      windows.push_back(sequences[i]);
    }else{
      for(unsigned int j = 0; j < sequences[i].size()-len+1;j++){
	vector<float>::iterator itbegin = sequences[i].begin() + j;
	vector<float>::iterator itend = sequences[i].begin()+ j + len;
	vector<float> tmp(itbegin,itend);
	windows.push_back(vector<float> (itbegin,itend));
      }
    }
  }
}

void getWindows_hop(vector<vector<int> > &sequences,vector<vector<int> > &windows, int len, int hop,vector<unsigned int> &lens){
  windows.resize(0);
  for(unsigned int i = 0; i < sequences.size(); i++){
    if(sequences[i].size() < (unsigned int) len){
      windows.push_back(sequences[i]);
      lens.push_back(1);
    }else{
      unsigned int n = 0;
      for(unsigned int j = 0; j < sequences[i].size()-len+1;j+=hop){
	vector<int>::iterator itbegin = sequences[i].begin() + j;
	vector<int>::iterator itend = sequences[i].begin()+ j + len;
	vector<int> tmp(itbegin,itend);
	windows.push_back(vector<int> (itbegin,itend));
	n++;
      }
      lens.push_back(n);
    }
  }
}

void getWindows_hop(vector<vector<float> > &sequences,vector<vector<float> > &windows, int len, int hop,vector<unsigned int> &lens){
  windows.resize(0);
  for(unsigned int i = 0; i < sequences.size(); i++){
    if(sequences[i].size() < (unsigned int) len){
      windows.push_back(sequences[i]);
      lens.push_back(1);
    }else{
      unsigned int n = 0;
      for(unsigned int j = 0; j < sequences[i].size()-len+1;j+=hop){
	vector<float>::iterator itbegin = sequences[i].begin() + j;
	vector<float>::iterator itend = sequences[i].begin()+ j + len;
	vector<float> tmp(itbegin,itend);
	windows.push_back(vector<float> (itbegin,itend));
	n++;
      }
      lens.push_back(n);
    }
  }
}

bool FileExists(string strFilename) {
  struct stat stFileInfo;
  if(stat(strFilename.c_str(),&stFileInfo) == 0) return true;
  else return false;
}
