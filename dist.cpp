#include <vector>
#include <iostream>
#include <queue>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "tools.h"

#define SC_COEFF 0.10
using namespace std;

float minimum(float a,float b,float c){
  float min = a;
  if(min > b) min = b;
  if(min > c) min = c;
  return min;
}

float DIST(vector<float> &a,vector<float> &b,int measure){
  switch (measure){
  case 1 : 
    return EUC(a,b);
    break;
  case 2:
    return DTW(a,b);
    break;
  case 3:
    return DTW_SC(a,b);
    break;
  case 4:
    return LB_KEOGH(a,b);
    break;
  default:
    cerr<<"Unsupported measure\n";
    return -1;
    break;
  }
  return 0;
}

/* This function returns the square of euclidean distance */
float EUC(vector<float> &a,vector<float> &b){
  if(a.size() != b.size()){
    fprintf(stderr,"Error : Sequences should be of same length.\n");
    return 0;
  }
  float dist = 0.0;
  for(unsigned int i = 0; i < a.size(); i++){
    dist += (a[i]-b[i])*(a[i]-b[i]);
  }
  return sqrt(dist);
}

/* This function returns the DTW distance */
float DTW(vector<float> &a,vector<float> &b){
  vector<vector<float> > dtw(a.size()+1,vector<float>(b.size()+1,0.0));
  for(unsigned int i = 1; i < b.size()+1; i++)
    dtw[0][i] = INT_MAX;
  for(unsigned int i = 1; i < a.size()+1; i++)
    dtw[i][0] = INT_MAX;
  dtw[0][0] = 0;
  for(unsigned int i = 1; i < a.size()+1; i++){
    for(unsigned int j = 1; j < b.size()+1; j++){
      float cost = (a[i-1]-b[j-1])*(a[i-1]-b[j-1]);
      dtw[i][j] = minimum(dtw[i-1][j]+cost,dtw[i][j-1]+cost,dtw[i-1][j-1]+cost);
    }
  }
  float ret = dtw[a.size()][b.size()];
  dtw.clear();
  return ret;
}

/* This function returns the DTW (with Sakoe and Chiba constraint) distance.
   The warp path is constrained to be within 10% of the main diagonal.
*/
float DTW_SC(vector<float> &a,vector<float> &b){
  int minlen = (unsigned int) a.size();
  if(a.size() > b.size()) minlen = (unsigned int) b.size();
  int constraint = (int) floor(SC_COEFF*minlen);
  vector<vector<float> > dtw(a.size()+1,vector<float>(b.size()+1,0.0));
  for(unsigned int i = 1; i < b.size()+1; i++)
    dtw[0][i] = INT_MAX;
  for(unsigned int i = 1; i < a.size()+1; i++)
    dtw[i][0] = INT_MAX;
  dtw[0][0] = 0;
  for(unsigned int i = 1; i < a.size()+1; i++){
    for(unsigned int j = 1; j < b.size()+1; j++){
      float cost = INT_MAX;
      if(abs((int)(i-j)) <= constraint) cost = (a[i-1]-b[j-1])*(a[i-1]-b[j-1]);
      dtw[i][j] = minimum(dtw[i-1][j]+cost,dtw[i][j-1]+cost,dtw[i-1][j-1]+cost);
    }
  }
  return dtw[a.size()][b.size()];
}


float findMin(vector<float> &vec, int start, int end){
  if(start < 0) start = 0;if(end > (int) vec.size()) end = (int) vec.size();
  float val = vec[start];
  for(int i = start; i < end; i++){
    if(vec[i] < val) val = vec[i];
  }
  return val;
}

float findMax(vector<float> &vec, int start, int end){
  if(start < 0) start = 0;if(end > (int) vec.size()) end = (int) vec.size();
  float val = vec[start];
  for(int i = start; i < end; i++){
    if(vec[i] > val) val = vec[i];
  }
  return val;
}

float LB_KEOGH(vector<float> &a, vector<float> &b){
  if(a.size() != b.size()){
    cerr<<"This measure supported for equal length series only\n";
    return 0.0;
  }
  int constraint = (int) floor(SC_COEFF*a.size());
  float lb_dist = 0;
  for(int i = 0; i < (int) a.size(); i++){
    float Li = findMin(a,i-constraint,i+constraint);
    float Ui = findMax(a,i-constraint,i+constraint);
    if(b[i] > Ui) lb_dist += pow((b[i]-Ui),2);
    if(b[i] < Li) lb_dist += pow((b[i]-Li),2);
  }
  return lb_dist;
}

void pairwiseDist(vector<vector<float> > &sequences1,vector<vector<float> > &sequences2,vector<vector<float> > &pairwise, int msr){
  for(unsigned int i = 0; i < sequences1.size();i++){
    vector<float> tmp(sequences2.size(),0.0);
    pairwise.push_back(tmp);
  }
  for(unsigned int i = 0; i < sequences1.size();i++){
    for(unsigned int j = 0; j < sequences2.size();j++){
      if(j < i){
	if(i < sequences2.size())
	  pairwise[i][j] = pairwise[j][i];
	else
	  pairwise[i][j] = DIST(sequences1[i],sequences2[j],msr);
      }else{
	pairwise[i][j] = DIST(sequences1[i],sequences2[j],msr);
      }
    }
  }
}

/* Finds the distance of query time series to the nn^{th} nearest neighbor in DB using the LBKeogh lower bound as proposed in Ratanamahantana and Keogh 2004.Returnsthe distance to the nn^th nearest neighbor and index stores its id in DB.
*/

class compare_A{
public:
  bool operator() (pair<float,unsigned int> a, pair<float,unsigned int>  b){
    return a.first > b.first;
  }
};

float queryNN(vector<float> &query, vector<vector<float> > &DB, int nn, unsigned int *index){
  float true_dist;
  priority_queue<pair<float,unsigned int>,vector<pair<float,unsigned int> >, compare_A > q;
  for(unsigned int i = 0; i < DB.size(); i++){
    if(q.size() < (unsigned int) nn){
      true_dist = DTW_SC(query,DB[i]);
      q.push(pair<float,unsigned int>(true_dist,i));
    }else{
      float LB_dist = LB_KEOGH(query,DB[i]);
      if(LB_dist < q.top().first){
	true_dist = LB_KEOGH(query,DB[i]);
	if(true_dist < q.top().first){
	  q.pop();
	  q.push(pair<float,unsigned int>(true_dist,i));
	}
      }
    }
  }
  *index = q.top().second;
  return q.top().first;
}
