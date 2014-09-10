#include <vector>
#include <map>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "tools.h"

using namespace std;
#define NEITHER       0
#define UP            1
#define LEFT          2
#define UP_AND_LEFT   3

float SMC(vector<int>& a, vector<int>& b){
  float sim = 0;unsigned int len;
  if(a.size() < b.size()) len = a.size();
  else len = b.size();
  if(len == 0) return 0.0;
  for(unsigned int i =0; i < len; i++){
    if(a[i] == b[i]){
      sim += 1;
    }
  }
  return sim/len;
}

//this function assumes that a and b are of equal length else the longer one is truncated.
float WSMC(vector<int>& a,vector<int>& b,map<pair<int,int>,float>& alphaMap){
  unsigned int len;
  if(a.size() < b.size())  len = a.size();
  else len = b.size();
  float sim = 0.0;
  for(unsigned int i = 0; i < len; i++){
    if(alphaMap.find(pair<int,int>(a[i],b[i])) != alphaMap.end())
      sim += alphaMap[pair<int,int>(a[i],b[i])];
    else
      sim += 1;
  }
  return (float) sim/len;
}

float LCS_DP(vector<int> &a, vector<int> &b){
  int len1 = a.size(); int len2 = b.size();
  int i,j;
  int *prev_row = (int*)malloc((len2+1) * sizeof(int));
  int *cur_row = (int*)malloc((len2+1) * sizeof(int));
  int *temp_ptr;
  for(i=0;i<=len2;i++){
    prev_row[i] = 0;
    cur_row[i] = 0;
  }
  for(i=1;i<=len1;i++){
    for(j=1;j<=len2;j++){
      if(a[i-1]==b[j-1])
	cur_row[j] = prev_row[j-1]+1;
      else{
	if(cur_row[j-1] > prev_row[j])
	  cur_row[j] = cur_row[j-1];
	else
	  cur_row[j] = prev_row[j];
      }
    }
    temp_ptr = prev_row;
    prev_row = cur_row;
    cur_row = temp_ptr;
  }  
  return prev_row[len2]/(sqrt(len1*len2));
}

void set_array_locns(int *row2,int mask_value,int range_start,int range_end){
  int range = range_end - range_start + 1;
  int mask_size = range/2;
  int mask[mask_size];
  for(int x=0;x<mask_size;x++)  mask[x] = mask_value;
  int locn_written_till = range_start-1;
  while(locn_written_till < range_end - mask_size){
    memcpy(row2+locn_written_till+1,mask,mask_size *sizeof(int));
    locn_written_till = locn_written_till + mask_size;
  }
  memcpy(row2+locn_written_till+1,mask,(range_end - locn_written_till -1) *sizeof(int));
}

int binary_search_for_location_hy(int *row,int low,int high,int value_to_search){
  int mid;
  while(high-low > 1){
    mid = (high + low) / 2;
    if(row[mid] > value_to_search) high = mid;
    else low = mid;
  }
  if(row[low]!=value_to_search) low = -1;
  return low;
}

float LCS_HY(vector<int> &a,vector<int> &b,int max_symbol,int max_occ_any_symbol){
  int len1 = a.size();int len2 = b.size();
  int i,symbol,*row = (int*)malloc((len2+1) * sizeof(int)),*matches_list;
  for(i=0;i<=len2;i++)  row[i] = 0;
  //if max_symbol and max_occ_any_symbol has not been computed, compute them now
  if((max_symbol == -1) || (max_occ_any_symbol == -1)){
    map<int,int> alphaMap;
    for(int i = 0; i < len1; i++){
      if(alphaMap.find(a[i]) == alphaMap.end())	alphaMap[a[i]] = 1;
      else alphaMap[a[i]] += 1;
    }
    for(int i = 0; i < len2; i++){
      if(alphaMap.find(b[i]) == alphaMap.end())	alphaMap[b[i]] = 1;
      else alphaMap[b[i]] += 1;
    }
    for(map<int,int>::iterator it = alphaMap.begin(); it != alphaMap.end(); it++){
      if(it->second > max_occ_any_symbol) max_occ_any_symbol = it->second;
      if(it->first > max_symbol) max_symbol = it->first;
    }
  }
  int index[max_symbol][max_occ_any_symbol+2];
  // Initialize index.
  for (int sym_cnt=0;sym_cnt<max_symbol;sym_cnt++){
    index[sym_cnt][0] = 1;
    index[sym_cnt][1] = len2+1;
  }
  // Create index.
  int prev_written_locn;
  for (i=len2;i>=1;i--){
    symbol = b[i-1];
    prev_written_locn = index[symbol-1][0];
    index[symbol-1][prev_written_locn+1] = i;
    index[symbol-1][0] = prev_written_locn+1;
  }
  // Update rows to get lcs.
  int matches_num,match_locn;
  for (i=1;i<=len1;i++){
    symbol = a[i-1];
    matches_list = &(index[symbol-1][0]);
    matches_num = matches_list[0];
    if (matches_num==1) continue;
    int new_value,value_to_search,match_point,start_point,end_point;
    for(int matches_cnt=2;matches_cnt<=matches_num;matches_cnt++){
      match_locn = matches_list[matches_cnt];
      new_value = row[match_locn-1] + 1;
      row[match_locn] = new_value;
      if (new_value > row[match_locn+1]){
	value_to_search = new_value;
	start_point = match_locn + 1;
	end_point = matches_list[matches_cnt-1];
	if(new_value > row[end_point-1]) 
	  match_point = end_point;
	else{
	  match_point = binary_search_for_location_hy(row,start_point,end_point,value_to_search);
	  if(match_point==-1) match_point = end_point;
	}
	int range = match_point - match_locn - 1;
	if (range <=4000){
	  for(int tmp=match_locn+1;tmp<=match_point-1;tmp++)
	    row[tmp] = new_value;
	}else{
	  set_array_locns(row,new_value,match_locn+1,match_point-1);
	}	
      }
    }
  }
  return row[len2]/sqrt(len1*len2);
}

/* This function converts two sequences into bitmaps and compares
   them to compute a distance between the two sequences. Originally 
   proposed in Kumar et al 2005 and used for anomaly detection in
   Wei et al 2005.
*/
float BITMAP(vector<int> &a,vector<int> &b, int level){
  map<vector<int>, float> ma;
  float maxcounta = 0;
  //slide a level length window along a and store counts of windows
  for(vector<int>::iterator it = a.begin(); it != a.end()-level;it += 1){
      vector<int> win(it,it+(level-1));
      if(ma.find(win) == ma.end()) ma[win] = 1;
      else ma[win] = ma[win]+1;      
      if(ma[win] > maxcounta) maxcounta = ma[win];
  }
  //normalize the counts by the maximum value
  for(map<vector<int>, float>::iterator it = ma.begin(); it != ma.end(); it++)
    it->second = it->second/maxcounta;
  map<vector<int>, float> mb;
  float maxcountb = 0;
  //slide a level length window along b and store counts of windows
  for(vector<int>::iterator it = b.begin(); it != b.end()-level;it += 1){
      vector<int> win(it,it+(level-1));
      if(mb.find(win) == mb.end()) mb[win] = 1;
      else mb[win] = mb[win]+1;      
      if(mb[win] > maxcountb) maxcountb = mb[win];
  }
  //normalize the counts by the maximum value
  for(map<vector<int>, float>::iterator it = mb.begin(); it != mb.end(); it++)
    it->second = it->second/maxcountb;
  //compute distance between the two bitmaps
  float dist = 0;
  for(map<vector<int>, float>::iterator it = ma.begin(); it != ma.end(); it++){
    if(mb.find(it->first) != mb.end())
      dist += (it->second - mb[it->first])* (it->second - mb[it->first]);
    else
      dist += it->second*it->second;
  }
  for(map<vector<int>, float>::iterator it = mb.begin(); it != mb.end(); it++){
    if(ma.find(it->first) == ma.end())
      dist += it->second*it->second;
  }
  return dist;
}
