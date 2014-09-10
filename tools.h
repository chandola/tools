#include <fstream>
#include <vector>
#include <map>
#include <string>
#include "fftw3.h"

using namespace std;
/* I/O functions for discrete sequences and time series */
void Tokenize(const string& str,vector<string> &tokens,const string& delimiters = " ");
void Tokenize(const string& str,vector<int> &tokens,const string& delimiters = " ");
void Tokenize(const string& str,vector<float> &tokens,const string& delimiters = " ");

int readSequences(vector<vector<int> > &sequences,ifstream &in);
int readSequences(vector<vector<float> > &sequences,ifstream &in);

void printSequences(vector<vector<int> > &pairwise, ofstream &out);
void printSequences(vector<vector<float> > &pairwise, ofstream &out);

void printSequence(vector<int> &sequence, ofstream &out);
void printSequence(vector<float> &sequence, ofstream &out);

void loadMap(map<pair<int,int>,float > &alphaMap, ifstream &in);
void printMap(map<pair<int,int>,float > &alphaMap, ofstream &out);
void findMax(vector<vector<int> > &sequences, int *max_symbol, int* max_occ_any_symbol);

void getWindows(vector<vector<int> > &sequences,vector<vector<int> > &windows, int len);
void getWindows(vector<vector<float> > &sequences,vector<vector<float> > &windows, int len);

void getWindows_hop(vector<vector<int> > &sequences,vector<vector<int> > &windows, int len, int hop,vector<unsigned int> &lens);
void getWindows_hop(vector<vector<float> > &sequences,vector<vector<float> > &windows, int len, int hop,vector<unsigned int> &lens);

bool FileExists(string strFilename);
#include <vector>

/* Functions to compute distance between time series */
/* DIST - Wrapper function for calling different distance measures
   1 - Euclidean
   2 - Dynamic Time Warp (DTW)
   3 - DTW with Sakuro Chiba Band
   4 - LB_KEOGH (Keogh, E. (2002).  Exact indexing of dynamic time warping)
*/
float DIST(vector<float> &a,vector<float> &b,int measure);

float EUC(vector<float> &a,vector<float> &b);
float DTW(vector<float> &a,vector<float> &b);
float DTW_SC(vector<float> &a,vector<float> &b);
float LB_KEOGH(vector<float> &a, vector<float> &b);

// Similarity between time series using cross correlation
float CROSSCORR(vector<float> &a, vector<float> &b,fftw_plan plan_forward, fftw_plan plan_backward, double *fin, fftw_complex *fout,fftw_complex *bin, double *bout);

/* pairwiseDist - Calculate pairwise distance between two sets of time series
 */
void pairwiseDist(vector<vector<float> > &sequences1,vector<vector<float> > &sequences2,vector<vector<float> > &pairwise, int msr);

/* queryNN - Finds the distance of query time series to the nn^{th} nearest neighbor in DB using the LBKeogh lower bound as proposed in Ratanamahantana and Keogh 2004. Returns the distance to the nn^th nearest neighbor and index stores its id in DB.
*/
float queryNN(vector<float> &query, vector<vector<float> > &DB, int nn, unsigned int *index);

/* Functions to compute distance/similarity for discrete (symbolic) sequences*/
//1. Simple Matching Coefficient
float SMC(vector<int> &a, vector<int> &b);
//2. Weighted Simple Matching Coefficient
float WSMC(vector<int> &a,vector<int> &b,map<pair<int,int>,float>&alphaMap);
//3. Normalized Longest Common Subsequence using Dynamic Programming (O(mn))
float LCS_DP(vector<int> &a, vector<int> &b);
//4. Normalized Longest Common Subsequence using Hybrid of Hunt Szymanski and Dynamic Programming (Budalakoti et al, 2007).
float LCS_HY(vector<int> &a, vector<int> &b,int max_symbol=-1,int max_occ_any_symbol=-1);
//5. Distance between bitmaps (Kumar et al 2005)
float BITMAP(vector<int> &a,vector<int> &b, int level);
