#include <vector>
#include <iostream>
#include "tools.h"
#include "fftw3.h"

/* internal function required for CROSSCORR */
void fftw_copy(fftw_complex *orig, fftw_complex *copy, int N){
  for(int j =0; j <= N/2; j++){
    copy[j][0] = orig[j][0];
    copy[j][1] = orig[j][1];
  }
  for(int j = N/2+1; j < N; j++){
    copy[j][0] = orig[N-j][0];
    copy[j][1] = -1*orig[N-j][1];
  }
}

/*
 This function computes the similarity as the maximum correlation between vector a and phase shifted version of vector b.
 The function assumes that the forward and backward plans for FFT have already been created using fin, fout and bin, bout, respectively.
 All arrays needs to be pre allocated with length equal to the length of vector a and b.  
*/

float CROSSCORR(std::vector<float> &a, std::vector<float> &b,fftw_plan plan_forward, fftw_plan plan_backward, double *fin, fftw_complex *fout,fftw_complex *bin, double *bout){
  if(a.size() != b.size()){
    std::cerr<<"The vectors a and b should of equal lengths\n";
    return 0.0;
  }
  int N = (int) a.size();
  fftw_complex *fta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftw_complex *ftb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  for(int j =0; j < N; j++){fin[j] = a[j];}
  fftw_execute(plan_forward);
  fftw_copy(fout,fta,N);
  for(int j =0; j < N; j++){fin[j] = b[j];}
  fftw_execute(plan_forward);
  fftw_copy(fout,ftb,N);
  //multiply fta and conj(ftb)
  for(int j = 0; j < N; j++){
    bin[j][0] = fta[j][0]*ftb[j][0] + fta[j][1]*ftb[j][1];
    bin[j][1] = -1*fta[j][0]*ftb[j][1] + fta[j][1]*ftb[j][0];
  }
  fftw_free(fta);fftw_free(ftb);
  fftw_execute(plan_backward);
  double max=0.0;
  for(int j = 0; j < N; j++)
    if(bout[j] > max) max = bout[j];
  return (float) max/N;
}
