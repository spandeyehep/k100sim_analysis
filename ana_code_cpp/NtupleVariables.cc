#define NtupleVariables_cxx
#include "NtupleVariables.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

bool NtupleVariables::IsUnique(vector<long> tracks, long thisTrack) {
	int occurance = 0;
	for(int i = 0; i < (int)tracks.size(); i++) {
		if(thisTrack == tracks.at(i)) occurance++;
		if(occurance > 1) return false;
	}
	return true;
}

vector<double> NtupleVariables::sort(vector<double> vec) {
	double temp;
	for(int i=1;i<vec.size();i++){
	    for(int j=i;j<vec.size();j++){
	      if((vec)[i-1] < (vec)[j] ){
			temp = (vec)[i-1];
			(vec)[i-1] = (vec)[j];
			(vec)[j] = temp;
	      }
	    }
	}
	return vec;
}