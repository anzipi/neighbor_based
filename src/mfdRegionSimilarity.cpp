
#include "mfdRegionSimilarity.h"
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#ifdef SUPPORT_OMP
#include <omp.h>
#endif
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

mfdRegionSimilarity::mfdRegionSimilarity(){
}

mfdRegionSimilarity::~mfdRegionSimilarity(){
}

int mfdRegionSimilarity::getNeighborSize(int rowSample, int colSample, double ** & dem, double ** & flag, 
	int totalRows, int totalCols, double noData, vector<int> & rNeighbor, vector<int> & cNeighbor){
	flag[rowSample][colSample] = 1;
	int neighborSize = 0;
	int cnt = 1;	
	vector<int> rBoundary;
	vector<int> cBoundary;
	while(cnt){			
		cnt = 0;		
		getBoundaryCord(neighborSize, rBoundary, cBoundary);
		for(int m = 0; m < rBoundary.size(); m++){
			int rowInterest = rowSample + rBoundary[m];
			int colInterest = colSample + cBoundary[m];
			if(rowInterest >= 0 && rowInterest < totalRows && colInterest >=0 && colInterest < totalCols){
				double interestValue = dem[rowInterest][colInterest];			
				if(abs(interestValue - noData) > VERY_SMALL && flag[rowInterest][colInterest] == 1){
					for(int n = 0; n < 8; n++){
						int r = rowInterest + rNeighbor[n];
						int c = colInterest + cNeighbor[n];
						if(r >= 0 && r < totalRows && c >=0 && c < totalCols){	
							if(abs(dem[r][c] - noData) > VERY_SMALL){
								if(dem[r][c] > interestValue && flag[r][c] != 1){
									flag[r][c] = 1;
									cnt++;
								}
							}
						}
					}
				}
			}
		}	
		if(cnt > 0){
			neighborSize++;	
			//cout << neighborSize << endl;
		}		
	}
	return neighborSize;
}


void mfdRegionSimilarity::getBoundaryCord(int neighborSize, vector<int> &rBoundary, vector<int> &cBoundary){
	if(neighborSize == 0){
		rBoundary.push_back(0);
		cBoundary.push_back(0);
	}else{
		//toprow
		for(int col = -neighborSize; col <= neighborSize; col++){
			rBoundary.push_back(-neighborSize);
			cBoundary.push_back(col);
		}
		//right col
		for(int row = -neighborSize + 1; row <= neighborSize - 1; row++){
			rBoundary.push_back(row);
			cBoundary.push_back(neighborSize);
		}
		//downrow
		for(int col = neighborSize; col >= -neighborSize; col--){
			rBoundary.push_back(neighborSize);
			cBoundary.push_back(col);
		}
		//leftrow
		for(int row = neighborSize - 1; row >= -neighborSize + 1; row--){
			rBoundary.push_back(row);
			cBoundary.push_back(-neighborSize);
		}
	}
	
	
}

void mfdRegionSimilarity::getNeighborEnv(int rowSample, int colSample, int totalRows, int totalCols, int neighborSize,
	double ** & dem, double ** & flag, double noData, vector<double> &neighborEnv){
	for(int row = -neighborSize; row <= neighborSize; row++){
		if(rowSample + row >= 0 && rowSample + row < totalRows){
			for(int col = -neighborSize; col <= neighborSize; col++){
				if(colSample + col >=0 && colSample + col < totalCols){
					if(flag[rowSample + row][colSample + col] == 1){						
						neighborEnv.push_back(dem[rowSample + row][colSample + col]);
					}
				}
			}
		}
	}
	//cout << rCord.size() << endl;
}

void mfdRegionSimilarity::frequencySta(double minEnv, double maxEnv, int number, vector<double> & neighborEnv, 
	double * &frequency){
	int * count = new int[number];
	for(int n = 0; n < number; n++){
		count[n] = 0;
	}
	for(int n = 0; n < neighborEnv.size(); n++){
		int index = int((neighborEnv[n] - minEnv) * number / (maxEnv - minEnv));		
		if(index == number){
			index = number - 1;
		}
		count[index] += 1;
	}
	for(int n = 0; n < number; n++){
		frequency[n] = count[n] * 1.0 / neighborEnv.size();		
	}
	delete [] count;
}

double mfdRegionSimilarity::histSimilarity(double * & frequencyPixel, double * & frequencySample, int number){
	double intersectValue = 0;
	double unionValue = 0;
	double similarityH;
	for(int m = 0; m < number; m++){
		intersectValue += min(frequencyPixel[m], frequencySample[m]);
		unionValue += max(frequencyPixel[m], frequencySample[m]);
	}
	//cout << intersectValue << ", " << unionValue << endl;
	if(unionValue == 0){
		similarityH = 0;
	}else{
		similarityH = intersectValue * 1.0 / unionValue;
	}
	//cout << similarityH << endl;
	return similarityH;
}
