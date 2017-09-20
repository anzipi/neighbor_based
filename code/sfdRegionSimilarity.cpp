#include <sstream>
#include "sfdRegionSimilarity.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

sfdRegionSimilarity::sfdRegionSimilarity(){
}

sfdRegionSimilarity::~sfdRegionSimilarity(){
}

//get the distance between the central pixel to the farthest pixel that flow into the central pixel directly or indirectly
// and sign the characteristic areas as 1 in the array "flag"
int sfdRegionSimilarity::getCharacterNeighbor(int rowSample, int colSample, double ** &dem, double ** &sfd,
	int totalRows, int totalCols, double noData, vector<int> & rNeighbor, vector<int> & cNeighbor,
	vector<int> & sfdNeighbor, double ** &flat){
	flat[rowSample][colSample] = 1;
	//初始信息
	int neighborSize = 0;
	int cnt = 0;
	int mark = 1;
	vector<int> rBoundary;
	vector<int> cBoundary;
	//根据流向确定邻域
	while(mark){		
		cnt = 0;
		getBoundaryCord(neighborSize, rBoundary, cBoundary);		
		for(int m = 0; m < rBoundary.size(); m++){
			//cout << neighborSize << ", " << rBoundary[m] << ", " << cBoundary[m] << endl;	
			int rowInterest = rowSample + rBoundary[m];
			int colInterest = colSample + cBoundary[m];
			if(rowInterest >= 0 && rowInterest < totalRows && colInterest >=0 && colInterest < totalCols){
				double interestValue = dem[rowInterest][colInterest];
				if(abs(interestValue - noData) > VERY_SMALL && flat[rowInterest][colInterest] == 1){								
					for(int n = 0; n < 8; n++){
						int r = rowInterest + rNeighbor[n];
						int c = colInterest + cNeighbor[n];
						if(r >= 0 && r < totalRows && c >=0 && c < totalCols){
							if(abs(sfd[r][c] - noData) > VERY_SMALL){							
								if(sfd[r][c] == sfdNeighbor[n] && flat[r][c] != 1){
									flat[r][c] = 1;
									//cout << n << endl;
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
			
		}else{
			mark = 0;
		}
		rBoundary.clear();
		cBoundary.clear();
	}
	//cout << neighborSize << endl;
	vector<int>().swap(rBoundary);
	vector<int>().swap(cBoundary);
	return neighborSize;
	//cout << "OK" << endl;
}

//collect the relative coordinates of the boundary pixels to the central pixel
void sfdRegionSimilarity::getBoundaryCord(int neighborSize, vector<int> &rBoundary, vector<int> &cBoundary){
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

//record the envValue of the characteristic neighborhood into a vector accroding to the flag array and neighborSize
void sfdRegionSimilarity::getNeighborEnv(int rowSample, int colSample, int totalRows, int totalCols, int neighborSize,
	double ** & env, double ** & flag, double noData, vector<double> &neighborEnv){//目前，没有提取邻域相对坐标
	for(int row = -neighborSize; row <= neighborSize; row++){
		if(rowSample + row >= 0 && rowSample + row < totalRows){
			for(int col = -neighborSize; col <= neighborSize; col++){
				if(colSample + col >= 0 && colSample + col < totalCols){
					if(flag[rowSample + row][colSample + col] == 1){						
						neighborEnv.push_back(env[rowSample + row][colSample + col]);
					}
				}
			}
		}
	}	
}

//statistic the frequency of envValue of each bin in the characteristic neighborhood 
void sfdRegionSimilarity::frequencySta(double minEnv, double maxEnv, int number, vector<double> & neighborEnv, 
	double * & frequency){
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

double sfdRegionSimilarity::histSimilarity(double * & frequencyPixel, double * & frequencySample, int number){
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

