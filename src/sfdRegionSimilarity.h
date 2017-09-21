#ifndef SFDREGIONSIMILARITY_H
#define SFDREGIONSIMILARITY_H
#include <omp.h>
#include "AscGrid.h"
#define VERY_SMALL 0.000001
using namespace std;

class sfdRegionSimilarity{

private:
	

public:
	sfdRegionSimilarity();
	~sfdRegionSimilarity();
	int getCharacterNeighbor(int rowSample, int colSample, double ** &dem, double ** &sfd, 
		int totalRows, int totalCols, double noData, vector<int> & rNeighbor, 
		vector<int> & cNeighbor, vector<int> & sfdNeighbor, double ** &flat);
	void getBoundaryCord(int neighborSize, vector<int> &rBoundary, vector<int> &cBoundary);
	void getNeighborEnv(int rowSample, int colSample, int totalRows, int totalCols, int neighborSize,
		double ** & env, double ** & flag, double noData, vector<double> & neighborEnv);
	void frequencySta(double minEnv, double maxEnv, int number, vector<double> & neighborEnv, 
		double * &frequency);
	double histSimilarity(double * & frequencyPixel, double * & frequencySample, int number);		
};

#endif
