#ifndef MFDREGIONSIMILARITY_H
#define MFDREGIONSIMILARITY_H
//#include <omp.h>
#include "AscGrid.h"
#define VERY_SMALL 0.000001
using namespace std;

class mfdRegionSimilarity{

private:
	

public:
	mfdRegionSimilarity();
	~mfdRegionSimilarity();
	int getNeighborSize(int rowSample, int colSample, double ** &dem, double ** &flag, int totalRows, 
		int totalCols, double noData, vector<int> & rNeighbor, vector<int> & cNeighbor);
	void getBoundaryCord(int neighborSize, vector<int> &rBoundary, vector<int> &cBoundary);
	void getNeighborEnv(int rowSample, int colSample, int totalRows, int totalCols, int neighborSize,
		double ** & dem, double ** & flag, double noData, vector<double> & neighborEnv);
	void frequencySta(double minEnv, double maxEnv, int number, vector<double> &neighborEnv, 
		double * &frequency);
	double histSimilarity(double * & frequencyPixel, double * & frequencySample, int number);
};

#endif
