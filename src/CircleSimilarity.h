#ifndef CIRCLESIMILARITY_H
#define CIRCLESIMILARITY_H
#include <omp.h>
#include "AscGrid.h"
#define VERY_SMALL 0.000001
using namespace std;

class CircleSimilarity{

private:
	double string_to_double( const std::string& s );

public:
	CircleSimilarity();
	~CircleSimilarity();	
	int getCharacteristicScale(int rowSample, int colSample, double envSample, double sampleVar, int totalRows, 
		int totalCols,	double ** &env, double noData, double std, int radiusMin, int radiusMax);		
	double getCircleSimilarity(int rowPixel, int colPixel, int rowSample, int colSample , double ** &env, 
		int totalRows, int totalCols, double Std, double varSample, double noData, int scaleIntegration, 
		double alpha);
	void getAnnulusCord(int radius, vector<int> & rCord, vector<int> & cCord);
	void getWeight(double * &weight, int neighborSize, double alpha);
	void arrayMoveRight(double *array, double * arrayMoved, int count);
	double getVectSimilarity(double *array1, double * array2, int length);
	double getMaxValue(double * array, int count);
};

#endif
