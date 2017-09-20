#include <sstream>
#include "CircleSimilarity.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

CircleSimilarity::CircleSimilarity(){
}

CircleSimilarity::~CircleSimilarity(){
}


int CircleSimilarity:: getCharacteristicScale(int rowSample, int colSample, double envSample, double sampleVar, int totalRows, 
	int totalCols,	double ** &env, double noData, double std, int radiusMin, int radiusMax){
	int flag = 1;
	int radius = 1;
	int indexAnnulus = 0;
	int indexDif = 0;
	vector<int> rCord;
	vector<int> cCord;
	double denominator = 2 * pow(std, 4) / sampleVar;
	int neighborSize = 0;
	double * annulusSimilarity = new double[radiusMax];//环形与中心栅格的相似度（环形内栅格与中心栅格相似度的平均值）
	double * similarityDif = new double[radiusMax - 1];//记录两相邻环形相似性的差
	while(flag == 1 && radius < radiusMax){
		getAnnulusCord(radius, rCord, cCord);
		double sumPointSimilarity = 0;
		int cnt = 0;
		for(int n = 0; n < rCord.size(); n++){
			int rowNeighbor = rowSample + rCord[n];
			if(rowNeighbor >= 0 && rowNeighbor < totalRows){
				int colNeighbor = colSample + cCord[n];
				if(colNeighbor >= 0 && colNeighbor < totalCols){
					double envNeighbor = env[rowNeighbor][colNeighbor];
					if(abs(envNeighbor - noData) > VERY_SMALL){						
						double numerator = pow((envNeighbor - envSample), 2);
						double pointSimilarity = exp(- numerator / denominator);
						sumPointSimilarity += pointSimilarity;
						cnt++;
					}
				}
			}			
		}
		rCord.clear();
		cCord.clear();
		if(cnt > 0){
			//cout << radius << ", " << sumPointSimilarity << ", " << cnt << endl;
			annulusSimilarity[indexAnnulus] = sumPointSimilarity / cnt;
			//cout << indexAnnulus << ", " << annulusSimilarity[indexAnnulus] << endl;
			if(indexAnnulus >= 1){
				indexDif = indexAnnulus - 1;
				similarityDif[indexDif] = annulusSimilarity[indexAnnulus - 1] - annulusSimilarity[indexAnnulus];				
				if(indexDif >= 2){
					if(similarityDif[indexDif - 1] > similarityDif[indexDif - 2] && similarityDif[indexDif - 1] > similarityDif[indexDif]){
						neighborSize = indexDif;						
						flag = 0;
					}
				}				
			}
			indexAnnulus++;	
			radius++;				
		}else{
			flag = 0;
		}
		
	}
	vector<int>().swap(rCord);
	vector<int>().swap(cCord);
	if(radius == radiusMax){
		double envMean;
		double envSum;
		double varSum;
		double difMean;
		int cnt = 0;
		for(int r = -radiusMax; r <= radiusMax; r++){
			int rowNeighbor = rowSample + r;
			if(rowNeighbor >= 0 && rowNeighbor < totalRows){
				for(int c = -radiusMax; c <= radiusMax; c++){
					int colNeighbor = colSample + c;
					if(colNeighbor >= 0 && colNeighbor < totalCols){
						double envNeighbor = env[rowNeighbor][colNeighbor];
						if(abs(envNeighbor - noData) > VERY_SMALL){	
							envSum += envNeighbor;
							cnt++;
						}
					}
				}
			}				
		}
		envMean = envSum / (cnt - 1);
		for(int r = -radiusMax; r <= radiusMax; r++){
			int rowNeighbor = rowSample + r;
			if(rowNeighbor >= 0 && rowNeighbor < totalRows){
				for(int c = -radiusMax; c <= radiusMax; c++){
					int colNeighbor = colSample + c;
					if(colNeighbor >= 0 && colNeighbor < totalCols){
						double envNeighbor = env[rowNeighbor][colNeighbor];
						if(abs(envNeighbor - noData) > VERY_SMALL){	
							varSum += pow(envNeighbor - envMean, 2);
						}
					}
				}
			}				
		}
		difMean = sqrt(varSum / (cnt - 1));
		if(difMean >= std){
			neighborSize = radiusMin;
		}else{
			neighborSize = radiusMax;
		}		
	}
	delete [] annulusSimilarity;
	delete [] similarityDif;
	return neighborSize;
}


double CircleSimilarity::getCircleSimilarity(int rowPixel, int colPixel, int rowSample, int colSample , double ** &env, 
	int totalRows, int totalCols, double Std, double varSample, double noData, int scaleIntegration, double alpha){
	double envPixel = env[rowPixel][colPixel];
	double envSample = env[rowSample][colSample];		
	double circleSimilarity = 0;
	
	double denominator = 2 * pow(Std,4) / varSample;
	double numerator = pow((envSample - envPixel),2);
	double * annulusSimilarity = new double [scaleIntegration + 1];
	annulusSimilarity[0] = exp(-numerator / denominator);
	vector<int> rCord;
	vector<int> cCord;
	for(int radius = 1; radius <= scaleIntegration; radius++){
		getAnnulusCord(radius, rCord, cCord);		
		double * envVectPixel = new double [rCord.size()];
		double * envVectSample = new double [rCord.size()];
		int count = 0;
		for(int n = 0; n < rCord.size(); n++){
			int rowNeighborPixel = rowPixel + rCord[n];
			int colNeighborPixel = colPixel + cCord[n];
			int rowNeighborSample = rowSample + rCord[n];
			int colNeighborSample = colSample + cCord[n];	
			if(rowNeighborPixel >= 0 && rowNeighborPixel < totalRows && rowNeighborSample >= 0 && rowNeighborSample < totalRows){
				if(colNeighborPixel >= 0 && colNeighborPixel < totalCols && colNeighborSample >= 0 && colNeighborSample < totalCols){
					double envNeighborPixel = env[rowNeighborPixel][colNeighborPixel];
					double envNeighborSample = env[rowNeighborSample][colNeighborSample];
					if(abs(envNeighborPixel - noData) > VERY_SMALL && abs(envNeighborSample - noData) > VERY_SMALL){
						envVectPixel[count] = envNeighborPixel;
						envVectSample[count] = envNeighborSample; 
						count ++;
					}
				}
			}								
		}
		rCord.clear();
		cCord.clear();		
		//用于存放上面得到的邻域向量的旋转结果			
		double * envVectPixelEnd = new double [count];
		double ** envVectSampleEnd  = new double *[count];
		double * vectSimilarity = new double[count];	
		for(int m = 0; m < count; m ++){
			envVectSampleEnd[m] = new double[count];
		}
		for(int m = 0; m < count; m ++){
			envVectPixelEnd[m] = envVectPixel[m];
			envVectSampleEnd[0][m] = envVectSample[m];				
		}

		for(int m = 1; m < count; m ++){
			arrayMoveRight(envVectSampleEnd[m - 1], envVectSampleEnd[m], count);
		}

		for(int m = 0; m < count; m ++){
			vectSimilarity[m] = getVectSimilarity(envVectPixelEnd, envVectSampleEnd[m], count);
		} 
		annulusSimilarity[radius] =  getMaxValue(vectSimilarity, count);
		for (int m = 0; m < count; m++){
			delete []  envVectSampleEnd[m];
		}
		delete[] envVectSampleEnd;
		delete[] envVectPixelEnd;
		delete[] envVectPixel;
		delete[] envVectSample;
		delete[] vectSimilarity;		
	}
	
	double * weight = new double[scaleIntegration + 1];
	getWeight(weight, scaleIntegration, alpha);
	for(int radius = 0; radius <= scaleIntegration; radius++){
		//cout << annulusSimilarity[radius] << ", " << weight[radius] << endl;
		circleSimilarity += annulusSimilarity[radius] * weight[radius];		
	}
	delete [] annulusSimilarity;
	delete [] weight;
	return circleSimilarity;	
}

void CircleSimilarity :: getAnnulusCord(int radius, vector<int> & rCord, vector<int> & cCord){
	int outer2 = radius * radius;//外环半径平方
	int inner2 = (radius - 1) * (radius - 1);//内环半径平方
	vector <int> rDist;//X方向上与中心栅格距离的绝对值
	vector <int> cDist;
	vector <int> rDistDescent;//rDist从小到大排序
	vector <int> cDistDescent;
	
	for(int r = radius; r >= 0; r--){
		for(int c = 0; c <= radius; c++){
			if( r * r + c * c  <= outer2 && r * r + c * c > inner2){
				rDist.push_back(r);
				cDist.push_back(c);
			}
		}
	}
	for(int v = 0; v < rDist.size(); v++){
		int r = rDist[rDist.size() - v - 1];
		int c = cDist[cDist.size() - v - 1];
		rDistDescent.push_back(r);
		cDistDescent.push_back(c);
	}
	for(int v = 0; v < rDist.size(); v++){
		int r = -rDist[v];
		int c = cDist[v];
		rCord.push_back(r);
		cCord.push_back(c);
	}
	for(int v = 0; v < rDist.size(); v++){
		int r = rDistDescent[v];
		int c = cDistDescent[v];
		if(r != rCord[rCord.size() - 1] || c != cCord[cCord.size() - 1]){
			rCord.push_back(r);
			cCord.push_back(c);
		}
	}
	for(int v = 0; v < rDist.size(); v++){
		int r = rDist[v];
		int c = -cDist[v];
		if(r != rCord[rCord.size() - 1] || c != cCord[cCord.size() - 1]){
			rCord.push_back(r);
			cCord.push_back(c);
		}
	}			
	for(int v = 0; v < rDist.size(); v++){
		int r = -rDistDescent[v];
		int c = -cDistDescent[v];
		if((r != rCord[rCord.size() - 1] || c != cCord[cCord.size() - 1]) && (r != rCord[0] || c != cCord[0])){
			rCord.push_back(r);
			cCord.push_back(c);
		}
	}
	vector<int>().swap(rDist);
	vector<int>().swap(cDist);
	vector<int>().swap(rDistDescent);
	vector<int>().swap(cDistDescent);

}

void CircleSimilarity:: getWeight(double * &weight, int neighborSize, double alpha){
	int length = neighborSize + 1;
	double * distance = new double[length];
	double sum = 0;	
	for(int d = 0; d < length; d++){
		distance[d] = pow(d + 0.5, alpha);	
		sum += 1 / distance[d];
	}
	for(int d = 0; d < length; d++){
		weight[d] = 1 / distance[d] / sum;
	}
	delete[] distance;
}

 void CircleSimilarity::arrayMoveRight(double * array, double * arrayMoved, int count){
	for(int l = 0; l < count; l++){
		arrayMoved[l] = array[l + 1];
	}
	arrayMoved[count - 1] = array[0];
}

double CircleSimilarity::getVectSimilarity(double *array1, double * array2, int length){
	double modeArray1, modeArray2, similar;
	double sum2Array1 = 0;
	double sum2Array2 = 0;
	double innerProduct = 0;
	for(int l = 0; l < length; l++){
		sum2Array1 += pow(array1[l], 2);
		sum2Array2 += pow(array2[l], 2);
		innerProduct += array1[l] * array2[l];
	}
	modeArray1 = sqrt(sum2Array1);
	modeArray2 = sqrt(sum2Array2);
	if(modeArray1 == 0 || modeArray2 == 0){
		similar = 0;
	}else{		
		similar = abs(innerProduct / (modeArray1 * modeArray2));
	}	
	return similar;	
}

double CircleSimilarity::getMaxValue(double * array, int count){
	double max = 0;
	for(int d = 0; d < count; d++){
		if(array[d] > max){
			max = array[d];
		}		
	}
	return max;
}

