﻿#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "AscGrid.h"
#include "AscGrid.cpp"
#include "mfdRegionSimilarity.h"
#include "mfdRegionSimilarity.cpp"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <ctime>
#define VERY_SMALL 0.000001
#define MEMORY_SIZE 5000
using namespace std;


double string_to_double( const std::string& s ){
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}
 
 //读取被特定符号分隔的字符
void parseStr(string str, char c, vector<string>& tokens){
	unsigned int posL = 0;
	unsigned int posR = 0;
	while(posR < str.length()-1){
		posR = str.find_first_of(c,posL);//字符c在字符串中第一次出现的位置by an
		string sub = str.substr(posL,posR-posL);//从字符串最左边开始到c字符出现之前的子串by an
		tokens.push_back(sub);//sub放入token by an
		posL = posR + 1;//posL改为c后第一位by an
	}
}

string& trimLineBreak(string &s){
	if (s.empty()){
		return s;
	}
	s.erase(0,s.find_first_not_of("\n"));//去除字符串前的回车，find_first_not_of("\n")查找第一个与回车不匹配的字符，返回它的位置
	s.erase(s.find_last_not_of("\n") + 1);//去除字符串后的回车，find_last_not_of("\n")查找不包含子串中的任何字符，返回最后一个位置 
	return s;
}

int getMin(int a, int b){
	if(a > b)
		return b;
	else
		return a;
}

double getSimilaityIntegration(double * similarity, int numOfLyr, double noData){
	double min = 1;
	for(int f = 0; f < numOfLyr; f++){
		if(abs(similarity[f] - noData) < VERY_SMALL){			
			min = noData;
			break;
		}else if(min > similarity[f]){
			min = similarity[f];
		}
	}
	return min;
}	

double getMaxValue(double* similarity, int length, double noData){
	double max = 0;
	int cnt = 0;
	for(int l = 0; l < length; l++){
		if(max < similarity[l]){
			max = similarity[l];
		}else if(abs(similarity[l] - noData) < VERY_SMALL){
			cnt++;
		}
	}
	if(cnt == length)
		return noData;
	else
		return max;
}

double getEnvRange(int totalRows, int totalCols, double**env, double noData){
	double max = 0;
	double min = 1000000;
	for(int m = 0; m < totalRows; m++){
		for(int n = 0; n < totalCols; n++){
			if(abs(env[m][n] - noData) > VERY_SMALL){
				 if(env[m][n] > max){
					 max = env[m][n];
				 }
				 if(env[m][n] < min){
					 min = env[m][n];
				 }		 
			}
		}
	}
	return max - min;
}


int main(int argc, char *argv[]){

	GDALAllRegister();	
    // prepare parameters begin
	char * demLyrsPath = argv[1];//高程数据，仅用来识别特征邻域 
	char * environLyrsPath = argv[2];//与邻域相似性计算有关的环境变量（很可能包括高程，此时需要重读一遍）
	char * trainSamplePath = argv[3];
	char * xName = argv[4];
	char * yName = argv[5];
	char * propertyName = argv[6];	
	char * propertyPredPath = argv[7];
	char * predUncertaintyPath = argv[8];
	char * segmentNumEnv = argv[9];	
	int envN = atoi(segmentNumEnv);
	time_t beginTime, endTime, usedTime;
	beginTime = clock();
	
	//读高程数据
	string demPath = demLyrsPath;	
	AscGrid demLyr;	
	demLyr.readAscGridGDAL(demPath);	
	int totalRows = demLyr.getNumOfRows(); 
	int totalCols = demLyr.getNumOfCols();
	double cellSize = demLyr.getCellSize();
	double lowerLeftX = demLyr.getXCor();
	double lowerLeftY = demLyr.getYCor();
	double noData = demLyr.getNodaVal();	
	double ** dem = new double * [totalRows];	
	for(int i = 0; i < totalRows; i++){
		dem[i] = new double[totalCols];	
	}
	dem = demLyr.values;
	cout << "read elevation data finished" << endl;	
	
	/*读取与邻域相似性计算有关的环境数据:1.获取环境变量文件名称；2.读环境变量数据*/
	vector<string> environLyrs;
	parseStr(string(environLyrsPath),'#',environLyrs);
	int numOfLyr = environLyrs.size();
	AscGrid * envLyr = new AscGrid[numOfLyr];	
	double *** envValue = new double **[numOfLyr];	
	double * minEnv = new double[numOfLyr];
	double * maxEnv = new double[numOfLyr];	
	for(int f = 0; f < numOfLyr; f++){
		envValue[f] = new double * [totalRows];		
		for(int i = 0; i < totalRows; i++){
			envValue[f][i] = new double [totalCols];			
		}
	}	
	for(int f = 0; f < numOfLyr; f++){
		envLyr[f].readAscGridGDAL(environLyrs[f]);
		envValue[f] = envLyr[f].values;
		minEnv[f] = envLyr[f].getMin();
		maxEnv[f] = envLyr[f].getMax();
		//cout << maxEnv[f] << ", " << minEnv[f] << endl;
	}
	delete [] envLyr;
	cout << "read env data finished" << endl;
	
	
	//读训练样点和验证样点数据
	vector<vector<string> > trainSampleList;
	vector<string> fields;
	FILE *pfin = NULL;
	pfin = fopen(trainSamplePath,"r");
	if(pfin == NULL)
		cerr<<"fail to open training sample file"<<endl;
	char row[500];
	while (fgets(row,500,pfin)!= NULL){
		string line = string(row);
		parseStr(line,',',fields);
		trainSampleList.push_back(fields);
		fields.clear();
	}
	fclose(pfin);	
	
	int numOfTrainSamples = trainSampleList.size() - 1; // the first row in training sampefile
	int columnsTrainSamples = trainSampleList[0].size();// 样点文件列数,sampleList[0]即样点文件中的第一行title
	int xIndexTrain = 0;
	int yIndexTrain = 0;
	int propertyIndexTrain = 0;	
	int* trainSampleRows = new int[numOfTrainSamples];
	int* trainSampleCols = new int[numOfTrainSamples];
	double* trainProperty = new double[numOfTrainSamples];
	for (int j = 0; j < columnsTrainSamples; j++){
		//注意去除换行符
		if (strcmp(xName, trimLineBreak(trainSampleList[0][j]).c_str()) == 0){
			xIndexTrain = j;
		}
		if (strcmp(yName, trimLineBreak(trainSampleList[0][j]).c_str()) == 0){
			yIndexTrain = j;
		}
		if (strncmp(propertyName, trimLineBreak(trainSampleList[0][j]).c_str(), strlen(propertyName)) == 0){
			propertyIndexTrain = j;
		}
	}
	//cout << xIndexTrain << ", " << yIndexTrain << ", " << propertyIndexTrain << endl;
	for(int i = 1; i <= numOfTrainSamples;i++){
		double curXTrain = string_to_double(trainSampleList[i][xIndexTrain]);//x属性列的值
		double curYTrain = string_to_double(trainSampleList[i][yIndexTrain]);
		trainSampleRows[i - 1] = totalRows - (int)((curYTrain - lowerLeftY) / cellSize) - 1;
		trainSampleCols[i-1] = (int)((curXTrain	- lowerLeftX) / cellSize);
		trainProperty[i - 1] = string_to_double(trainSampleList[i][propertyIndexTrain]);
	}	
	cout << "read training sampefiles finished" << endl;	
	
	//定义三个vector，分别用来存放邻域8个位置的相对坐标和相应栅格流入中心栅格时的流向值
	int array1[] = {0, -1, -1, -1, 0, 1, 1, 1};
	int array2[] = {-1, -1, 0, 1, 1, 1, 0, -1};
	int array3[] = {1, 2, 4, 8, 16, 32, 64, 128};
	vector <int> rNeighbor;
	vector <int> cNeighbor;
	vector <int> sdfNeighbor;	
	for(int i = 0; i < 8; i++){		
		rNeighbor.push_back(array1[i]);
		cNeighbor.push_back(array2[i]);
		sdfNeighbor.push_back(array3[i]);
	}
	
	mfdRegionSimilarity neighbor;
	vector<double> neighborEnv;	
	double ** flag = new double *[totalRows];
	double *** frequencyTrain = new double **[numOfTrainSamples];
	for(int row = 0; row < totalRows; row++){			
		flag[row] = new double[totalCols];		
	}	
	for(int i = 0; i < numOfTrainSamples; i++){
		frequencyTrain[i] = new double * [numOfLyr];
		for(int f = 0; f < numOfLyr; f++){
			frequencyTrain[i][f] = new double[envN];
			for(int j = 0; j < envN; j++){
				frequencyTrain[i][f][j] = 0;
			}
		}
	}	
	for(int i = 0; i < numOfTrainSamples; i++){
		int rowTrain = trainSampleRows[i];
		int colTrain = trainSampleCols[i];	
		double demTrain = dem[rowTrain][colTrain];
		if(abs(demTrain - noData) > VERY_SMALL){		
			for(int row = 0; row < totalRows; row++){		
				for(int col = 0; col < totalCols; col++){
					flag[row][col] = noData;
				}
			}		
			int neighborSize = neighbor.getNeighborSize(rowTrain, colTrain, dem, flag, totalRows, 
				totalCols, noData, rNeighbor, cNeighbor);						
			for(int f = 0; f < numOfLyr; f++){	
				//cout<<"f = " << f << endl;
				neighbor.getNeighborEnv(rowTrain, colTrain, totalRows, totalCols, neighborSize, 
					envValue[f], flag, noData, neighborEnv);		
				neighbor.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequencyTrain[i][f]);
				neighborEnv.clear();
			}
		}			
	}	
	cout << "calculating frequencyTrain finished" << endl;
	vector<double>().swap(neighborEnv);
//sleep(1000);

	
	double ** propertyPredValue = new double* [totalRows];//区域推测值和推测不确定性
	double ** predUncertainty = new double * [totalRows];
	
	for(int i = 0; i < totalRows; i++){
		propertyPredValue[i] = new double[totalCols];
		predUncertainty[i] = new double[totalCols];
	}
	omp_set_num_threads(8);		
	#pragma omp parallel for	
	//计算验证样点和每个训练样点的相似性
	for(int i = 0; i < totalRows; i++){
		for(int j = 0; j < totalCols; j++){
			int demPixel = dem[i][j];		
			if(abs(demPixel - noData) > VERY_SMALL){
				double * sampleSimilarity = new double [numOfTrainSamples];
				for(int train = 0; train < numOfTrainSamples; train++){
					sampleSimilarity[train] = noData;
				}
				for(int row = 0; row < totalRows; row++){		
					for(int col = 0; col < totalCols; col++){
						flag[row][col] = noData;
					}
				}
				mfdRegionSimilarity pixel;
				int neighborSize = pixel.getNeighborSize(i, j, dem, flag, totalRows, totalCols, 
					noData, rNeighbor, cNeighbor);
				double ** frequency = new double * [numOfLyr];
				vector<double> neighborEnv;
				for(int f = 0; f < numOfLyr; f++){
					frequency[f] = new double [envN];
					for(int k = 0; k < envN;k++){
						frequency[f][k] = 0;
					}
					pixel.getNeighborEnv(i, j, totalRows, totalCols, neighborSize, envValue[f],
						flag, noData, neighborEnv);
					pixel.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequency[f]);
					neighborEnv.clear();				
				}
				//cout << "neighborSize = " << neighborSize << endl;
				vector<double>().swap(neighborEnv);
				for(int train = 0; train < numOfTrainSamples; train++){		
					int rowTrain = trainSampleRows[train];
					int colTrain = trainSampleCols[train];	
					double demTrain = dem[rowTrain][colTrain];
					if(abs(demTrain - noData) > VERY_SMALL){	
						double * similarityH = new double [numOfLyr];
						for(int f = 0; f < numOfLyr; f++){					
							similarityH[f] = pixel.histSimilarity(frequency[f], frequencyTrain[train][f], envN);			
							
						}
						double similarityEnv = getSimilaityIntegration(similarityH, numOfLyr, noData);
						sampleSimilarity[train]= similarityEnv;
						//cout << i << ", " << j << ", " << similarityEnv << endl;
						delete [] similarityH;	
					}
				}				
				double sampleSimilarityMax = getMaxValue(sampleSimilarity, numOfTrainSamples, noData);
				if(abs(sampleSimilarityMax - noData) < VERY_SMALL){
					predUncertainty[i][j] = noData;
					propertyPredValue[i][j] = noData;
				}else if(sampleSimilarityMax == 0){
					predUncertainty[i][j] = 1;
					propertyPredValue[i][j] = noData;
				}else{
					predUncertainty[i][j] = 1 - sampleSimilarityMax;
					double sumProperty = 0;
					double sumSimilarity = 0;				
					for(int train = 0; train < numOfTrainSamples; train++){	
						if(abs(sampleSimilarity[train] - noData) > VERY_SMALL){														
							sumProperty += trainProperty[train] * sampleSimilarity[train];
							sumSimilarity += sampleSimilarity[train];
						}
					}
					propertyPredValue[i][j] = sumProperty / sumSimilarity;
					cout << i << ", " << j << ", " << propertyPredValue[i][j] << endl;
				}
			}
		}			
	}	
	delete [] envValue;
	delete [] minEnv;
	delete [] maxEnv;
	
	//验证样点属性推理
	
	cout << "predicting area finished" << endl;
	//将预测结果输出
	AscGrid pred(totalCols, totalRows, lowerLeftX, lowerLeftY, cellSize, noData, propertyPredValue);
	string propertyOut = propertyPredPath;
	pred.createAscGridGADL(environLyrs[0], propertyOut);
	pred.writeAscGridGDAL(propertyOut);
	
	AscGrid uncertainty(totalCols, totalRows, lowerLeftX, lowerLeftY, cellSize, noData, predUncertainty);
	string uncertaintyOut = predUncertaintyPath;
	uncertainty.createAscGridGADL(environLyrs[0], uncertaintyOut);
	uncertainty.writeAscGridGDAL(uncertaintyOut);
	cout << "output result finished" << endl;
	
	for(int i = 0; i < totalRows; i++){
		delete [] propertyPredValue[i];
		delete [] predUncertainty[i];
	}
				
	delete [] propertyPredValue;
	delete [] predUncertainty;
	
	endTime = clock();
	usedTime = double(endTime - beginTime)/CLOCKS_PER_SEC;
	cout << "time used is " << usedTime << ", seconds" << endl;
	delete[] flag;
	delete[] dem;
	delete[] trainSampleRows;
	delete[] trainSampleCols;
	delete[] trainProperty;	
	cout << "OK!" << endl;	
	return 0;
	
}

