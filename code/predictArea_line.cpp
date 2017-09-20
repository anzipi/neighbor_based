#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "AscGrid.h"
#include "AscGrid.cpp"
#include "LineSimilarity.h"
#include "LineSimilarity.cpp"
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
		}else{
			continue;
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
	char * demLyrsPath = argv[1];//环境变量文件，用#隔开	
	char * environLyrsPath = argv[2];//与邻域相似性计算有关的环境变量（很可能包括高程，此时需要重读一遍）
	char * trainSamplePath = argv[3];
	char * xName = argv[4];
	char * yName = argv[5];
	char * propertyName = argv[6];
	char * propertyPredPath = argv[7];
	char * predUncertaintyPath = argv[8];
	char * breakNumber = argv[9];
	
	int envN = atoi(breakNumber);
	
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
	
	//读训练样点数据
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
		vector<string>().swap(fields);
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
	double *curXTrain = new double[numOfTrainSamples];
	double *curYTrain = new double[numOfTrainSamples];
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
		curXTrain[i - 1] = string_to_double(trainSampleList[i][xIndexTrain]);//x属性列的值
		curYTrain[i - 1] = string_to_double(trainSampleList[i][yIndexTrain]);
		trainSampleRows[i - 1] = totalRows - (int)((curYTrain[i-1] - lowerLeftY) / cellSize) - 1;
		trainSampleCols[i-1] = (int)((curXTrain[i-1] - lowerLeftX) / cellSize);
		trainProperty[i - 1] = string_to_double(trainSampleList[i][propertyIndexTrain]);
	}
	cout << "read training samplefiles finished" << endl;	
	vector<vector<string> >().swap(trainSampleList);
	
	//计算训练样点的流向长度和邻域特征
	double * flowLengthTrain = new double[numOfTrainSamples];
	double *** frequencyTrain = new double **[numOfTrainSamples];
	for(int i = 0; i < numOfTrainSamples; i++){
		frequencyTrain[i] = new double * [numOfLyr];
		for(int f = 0; f < numOfLyr; f++){
			frequencyTrain[i][f] = new double[envN];
			for(int j = 0; j < envN; j++){
				frequencyTrain[i][f][j] = 0;
			}
		}
	}
	LineSimilarity neighbor;
	vector<int> rCord;
	vector<int> cCord;	
	vector<double> neighborDEM;
	vector<double> neighborEnv;
	for(int i = 0; i < numOfTrainSamples; i++){		
		flowLengthTrain[i] = 0;
		int rowTrain = trainSampleRows[i];
		int colTrain = trainSampleCols[i];
		neighbor.getCharacterNeighbor(rowTrain, colTrain, dem, totalRows, totalCols, noData, 
			rCord, cCord, neighborDEM, flowLengthTrain[i]);
		//cout << flowLengthTrain[i] << endl;
		for(int f = 0; f < numOfLyr; f++){			
			neighbor.getNeighborEnv(rowTrain, colTrain, rCord, cCord, envValue[f], 
				neighborEnv, noData);
			//cout << neighborEnv.size() << endl;
			neighbor.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequencyTrain[i][f]);
			neighborEnv.clear();	
		}
		rCord.clear();
		cCord.clear();
		neighborDEM.clear();		
	}
	cout << "calculation features of training samples finished" << endl;
	vector<int>().swap(rCord);
	vector<int>().swap(cCord);		
	vector<double>().swap(neighborDEM);		
	vector<double>().swap(neighborEnv);
	
	double ** propertyPredValue = new double*[totalRows];//区域推测值和推测不确定性
	double ** predUncertainty = new double*[totalRows];

	for(int i = 0; i < totalRows; i++){
		propertyPredValue[i] = new double[totalCols];
		predUncertainty[i] = new double[totalCols];		
		for(int j = 0; j < totalCols; j++){
			propertyPredValue[i][j] = noData;
			predUncertainty[i][j] = noData;			
		}
	}
		
	//omp_set_num_threads(4);		
	//#pragma omp parallel for		
	for(int i = 20; i < totalRows - 20; i++){
		for(int j = 20; j < totalCols - 20; j++){									
			//计算全区的特征邻域大小		
			if(abs(dem[i][j] - noData) > VERY_SMALL){
				vector<int> rCord;
				vector<int> cCord;	
				vector<double> neighborDEM;
				vector<double> neighborEnv;
				double flowLength = 0;
				double * sampleSimilarity = new double[numOfTrainSamples];
				neighbor.getCharacterNeighbor(i, j, dem, totalRows, totalCols, noData, rCord, cCord, 
					neighborDEM, flowLength);
				int neighborSize = rCord.size() - 1;
				double ** frequency = new double * [numOfLyr];
				for(int f = 0; f < numOfLyr; f++){
					frequency[f] = new double[envN];
					for(int k = 0; k < envN;k++){
						frequency[f][k] = 0;
					}
					neighbor.getNeighborEnv(i, j, rCord, cCord, envValue[f], neighborEnv, noData);
					neighbor.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequency[f]);			
				//cout << similarityH[f] << endl;
					neighborEnv.clear();	
				}
				for(int train = 0; train < numOfTrainSamples; train++){
					double similarityL = neighbor.getLengthSimilarity(flowLength, flowLengthTrain[train]);
					double * similarityH = new double [numOfLyr];
					for(int f = 0; f < numOfLyr; f++){						
						similarityH[f] = neighbor.histSimilarity(frequency[f], frequencyTrain[train][f], envN);						
					}				
					double similarityEnv = getSimilaityIntegration(similarityH, numOfLyr, noData);					
					sampleSimilarity[train] = similarityEnv;
					//cout << i << ", " << j << ", " << train << ", " << sampleSimilarity[train] << endl;
					delete [] similarityH;
					
				}
				for(int f = 0; f < numOfLyr; f++){
					delete [] frequency[f];
				}
				delete [] frequency;
				vector<int>().swap(rCord);
				vector<int>().swap(cCord);		
				vector<double>().swap(neighborDEM);		
				vector<double>().swap(neighborEnv);
				
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
					double maxSimilarityProperty = 0;
					double maxSimilarity = 0;
					for(int train = 0; train < numOfTrainSamples; train++){	
						if(abs(sampleSimilarity[train] - noData) > VERY_SMALL){														
							sumProperty += trainProperty[train] * sampleSimilarity[train];
							sumSimilarity += sampleSimilarity[train];
							/*if(sampleSimilarity[train] > maxSimilarity){
								maxSimilarity = sampleSimilarity[train];
								maxSimilarityProperty = maxSimilarity * trainProperty[train];
							}*/
							
						}
					}
					propertyPredValue[i][j] = sumProperty / sumSimilarity;
					//propertyPredValue[i] = maxSimilarityProperty + (1 - maxSimilarity) * (sumProperty - maxSimilarityProperty) / (sumSimilarity - maxSimilarity);

				}
				delete [] sampleSimilarity;
			}			
			cout << i << ", " << j << ", " << propertyPredValue[i][j] << ", " << predUncertainty[i][j] << endl;
		}
		
	}
	cout << "calculation finished" << endl;
	delete [] dem;
	delete [] envValue;
	delete [] frequencyTrain;
	delete [] flowLengthTrain;
	delete [] trainSampleRows;
	delete [] trainSampleCols;
	delete [] trainProperty;
	delete [] curXTrain;
	delete [] curYTrain;
	delete [] minEnv;
	delete [] maxEnv;
	
	endTime = clock();
	usedTime = double(endTime - beginTime)/CLOCKS_PER_SEC;
	cout << "time used is " << usedTime << ", seconds" << endl;

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
	cout << "OK!" << endl;		
	return 0;
	
}

