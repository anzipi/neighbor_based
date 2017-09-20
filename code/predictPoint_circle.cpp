#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "AscGrid.h"
#include "AscGrid.cpp"
#include "CircleSimilarity.h"
#include "CircleSimilarity.cpp"
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


int main(int argc, char *argv[]){

	GDALAllRegister();	
    // prepare parameters begin
	char * environLyrsPath = argv[1];//环境变量文件，用#隔开
	char * varLyrsPath = argv[2];
    char * trainSamplePath = argv[3];
	char * testSamplePath = argv[4];
	char * xName = argv[5];
	char * yName = argv[6];
	char * propertyName = argv[7];
	char * propertyPredPath = argv[8];
	char * neighborMin = argv[9];
	char * neighborMax = argv[10];
	char * attenuationPower = argv[11];

	int radiusMin = atoi(neighborMin);
	int radiusMax = atoi(neighborMax);
	double alpha = atof(attenuationPower);
	
	time_t beginTime, endTime, usedTime;
	beginTime = clock();
		
	//获取环境变量文件名称
	vector<string> environLyrs;
	parseStr(string(environLyrsPath),'#', environLyrs);
	int numOfLyr = environLyrs.size();
	
	//获取环境变量var文件名称
	vector<string> envVarLyrs;
	parseStr(string(varLyrsPath),'#', envVarLyrs);	
	
	//读环境变量值和var值
	AscGrid * envLyr = new AscGrid[numOfLyr];
	AscGrid * varLyr = new AscGrid[numOfLyr];
	double *** envValue = new double **[numOfLyr];
	double *** varValue = new double **[numOfLyr];
	double * envStd = new double[numOfLyr];
	envLyr[0].readAscGridGDAL(environLyrs[0]);
	varLyr[0].readAscGridGDAL(envVarLyrs[0]);
	int totalRows = envLyr[0].getNumOfRows(); 
	int totalCols = envLyr[0].getNumOfCols();
	double cellSize = envLyr[0].getCellSize();
	double lowerLeftX = envLyr[0].getXCor();
	double lowerLeftY = envLyr[0].getYCor();
	double noData = envLyr[0].getNodaVal();	
	for(int f = 0; f < numOfLyr; f++){
		envValue[f] = new double * [totalRows];	
		varValue[f] = new double *[totalRows];
		for(int i = 0; i < totalRows; i++){
			envValue[f][i] = new double [totalCols];	
			varValue[f][i] = new double [totalCols];
		}	
	}	
	for(int f = 1; f < numOfLyr; f++){
		envLyr[f].readAscGridGDAL(environLyrs[f]);	
		varLyr[f].readAscGridGDAL(envVarLyrs[f]);
	}
	
	for(int f = 0; f < numOfLyr; f++){		
		envValue[f] = envLyr[f].values;
		envStd[f] = envLyr[f].getStdValue();
		varValue[f] = varLyr[f].values;
	}
	cout << "read env data finished" << endl;
	delete[] envLyr;	
	delete[] varLyr;
	
	
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
	
	vector<vector<string> > testSampleList;
	vector<string> field;
	FILE *pfins = NULL;
	pfins = fopen(testSamplePath,"r");
	if(pfins == NULL)
		cerr<<"fail to open testing sample file"<<endl;
	char rows[500];
	while (fgets(rows,500,pfins)!= NULL){
		string line = string(rows);
		parseStr(line,',',field);
		testSampleList.push_back(field);
		field.clear();
	}
	fclose(pfins);
	
	int numOfTrainSamples = trainSampleList.size() - 1; // the first row in training sampefile
	int columnsTrainSamples = trainSampleList[0].size();// 样点文件列数,sampleList[0]即样点文件中的第一行title
	int xIndexTrain = 0;
	int yIndexTrain = 0;
	int propertyIndexTrain = 0;
	
	int numOfTestSamples = testSampleList.size() - 1; // the firt row in Testing sampefile
	int columnsTestSamples = testSampleList[0].size();// 样点文件列数,sampleList[0]即样点文件中的第一行title
	int xIndexTest = 0;
	int yIndexTest = 0;
	int propertyIndexTest = 0;	
	
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
	int* testSampleRows = new int[numOfTestSamples];
	int* testSampleCols = new int[numOfTestSamples];
	double* testProperty = new double[numOfTestSamples];
	
	for (int j = 0; j < columnsTestSamples; j++){
		//注意去除换行符
		if (strcmp(xName, trimLineBreak(testSampleList[0][j]).c_str()) == 0){
			xIndexTest = j;
		}
		if (strcmp(yName, trimLineBreak(testSampleList[0][j]).c_str()) == 0){
			yIndexTest = j;
		}
		if (strncmp(propertyName, trimLineBreak(testSampleList[0][j]).c_str(), strlen(propertyName)) == 0){
			propertyIndexTest = j;
		}
	}
	//cout << xIndexTest << ", " << yIndexTest << ", " << propertyIndexTest << endl;
	for(int i = 1; i <= numOfTestSamples;i++){
		double curXTest = string_to_double(testSampleList[i][xIndexTest]);//x属性列的值
		double curYTest = string_to_double(testSampleList[i][yIndexTest]);		
		testSampleRows[i - 1] = totalRows - (int)((curYTest - lowerLeftY) / cellSize) - 1;
		testSampleCols[i - 1] = (int)((curXTest	- lowerLeftX) / cellSize);		
		testProperty[i - 1] = string_to_double(testSampleList[i][propertyIndexTest]);
	}
	cout << "read sampefiles finished" << endl;
	vector<vector<string> >().swap(trainSampleList);
	vector<vector<string> >().swap(testSampleList);
	
	//读取训练样点的对全区的Var，识别其特征邻域
	CircleSimilarity train;
	double ** varTrain = new double*[numOfTrainSamples];
	int ** scaleTrain = new int *[numOfTrainSamples];
	for(int i = 0; i < numOfTrainSamples; i++){
		varTrain[i] = new double[numOfLyr];
		scaleTrain[i] = new int[numOfLyr];		
	}
	for(int i = 0; i < numOfTrainSamples; i++){		
		int rowTrain = trainSampleRows[i];
		int colTrain = trainSampleCols[i];
		for(int f = 0; f < numOfLyr; f++){
			double envTrain = envValue[f][rowTrain][colTrain];
			//cout << envValue[f][rowTrain][colTrain]<< ", " << envStd[f]<< endl;
			if(abs(envTrain - noData) > VERY_SMALL){
				varTrain[i][f] = varValue[f][rowTrain][colTrain];				
				scaleTrain[i][f] = train.getCharacteristicScale(rowTrain, colTrain, envTrain, varTrain[i][f], totalRows, 
					totalCols,	envValue[f], noData, envStd[f], radiusMin, radiusMax);	
			//cout << i << "," << scaleTrain[i][f] << endl;
			}else{
				varTrain[i][f] = noData;
				scaleTrain[i][f] = noData;
			}
		}		
	}
	
	/*
	for(int i = 0; i < numOfTrainSamples; i++){
		for(int f = 0; f < numOfLyr; f++){
			cout << i << ", " << varTrain[i][f] << ", "<< scaleTrain[i][f] << endl;
		}
	}
	sleep(1000);*/
	cout << "calculate var and scale of training samples finished" << endl;
	double * propertyPredValue = new double[numOfTestSamples];//区域推测值和推测不确定性
	double * predUncertainty = new double[numOfTestSamples];
	
	//omp_set_num_threads(8);		
	//#pragma omp parallel for	
	for(int i = 0; i < numOfTestSamples; i++){	
		int rowTest = testSampleRows[i];
		int colTest = testSampleCols[i];	
		double * sampleSimilarity = new double[numOfTrainSamples];//待推测点与每个训练样点在样点级别的相似性						
		//待推测点和样点的相似性计算
		for(int j = 0; j < numOfTrainSamples; j++){	
			int rowTrain = trainSampleRows[j];
			int colTrain = trainSampleCols[j];
			double * envSimilarity = new double[numOfLyr];//变量级别的相似性
			for(int f = 0; f < numOfLyr; f++){
				envSimilarity[f] = noData;
			}									
			//待推测点和样点在样点级别的相似性
			for(int f = 0; f < numOfLyr; f++){
				double envTest = envValue[f][rowTest][colTest];
				double envTrain = envValue[f][rowTrain][colTrain];				
				if(abs(envTest - noData) > VERY_SMALL && abs(envTrain - noData) > VERY_SMALL){
					CircleSimilarity test;
					double varTest = varValue[f][rowTest][colTest];
					int scaleTest = test.getCharacteristicScale(rowTest, colTest, envTest, varTest, totalRows, 
						totalCols,	envValue[f], noData, envStd[f], radiusMin, radiusMax);
					int scaleIntegration = getMin(scaleTest, scaleTrain[j][f]);
					//cout << scaleTest << ", " << scaleTrain[j][f] << ", " << scaleIntegration << endl;
					envSimilarity[f] = test.getCircleSimilarity(rowTest, colTest, rowTrain, colTrain, envValue[f], totalRows,
						totalCols, envStd[f], varTrain[j][f], noData, scaleIntegration, alpha);	
				}					
			}
			sampleSimilarity[j] = getSimilaityIntegration(envSimilarity, numOfLyr, noData);
			//cout << i << "," << j << ", " << sampleSimilarity[j] << endl;
			delete[] envSimilarity;					
		}
		//待推测点土壤属性值和不确定性计算
		double sampleSimilarityMax = getMaxValue(sampleSimilarity, numOfTrainSamples, noData);		
		if(abs(sampleSimilarityMax - noData) < VERY_SMALL){
			predUncertainty[i] = noData;
			propertyPredValue[i] = noData;
		}else if(sampleSimilarityMax == 0){
			predUncertainty[i] = 1;
			propertyPredValue[i] = noData;
		}else{
			predUncertainty[i] = 1 - sampleSimilarityMax;
			double sumProperty = 0;
			double sumSimilarity = 0;	
			for(int j = 0; j < numOfTrainSamples; j++){	
				if(abs(sampleSimilarity[j] - noData) > VERY_SMALL){														
					sumProperty += trainProperty[j] * sampleSimilarity[j];
					sumSimilarity += sampleSimilarity[j];
				}
			}
			propertyPredValue[i] = sumProperty / sumSimilarity;
				
		}
		delete[] sampleSimilarity;
		//cout << i << propertyPredValue[i] << endl;				
	}		
	cout << "soil property prediction of the testing samples finished" << endl;	

	//将预测结果输出
	ofstream fin(propertyPredPath); 
	if(!fin){
		cout << "Unable to open outfile"<<endl;// terminate with error
	}
	fin << "ID,testProperty,pred_property,Uncertainty" << endl;
	for(int i = 0; i < numOfTestSamples; i++){
		fin << i<< "," << testProperty[i] <<","<< propertyPredValue[i] << "," << predUncertainty[i]<< endl;
	}		
	fin.close();
	cout << "output result finished" << endl;

	vector<string>().swap(environLyrs);
	delete[] envValue;
	delete[] propertyPredValue;
	delete[] predUncertainty;
	delete[] trainSampleRows;
	delete[] trainSampleCols;
	delete[] trainProperty;
	delete[] envStd;
	delete[] varTrain;	
	delete[] scaleTrain;
	endTime = clock();
	usedTime = double(endTime - beginTime)/CLOCKS_PER_SEC;
	cout << "time used is " << usedTime << ", seconds" << endl;	
	cout << "OK!" << endl;	
	return 0;
	
}

