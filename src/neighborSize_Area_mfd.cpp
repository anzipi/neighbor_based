#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "AscGrid.h"
#include "AscGrid.cpp"
#include "PointSimilarity.h"
#include "PointSimilarity.cpp"
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
	char * demLyrsPath = argv[1];//环境变量文件，用#隔开	
	char * trainSamplePath = argv[2];
	char * xName = argv[3];
	char * yName = argv[4];
	char * propertyName = argv[5];
	char * propertyPredPath = argv[6];
	char * predUncertaintyPath = argv[7];

	
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
	
	//定义三个vector，分别用来存放邻域8个位置的相对坐标和相应栅格流入中心栅格时的流向值
	int array1[] = {0, -1, -1, -1, 0, 1, 1, 1};
	int array2[] = {-1, -1, 0, 1, 1, 1, 0, -1};
	int array3[] = {1, 2, 4, 8, 16, 32, 64, 128};
	vector <int> rNeighbor(8);
	vector <int> cNeighbor(8);
	vector <int> sdfNeighbor(8);
	for(int i = 0; i < 8; i++){
		rNeighbor[i] = array1[i];
		cNeighbor[i] = array2[i];
		sdfNeighbor[i] = array3[i];
	}
	delete [] array1;
	delete [] array2;
	delete [] array3;
	
	mfdRegionSimilarity train;	
	vector<double> neighborEnv;
	//int * areaNeighborTrain = new int [numOfTrainSamples];
	//double *** featureTrain = new double **[numOfTrainSamples * numOfLyr];
	double ** flag = new double *[totalRows];
	double *** frequencyTrain = new double **[numOfLyr];
	for(int row = 0; row < totalRows; row++){			
		flag[row] = new double[totalCols];		
	}
	for(int f = 0; f < numOfLyr; f++){
		frequencyTrain[f] = new double * [numOfTrainSamples];
		for(int i = 0; i < numOfTrainSamples; i++){
			frequencyTrain[f][i] = new double[envN];
			for(int j = 0; j < envN; j++){
				frequencyTrain[f][i][j] = 0;
			}
		}
	}	
	for(int i = 0; i < numOfTrainSamples; i++){
		int rowTrain = trainSampleRows[i];
		int colTrain = trainSampleCols[i];			
		for(int row = 0; row < totalRows; row++){		
			for(int col = 0; col < totalCols; col++){
				flag[row][col] = noData;
			}
		}		
		int neighborSize = train.getNeighborSize(rowTrain, colTrain, dem, flag, totalRows, 
			totalCols, noData, rNeighbor, cNeighbor);				
		for(int f = 0; f < numOfLyr; f++){	
			//cout<<"f = " << f << endl;
			train.getNeighborEnv(rowTrain, colTrain, totalRows, totalCols, neighborSize, 
				envValue[f], flag, noData, neighborEnv);		
			train.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequencyTrain[f][i]);
			neighborEnv.clear();
		}	
	}	
	cout << "calculating frequencyTrain finished" << endl;
	vector<double>.swap(neighborEnv);
	
	
					
	double ** neighborSize = new double *[totalRows];	
	for(int i = 0; i < totalRows; i++){
		propertyPredValue[i] = new double[totalCols];
		predUncertainty[i] = new double[totalCols];
		neighborSize[i] = new double[totalCols];
		flag[i] = new int[totalCols];	
		for(int j = 0; j < totalCols; j++){
			propertyPredValue[i][j] = noData;
			predUncertainty[i][j] = noData;
			neighborSize[i][j] = noData;
			flag[i][j] = noData;			
		}
	}
	
	//omp_set_num_threads(4);		
	//#pragma omp parallel for		
	for(int i = 0; i < totalRows; i++){
		for(int j = 0; j < totalCols; j++){				
			if(abs(dem[i][j] - noData) > VERY_SMALL){
				for(int row = 0; row < totalRows; row++){		
					for(int col = 0; col < totalCols; col++){
						flag[row][col] = noData;
					}
				}
				mfdRegionSimilarity pixel;
				neighborSize[i][j] = pixel.getNeighborSize(i, j, dem, flag, totalRows, totalCols,
					noData,rNeighbor, cNeighbor);				
				cout << i << ", " << j << ", " << neighborSize[i][j] << ", " << endl;
			}
			
		}	
	}
	
	string fileName = "./heshan/result/maxDist2.tif";
	AscGrid size(totalCols, totalRows, lowerLeftX, lowerLeftY, cellSize, noData, dist);
	size.createAscGridGADL(demPath, fileName);		
	size.writeAscGridGDAL(fileName);
	cout << "soil property prediction of the whole area finished" << endl;

	
	cout << "calculation finished" << endl;
	endTime = clock();
	usedTime = double(endTime - beginTime)/CLOCKS_PER_SEC;
	cout << "time used is " << usedTime << ", seconds" << endl;

	//将预测结果输出
	/*AscGrid pred(totalCols, totalRows, lowerLeftX, lowerLeftY, cellSize, noData, propertyPredValue);
	pred.createAscGridGADL(environLyrs[0], propertyPredPath);
	pred.writeAscGridGDAL(propertyPredPath);
	
	AscGrid uncertainty(totalCols, totalRows, lowerLeftX, lowerLeftY, cellSize, noData, predUncertainty);
	uncertainty.createAscGridGADL(environLyrs[0], predUncertaintyPath);
	uncertainty.writeAscGridGDAL(predUncertaintyPath);
	cout << "output result finished" << endl;*/
	
	delete[] flag;
	delete [] dem;
	delete [] neighborSize;
	delete[] trainSampleRows;
	delete[] trainSampleCols;
	delete[] trainProperty;
	delete[] curXTrain;
	delete[] curYTrain;
		
	cout << "OK!" << endl;	
	return 0;
	
}

