#if (defined _DEBUG) && (defined MSVC) && (defined VLD)
#include "vld.h"
#endif /* Run Visual Leak Detector during Debug */
// 1. Include the header of this cpp if existed, e.g., predictPoint_sfd.h
// 2. Include C/C++ system files.
#ifdef SUPPORT_OMP
#include <omp.h>
#endif
// 3. Include Other libraries' .h files.
#include "utilities.h" // from UtilsClass of Lries-2415
#include "AscGrid.h"
// 4. Include other header files of current project if existed.
#include "common_func.h"
#include "mfdRegionSimilarity.h"

//#define MEMORY_SIZE 5000
using namespace std;

int main(int argc, char *argv[]){

	GDALAllRegister();
    // prepare parameters begin
	char * demLyrsPath = argv[1];//高程数据，仅用来识别特征邻域 
	char * environLyrsPath = argv[2];//与邻域相似性计算有关的环境变量（很可能包括高程，此时需要重读一遍）
	char * trainSamplePath = argv[3];
	char * testSamplePath = argv[4];
	char * xName = argv[5];
	char * yName = argv[6];
	char * propertyName = argv[7];	
	char * propertyPredPath = argv[8];
	char * segmentNumEnv = argv[9];	
    
    double **flag = NULL;
    double ** sampleSimilarity = NULL;
    double * propertyPredValue = NULL;
    double * predUncertainty = NULL;
    int* trainSampleRows = NULL;
    int* trainSampleCols = NULL;
    double* trainProperty = NULL;

	int envN = atoi(segmentNumEnv);
	time_t beginTime, endTime, usedTime;
	beginTime = clock();

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
	
	trainSampleRows = new int[numOfTrainSamples];
	trainSampleCols = new int[numOfTrainSamples];
	trainProperty = new double[numOfTrainSamples];
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
	flag = new double *[totalRows];
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
			//cout << neighborSize << endl;
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
	
	//sleep(1000);

	
	propertyPredValue = new double[numOfTestSamples];//区域推测值和推测不确定性
	predUncertainty = new double[numOfTestSamples];
	sampleSimilarity = new double *[numOfTestSamples];	
	for(int i = 0; i < numOfTestSamples; i++){
		sampleSimilarity[i] = new double[numOfTrainSamples];
		for(int j = 0; j < numOfTrainSamples; j++){
			sampleSimilarity[i][j] = noData;
		}
	}	 
	//计算验证样点和每个训练样点的相似性
	for(int i = 0; i < numOfTestSamples; i++){
		int rowTest = testSampleRows[i];
		int colTest = testSampleCols[i];
		int demTest = dem[rowTest][colTest];
		if(abs(demTest - noData) > VERY_SMALL){
			for(int row = 0; row < totalRows; row++){		
				for(int col = 0; col < totalCols; col++){
					flag[row][col] = noData;
				}
			}		
			int neighborSize = neighbor.getNeighborSize(rowTest, colTest, dem, flag, totalRows, 
				totalCols, noData, rNeighbor, cNeighbor);
			double ** frequency = new double * [numOfLyr];
			for(int f = 0; f < numOfLyr; f++){
				frequency[f] = new double [envN];
				for(int k = 0; k < envN;k++){
					frequency[f][k] = 0;
				}
				neighbor.getNeighborEnv(rowTest, colTest, totalRows, totalCols, neighborSize, 
					envValue[f], flag, noData,neighborEnv);
				neighbor.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequency[f]);
				neighborEnv.clear();				
			}
			//cout << "neighborSize = " << neighborSize << endl;
			
			for(int j = 0; j < numOfTrainSamples; j++){	
				int rowTrain = trainSampleRows[j];
				int colTrain = trainSampleCols[j];
				double demTrain = dem[rowTrain][colTrain];
				if(abs(demTrain - noData) > VERY_SMALL){
					double * similarityH = new double [numOfLyr];
					for(int f = 0; f < numOfLyr; f++){					
						similarityH[f] = neighbor.histSimilarity(frequency[f], frequencyTrain[j][f], envN);			
						
					}
					double similarityEnv = getSimilaityIntegration(similarityH, numOfLyr, noData);
					sampleSimilarity[i][j] = similarityEnv;
					//cout << i << ", " << j << ", " << similarityEnv << endl;
					delete [] similarityH;						
				}									
			}				
		}			
	}
	vector<double>().swap(neighborEnv);
	delete [] envValue;
	delete [] minEnv;
	delete [] maxEnv;
	
	//验证样点属性推理
	for(int i = 0; i < numOfTestSamples; i++){
		double sampleSimilarityMax = getMaxValue(sampleSimilarity[i], numOfTrainSamples, noData);
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
				if(abs(sampleSimilarity[i][j] - noData) > VERY_SMALL){														
					sumProperty += trainProperty[j] * sampleSimilarity[i][j];
					sumSimilarity += sampleSimilarity[i][j];
				}
			}
			propertyPredValue[i] = sumProperty / sumSimilarity;
		}
	}
	cout << "predicting testing samples finished" << endl;
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
	
	/*
	ofstream out(similarityPath); 
	if(!out){
		cout << "Unable to open outfile"<<endl;// terminate with error
	}
	out << "trainID,testID,similarity" << endl;
	for(int i = 0; i < numOfTestSamples; i++){
		for(int j = 0; j < numOfTrainSamples; j++){
			out << i<< "," << j <<","<< sampleSimilarity[i][j]<< endl;
		}		
	}		
	out.close();*/
	
	
	endTime = clock();
	usedTime = double(endTime - beginTime)/CLOCKS_PER_SEC;
	cout << "time used is " << usedTime << ", seconds" << endl;
	delete[] flag;
	delete[] dem;
	delete[] sampleSimilarity;
	delete[] propertyPredValue;
	delete[] predUncertainty;
	delete[] trainSampleRows;
	delete[] trainSampleCols;
	delete[] trainProperty;	
	cout << "OK!" << endl;	
	return 0;
	
}

