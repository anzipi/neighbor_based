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
#include "sfdRegionSimilarity.h"


#define MEMORY_SIZE 5000

using namespace std;


int main(int argc, char *argv[]){

	GDALAllRegister();	
    // prepare parameters begin
	char * demLyrsPath = argv[1];//filled dem data  
	char * flowDirectionPath = argv[2];//流flow direction file
	char * environLyrsPath = argv[3];//covariate data
	char * trainSamplePath = argv[4];
	char * testSamplePath = argv[5];
	char * xName = argv[6];
	char * yName = argv[7];
	char * propertyName = argv[8];
	char * propertyPredPath = argv[9];
	char * intervalNumber = argv[10];
	
	time_t beginTime, endTime, usedTime;
	beginTime = clock();
	int envN = atoi(intervalNumber);
	
	//读高程数据
	string demPath = demLyrsPath;
	string sfdPath = flowDirectionPath;
	AscGrid demLyr, sfdLyr;	
	demLyr.readAscGridGDAL(demPath);	
	sfdLyr.readAscGridGDAL(sfdPath);
	int totalRows = demLyr.getNumOfRows(); 
	int totalCols = demLyr.getNumOfCols();
	double cellSize = demLyr.getCellSize();
	double lowerLeftX = demLyr.getXCor();
	double lowerLeftY = demLyr.getYCor();
	double noData = demLyr.getNodaVal();	
	double ** dem = new double * [totalRows];	
	double ** sfd = new double *[totalRows];
	for(int i = 0; i < totalRows; i++){
		dem[i] = new double[totalCols];
		sfd[i] = new double [totalCols];
	}
	dem = demLyr.values;
	sfd = sfdLyr.values;
	cout << "read elevation and sfd data finished" << endl;

	/*read env data:1.get the file names of env layers, 2.read env layers*/
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
	
	//read training samples and testing samples
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
	
	//define 3 vectors to record the relative coordinates of 8 pixels in the 3*3 neighborhood in rectangle shape
	// and the value of flow direction if the pixel flow into the central location
	int array1[] = {0, -1, -1, -1, 0, 1, 1, 1};
	int array2[] = {-1, -1, 0, 1, 1, 1, 0, -1};
	int array3[] = {1, 2, 4, 8, 16, 32, 64, 128};
	vector <int> rNeighbor;
	vector <int> cNeighbor;
	vector <int> sfdNeighbor;
	for(int i = 0; i < 8; i++){
		rNeighbor.push_back(array1[i]);
		cNeighbor.push_back(array2[i]);
		sfdNeighbor.push_back(array3[i]);
	}
	
	//determin the neighborhood size of training samples
	// and calculate the frequency of env values from the histogram of the neighborEnv
	double ** flag = new double * [totalRows];
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
	for(int i = 0; i < totalRows; i++){
		flag[i] = new double [totalCols];		
	}
	sfdRegionSimilarity train;	
	vector<double> neighborEnv;
	for(int i = 0; i < numOfTrainSamples; i++){		
		for(int row = 0; row < totalRows; row++){			
			for(int col = 0; col < totalCols; col++){
				flag[row][col] = noData;
			}
		}
		int rowTrain = trainSampleRows[i];
		int colTrain = trainSampleCols[i];
		double demTrain = dem[rowTrain][colTrain];
		if(abs(demTrain - noData) > VERY_SMALL){
			int neighborSize = train.getCharacterNeighbor(rowTrain, colTrain, dem, sfd, totalRows, 
				totalCols, noData, rNeighbor, cNeighbor, sfdNeighbor, flag);
			//cout << neighborSize << endl;
			
			
			for(int f = 0; f < numOfLyr; f++){			
				train.getNeighborEnv(rowTrain, colTrain, totalRows, totalCols, neighborSize, envValue[f], 
					flag, noData, neighborEnv);				
				//cout << neighborEnv.size() << endl;
				train.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequencyTrain[i][f]);
				neighborEnv.clear();	
			}
		}			
		//cout << rCord.size() << endl;			
	}	
	cout << "calculation neighborSize of training samples finished" << endl;
	//sleep(1000);
	
	double * propertyPredValue = new double[numOfTestSamples];
	double * predUncertainty = new double[numOfTestSamples];
	double ** sampleSimilarity = new double*[numOfTestSamples];
	for(int i = 0; i < numOfTestSamples; i++){
		sampleSimilarity[i] = new double [numOfTrainSamples];	
		for(int j = 0; j < numOfTrainSamples; j++){
			sampleSimilarity[i][j] = noData;
		}
	}
	sfdRegionSimilarity test;	
	//calculate the env similarity over the spatial neighborhood between training and testing samples
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
			int neighborSize = test.getCharacterNeighbor(rowTest, colTest, dem, sfd, totalRows, totalCols,
				noData,	rNeighbor, cNeighbor, sfdNeighbor, flag);
			double ** frequency = new double * [numOfLyr];
			for(int f = 0; f < numOfLyr; f++){
				frequency[f] = new double [envN];
				for(int k = 0; k < envN;k++){
					frequency[f][k] = 0;
				}
				test.getNeighborEnv(rowTest, colTest, totalRows, totalCols, neighborSize, envValue[f], 
					flag, noData, neighborEnv);
				test.frequencySta(minEnv[f], maxEnv[f], envN, neighborEnv, frequency[f]);
				neighborEnv.clear();				
			}
			for(int j = 0; j < numOfTrainSamples; j++){
				int rowTrain = trainSampleRows[j];
				int colTrain = trainSampleCols[j];
				double demTrain = dem[rowTrain][colTrain];
				if(abs(demTrain - noData) > VERY_SMALL){
					double * similarityH = new double [numOfLyr];
					for(int f = 0; f < numOfLyr; f++){					
						similarityH[f] = test.histSimilarity(frequency[f], frequencyTrain[j][f], envN);			
						//cout << similarityH[f] << endl;					
					}
					double similarityEnv = getSimilaityIntegration(similarityH, numOfLyr, noData);
					sampleSimilarity[i][j] = similarityEnv;
					//cout << i << ", " << j << ", " << similarityEnv << endl;
					delete [] similarityH;
				}					
			}
			for(int f = 0; f < numOfLyr; f++){
				delete [] frequency[f];
			}
			delete [] frequency;				
			
		}			
		//cout << neighborSize << endl;		
	}
	vector<double>().swap(neighborEnv);	
	delete [] frequencyTrain;
	delete[] trainSampleRows;
	delete[] trainSampleCols;
	delete [] dem;
	delete [] envValue;	
	
	//evaluate the soil property value and calculate the predition uncertainty of testing samples
	for(int i = 0; i < numOfTestSamples; i++){
		double sampleSimilarityMax = getMaxValue(sampleSimilarity[i], numOfTrainSamples, noData);
		//cout << sampleSimilarityMax << endl;
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
			double maxSimilarityProperty = 0;
			double maxSimilarity = 0;
			for(int j = 0; j < numOfTrainSamples; j++){	
				if(abs(sampleSimilarity[i][j] - noData) > VERY_SMALL){														
					sumProperty += trainProperty[j] * sampleSimilarity[i][j];
					sumSimilarity += sampleSimilarity[i][j];
				}
				/*if(abs(sampleSimilarity[i][j] - noData) > VERY_SMALL){														
					sumProperty += trainProperty[j] * sampleSimilarity[i][j];
					sumSimilarity += sampleSimilarity[i][j];
					if(sampleSimilarity[i][j] > maxSimilarity){
						maxSimilarity = sampleSimilarity[i][j];
						maxSimilarityProperty = maxSimilarity * trainProperty[j];
					}
					
				}*/
			}
			propertyPredValue[i] = sumProperty / sumSimilarity;
			//propertyPredValue[i] = maxSimilarityProperty + (1 - maxSimilarity) * (sumProperty - maxSimilarityProperty) / (sumSimilarity - maxSimilarity);

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
	
	
	/*ofstream out(similarityPath); //输出样点相似性
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

	
	delete[] propertyPredValue;
	delete[] predUncertainty;
	delete[] trainProperty;
	
	cout << "OK!" << endl;	
	return 0;
	
}

