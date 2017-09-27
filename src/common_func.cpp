#include "common_func.h"
#include "utilities.h"
#include <vector>
#include <map>

char* string_to_char(const string s) {
    char* c;
    const int len = s.length();
    c = new char[len + 1];
    strcpy(c, s.c_str());
    return c;
}

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

float getSimilaityIntegration(float * similarity, int numOfLyr, double noData){
	float min = 1;
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

float getMaxValue(float* similarity, int length, float noData){
	float max = 0;
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


void readSamples(char* samplePath, string xName, string yName, string propertyName,
	vector<double> & xSamples, vector<double> & ySamples, vector<double> & attrSample){
		vector<vector<string> > sampleList;
		vector<string> fields;
		FILE *pfin = NULL;
		pfin = fopen(samplePath,"r");
		if(pfin == NULL)
			cerr<<"fail to open sample file"<<endl;
		char row[500];
		while (fgets(row,500,pfin)!= NULL){
			string line = string(row);
			parseStr(line,',',fields);
			sampleList.push_back(fields);
			fields.clear();
		}
		fclose(pfin);	

		int numOfSamples = sampleList.size() - 1; // the first row in sampefile
		int columnsSamples = sampleList[0].size();// 样点文件列数,sampleList[0]即样点文件中的第一行title
		int xIndex = 0;
		int yIndex = 0;
		int propertyIndex = 0;	
		for (int j = 0; j < columnsSamples; j++){
			//注意去除换行符
			if (xName == sampleList[0][j]){
				xIndex = j;
			}
			if (yName == sampleList[0][j]){
				yIndex = j;
			}
			if (propertyName == sampleList[0][j]){
				propertyIndex = j;
			}
		}
		//cout << xIndex << ", " << yIndex << ", " << propertyIndex << endl;
		for(int i = 1; i <= numOfSamples;i++){
			double curX = string_to_double(sampleList[i][xIndex]);//x属性列的值
			double curY = string_to_double(sampleList[i][yIndex]);
			xSamples.push_back(curX);
			ySamples.push_back(curY);
			double propertyValue = string_to_double(sampleList[i][propertyIndex]);
			attrSample.push_back(propertyValue);
		}	
		//cout << "read sampefiles finished" << endl;	
}

void getExtremeEnv(vector<float> env, float &minValue, float &maxValue){
	int neighborSize = env.size();
	for(int n = 0; n < neighborSize; n++){
		if(minValue > env[n]){
			minValue = env[n];
		}
		if(maxValue < env[n]){
			maxValue = env[n];
		}

	}
}

void frequencySta(float minValue, float maxValue, int freqnum, vector<float> neighborEnv, float * &frequency){
	int neighborSize = neighborEnv.size();
	for(int n = 0; n < neighborSize; n++){
		int index = int((neighborEnv[n] - minValue) * freqnum / (maxValue - minValue));		
		if(index == freqnum){
			index = freqnum - 1;
		}
		frequency[index] += 1;
	}
	for(int n = 0; n < freqnum; n++){
		frequency[n] = frequency[n]  / neighborSize;		
	}	
}

float histSimilarity(float * & frequencyPixel, float * & frequencySample, int number){
	float intersectValue = 0;
	float unionValue = 0;
	float similarityH;
	for(int m = 0; m < number; m++){
		intersectValue += min(frequencyPixel[m], frequencySample[m]);
		unionValue += max(frequencyPixel[m], frequencySample[m]);
	}
	//cout << intersectValue << ", " << unionValue << endl;
	if(unionValue == 0){
		similarityH = 0;
	}else{
		similarityH = intersectValue / unionValue;
	}
	//cout << similarityH << endl;
	return similarityH;
}

void predictProperty(tdpartition * pred_part, tdpartition *uncer_part, int nx, int ny, 
	map<int, Cell*> cells_map, float *** frequencyTrain, float *minEnv, float *maxEnv, 
	vector<double> attrs_train, int env_num, int freqnum, int num_train){
	for(int j = 0; j < ny; j++){
		for(int i = 0; i < nx; i++){
			int pixel_idx = j * nx + i;
			map<int ,Cell*>::iterator iter;
			iter = cells_map.find(pixel_idx);
			if(iter != cells_map.end()){
				Cell* pixel_cell = cells_map.at(pixel_idx);
				vector<int> pixel_upIdxes;
				pixel_cell->getUpCellIndexes(pixel_upIdxes);
				map<int, vector<float> > pixel_upEnvValues;
				for (int kk = 0; kk < env_num; kk++) {
					vector<float> tmpvalues;
					pixel_upEnvValues.insert(make_pair(kk, tmpvalues));
				}
				for (vector<int>::iterator it = pixel_upIdxes.begin(); it != pixel_upIdxes.end(); it++) {
					vector<float> tmpEnvvs = cells_map.at(*it)->getEnvValue();
					if (tmpEnvvs.size() != env_num) {
						continue; // although this may not happen, just check.
					}
					for (int kk = 0; kk < env_num; kk++) {
						pixel_upEnvValues.at(kk).push_back(tmpEnvvs[kk]);						
					}
				}			
				
				float ** frequency = new float *[env_num];
				for(int kk = 0; kk < env_num; kk++){
					frequency[kk] = new float[freqnum];
					for(int n = 0; n < freqnum; n++){
						frequency[kk][n] = 0;
					}
				}
				for(int kk = 0; kk < env_num; kk++){
					frequency[kk] = new float[freqnum];
					frequencySta(minEnv[kk], maxEnv[kk], freqnum, pixel_upEnvValues.at(kk), frequency[kk]);					
				}
				float *sampleSimilarity = new float[num_train];
				for(int train = 0; train < num_train; train ++){
					float * envSimilarity = new float[env_num];
					for(int kk = 0; kk < env_num; kk++){
						envSimilarity[kk] = histSimilarity(frequency[kk], frequencyTrain[train][kk], freqnum);
					}
					sampleSimilarity[train] = getSimilaityIntegration(envSimilarity, env_num, NODATA_VALUE);
				}
				float maxSimilarity = getMaxValue(sampleSimilarity, num_train, NODATA_VALUE);
				float uncertainty = 1 -maxSimilarity;
				float sumPredSimilarity = 0, sumSimilarity = 0;
				for(int train = 0; train < num_train; train ++){
					sumPredSimilarity += attrs_train[train] * sampleSimilarity[train];
					sumSimilarity += sampleSimilarity[train];
				}
				float pred = sumPredSimilarity / sumSimilarity;
				pred_part->setData(i, j, pred);
				uncer_part->setData(i, j, uncertainty);
			}else{
				pred_part->setData(i, j, NODATA_VALUE);
				uncer_part->setData(i, j, NODATA_VALUE);
			}

			
		}

	}
}
