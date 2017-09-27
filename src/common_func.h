#ifndef H_COMMON_FUNC
#define H_COMMON_FUNC
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include "commonLib.h"
#include "flowInCells.h"

using namespace::std;
#define VERY_SMALL 0.000001

char* string_to_char(const string s);
double string_to_double( const std::string& s );
void parseStr(string str, char c, vector<string>& tokens);
string& trimLineBreak(string &s);
int getMin(int a, int b);
double getSimilaityIntegration(double * similarity, int numOfLyr, double noData);
double getMaxValue(double* similarity, int length, double noData);
void readSamples(char* samplePath, string xName, string yName, string propertyName,
	vector<double> & xSamples, vector<double> & ySamples, vector<double> & attrSample);
void getExtremeEnv(vector<float> env, float &minValue, float &maxValue);
void frequencySta(float minValue, float maxValue, int freqnum, vector<float> neighborEnv,
	float * &frequency);
float histSimilarity(float * & frequencyPixel, float * & frequencySample, int number);
void predictProperty(tdpartition * pred_part, tdpartition *uncer_part, int nx, int ny, 
	map<int, Cell*> cells_map, float *** frequencyTrain, float *minEnv, float *maxEnv, 
	vector<double> attrs_train, int env_num, int freqnum, int num_train);
#endif /* H_COMMON_FUNC */
