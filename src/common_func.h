#ifndef H_COMMON_FUNC
#define H_COMMON_FUNC
#include <string>
#include <vector>
#include <sstream>
using namespace::std;
#define VERY_SMALL 0.000001

char* string_to_char(const string s);
double string_to_double( const std::string& s );
void parseStr(string str, char c, vector<string>& tokens);
string& trimLineBreak(string &s);
int getMin(int a, int b);
double getSimilaityIntegration(double * similarity, int numOfLyr, double noData);
double getMaxValue(double* similarity, int length, double noData);

#endif /* H_COMMON_FUNC */