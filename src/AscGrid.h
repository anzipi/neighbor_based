//for environmental data input and
//inference result output
//using ascgrid file format

#ifndef ASCGRID_H
#define ASCGRID_H
#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <vector>
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_string.h"
using namespace std;
class AscGrid
{
	public:
		AscGrid();
		AscGrid(int col, int row, double xCor,
			double yCor, double cell, double noData, double** values);
		~AscGrid();

		void setHeadInfo(int col, int row, double xCor,
			double yCor, double cell, double noData);
		void setValues(double** value);
		void setCellValue(int row, int col, double value);
		void setOutputFileName(string filename);
        //void setStdValue(double stdv);//wtf
		double getValue(int row, int col);
		double** getValMatrix();//anym
		
		double getMax();
		double getMin();
		int getNumOfRows();
		int getNumOfCols();
		double getXCor();
		double getYCor();
		double getCellSize();
		double getNodaVal();
		double getStdValue();//wtf

		//bool nearNoDataCell(int col, int row);
		bool isValidCell(int col, int row);
		//bool isValidCell(int col, int row, double value);

		double * readAscGridGDAL_Block(string filename, int offsetX, int offsetY, int totalRows, int totalCols);
		
		void readAscGridGDAL(string filename); //added by jiangjc
		void createAscGridGADL(string srcfilename,string filename);//added by jiangjc
		void writeAscGridGDAL(string filename);//added by jiangjc
		void readAscGrid(string filename);
		void readAscGrid(string filename, double x_min, double y_min, double x_max, double y_max);
		void writeAscGrid(string filename);
		void freeValues();
		double** values;
		double* valuesStorage;

	private:
		int totalCols;
		int totalRows;
		int subRows;
		int subCols;
		int startRow, startCol;
		double xCor, yCor;
		double xCor_sub,yCor_sub;
		double cellSize;
		double noDataVal;
		double stdvalue;//wtf

		string projection;//added by jiangjc
		double* pTransform; //added by jiangjc
		//for reading and writing
		string inAscFileName;
		string outAscFileName;

        //这个grid数据的最大值和最小值
		double m_dMax, m_dMin;
		void initGridStatistics();
};

#endif
