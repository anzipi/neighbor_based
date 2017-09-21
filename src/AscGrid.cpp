#include "AscGrid.h"
#include <math.h>
#include <iostream>
#include<cmath>
#define VERY_SMALL 0.000001
using namespace std;

AscGrid::AscGrid(){

	totalCols = 0;
    totalRows = 0;
	subCols = -9999;
	subRows = -9999;
    xCor = 0; //top left x
    yCor = 0; 
	xCor_sub = -9999;
	yCor_sub = -9999;
    cellSize = 1;
    noDataVal = -9999;
	values = NULL;
	this->pTransform = new double[6];
    initGridStatistics();
}

AscGrid::~AscGrid(){

}

void AscGrid::freeValues(){
	if(values != NULL){
		delete[] values;
	}
	if(valuesStorage != NULL)
		delete[] valuesStorage;
}

AscGrid::AscGrid(int col, int row, double xCor,
	double yCor, double cell, double noData,double** values){
		this->totalCols = col;
		this->totalRows = row;
		subCols = -9999;
		subRows = -9999;
		xCor_sub = -9999;
		yCor_sub = -9999;
		this->xCor = xCor;
		this->yCor = yCor;
		this->cellSize = cell;
		this->noDataVal = noData;
		this->values = values;
		this->pTransform = new double[6];
		initGridStatistics();
}


void AscGrid::setHeadInfo(int col, int row, double xCor,
            double yCor, double cell, double noData){
                this->totalCols = col;
                this->totalRows = row;
                this->xCor = xCor;
                this->yCor = yCor;
                this->cellSize = cell;
                this->noDataVal = noData;
}

void AscGrid::setValues(double** newValues){
	this->values = newValues;
}


void AscGrid::setCellValue(int row, int col, double value){
	values[row][col] = value;
}

void AscGrid::setOutputFileName(string filename){
	outAscFileName = filename;
}

double AscGrid::getValue(int row, int col){
	return values[row][col];
}

double** AscGrid::getValMatrix(){
	return values;
}

double AscGrid::getMax(){
	m_dMax = numeric_limits<double>::min();
	if (subRows != -9999 && subCols != -9999){
		for (int i = 0; i < subRows; i++)
		{
			 for (int j = 0; j < subCols; j++)
			 {
				 if (values[i][j] > m_dMax && abs(values[i][j] - noDataVal)>VERY_SMALL){
					 m_dMax = values[i][j];
				 }
			  }
		 }
	}else{
		for (int i = 0; i < totalRows; i++)
		{
			 for (int j = 0; j < totalCols; j++)
			 {
				 if (values[i][j] > m_dMax && abs(values[i][j] - noDataVal)>VERY_SMALL){
					 m_dMax = values[i][j];
				 }
			  }
		 }
	}
	return m_dMax;
}

double AscGrid::getMin()
{
	m_dMin = numeric_limits<double>::max();
	if (subRows != -9999 && subCols != -9999){
		for (int i = 0; i < subRows; i++)
		{
			for (int j = 0; j < subCols; j++)
			{
				if (values[i][j] < m_dMin && abs(values[i][j] - noDataVal) > VERY_SMALL)
				{
					m_dMin = values[i][j];
				}
			}
		}
	}else{
		for (int i = 0; i < totalRows; i++)
		{
			for (int j = 0; j < totalCols; j++)
			{
				if (values[i][j] < m_dMin && abs(values[i][j] - noDataVal) > VERY_SMALL)
				{
					m_dMin = values[i][j];
				}
			}
		}
	}
    //cout<<"min"<<m_dMin<<endl;
    return m_dMin;
}

bool AscGrid::isValidCell(int Row, int Col)
{
    if (Col < 0 || Col >= totalCols || Row < 0 || Row >= totalRows)
        return false;
    return true;
}

double * AscGrid::readAscGridGDAL_Block(string filename, int offsetX, int offsetY, int totalRows, int totalCols){   // added by jiangjc

	if (filename != "")
	{

		GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly) ;
		if( poDatasetsrc == NULL )
		{
			cout<<"[ERROR] data file is not open correct"<<endl;
			exit(1);
		}
		GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand(1);


		valuesStorage = new double[totalRows*totalCols];

		poBandsrc->RasterIO( GF_Read, offsetX, offsetY, totalCols, totalRows, valuesStorage, totalCols, totalRows, GDT_Float64, 0, 0 );//gdt_float64

        //cout<<"valuesStorage[0]"<<valuesStorage[0]<<endl;
		//close gdal
		if (poDatasetsrc != NULL)
		{
			GDALClose((GDALDatasetH)poDatasetsrc);
			poDatasetsrc = NULL;
		}

	}
	else
	{
		cerr<<"need to specify the file name"<<endl;
	}

	return valuesStorage;
}



void AscGrid::readAscGridGDAL(string filename){// added by jiangjc
	if (filename != "")
	{

		GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly) ;
		if( poDatasetsrc == NULL )
		{
			cout<<"[ERROR] data file is not open correct"<<endl;
			exit(1);
		}

		poDatasetsrc->GetGeoTransform(pTransform);
		projection = poDatasetsrc->GetProjectionRef();
		GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand(1);
		double *pdfMin=new double,*pdfMax=new double,*pdfMean=new double,*pdfStdDev=new double;
        poBandsrc->ComputeStatistics(false,pdfMin,pdfMax,pdfMean,pdfStdDev,NULL,NULL);
        stdvalue = pdfStdDev[0];//wtf
		//read in the file header information
		totalRows = poBandsrc->GetYSize();
		totalCols = poBandsrc->GetXSize();
		xCor = pTransform[0];
		yCor = pTransform[3] - pTransform[1]*totalRows;// 由左上角变为左下角
        cellSize = pTransform[1];// w-e pixel resolution
		noDataVal = poBandsrc->GetNoDataValue();
		//read in the body
		values = new double*[totalRows];
		valuesStorage = new double[totalRows*totalCols];
		for (int rowId = 0; rowId < totalRows; rowId++)
		{
			values[rowId] = &valuesStorage[rowId*totalCols];
		}

		poBandsrc->RasterIO( GF_Read, 0, 0, totalCols, totalRows, valuesStorage, totalCols, totalRows, GDT_Float64, 0, 0 );//wtf   GDT_Float64
		//close gdal
		if (poDatasetsrc != NULL)
		{
			GDALClose((GDALDatasetH)poDatasetsrc);
			poDatasetsrc = NULL;
		}
	}
	else{
		cerr<<"need to specify the file name"<<endl;
	}

}
void AscGrid::readAscGrid(string filename){
	ifstream AscReader; // ifstream: 读操作（输入）的文件类
    if (filename != ""){
        inAscFileName = filename;
		AscReader.open(inAscFileName.c_str(),ios::in);
		if (!AscReader.is_open()){
			cerr<<"cannot open file "<< filename << endl;
		}else{
			//read in the file header information
			string tmp;
			for (int i = 0; i < 6; i++)
			{
				getline(AscReader,tmp);//按行读取数据，AscReader表示一个输入流，此语句表示将从输入流中读取的数据放入tmp中
				int j = tmp.find_last_of(' ');//查找字符串中最后一个出现空格的位置
				int end = tmp.find('\n');// substr主要功能是复制子字符串，要求从指定位置开始，并具有指定的长度
				string extractNum = (tmp.substr(j+1,end)).c_str();//c语言中没有string类型，必须通过string类对象的成员函数c_str()把string 对象转换成c中的字符串样式
				switch (i){
					case 0:
						totalCols = atoi(extractNum.c_str());// atoi字符串转整型
						break;
					case 1:
						totalRows = atoi(extractNum.c_str());
						break;
					case 2:
						xCor = atof(extractNum.c_str()); //字符串转float
						break;
					case 3:
						yCor = atof(extractNum.c_str());
						break;
					case 4:
						cellSize = atof(extractNum.c_str());
						break;
					case 5:
						noDataVal = atof(extractNum.c_str());
						break;
				}
				startRow = 0;
				startCol = 0;
			}// end of reading the header
			//readin the body
			if (totalRows != 0){
				values = new double*[totalRows];
				valuesStorage = new double[totalRows * totalCols];
				for (int rowId = 0; rowId < totalRows; rowId++)
				{
					values[rowId] = &valuesStorage[rowId * totalCols];
				}
				for (int rowId = 0; rowId < totalRows; rowId++)
				{
					string tmp;
					vector<string> stringVal;//to store values in a row
					getline(AscReader,tmp);
					string::size_type pos1, pos2;
					pos1 = 0;
					pos2 = tmp.find_first_of(' ',pos1); //从pos1开始查找在字符串中第一个出现空格的位置
					while(string::npos > pos2){
						stringVal.push_back(tmp.substr(pos1,pos2-pos1));
						pos1 = pos2 + 1;
						pos2 = tmp.find_first_of(' ',pos1);
					}
					stringVal.push_back(tmp.substr(pos1));
					for (int j = 0; j < totalCols; j++)
					{
						this->values[rowId][j] = atof(stringVal[j].c_str());
						double aa = values[rowId][j];
					}
				}

			}
			else
			{
				this->values = NULL;
			}
			AscReader.close();
		}
    }
	else{
		cerr<<"need to specify the file name"<<endl;
	}
}

void AscGrid::readAscGrid(string filename, double x_min, double y_min, double x_max, double y_max){
	ifstream AscReader;
	if (filename != "")
    {
        inAscFileName = filename;
		AscReader.open(inAscFileName.c_str(),ios::in);

		if (!AscReader.is_open()){
			cerr<<"cannot open file "<< filename << endl;
		}else{
			//read in the file header information
			char tmp[100];
			for (int ih = 0; ih < 6; ih++)
			{
				AscReader.getline(tmp,100);
				int j = string(tmp).find_last_of(' ');
				int end = string(tmp).length();
				string extractNum = string(tmp).substr(j+1,end-j-1);
				switch (ih)
				{
					case 0:
						totalCols = atoi(extractNum.c_str());
						break;
					case 1:
						totalRows = atoi(extractNum.c_str());
						break;
					case 2:
						xCor = atof(extractNum.c_str());
						break;
					case 3:
						yCor = atof(extractNum.c_str());
						break;
					case 4:
						cellSize = atof(extractNum.c_str());
						break;
					case 5:
						noDataVal = atof(extractNum.c_str());
						break;
				}

			}// end of reading the header
			subCols = (int) ((x_max - x_min)/cellSize);
			subRows = (int) ((y_max - y_min)/cellSize);
			startRow = totalRows - (int) ((y_max - yCor)/cellSize);
			startCol = (int) ((x_min - xCor)/cellSize);
			xCor_sub = x_min;
			yCor_sub = y_min;


			//readin the body
			if (subRows != 0){
				values = new double*[subRows];
				for (int ir = 0; ir < subRows; ir++)
				{
					values[ir] = new double[subCols];
				}
				//skip the first subRows lines
				for (int ix = 0; ix < startRow; ix++)
				{
					string skip="";
					getline(AscReader,skip);
				}

				for(int iii = startRow; iii < totalRows; iii ++){
					string tmp1;
					vector<string> stringVal;//to store values in a row
					getline(AscReader,tmp1);
					string::size_type pos1, pos2;
					pos1 = 0;
					pos2 = tmp1.find_first_of(' ',pos1);
					while(string::npos > pos2){
						stringVal.push_back(tmp1.substr(pos1,pos2-pos1));
						pos1 = pos2 + 1;
						pos2 = tmp1.find_first_of(' ',pos1);
					}
					//stringVal.push_back(tmp1.substr(pos1));
					for (int jjj = startCol; jjj < totalCols; jjj++)
					{
						values[iii-startRow][jjj-startCol] = atof(stringVal[jjj].c_str());
					}
				}
			}else{
				this->values = NULL;
			}
			AscReader.close();
		}

    }else{
		cerr<<"need to specify the file name"<<endl;
	}
}

void AscGrid::createAscGridGADL(string srcfilename,string filename){
    if(filename !=""){
		// 读取投影信息
		GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen( srcfilename.c_str(), GA_ReadOnly );
		if( poDatasetsrc == NULL)
		{
			//do something
			cout<<"[ERROR] data file is not open correct"<<endl;
			exit(1);
		}
		poDatasetsrc->GetGeoTransform(pTransform);
		projection = poDatasetsrc->GetProjectionRef();
		if (poDatasetsrc != NULL)
		{
			//
			//poDatasetsrc->~GDALDataset();
			GDALClose((GDALDatasetH)poDatasetsrc);
			poDatasetsrc = NULL;
		}
		//新建tiff文件
		GDALDriver* poDriver = NULL;
		string format = "GTiff";
		poDriver = GetGDALDriverManager()->GetDriverByName(format.c_str());
		if (poDriver == NULL)
		{
			cout<<"[ERROR] poDriver is NULL."<<endl;
			exit(1);
		}
		GDALDataset* poDataset = poDriver->Create(filename.c_str(), totalCols, totalRows, 1, GDT_Float64, NULL);//wtf   64
		poDataset->SetGeoTransform( pTransform );
		poDataset->SetProjection(projection.c_str());
		delete []this->pTransform;

		if (poDataset == NULL)
	    {
			 //do something
			 cout<<"[ERROR] poDatasetdest is NULL"<<endl;
			 exit(1);
	    }
		if (poDataset != NULL)
	    {
			//
			poDataset->FlushCache();
			//poDataset->~GDALDataset();
			GDALClose((GDALDatasetH)poDataset);
			poDataset = NULL;
	    }
	}else{
		cerr<<"need to specify the output file name"<<endl;
	}
}
void AscGrid::writeAscGridGDAL(string filename){
    GDALDataset* poDataset = NULL;
	poDataset = (GDALDataset *) GDALOpen( filename.c_str(), GA_Update );
	if( poDataset == NULL /*检查是否正常打开文件*/)
	{
		//do something
		cout<<"[ERROR] data file is not open correct"<<endl;
		exit(1);
	}
	GDALRasterBand*	poBanddest = poDataset->GetRasterBand(1);
	if (poBanddest == NULL)
	{
		//do something
		cout<<"[ERROR] poBanddest is NULL"<<endl;
		exit(1);
	}

	poBanddest->SetNoDataValue(noDataVal);

	valuesStorage = new double[totalRows*totalCols];
	for (int i = 0;i<totalRows;i++ )
	{
		for (int j=0;j<totalCols; j++)
		{
			valuesStorage[i*totalCols + j] = values[i][j];
		}
	}
	poBanddest->RasterIO(GF_Write, 0,0, totalCols, totalRows, valuesStorage, totalCols, totalRows, GDT_Float64, 0, 0);//wtf 64
	poDataset->FlushCache();

	if (poDataset != NULL)
	{
		//poDataset->~GDALDataset();
		GDALClose((GDALDatasetH)poDataset);
		poDataset = NULL;
	}

}
void AscGrid::writeAscGrid(string filename){
	ofstream AscWriter;
	if (filename != "")
    {
        outAscFileName = filename;
		AscWriter.open(outAscFileName.c_str(),ios::out);
		if(!AscWriter.is_open()){
			cerr<<"fail to open file "<< filename <<endl;
		}else{
			if (subCols != -9999)
				totalCols = subCols;
			if (subRows != -9999)
				totalRows = subRows;
			if (xCor_sub != -9999)
				xCor = xCor_sub;
			if (yCor_sub != -9999)
				yCor = yCor_sub;
			AscWriter<<"ncols "<<this->totalCols<<endl;
			AscWriter<<"nrows "<< this->totalRows<<endl;
			AscWriter<<"xllcorner "<< this->xCor<<endl;
			AscWriter<<"yllcorner "<<this->yCor<<endl;
			AscWriter<<"cellsize "<<this->cellSize<<endl;
			AscWriter<<"NODATA_value "<<this->noDataVal<<endl;
			for (int i = 0; i < totalRows; i++)
			{
				for (int j = 0; j < totalCols - 1; j++)
				{
					AscWriter<<this->values[i][j];
					AscWriter<<' ';
				}
				AscWriter<<this->values[i][totalCols - 1];
				AscWriter<<'\n';
			}
			AscWriter.close();
		}
    }else{
		cerr<<"need to specify the output file name"<<endl;
	}
}

int AscGrid::getNumOfCols(){
	if (subCols != -9999)
		return subCols;
	return totalCols;
}

int AscGrid::getNumOfRows(){
	if (subRows != -9999)
		return subRows;
	return totalRows;
}

double AscGrid::getXCor(){
	if (xCor_sub != -9999)
		return xCor_sub;
	return xCor;
}

double AscGrid::getYCor(){
	if(yCor_sub != -9999)
		return yCor_sub;
	return yCor;
}

double AscGrid::getCellSize(){
	return cellSize;
}

double AscGrid::getNodaVal(){
	return noDataVal;
}

void AscGrid::initGridStatistics()
{
    m_dMax = noDataVal;
    m_dMin = noDataVal;
	
}
double AscGrid::getStdValue()//wtf
{
    return stdvalue;
}