#if (defined _DEBUG) && (defined MSVC) && (defined VLD)
#include "vld.h"
#endif /* Run Visual Leak Detector during Debug */
// 1. Include the header of this cpp if existed, e.g., predictPoint_sfd.h
// 2. Include C/C++ system files.
// 3. Include Other libraries' .h files.
#include <mpi.h>  // mpich or openmpi
#ifdef SUPPORT_OMP
#include <omp.h>
#endif
#include "commonLib.h" // from TauDEM
#include "linearpart.h"
#include "createpart.h"
#include "initneighbor.h"
#include "utilities.h" // from UtilsClass of Lries-2415

// 4. Include other header files of current project if existed.
#include "AscGrid.h"
#include "common_func.h"
#include "sfdRegionSimilarity.h"


//#define MEMORY_SIZE 5000

using namespace std;

/*!
 * \brief Prediction based on neighbor defined by single flow direction model.
 *        Parallized based TauDEM parallization framework. Modified by Liangjun.
 * \param[in] fdmodel Flow direction model. 0 - D8, 1 - Dinf.
 * \param[in] flowdirf Full path of single flow direction.
 * \param[in] envfs Vector of full paths of environmental variables.
 * \param[in] xname Coordinate x
 * \param[in] yname Coordinate y
 * \param[in] attrname Property name
 * \param[in] freqnum Number of section for frequency calculation.
 * \param[in] outtype Output method. 0 - predict the whole study area
 *                                   1 - predict the validation point only
 *                                   2 - both 0 and 1
 * \param[in] validf Validation sample path.
 * \param[in] outf Output prediction raster.
 */
int predict_point_sfd(int fdmodel, char * flowdirf, vector<string> envfs, char * trainf,
                      string xname, string yname, string attrname, int freqnum, int outtype,
                      char *validf=NULL, char *outf=NULL)
{
    MPI_Init(NULL, NULL);
    {
        int rank, size;
        MPI_Comm_rank(MCW, &rank);
        MPI_Comm_size(MCW, &size);
        if (rank == 0) {
            cout << endl;
            cout << "Neighbor-based prediction method by flow direction model. Version 1.1" << endl;
            cout << endl;
            fflush(stdout);
            if (size > 1) {
                cout << "At the moment, you can only use '-n 1' or not use mpiexec!" << endl;
                MPI_Abort(MCW, 1);
                return 1;
            }
        }
        // begin timer
        double begint = MPI_Wtime();
        // read tiff header information using tiffIO from flow direction
        tiffIO *fdir_rst = NULL;
        if (fdmodel == 0) {
            fdir_rst = new tiffIO(flowdirf, SHORT_TYPE);
        } else if (fdmodel == 1) {
            fdir_rst = new tiffIO(flowdirf, FLOAT_TYPE);
        }
        long totalX = fdir_rst->getTotalX();
        long totalY = fdir_rst->getTotalY();
        double dx = fdir_rst->getdxA();
        double dy = fdir_rst->getdyA();
        if (rank == 0) {
            cout << "Row: " << totalY << ", Col: " << totalX << ", cellsize: " << dx << endl;
        }

        // read DEM data into partition
        tdpartition *fdir_part = NULL;
        fdir_part = CreateNewPartition(fdir_rst->getDatatype(), totalX, totalY, dx, dy, fdir_rst->getNodata());
        // get the size of current partition
        int nx = fdir_part->getnx();
        int ny = fdir_part->getny();
        int xstart, ystart;
        fdir_part->localToGlobal(0, 0, xstart, ystart); // calculate current partition's first cell's position
        fdir_part->savedxdyc(*fdir_rst);
        fdir_rst->read(xstart, ystart, ny, nx, fdir_part->getGridPointer()); // get the current partition's pointer

        // read parameters data into *partition
        int env_num = envfs.size();
        linearpart<float> *env_parts = new linearpart<float>[env_num];
        for (int num = 0; num < env_num; num++)
        {
            tiffIO paramsf(string_to_char(envfs[num]), FLOAT_TYPE);
            if (!fdir_rst->compareTiff(paramsf))
            {
                printf("File size do not match\n%s\n", envfs[num]);
                MPI_Abort(MCW, 5);
                return 1;
            }
            env_parts[num].init(totalX, totalY, dx, dy, MPI_FLOAT, *((float *)paramsf.getNodata()));
            paramsf.read(xstart, ystart, ny, nx, env_parts[num].getGridPointer());
        }
        double readt = MPI_Wtime(); // record reading time

        tdpartition *pred_part = NULL;
        if (outtype != 1) {
            // create empty partition to store new result
            pred_part = CreateNewPartition(FLOAT_TYPE, totalX, totalY, dx, dy, NODATA_VALUE);
        }

        // COMPUTING CODE BLOCK
        if (rank == 0) cout << "=============  Initial cells without inflow  =============" << endl;
        tdpartition *neighbor;
        neighbor = CreateNewPartition(SHORT_TYPE, totalX, totalY, dx, dy, NODATA_VALUE);
        // share information between ranks and set borders to zero
        fdir_part->share();
        neighbor->clearBorders();
        node temp;
        queue<node> que;
        int useOutlets = 0;
        long numOutlets = 0;
        int *outletsX = 0, *outletsY = 0;
        //Count the flow receiving neighbors and put node with no contributing neighbors on que
        if (fdmodel == 0) {
            initNeighborD8up(neighbor, fdir_part, &que, nx, ny, useOutlets, outletsX, outletsY, numOutlets);
        } else if (fdmodel == 1) {
            initNeighborDinfup(neighbor, fdir_part, &que, nx, ny, useOutlets, outletsX, outletsY, numOutlets);
        }
        int noflowincell_num_part = que.size();
        int noflowincell_num = 0;
        MPI_Allreduce(&noflowincell_num_part, &noflowincell_num, 1, MPI_INT, MPI_SUM, MCW);
        if (rank == 0) {
            cout << "Initial cells without inflow done, with a count of " << noflowincell_num << endl;
            cout << "=============  Trace downslope and store upstream cells and variables "
                "for each cell =============" << endl;
        }
        /*!
         * Map to store the nodes of upstream inflow cells
         * key is index of cell, i.e., row * nRows + col
         */
        map<int, vector<node> > upCellRC_map;
        bool finished_flag = false;
        int i, j, in, jn, k, cellidx, cellidx_n;
        int ig, jg, ing, jng; // global row and col
        short tempNeighbor;
        short fd8; // d8 flow direction
        float angle; // dinf flow direction
        while (!finished_flag) {
            while (!que.empty())
            {
                //Takes next node with no contributing neighbors
                temp = que.front();
                que.pop();
                i = temp.x;  // col
                j = temp.y;  // row
                cout << "remain que num: " << que.size() << ", current x: " << i << ", y: " << j << 
                    ", upCellRC_map: " << upCellRC_map.size() << endl;
                fdir_part->localToGlobal(i, j, ig, jg);
                cellidx = ig + jg * totalY;
                if (!fdir_part->hasAccess(i, j) || fdir_part->isNodata(i, j)) {
                    continue;
                }
                // Evaluate up slope inflow cells
                for (k = 1; k <= 8; k ++) {
                    in = i + d1[k];
                    jn = j + d2[k];
                    fdir_part->localToGlobal(in, jn, ing, jng);
                    cellidx_n = ing + jng * totalY;
                    if (!fdir_part->hasAccess(in, jn) || fdir_part->isNodata(in, jn)) {
                        continue;
                    }
                    if (fdmodel == 0) {
                        fdir_part->getData(in, jn, fd8);
                        if (abs(fd8 - k) != 4) {
                            continue;
                        }
                    }
                    else if (fdmodel == 1) {
                        fdir_part->getData(in, jn, angle);
                        double p = prop(angle, (k + 4) % 8, dx, dy);
                        /// this can be a parameter, which means how many inflow water
                        /// can be seen as the upstream of the current cell.
                        double thresh = 0.1;
                        if (p < thresh) {
                            continue;
                        }
                    }

                    vector<node> tmp_node_vec;
                    node node_n;
                    node_n.x = in;
                    node_n.y = jn;
                    if (upCellRC_map.find(cellidx_n) == upCellRC_map.end()) {
                        tmp_node_vec.push_back(node_n);
                        upCellRC_map.insert(make_pair(cellidx_n, tmp_node_vec));
                    }
                    else {
                        upCellRC_map.at(cellidx_n).push_back(node_n);
                    }
                    if (upCellRC_map.find(cellidx) == upCellRC_map.end()) {
                        vector<node> tmp_nodes;
                        upCellRC_map.insert(make_pair(cellidx, tmp_nodes));
                    }
                    vector<node> upnodes = upCellRC_map.at(cellidx_n);
                    for (vector<node>::iterator it = upnodes.begin(); it != upnodes.end(); it++) {
                        upCellRC_map.at(cellidx).push_back(*it);
                    }
                }
                // End evaluate up slope inflow cells
                // check and push to que the newest cells without inflow
                if (fdmodel == 0) {
                    fdir_part->getData(i, j, fd8);
                    in = i + d1[fd8];
                    jn = j + d2[fd8];
                    neighbor->addToData(in, jn, (short)-1);
                    if (fdir_part->isInPartition(in, jn) && neighbor->getData(in, jn, tempNeighbor) == 0) {
                        temp.x = in;
                        temp.y = jn;
                        que.push(temp);
                    }
                } 
                else {
                    fdir_part->getData(i, j, angle);
                    for (k = 1; k <= 8; k++) {
                        double p = prop(angle, k, dx, dy);
                        if (p < 0.) continue;
                        in = i + d1[k];
                        jn = j + d2[k];
                        neighbor->addToData(in, jn, (short)-1);
                        if (fdir_part->isInPartition(in, jn) && neighbor->getData(in, jn, tempNeighbor) == 0) {
                            temp.x = in;
                            temp.y = jn;
                            que.push(temp);
                        }
                    }
                }
                neighbor->addBorders();
                for (i = 0; i < nx; i++) {
                    if (neighbor->getData(i, -1, tempNeighbor) != 0 && 
                        neighbor->getData(i, 0, tempNeighbor) == 0) {
                        temp.x = i;
                        temp.y = 0;
                        que.push(temp);
                    }
                    if (neighbor->getData(i, ny, tempNeighbor) != 0 && 
                        neighbor->getData(i, ny - 1, tempNeighbor) == 0) {
                        temp.x = i;
                        temp.y = ny - 1;
                        que.push(temp);
                    }
                }
                neighbor->clearBorders();
                finished_flag = que.empty();
            }
        }
        int maxidx = 0;
        int maxcount = 0;
        for (map<int, vector<node> >::iterator it = upCellRC_map.begin(); it != upCellRC_map.end(); it++) {
            node cur_node;
            cur_node.x = it->first % totalY;
            cur_node.y = it->first / totalY;
            if (find(it->second.begin(), it->second.end(), cur_node) == it->second.end()) {
                it->second.push_back(cur_node);
            }
            if (it->second.size() > maxcount) {
                maxcount = it->second.size();
                maxidx = it->first;
            }
        }
        vector<node> maxupcells = upCellRC_map.at(maxidx);
        cout << "Total cells with inflow cells: " << upCellRC_map.size() << endl;
        cout << "The maximum inflow cells number is: " << maxcount << endl;
        for (vector<node>::iterator it = maxupcells.begin(); it != maxupcells.end(); it++) {
            pred_part->setData(it->x, it->y, 1.f);
        }
        for (j = 0; j < ny; j++) // rows
        {
            for (i = 0; i < nx; i++) // cols
            {
                ;
            }
        }
        // END COMPUTING CODE BLOCK
        double computet = MPI_Wtime(); // record computing time
        if (outtype != 1) {
            // create and write TIFF file
            float nodata = NODATA_VALUE;
            tiffIO destTIFF(outf, FLOAT_TYPE, &nodata, *fdir_rst);
            destTIFF.write(xstart, ystart, ny, nx, pred_part->getGridPointer());
        }
        double writet = MPI_Wtime(); // record writing time

        double dataRead, compute, write, total, tempd;
        dataRead = readt - begint;
        compute = computet - readt;
        write = writet - computet;
        total = writet - begint;

        MPI_Allreduce(&dataRead, &tempd, 1, MPI_DOUBLE, MPI_MAX, MCW);
        dataRead = tempd;
        MPI_Allreduce(&compute, &tempd, 1, MPI_DOUBLE, MPI_MAX, MCW);
        compute = tempd;
        MPI_Allreduce(&write, &tempd, 1, MPI_DOUBLE, MPI_MAX, MCW);
        write = tempd;
        MPI_Allreduce(&total, &tempd, 1, MPI_DOUBLE, MPI_MAX, MCW);
        total = tempd;

        if (rank == 0)
            printf("Processors: %d\nRead time: %f\nCompute time: %f\nWrite time: %f\nTotal time: %f\n",
            size, dataRead, compute, write, total);

        /// free memory
        delete fdir_part;
        delete[] env_parts;
        if (NULL != pred_part) delete pred_part;
    }
    MPI_Finalize();
    return 0;
}

/*!
 * \brief Original serial version by an.
 * \deprecated Replaced by parallized version.
 */
int predict_point_sfd_serial(string demPath, string sfdPath, vector<string> environLyrs, string trainSamplePath,
    string testSamplePath, string xName, string yName, string propertyName, string propertyPredPath, int envN)
{
    time_t beginTime, endTime, usedTime;
	beginTime = clock();

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
    /*vector<string> environLyrs;
    parseStr(string(environLyrsPath),'#',environLyrs);*/
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
	pfin = fopen(trainSamplePath.c_str(),"r");
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
    pfins = fopen(testSamplePath.c_str(), "r");
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
	
	int numOfTrainSamples = trainSampleList.size() - 1; // the first row in training sample file
	int columnsTrainSamples = trainSampleList[0].size();// 样点文件列数,sampleList[0]即样点文件中的第一行title
	int xIndexTrain = 0;
	int yIndexTrain = 0;
	int propertyIndexTrain = 0;
	
	int numOfTestSamples = testSampleList.size() - 1; // the first row in Testing sample file
	int columnsTestSamples = testSampleList[0].size();// 样点文件列数,sampleList[0]即样点文件中的第一行title
	int xIndexTest = 0;
	int yIndexTest = 0;
	int propertyIndexTest = 0;	
	
	int* trainSampleRows = new int[numOfTrainSamples];
	int* trainSampleCols = new int[numOfTrainSamples];
	double* trainProperty = new double[numOfTrainSamples];
	for (int j = 0; j < columnsTrainSamples; j++){
		//注意去除换行符
        if (strcmp(xName.c_str(), trimLineBreak(trainSampleList[0][j]).c_str()) == 0){
			xIndexTrain = j;
		}
        if (strcmp(yName.c_str(), trimLineBreak(trainSampleList[0][j]).c_str()) == 0){
			yIndexTrain = j;
		}
        if (strncmp(propertyName.c_str(), trimLineBreak(trainSampleList[0][j]).c_str(), strlen(propertyName.c_str())) == 0){
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
		if (strcmp(xName.c_str(), trimLineBreak(testSampleList[0][j]).c_str()) == 0){
			xIndexTest = j;
		}
		if (strcmp(yName.c_str(), trimLineBreak(testSampleList[0][j]).c_str()) == 0){
			yIndexTest = j;
		}
		if (strncmp(propertyName.c_str(), trimLineBreak(testSampleList[0][j]).c_str(), strlen(propertyName.c_str())) == 0){
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
	cout << "read sample files finished" << endl;
	
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
	
	// determine the neighborhood size of training samples
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
	
	//evaluate the soil property value and calculate the prediction uncertainty of testing samples
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


int main(int argc, char *argv[]){
    GDALAllRegister();
    if (argc < 10) {
        cout << "Input error, at least 9 parameters are required!" << endl;
    }
    int flowdir_model = atoi(argv[1]);
    char *flowDirectionPath = argv[2];
    char *environLyrsPath = argv[3];
    char *trainSamplePath = argv[4];
    char *xName = argv[5];
    char *yName = argv[6];
    char *propertyName = argv[7];
    int envN = atoi(argv[8]);
    /*
     * 0 - predict the whole study area
     * 1 - predict the validation point only
     * 2 - both 0 and 1
     */
    int outType = atoi(argv[9]);
    char *validationSamplePath = NULL;
    char *outPath = NULL;
    if (outType == 0) {
        outPath = argv[10];
    }
    else if (outType == 1) {
        validationSamplePath = argv[10];
    }
    else if (outType == 2) {
        if (argc < 11) {
            cout << "Validation sample file and the output path "
                "are both required when outType is 2." << endl;
        }
        validationSamplePath = argv[10];
        outPath = argv[11];
    }
    else {
        cout << "outType must be 0, 1, or 2!" << endl;
    }

    vector<string> environLyrs = SplitString(environLyrsPath, '#');
    cout << endl << "INPUT SUMMARY: " << endl;
    if (flowdir_model == 0) cout << "D8 ";
    else if (flowdir_model == 1) cout << "Dinf ";
    cout << "Flow Direction: " << flowDirectionPath << endl;
    cout << "Environment variables: " << endl;
    for (vector<string>::iterator it = environLyrs.begin(); it != environLyrs.end(); it++) {
        cout << "    " << *it << endl;
    }
    cout << "Training samples path: " << trainSamplePath << endl;
    cout << "XName: " << xName << ", YName: " << yName << ", Property: " << propertyName << endl;
    if (NULL != validationSamplePath)
        cout << "Validation samples path: " << validationSamplePath << endl;
    if (NULL != outPath)
        cout << "Output prediction: " << outPath << endl;

    predict_point_sfd(flowdir_model, flowDirectionPath, environLyrs, trainSamplePath, xName, yName,
        propertyName, envN, outType, validationSamplePath, outPath);
}