#if (defined _DEBUG) && (defined MSVC) && (defined VLD)
#include "vld.h"
#endif /* Run Visual Leak Detector during Debug */
// 1. Include the header of this cpp if existed, e.g., predictPoint_sfd.h
// 2. Include C/C++ system files.
// 3. Include Other libraries' .h files.
#include <mpi.h>  // mpich or openmpi
#ifdef SUPPORT_OMP
#include <omp.h>
#include <typeinfo.h>//an
#endif
#include "commonLib.h" // from TauDEM
#include "linearpart.h"
#include "createpart.h"
#include "initneighbor.h"
#include "utilities.h" // from UtilsClass of Lries-2415

// 4. Include other header files of current project if existed.
//#include "flowInCells.h"
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
 * \param[in] predf Output prediction raster.
 *\param[in] uncerf uncertainty raster.
 */

int predict_point_sfd(int fdmodel, char *flowdirf, char *streamf, vector<string> envfs,
    char *trainf,string xname, string yname, string attrname, int freqnum, int outtype,
    char *validf= NULL, char *predf = NULL, char *uncerf = NULL)
{
    MPI_Init(NULL, NULL);
    {
        int rank, size;
        MPI_Comm_rank(MCW, &rank);
        MPI_Comm_size(MCW, &size);
        if (rank == 0) {
            cout << endl;
            cout << "Neighbor-based prediction method by flow direction model. Version 1.2" << endl;
            cout << endl << "INPUT SUMMARY: " << endl;
            if (fdmodel == 0) cout << "D8 ";
            else if (fdmodel == 1) cout << "Dinf ";
            cout << "Flow Direction: " << flowdirf << endl;
            cout << "Stream: " << streamf << endl;
            cout << "Environment variables: " << endl;
            for (vector<string>::iterator it = envfs.begin(); it != envfs.end(); it++) {
                cout << "    " << *it << endl;
            }
            cout << "Training samples path: " << trainf << endl;
            cout << "XName: " << xname << ", YName: " << yname << ", Property: " << attrname << endl;
            if (NULL != validf)
                cout << "Validation samples path: " << validf << endl;
            if (NULL != predf)
                cout << "Output prediction: " << predf << endl;
            if (NULL != uncerf)
                cout << "Output uncertainty: " << uncerf << endl;
            fflush(stdout);
            //if (size > 1) {
            //    cout << "At the moment, you can only use '-n 1' or not use mpiexec!" << endl;
            //    MPI_Abort(MCW, 1);
            //    return 1;
            //}
        }
        // begin timer
        double begint = MPI_Wtime();
        // read tiff header information using tiffIO from flow direction
        tiffIO *fdir_rst = NULL;
        if (!FileExists(flowdirf)){
            printf("File not existed!\n%s\n", flowdirf);
            MPI_Abort(MCW, 5);
            return 1;
        }
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

        // read flow direction data into partition
        tdpartition *fdir_part = NULL;
        fdir_part = CreateNewPartition(fdir_rst->getDatatype(), totalX, totalY, dx, dy, fdir_rst->getNodata());
        // get the size of current partition
        int nx = fdir_part->getnx();
        int ny = fdir_part->getny();
        int xstart, ystart;
        fdir_part->localToGlobal(0, 0, xstart, ystart); // calculate current partition's first cell's position
        fdir_part->savedxdyc(*fdir_rst);
        fdir_rst->read(xstart, ystart, ny, nx, fdir_part->getGridPointer()); // get the current partition's pointer

        // read stream data into partition
        linearpart<short> *stream_part = new linearpart<short>;
        if (!FileExists(streamf)){
            printf("File not existed!\n%s\n", streamf);
            MPI_Abort(MCW, 5);
            return 1;
        }
        tiffIO stream_rst(streamf, SHORT_TYPE);
        if (!fdir_rst->compareTiff(stream_rst))
        {
            printf("File size do not match\n%s\n", streamf);
            MPI_Abort(MCW, 5);
            return 1;
        }
        stream_part->init(totalX, totalY, dx, dy, MPI_SHORT, *((short *)stream_rst.getNodata()));
        stream_rst.read(xstart, ystart, ny, nx, stream_part->getGridPointer());
        // read parameters data into partition
        int env_num = envfs.size();
        linearpart<float> *env_parts = new linearpart<float>[env_num];
        for(int num = 0; num < env_num; num++)
        {
            if (!FileExists(envfs[num])){
                cout << "File not existed: " << envfs[num] << endl;
                MPI_Abort(MCW, 5);
                return 1;
            }
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

        // COMPUTING CODE BLOCK
        if (rank == 0) cout << "=============  Construct in flow relationships  =============" << endl;
        /*!
         * Map of the ``Cell`` instance of each cell in current partition
         * key is global index of cell, i.e., row * nCols + col => jglobal * totalX + iglobal 
         * value is the ``Cell`` instance
        */
        map<int, Cell*> cells_map;
        /*!
         * Map of the cell index to position index of global data
         * key is global index of cell
         * value is the position index excluding NoData
         */
        map<int, int> position_map;
        int i, j, k, ig, jg;
        int in, jn;
        short tmpstream;
        short fd8; // d8 flow direction
        float angle; // dinf flow direction

        vector<int> globIndex_rank;
        /// get global index of valid cells
        for (int j = 0; j < ny; j++) { // row
            for (int i = 0; i < nx; i++) { // col
                fdir_part->localToGlobal(i, j, ig, jg);
                int ginx = jg * totalX + ig;
                // extract valid cells
                if (!fdir_part->hasAccess(i, j) || fdir_part->isNodata(i, j)) continue;
                if (!stream_part->hasAccess(i, j) || stream_part->getData(i, j, tmpstream) > 0) continue;
                bool env_nodata = false;
                for (k = 0; k < env_num; k++) {
                    if (!env_parts[k].hasAccess(i, j) || env_parts[k].isNodata(i, j)) {
                        env_nodata = true;
                        break;
                    }
                }
                if (env_nodata) continue;
                globIndex_rank.push_back(ginx);
            }
        }
        vector<int>(globIndex_rank).swap(globIndex_rank);
        int validCellNum_rank = globIndex_rank.size();

        // create arrays need to be sent and gathered.
        int *globIndex_ra = new int[validCellNum_rank];
        float *flowdir_ra = new float[validCellNum_rank];
        float **env_ra = new float*[env_num];
        for (k = 0; k < env_num; k++) {
            env_ra[k] = new float[validCellNum_rank];
        }
        float * minEnv = new float[env_num];
        float * maxEnv = new float[env_num];
        for (k = 0; k < env_num; k++){
            minEnv[k] = MAXFLOAT;
            maxEnv[k] = MINFLOAT;
        }
        for (i = 0; i < validCellNum_rank; i++) {
            globIndex_ra[i] = globIndex_rank[i];
            int cur_ig = globIndex_ra[i] % totalX;
            int cur_jg = globIndex_ra[i] / totalX;
            int cur_il, cur_jl;
            fdir_part->globalToLocal(cur_ig, cur_jg, cur_il, cur_jl);
            if (fdmodel == 0) {
                fdir_part->getData(cur_il, cur_jl, fd8);
                flowdir_ra[i] = (float)fd8;
            } 
            else if (fdmodel == 1) {
                fdir_part->getData(cur_il, cur_jl, flowdir_ra[i]);
            }
            for (k = 0; k < env_num; k++) {
                env_parts[k].getData(cur_il, cur_jl, env_ra[k][i]);
                if (minEnv[k] > env_ra[k][i]){
                    minEnv[k] = env_ra[k][i];
                }
                if (maxEnv[k] < env_ra[k][i]){
                    maxEnv[k] = env_ra[k][i];
                }
            }
        }
        for (k = 0; k < env_num; k++){
            MPI_Allreduce(MPI_IN_PLACE, &minEnv[k], 1, MPI_FLOAT, MPI_MIN, MCW);
            MPI_Allreduce(MPI_IN_PLACE, &maxEnv[k], 1, MPI_FLOAT, MPI_MAX, MCW);
            //cout << "Env No:" << k << ", " << minEnv[k] << ", " << maxEnv[k] << endl;
        }

        int *localCellCount = new int[size]; /// cell numbers in each partition
        int validCellNum_all;
        MPI_Allreduce(&validCellNum_rank, &validCellNum_all, 1, MPI_INT, MPI_SUM, MCW);
        MPI_Allgather(&validCellNum_rank, 1, MPI_INT, localCellCount, 1, MPI_INT, MCW);
        //cout << "rank: " << rank << ", valid cell: " << validCellNum_rank << ", all: " << validCellNum_all << endl;
        int *displs = new int[size];
        displs[0] = 0;
        for (i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + localCellCount[i - 1];
        }

        // create arrays to store data of entire area.
        int *validCellIndexes = new int[validCellNum_all];
        int *globIndex_all = new int[validCellNum_all];
        float *flowdir_all = new float[validCellNum_all];
        float **env_all = new float*[env_num];
        for (k = 0; k < env_num; k++) {
            env_all[k] = new float[validCellNum_all];
        }

        MPI_Allgatherv(globIndex_ra, validCellNum_rank, MPI_INT, validCellIndexes, localCellCount, displs, MPI_INT, MCW);
        MPI_Allgatherv(flowdir_ra, validCellNum_rank, MPI_FLOAT, flowdir_all, localCellCount, displs, MPI_FLOAT, MCW);
        for (k = 0; k < env_num; k++) {
            MPI_Allgatherv(env_ra[k], validCellNum_rank, MPI_FLOAT, env_all[k], localCellCount, displs, MPI_FLOAT, MCW);
        }
        // release memory on rank
        globIndex_rank.clear();
        Release1DArray(globIndex_ra);
        Release1DArray(flowdir_ra);
        Release2DArray(env_num, env_ra);

        // Construct cells_map and position_map
        for (i = 0; i < validCellNum_all; i++) {
            position_map.insert(make_pair(validCellIndexes[i], i));
            Cell* curcell = new Cell(validCellIndexes[i], flowdir_all[i]);
            for (k = 0; k < env_num; k++) {
                curcell->addEnvValue(env_all[k][i]);
            }
            cells_map.insert(make_pair(validCellIndexes[i], curcell));
        }
        // Add in flow cells for each cell
        for (map<int, Cell*>::iterator it = cells_map.begin(); it != cells_map.end(); it++) {
            int cur_ig = it->first % totalX;
            int cur_jg = it->first / totalX;
            if (fdmodel == 0) {
                int curd8 = (int)it->second->flowDirection();
                in = cur_ig + d1[curd8];
                jn = cur_jg + d2[curd8];
                int cur_idx = jn * totalX + in;
                if (cells_map.find(cur_idx) == cells_map.end()) continue;
                cells_map.at(cur_idx)->addUpCell(it->second);
            } 
            else if (fdmodel == 1) {
                angle = it->second->flowDirection();
                for (k = 1; k <= 8; k++) {
                    double p = prop(angle, k, dx, dy);
                    double thresh = 0.01;
                    if (p < thresh) continue;
                    in = cur_ig + d1[k];
                    jn = cur_jg + d2[k];
                    int cur_idx = jn * totalX + in;
                    if (cells_map.find(cur_idx) == cells_map.end()) continue;
                    cells_map.at(cur_idx)->addUpCell(it->second);
                }
            }
        }
        // release memory on rank
        Release1DArray(globIndex_all);
        Release1DArray(flowdir_all);
        Release2DArray(env_num, env_all);

        // Read training samples
        vector<double> attrs_train;	
        int num_train = 0;
        int *train_idx;
        double *train_attrs;
        if (rank == 0) {
            cout << "============= Read training samples  =============" << endl;
            vector<double> xTrain;
            vector<double> yTrain;
            if (!FileExists(trainf)){
                cout << "File not existed: " << trainf << endl;
                MPI_Abort(MCW, 5);
                return 1;
            }
            readSamples(trainf, xname, yname, attrname, xTrain, yTrain, attrs_train);
            num_train = xTrain.size();
            train_idx = new int[num_train];
            train_attrs = new double[num_train];
            for (int train = 0; train < num_train; train++){
                int row, col;
                fdir_rst->geoToGlobalXY(xTrain[train], yTrain[train], col, row);
                int tmpidx = row * totalX + col;
                if (cells_map.find(tmpidx) == cells_map.end()) {
                    continue;
                }
                train_idx[train] = tmpidx;
                train_attrs[train] = attrs_train[train];
            }
            cout << "============= Read training samples finished =============" << endl;
        }
        MPI_Bcast(&num_train, 1, MPI_INT, 0, MCW);
        if (rank != 0) {
            train_idx = new int[num_train];
            train_attrs = new double[num_train];
        }
        MPI_Bcast(train_idx, num_train, MPI_INT, 0, MCW);
        MPI_Bcast(train_attrs, num_train, MPI_DOUBLE, 0, MCW);

        // NOW, each rank (processor) has a some cells_map!
        if (rank == 0) {
            cout << "============= calculate neighborhood feature for training samples  =============" << endl;
        }
        for (i = 0; i < num_train; i++) {
            Cell* cur_cell = cells_map.at(train_idx[i]);
            cur_cell->setAttribute(train_attrs[i]);
            cur_cell->setUnCertainty(0.);
            if (!cur_cell->frequencyStats(cells_map, env_num, minEnv, maxEnv, freqnum)) {
                cout << "Calculate frequency of environment variables failed on " << train_idx[i] << endl;
                MPI_Abort(MCW, 6);
                return 1;
            }
            //cur_cell->printFrequency();
        }
        //for (i = 0; i < num_train; i++) {
        //    cout << train_idx[i] << ": " << cells_map.at(train_idx[i])->getAttribute() << endl;
        //}

        tdpartition *pred_part = NULL;
        tdpartition *uncer_part = NULL;
        if (outtype != 1) {
            if (rank == 0) {
                cout << "============= Predict on the whole area  =============" << endl;
            }
            //calculate predictive value for pred_part
            pred_part = CreateNewPartition(FLOAT_TYPE, totalX, totalY, dx, dy, NODATA_VALUE);
            uncer_part = CreateNewPartition(FLOAT_TYPE, totalX, totalY, dx, dy, NODATA_VALUE);

            for (int j = 0; j < ny; j++){ // row
                for (int i = 0; i < nx; i++){ // col
                    int iglob, jglob;
                    pred_part->localToGlobal(i, j, iglob, jglob);
                    int pixel_idx = jglob * totalX + iglob;
                    if (cells_map.find(pixel_idx) == cells_map.end()) {
                        continue;
                    }
                    Cell *cur_cell = cells_map.at(pixel_idx);
                    //cout << rank << ": " << pixel_idx << ": ";
                    bool istrain = false;
                    for (k = 0; k < num_train; k++) {
                        if (train_idx[k] == pixel_idx) {
                            istrain = true;
                            break;
                        }
                    }
                    if (!istrain) {
                        cur_cell->frequencyStats(cells_map, env_num, minEnv, maxEnv, freqnum);
                        cur_cell->predictProperty(cells_map, env_num, minEnv, maxEnv, num_train, train_idx);
                    }
                    pred_part->setData(i, j, cur_cell->getAttribute());
                    uncer_part->setData(i, j, cur_cell->getUnCertainty());
                    //cout << cur_cell->getAttribute() << ", " << cur_cell->getUnCertainty() << endl;
                }
            }
        }
        // END COMPUTING CODE BLOCK

        double computet = MPI_Wtime(); // record computing time
        if (outtype != 1) {
            // create and write TIFF file to output predictive soil property
            float nodata = NODATA_VALUE;
            tiffIO predTIFF(predf, FLOAT_TYPE, &nodata, *fdir_rst);
            predTIFF.write(xstart, ystart, ny, nx, pred_part->getGridPointer());

			// create and write TIFF file to output prediction uncertainty			
			tiffIO uncerTIFF(uncerf, FLOAT_TYPE, &nodata, *fdir_rst);
			uncerTIFF.write(xstart, ystart, ny, nx, uncer_part->getGridPointer());
			//cout << "finished" << endl;
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
        delete stream_part;
        delete[] env_parts;
        if (NULL != pred_part) delete pred_part;
    }
    MPI_Finalize();
    return 0;
}

int main(int argc, char *argv[]){
    GDALAllRegister();
    if (argc < 11) {
        cout << "Input error, at least 10 parameters are required!" << endl;
    }
    int flowdir_model = atoi(argv[1]);
    char *flowDirectionPath = argv[2];
    char *streamPath = argv[3];
    char *environLyrsPath = argv[4];
    char *trainSamplePath = argv[5];
    char *xName = argv[6];
    char *yName = argv[7];
    char *propertyName = argv[8];
    int envN = atoi(argv[9]);
    /*
     * 0 - predict the whole study area
     * 1 - predict the validation point only
     * 2 - both 0 and 1
     */
    int outType = atoi(argv[10]);
    char *validationSamplePath = NULL;
    char *predPath = NULL;
	char *uncerPath = NULL;
    if (outType == 0) {
        predPath = argv[11];
		uncerPath = argv[12];
    }
    else if (outType == 1) {
        validationSamplePath = argv[11];
    }
    else if (outType == 2) {
        if (argc < 11) {
            cout << "Validation sample file and the output path "
                "are both required when outType is 2." << endl;
        }
        validationSamplePath = argv[11];
        predPath = argv[12];
		uncerPath = argv[13];
    }
    else {
        cout << "outType must be 0, 1, or 2!" << endl;
    }

    vector<string> environLyrs = SplitString(environLyrsPath, '#');
    
    predict_point_sfd(flowdir_model, flowDirectionPath, streamPath, environLyrs,
        trainSamplePath, xName, yName,propertyName, envN, outType,
        validationSamplePath, predPath, uncerPath);
	//system("pause");
	return 0;

}