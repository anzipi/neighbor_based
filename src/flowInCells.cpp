#include <queue>
#include <deque>
#include <algorithm>
#include "flowInCells.h"
#include <iostream>
#include "utilities.h"
#include "common_func.h"
using namespace std;

Cell::Cell(int idx, float flowdir) : 
m_globalIndex(idx), m_upCellsNum(0), m_flowdir(flowdir), m_envNum(-1), m_intervalNum(-1),
m_predictAttr(NODATA_VALUE), m_uncertainty(NODATA_VALUE), m_envFrequency(NULL){
    
}
Cell::~Cell() {

}
void Cell::addUpCell(Cell* up_cell) {
    m_upCells.push_back(up_cell);
    m_upCellsNum += (up_cell->upCellsCount() + 1);
}
void Cell::addEnvValue(float envv) {
    m_envValues.push_back(envv);
    m_envNum = m_envValues.size();
}
void Cell::addEnvValue(vector<float>& envvs) {
    for (vector<float>::iterator it = envvs.begin(); it != envvs.end(); it++) {
        m_envValues.push_back(*it);
    }
    m_envNum = m_envValues.size();
}
void Cell::getUpCellIndexes(vector<int>& up_indexes) {
    up_indexes.push_back(m_globalIndex);
    deque<Cell*> que;
    if (!m_upCells.empty())
    {
        for (vector<Cell*>::iterator it = m_upCells.begin(); it != m_upCells.end(); it++) {
            if ((*it)->upCellsCount() == 0) {
                up_indexes.push_back((*it)->globalIndex());
            }
            else {
                que.push_back(*it);
            }
        }
    }
    while (!que.empty()) {
        Cell *temp = que.front();
        que.pop_front();
        if (find(up_indexes.begin(), up_indexes.end(), temp->globalIndex()) == up_indexes.end()) {
            up_indexes.push_back(temp->globalIndex());
        }
        vector<Cell*> tmp_ups = temp->upCells();
        for (vector<Cell*>::iterator it = tmp_ups.begin(); it != tmp_ups.end(); it++) {
            //cout << (*it)->globalIndex() << endl;
            if ((*it)->upCellsCount() == 0) {
                if (find(up_indexes.begin(), up_indexes.end(), (*it)->globalIndex()) == up_indexes.end()) {
                    up_indexes.push_back((*it)->globalIndex());
                }
            }
            else {
                if (find(que.begin(), que.end(), *it) == que.end()) {
                    que.push_back(*it);
                }
            }
        }
    }
}

bool Cell::frequencyStats(map<int, Cell*>& cellmap, int envnum, float* minv, float* maxv, int interval_num) {
    vector<int> cur_upIdxes;
    getUpCellIndexes(cur_upIdxes);
    cout << "upflow num: " << cur_upIdxes.size() << endl;
    m_intervalNum = interval_num;
    if (m_envNum != envnum) {
        cout << "The input env num " << envnum << " should be " << m_envNum << endl;
        return false;
    }
    m_envFrequency = new float*[envnum];
    for (int kk = 0; kk < envnum; kk++){
        m_envFrequency[kk] = new float[interval_num];
        for (int j = 0; j < interval_num; j++) {
            m_envFrequency[kk][j] = 0.f;
        }
    }
    for (vector<int>::iterator it = cur_upIdxes.begin(); it != cur_upIdxes.end(); it++) {
        vector<float> tmpEnvvs = cellmap.at(*it)->getEnvValue();
        if (tmpEnvvs.size() != envnum) {
            continue; // although this may not happen, just check.
        }
        for (int kk = 0; kk < envnum; kk++) {
            int index = int(floor((tmpEnvvs[kk] - minv[kk]) * interval_num / (maxv[kk] - minv[kk])));
            m_envFrequency[kk][index] += 1;
        }
    }
    /// calculate percent
    for (int kk = 0; kk < envnum; kk++){
        for (int j = 0; j < interval_num; j++) {
            m_envFrequency[kk][j] /= cur_upIdxes.size();
        }
    }
    return true;
}

void Cell::printFrequency() {
    if (m_envNum < 0 || m_intervalNum < 0 || NULL == m_envFrequency) return;
    for (int k = 0; k < m_envNum; k++) {
        cout << "Cell index: " << m_globalIndex << ", env variable index: " << k << endl;
        for (int i = 0; i < m_intervalNum; i++) {
            cout << m_envFrequency[k][i] << ", ";
        }
        cout << endl;
    }
}

bool Cell::predictProperty(map<int, Cell*>& cells_map, int env_num, float* minv, float* maxv, int num_train, int* train_idx) {
    if (NULL == m_envFrequency || m_envNum != env_num || m_intervalNum < 0) {
        return false;
    }
    float *sampleSimilarity = new float[num_train];
    for (int i = 0; i < num_train; i++){
        Cell *cur_trainCell = cells_map.at(train_idx[i]);
        float * envSimilarity = new float[env_num];
        for (int kk = 0; kk < env_num; kk++){
            envSimilarity[kk] = histSimilarity(m_envFrequency[kk], cur_trainCell->getEnvFrequency()[kk], m_intervalNum);
        }
        sampleSimilarity[i] = getSimilaityIntegration(envSimilarity, env_num, NODATA_VALUE);
    }
    float maxSimilarity = getMaxValue(sampleSimilarity, num_train, NODATA_VALUE);
    float uncertainty = 1.f - maxSimilarity;
    float sumPredSimilarity = 0.f, sumSimilarity = 0.f;
    for (int i = 0; i < num_train; i++){
        Cell *cur_trainCell = cells_map.at(train_idx[i]);
        sumPredSimilarity += cur_trainCell->getAttribute() * sampleSimilarity[i];
        sumSimilarity += sampleSimilarity[i];
    }
    float pred = sumPredSimilarity / sumSimilarity;
    if (pred != pred) {
        pred = NODATA_VALUE;
    }
    m_uncertainty = uncertainty;
    m_predictAttr = pred;
}
