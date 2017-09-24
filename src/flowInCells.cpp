#include <queue>
#include <deque>
#include <algorithm>
#include "flowInCells.h"
#include <iostream>
using namespace std;
Cell::Cell(int idx) : m_globalIndex(idx), m_upCellsNum(0) {
    
}
Cell::~Cell() {

}
void Cell::addUpCell(Cell* up_cell) {
    m_upCells.push_back(up_cell);
    m_upCellsNum += (up_cell->upCellsCount() + 1);
}
void Cell::addEnvValue(float envv) {
    m_envValues.push_back(envv);
}
void Cell::addEnvValue(vector<float>& envvs) {
    for (vector<float>::iterator it = envvs.begin(); it != envvs.end(); it++) {
        m_envValues.push_back(*it);
    }
}
void Cell::getUpCellIndexes(vector<int>& up_indexes) {
    up_indexes.push_back(m_globalIndex);
    deque<Cell*> que;
    if (m_upCellsNum > 0)
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