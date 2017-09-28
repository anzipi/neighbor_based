#ifndef H_CELL_CLASS
#define H_CELL_CLASS

#include <vector>
#include <map>
using namespace std;

class Cell{
public:
    Cell(int idx, float flowdir);
    ~Cell();
    void addUpCell(Cell* up_cell);
    int upCellsCount() { return m_upCellsNum; }
    int globalIndex() { return m_globalIndex; }
    float flowDirection() { return m_flowdir; }
    vector<Cell*> upCells() { return m_upCells; }
    void getUpCellIndexes(vector<int>& up_indexes);
    void addEnvValue(float envv);
    void addEnvValue(vector<float>& envvs);
    vector<float>& getEnvValue() { return m_envValues; }
    // similarity calculation related
    void setAttribute(float attr) { m_predictAttr = attr; }
    float getAttribute() { return m_predictAttr; }
    void setUnCertainty(float unc) { m_uncertainty = unc; }
    float getUnCertainty() { return m_uncertainty; }
    bool frequencyStats(map<int, Cell*>& cellmap, int envnum, float* minv, float* maxv, int interval_num);
    float** getEnvFrequency() { return m_envFrequency; }
    void printFrequency();
    bool predictProperty(map<int, Cell*>& cells_map, int env_num, float* minv, float* maxv, int num_train, int* train_idx);
private:
    int            m_globalIndex;
    int            m_upCellsNum;
    int            m_envNum;
    int            m_intervalNum;
    float          m_flowdir;
    float          m_predictAttr;
    float          m_uncertainty;
    vector<float>  m_envValues;
    vector<Cell*>  m_upCells;
    float**        m_envFrequency;
};

#endif // !H_CELL_CLASS