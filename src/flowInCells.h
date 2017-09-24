#include <vector>
using namespace std;

class Cell{
public:
    Cell(int idx);
    ~Cell();
    void addUpCell(Cell* up_cell);
    int upCellsCount() { return m_upCellsNum; }
    int globalIndex() { return m_globalIndex; }
    vector<Cell*> upCells() { return m_upCells; }
    void getUpCellIndexes(vector<int>& up_indexes);
    void addEnvValue(float envv);
    void addEnvValue(vector<float>& envvs);
    vector<float>& getEnvValue() { return m_envValues; }
private:
    int            m_globalIndex;
    int            m_upCellsNum;
    vector<float>  m_envValues;
    vector<Cell*>  m_upCells;
};