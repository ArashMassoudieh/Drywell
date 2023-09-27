#ifndef CELL_H
#define CELL_H

#include <map>
#include <string>

using namespace std;

enum class _time {past, current, both};

class Grid;

class Cell
{
public:
    Cell();
    void InitiateQuans();
    void SetTheta(const double &val, const _time &t);
    double Theta(const _time &t) const;
    double getValue(const string quan) const;
private:
    map<string,double> quants;
    Grid *parent;
};

#endif // CELL_H
