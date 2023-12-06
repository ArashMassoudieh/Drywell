#ifndef CELL_H
#define CELL_H

#include <map>
#include <string>

using namespace std;

enum class edge {left, right, up, down};
enum class _time {past, current, both};
enum class boundaryType {none, gradient, fixedpressure, fixedmoisture, symmetry};

struct boundary
{
    boundaryType type;
    edge boundary_edge;
    double value;
};

class Grid;

class Cell
{
public:
    Cell();
    ~Cell();
    Cell(const Cell &RHS);
    Cell& operator=(const Cell &RHS);
    void InitiateQuans();
    void SetTheta(const double &val, const _time &t);
    double Theta(const _time &t) const;
    double getValue(const string &quan) const;
    void SetValue(const string &quan, const double &value);
    boundary Boundary;
    void SetBoundary(boundaryType typ, edge boundaryEdge, const double &value = 0);
    double H(_time t) const;
    void SetQOut(const double &val)
    {
        qOut = val;
    }
    double QOut()
    {
        return qOut;
    }
private:
    map<string,double> quants;
    Grid *parent;
    double qOut = 0;

};

#endif // CELL_H
