#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cell.h>
#include "Vector_arma.h"

using namespace std;

class Grid
{
public:
    Grid();
    Grid(int nz, int nr, double depth, double radius);
    CVector_arma Residual(const CVector_arma &X, const double &dt);
    void SetStateVariable(const CVector_arma &X);
    double getVal(int i, int j, const string &val, const edge &ej) const;
private:
    vector<vector<Cell>> cells;
    unsigned int nz;
    unsigned int nr;
    double dz;
    double dr;
    double K(int i,int j,const edge &ej) const;
    double D(int i,int j,const edge &ej) const;
    double invC(int i,int j,const edge &ej) const;
};

#endif // GRID_H
