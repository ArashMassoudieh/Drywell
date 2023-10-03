#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cell.h>
#include "Vector_arma.h"
#include "Matrix_arma.h"

using namespace std;

class Grid
{
public:
    Grid();
    Grid(int nz, int nr, double depth, double radius);
    CVector_arma Residual(const CVector_arma &X, const double &dt);
    void SetStateVariable(const CVector_arma &X);
    CMatrix_arma Jacobian(const CVector_arma &X, const double &dt);
    double getVal(int i, int j, const string &val, const edge &ej) const;
    bool OneStepSolve(const double &dt);
    Cell* cell(int i, int j)
    {
        return &cells[i][j];
    }
    void write_to_vtp(const string &name) const;
private:
    vector<vector<Cell>> cells;
    unsigned int nz;
    unsigned int nr;
    double dz;
    double dr;
    double K(int i,int j,const edge &ej) const;
    double D(int i,int j,const edge &ej) const;
    double invC(int i,int j,const edge &ej) const;
    Cell* Neighbour(int i, int j, const edge &ej);
    CVector_arma GetStateVariable() const;

};

#endif // GRID_H
