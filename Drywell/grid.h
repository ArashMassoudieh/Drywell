#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cell.h>
#include "Vector_arma.h"
#include "Matrix_arma.h"
#include "Matrix.h"
#include "Vector.h"

//#define CMatrix_arma CMatrix
//#define CVector_arma CVector

using namespace std;

struct _solution_state
{
    double t;
    double dt;
    unsigned int NI_max=20;
    unsigned int NI_min=5;
    double dt_scale_factor = 0.75;
    int number_of_iterations = 0;
    double dt_scale_factor_fail = 0.2;
    int max_iterations = 100;
};

class Grid
{
public:
    Grid();
    Grid(int nz, int nr, double depth, double radius);
    CVector_arma Residual(const CVector_arma &X, const double &dt);
    void SetStateVariable(const CVector_arma &X,const _time &t=_time::current);
    CMatrix_arma Jacobian(const CVector_arma &X, const double &dt);
    double getVal(int i, int j, const string &val, const edge &ej) const;
    bool OneStepSolve(const double &dt);
    bool Solve(const double &t0, const double &dt0, const double &t_end, const double &write_interval);
    Cell* cell(int i, int j)
    {
        if (i>=0 && j>=0 && i<nz & j<nr)
            return &cells[i][j];
        else
            return nullptr;
    }
    void write_to_vtp(const string &name) const;
    void write_to_vtp(const string &name,const CMatrix &res) const;
    void WriteResults(const string &filename);
    _solution_state Solution_State;
    CMatrix H();
    CMatrix Se();
    CMatrix Theta(_time t);
    vector<CMatrix> results;

private:
    vector<vector<Cell>> cells;
    unsigned int nz;
    unsigned int nr;
    double dz;
    double dr;
    double K(int i,int j,const edge &ej);
    double D(int i,int j,const edge &ej);
    double H(int i,int j,const edge &ej, _time t);
    double invC(int i,int j,const edge &ej);
    Cell* Neighbour(int i, int j, const edge &ej, bool op=false);
    CVector_arma GetStateVariable(const _time &t=_time::current) const;

};

#endif // GRID_H
