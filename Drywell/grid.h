#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cell.h>
#include "Vector_arma.h"
#include "Matrix_arma.h"
#include "Matrix.h"
#include "Vector.h"
#include "BTC.h"
#include "interface.h"

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
    ~Grid();
    Grid(const Grid &RHS);
    Grid& operator=(const Grid &RHS);
    Grid(int nz, int nr, double depth, double radius);
    CVector_arma Residual(const CVector_arma &X, const double &dt);
    void SetStateVariable(const CVector_arma &X,const _time &t=_time::current);
    CMatrix_arma Jacobian(const CVector_arma &X, const double &dt);
    double getVal(int i, int j, const string &val, const edge &ej);
    bool OneStepSolve(const double &dt);
    bool Solve(const double &t0, const double &dt0, const double &t_end, const double &write_interval);
    Cell* cell(int i, int j)
    {
        return &cells[i][j];
    }
    void write_to_vtp(const string &name) const;
    void write_to_vtp(const string &name,const CMatrix &res) const;
    void WriteResults(const string &filename);
    CTimeSeries<double> ExtractMoisture(int i, int j);
    _solution_state Solution_State;
    CMatrix H();
    CMatrix Se();
    CMatrix Theta(_time t);
    vector<CMatrix> results;
    bool SetProp(const string &propname, const string &value);
    CTimeSeries<double> &WaterDepth()
    {
        return Well_Water_Depth;
    }
private:
    vector<vector<Cell>> cells;
    vector<vector<Interface>> interfaces_r;
    vector<vector<Interface>> interfaces_z;
    unsigned int nz;
    unsigned int nr;
    double dz;
    double dr;
    double K(int i,int j,const edge &ej);
    double D(int i,int j,const edge &ej);
    double invC(int i,int j,const edge &ej);
    Cell* Neighbour(int i, int j, const edge &ej, bool op=false);
    CVector_arma GetStateVariable(const _time &t=_time::current) const;
    double well_H = 0;
    double well_H_old = 0;
    double r_w;
    double beta;
    double alpha;
    CTimeSeries<double> inflow;
    CTimeSeries<double> Well_Water_Depth;

};

#endif // GRID_H
