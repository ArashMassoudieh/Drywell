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
#include "propertygenerator.h"

//#define CMatrix_arma CMatrix
//#define CVector_arma CVector

//using namespace std;

struct _solution_state
{
    double t;
    double dt;
    unsigned int NI_max=30;
    unsigned int NI_min=10;
    double dt_scale_factor = 0.75;
    int number_of_iterations = 0;
    double dt_scale_factor_fail = 0.2;
    int max_iterations = 100;
    double err;
    double lambda=1;
    double lambda_reduction_factor=0.1;
    double min_time_step = 1e-12;
    unsigned int count_error_expanding_allowed = 10;
    enum class _solution_method {NR, LM, NLC} Solution_Method = _solution_method::NLC;
    double max_time_step = 0.1;
};

class Grid
{
public:
    Grid();
    ~Grid();
    Grid(const Grid &RHS);
    Grid& operator=(const Grid &RHS);
    Grid(int nz, int nr, double depth, double radius);
    CVector_arma Residual(const CVector_arma &X, const double &dt, bool resetstatevariables = false);
    CVector_arma Residual_TR(const CVector_arma &X, const double &dt, bool resetstatevariables = false);
    double CalcOutFlow();
    void SetStateVariable(const CVector_arma &X,const _time &t=_time::current);
    void SetStateVariable_TR(const CVector_arma &X,const _time &t=_time::current);
    CMatrix_arma Jacobian(const CVector_arma &X, const double &dt);
    double getVal(int i, int j, const std::string &val, const edge &ej);
    double getVal(int i, int j, prop val, const edge& ej);
    bool OneStepSolve(const double &dt);
    bool OneStepSolveLM(const double &dt);
    bool OneStepSolve_no_lamba_correction(const double &dt);
    bool Solve(const double &t0, const double &dt0, const double &t_end, const double &write_interval);
    Cell* cell(int i, int j)
    {
        return &cells[i][j];
    }
    void write_to_vtp(const std::string &name) const;
    void write_to_vtp(const std::string &name,const CMatrix &res,const std::string &quanname="Moisture Content", const double &scale = 1) const;
    void WriteResults(const std::string &filename);
    double Max(const std::string &quan);
    double Min(const std::string &quan);
    void WriteResults(const std::string &quan, const std::string &filename);
    CTimeSeries<double> ExtractMoisture(int i, int j);
    _solution_state Solution_State;
    CMatrix H();
    void UpdateH();
    CMatrix Se();
    double TotalWaterContent();
    double WellWaterContent();
    CMatrix QuanMatrix(const std::string &quan);
    CMatrix Theta(_time t);
    std::vector<CMatrix> results;
    bool SetProp(const std::string &propname, const std::string &value);
    CTimeSeries<double> &WaterDepth()
    {
        return Well_Water_Depth;
    }
    CTimeSeries<double> &TotalWaterInSoil()
    {
        return Total_Water_Content;
    }
    CTimeSeries<double> &TotalWaterInWell()
    {
        return Well_Water_Content;
    }
    CTimeSeries<double> &OutFlow()
    {
        return Outflow;
    }
    bool AssignProperty(PropertyGenerator *prop);
    void SetName(std::string _name) {name = _name;};
private:
    std::vector<std::vector<Cell>> cells;
    std::vector<std::vector<Interface>> interfaces_r;
    std::vector<std::vector<Interface>> interfaces_z;
    unsigned int nz;
    unsigned int nr;
    double dz;
    double dr;
    double K(int i,int j,const edge &ej);
    double D(int i,int j,const edge &ej);
    double invC(int i,int j,const edge &ej);
    Cell* Neighbour(int i, int j, const edge &ej, bool op=false);
    CVector_arma GetStateVariable(const _time &t=_time::current) const;
    CVector_arma GetStateVariable_TR(const _time &t=_time::current) const;
    double well_H = 0;
    double well_H_old = 0;
    double r_w;
    double beta;
    double alpha;
    CTimeSeries<double> inflow;
    CTimeSeries<double> Well_Water_Depth;
    CTimeSeries<double> Outflow;
    CTimeSeries<double> Total_Water_Content;
    CTimeSeries<double> Well_Water_Content;
    std::string name;

};

double upstream(const double &x1, const double &x2, const double &val1, const double &val2)
{
    if (x1>x2)
        return val1;
    else
        return val2;
}

#endif // GRID_H
