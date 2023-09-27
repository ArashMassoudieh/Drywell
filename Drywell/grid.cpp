#include "grid.h"

Grid::Grid()
{

}

Grid::Grid(int _nz, int _nr, double depth, double radius)
{
    nz = _nz;
    nr = _nr;
    dz = depth/nr;
    dr = radius/nz;
    cells.resize(nz);
    for (unsigned int i=0; i<nz; i++)
        cells[i].resize(nr);
}

CVector_arma Grid::Residual(const CVector_arma &X, const double &dt)
{
    CVector_arma Res(nr*nz);
    SetStateVariable(X);
    for (unsigned int i=0; i<nz; i++)
        for (unsigned int j=0; j<nr; j++)
        {

        }
    return Res;
}

void Grid::SetStateVariable(const CVector_arma &X)
{
    for (unsigned int i=0; i<nz; i++)
        for (unsigned int j=0; j<nr; j++)
            cells[i][j].SetTheta(X[j+nz*i],_time::current);

}

double Grid::K(int i,int j,const edge &ej) const
{
    double Se;
    if (ej == edge::left)
        Se = (max(cells[i][j].Theta(_time::current),cells[i][j-1].Theta(_time::current))-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i][j-1].getValue("theta_r"))/(0.5*cells[i][j].getValue("theta_s")+0.5*cells[i][j-1].getValue("theta_s")-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i][j-1].getValue("theta_r"));
    else if (ej == edge::right)
        Se = (max(cells[i][j].Theta(_time::current),cells[i][j-1].Theta(_time::current))-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i][j+1].getValue("theta_r"))/(0.5*cells[i][j].getValue("theta_s")+0.5*cells[i][j+1].getValue("theta_s")-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i][j+1].getValue("theta_r"));
    else if (ej == edge::down)
        Se = (max(cells[i][j].Theta(_time::current),cells[i+1][j].Theta(_time::current))-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i+1][j].getValue("theta_r"))/(0.5*cells[i][j].getValue("theta_s")+0.5*cells[i+1][j].getValue("theta_s")-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i+1][j].getValue("theta_r"));
    else if (ej == edge::up)
        Se = (max(cells[i][j].Theta(_time::current),cells[i-1][j].Theta(_time::current))-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i-1][j].getValue("theta_r"))/(0.5*cells[i][j].getValue("theta_s")+0.5*cells[i-1][j].getValue("theta_s")-0.5*cells[i][j].getValue("theta_r")-0.5*cells[i-1][j].getValue("theta_r"));


}

double Grid::invC(int i,int j,const edge &ej) const
{

}

double Grid::getVal(int i, int j, const string &val, const edge &ej)
{
    if (ej == edge::left)
        return 0.5*(cells[i][j].getValue(val)+cells[i][j-1].getValue(val));
    else if (ej == edge::right)
        return 0.5*(cells[i][j].getValue(val)+cells[i][j+1].getValue(val));
    else if (ej == edge::down)
        return 0.5*(cells[i][j].getValue(val)+cells[i+1][j].getValue(val));
    else
        return 0.5*(cells[i][j].getValue(val)+cells[i-1][j].getValue(val));

}
