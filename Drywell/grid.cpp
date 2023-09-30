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
            if (cells[i][j].Boundary.type==boundaryType::none)
            {   double r = (j+0.5)*dr;
                Res[j+nz*i] = (cells[i][j].Theta(_time::current)-cells[i][j].Theta(_time::past))/dt;
                Res[j+nz*i] += -1/r*pow(1.0/dr,2)*((r+dr/2)*D(i,j,edge::right)*(cells[i][j+1].Theta(_time::current)-cells[i][j].Theta(_time::current))-(r-dr/2)*D(i,j,edge::left)*(cells[i][j].Theta(_time::current)-cells[i][j-1].Theta(_time::current)));
                Res[j+nz*i] += -pow(1.0/dz,2)*(D(i,j,edge::down)*(cells[i+1][j].Theta(_time::current)-cells[i][j].Theta(_time::current))-D(i,j,edge::up)*(cells[i][j].Theta(_time::current)-cells[i-1][j].Theta(_time::current)));
                Res[j+nz*i] += -1/dz*(K(i,j,edge::up)-K(i,j,edge::down));
            }
            else if (cells[i][j].Boundary.type==boundaryType::fixedmoisture)
            {
                Res[j+nz*i] = cells[i][j].Theta(_time::current) - cells[i][j].Boundary.value;
            }
            else if (cells[i][j].Boundary.type==boundaryType::fixedpressure)
            {
                Res[j+nz*i] = cells[i][j].H(_time::current) - cells[i][j].Boundary.value;
            }
            else if (cells[i][j].Boundary.type==boundaryType::gradient)
            {
                if (cells[i][j].Boundary.boundary_edge == edge::down || cells[i][j].Boundary.boundary_edge == edge::up)
                    Res[i+nz*j] = (cells[i][j].Theta(_time::current)-Neighbour(i,j,cells[i][j].Boundary.boundary_edge)->Theta(_time::current))/dz;
                else
                    Res[i+nz*j] = (cells[i][j].Theta(_time::current)-Neighbour(i,j,cells[i][j].Boundary.boundary_edge)->Theta(_time::current))/dr;

            }

        }
    return Res;
}

void Grid::SetStateVariable(const CVector_arma &X)
{
    for (unsigned int i=0; i<nz; i++)
        for (unsigned int j=0; j<nr; j++)
            cells[i][j].SetTheta(X[j+nz*i],_time::current);

}

double Grid::D(int i,int j,const edge &ej) const
{
    return K(i,j,ej)*invC(i,j,ej);
}

double Grid::K(int i,int j,const edge &ej) const
{
    double Se = max((max(cells[i][j].Theta(_time::current),cells[i][j-1].Theta(_time::current))-getVal(i,j,"theta_r",ej))/(getVal(i,j,"theta_s",ej)-getVal(i,j,"theta_r",ej)),1e-6);
    double m = 1.0-1.0/getVal(i,j,"n",ej);
    double K = pow(Se,0.5)*getVal(i,j,"Ks",ej)*pow(1-pow(1-pow(Se,1.0/m),m),2);
    return K;

}

double Grid::invC(int i,int j,const edge &ej) const
{
    double C;
    double Se = max((max(cells[i][j].Theta(_time::current),cells[i][j-1].Theta(_time::current))-getVal(i,j,"theta_r",ej))/(getVal(i,j,"theta_s",ej)-getVal(i,j,"theta_r",ej)),1e-6);
    double m = 1.0-1.0/getVal(i,j,"n",ej);
    if (getVal(i,j,"theta",ej)>getVal(i,j,"theta_s",ej))
        C = 1.0/getVal(i,j,"epsilon",ej);
    else
        C = 1.0/getVal(i,j,"alpha",ej)/(getVal(i,j,"theta_s",ej)-getVal(i,j,"theta_r",ej))*pow(pow(Se,-m)-1,-m)*pow(Se,-m-1);
    return C;
}

double Grid::getVal(int i, int j, const string &val, const edge &ej) const
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

Cell* Grid::Neighbour(int i, int j, const edge &ej)
{
    if (ej == edge::down)
        return &cells[i-1][j];
    else if (ej == edge::up)
        return &cells[i+1][j];
    else if (ej == edge::left)
        return &cells[i][j+1];
    else if (ej == edge::right)
        return &cells[i][j-1];

}
