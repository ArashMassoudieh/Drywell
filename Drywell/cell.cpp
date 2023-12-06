#include "cell.h"
#include <cmath>

Cell::Cell()
{
    InitiateQuans();
}

void Cell::InitiateQuans()
{
    quants["theta"] = 0.2;
    quants["theta_past"] = 0.2;
    quants["theta_s"] = 0.4;
    quants["theta_r"] = 0.08;
    quants["alpha"] = 100;
    quants["n"] = 1.58;
    quants["Ks"] = 0.3;
    quants["epsilon"] = 0.01;
    Boundary.type = boundaryType::none;
}

Cell::~Cell()
{

}
Cell::Cell(const Cell &RHS)
{
    quants = RHS.quants;
    Boundary = RHS.Boundary;
}

Cell& Cell::operator=(const Cell &RHS)
{
    quants = RHS.quants;
    Boundary = RHS.Boundary;
    return *this;
}

double Cell::getValue(const string &quan) const
{
    if (quants.count(quan)>0)
        return quants.at(quan);
    else
        return 0;
}

void Cell::SetValue(const string &quan, const double &value)
{
    if (quants.count(quan)>0)
        quants.at(quan) = value;
}

void Cell::SetTheta(const double &val,const _time &t)
{
    if (t==_time::current)
        quants["theta"] = val;
    else if (t==_time::past)
        quants["theta_past"]=val;
    else
    {
        quants["theta"] = val;
        quants["theta_past"]=val;
    }
}

double Cell::Theta(const _time &t) const
{
    if (t==_time::current)
          return quants.at("theta");
    else if (t==_time::past)
        return quants.at("theta_past");
    else
    {
        return quants.at("theta");
    }
}

void Cell::SetBoundary(boundaryType typ, edge boundaryEdge, const double &value)
{
    Boundary.boundary_edge = boundaryEdge;
    Boundary.type = typ;
    Boundary.value = value;
}

double Cell::H(_time t) const
{
    double Se = (Theta(t)-quants.at("theta_r"))/(quants.at("theta_s")-quants.at("theta_r"));
    double H;
    H = -1.0/quants.at("alpha")*pow(pow(min(Se,0.999),-(quants.at("n")-1)/quants.at("n")-1),quants.at("n"));
    if (Se>0.999)
       H += (Se-0.999)/quants.at("epsilon");

    return H;
}


