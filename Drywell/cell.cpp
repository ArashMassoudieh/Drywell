#include "cell.h"

Cell::Cell()
{
    InitiateQuans();
}

void Cell::InitiateQuans()
{
    quants["theta"] = 1;
    quants["theta_past"] = 1;
    quants["theta_s"] = 0.4;
    quants["theta_r"] = 0.08;
    quants["alpha"] = 1;
    quants["n"] = 1.58;
    quants["Ks"] = 0.01;
    quants["epsilon"] = 0.01;
    Boundary.type = boundaryType::none;
}

double Cell::getValue(const string quan) const
{
    if (quants.count(quan)>0)
        return quants.at(quan);
    else
        return 0;
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
