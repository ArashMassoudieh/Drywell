#include "cell.h"
#include <cmath>

Cell::Cell()
{
    InitiateQuans();
}

void Cell::InitiateQuans()
{
    /*
    quants["theta"] = 0.25;
    quants["theta_past"] = 0.2;
    quants["theta_s"] = 0.4;
    quants["theta_r"] = 0.03;
    quants["alpha"] = 90;
    quants["n"] = 1.5;
    quants["Ks"] = 1;
    quants["epsilon"] = 0.01;
    */
    quants.theta = 0.25;
    quants.theta_past = 0.2;
    quants.theta_s = 0.4;
    quants.theta_r = 0.03;
    quants.alpha = 90;
    quants.n = 1.5;
    quants.Ks = 1;
    quants.epsilon = 0.01;
    Boundary.type = boundaryType::none;
}

Cell::~Cell()
{

}
Cell::Cell(const Cell &RHS)
{
    quants = RHS.quants;
    Boundary = RHS.Boundary;
    center_r = RHS.center_r;
    center_z = RHS.center_z;
}

Cell& Cell::operator=(const Cell &RHS)
{
    quants = RHS.quants;
    Boundary = RHS.Boundary;
    center_r = RHS.center_r;
    center_z = RHS.center_z;
    return *this;
}

double Cell::getValue(const string &quan) const
{
    /*if (quants.count(quan)>0)
        return quants.at(quan);
    else
        return 0;
    */
    
    if (quan == "theta")
        return quants.theta;
    else if (quan == "theta_past")
        return quants.theta_past;
    else if (quan == "theta_s")
        return quants.theta_s;
    else if (quan == "theta_r")
        return quants.theta_r;
    else if (quan == "alpha")
        return quants.alpha;
    else if (quan == "n")
        return quants.n;
    else if (quan == "Ks")
        return quants.Ks;
    else if (quan == "epsilon")
        return quants.epsilon;
    else
        return 0;
}

void Cell::SetValue(const string &quan, const double &value)
{
    /*if (quants.count(quan)>0)
        quants.at(quan) = value;
    */
    if (quan == "theta")
        quants.theta = value;
    else if (quan == "theta_past")
        quants.theta_past = value;
    else if (quan == "theta_s")
        quants.theta_s = value;
    else if (quan == "theta_r")
        quants.theta_r = value;
    else if (quan == "alpha")
        quants.alpha = value;
    else if (quan == "n")
        quants.n = value;
    else if (quan == "Ks")
        quants.Ks = value;
    else if (quan == "epsilon")
        quants.epsilon = value;
}

double Cell::getValue(prop quan) const
{
    switch (quan)
    {
        case prop::alpha:
            return quants.alpha;
            break;
        case prop::epsilon:
            return quants.epsilon;
            break;
        case prop::Ks:
            return quants.Ks;
            break;
        case prop::n:
            return quants.n;
            break;
        case prop::theta:
            return quants.theta;
            break;
        case prop::theta_past:
            return quants.theta_past;
            break;
        case prop::theta_r:
            return quants.theta_r;
            break;
        case prop::theta_s:
            return quants.theta_s;
            break;
    }
    
}

void Cell::SetValue(prop quan, const double& value)
{
    switch (quan)
    {
    case prop::alpha:
        quants.alpha = value;
        break;
    case prop::epsilon:
        quants.epsilon = value;
        break;
    case prop::Ks:
        quants.Ks = value;
        break;
    case prop::n:
        quants.n = value;
        break;
    case prop::theta:
        quants.theta = value;
        break;
    case prop::theta_past:
        quants.theta_past = value;
        break;
    case prop::theta_r:
        quants.theta_r = value;
        break;
    case prop::theta_s:
        quants.theta_s = value;
        break;
    }

}

void Cell::SetTheta(const double &val,const _time &t)
{
    if (t==_time::current)
        quants.theta = val;
    else if (t==_time::past)
        quants.theta_past=val;
    else
    {
        quants.theta = val;
        quants.theta_past=val;
    }
}

double Cell::Theta(const _time &t) const
{
    if (t==_time::current)
          return quants.theta;
    else if (t==_time::past)
        return quants.theta_past;
    else
    {
        return quants.theta;
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
    double Se = (Theta(t)-quants.theta_r)/(quants.theta_s-quants.theta_r);
    double H;
    H = -1.0/quants.alpha*pow(pow(max(min(Se,0.999),0.001),-quants.n/(quants.n-1))-1,1/quants.n);
    if (Se>0.999)
       H += (Se-0.999)/quants.epsilon;

    return H;
}


