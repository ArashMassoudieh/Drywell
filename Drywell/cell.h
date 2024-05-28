#ifndef CELL_H
#define CELL_H

#include <map>
#include <string>

//using namespace std;

enum class edge {left, right, up, down};
enum class _time {past, current, both};
enum class boundaryType {none, gradient, fixedpressure, fixedmoisture, symmetry};
enum class prop {theta, theta_past, theta_s, theta_r, alpha, n, Ks, epsilon, H};
struct boundary
{
    boundaryType type;
    edge boundary_edge;
    double value;
};

struct _quantypes
{
    double theta = 0.25;
    double theta_past = 0.2;
    double theta_s = 0.4;
    double theta_r = 0.03;
    double alpha = 90;
    double n = 1.5;
    double Ks = 1;
    double epsilon = 0.01;
    double H;

};

class Grid;

class Cell
{
public:
    Cell();
    ~Cell();
    Cell(const Cell &RHS);
    Cell& operator=(const Cell &RHS);
    void InitiateQuans();
    void SetTheta(const double &val, const _time &t);
    double Theta(const _time &t) const;
    double getValue(const std::string &quan) const;
    double getValue(prop property) const;
    void SetValue(const std::string &quan, const double &value);
    void SetValue(prop property, const double& value);
    boundary Boundary;
    void SetBoundary(boundaryType typ, edge boundaryEdge, const double &value = 0);
    double H(_time t, bool calc = true) const;
    void SetQOut(const double &val)
    {
        qOut = val;
    }
    double QOut()
    {
        return qOut;
    }
    double z() const {return center_z;}
    double r() const {return center_r;}
    void setz(const double &val) {center_z = val;}
    void setr(const double &val) {center_r=val;}
private:
    //std::map<string,double> quants;
    _quantypes quants;
    Grid *parent;
    double qOut = 0;
    double center_z;
    double center_r;

};

#endif // CELL_H
