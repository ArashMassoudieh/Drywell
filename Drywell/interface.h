#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>
#include <map>

//using namespace std;
class Grid;

class Interface
{
public:
    Interface();
    ~Interface();
    Interface(const Interface &RHS);
    Interface& operator=(const Interface &RHS);
    double getValue(const std::string quan) const;
private:
    std::map<std::string,double> quants;
    Grid *parent;
};

#endif // INTERFACE_H
