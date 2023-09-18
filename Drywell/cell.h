#ifndef CELL_H
#define CELL_H

#include <map>
#include <string>

using namespace std;

class Grid;

class Cell
{
public:
    Cell();

private:
    map<string,double> quants;
    Grid *parent;
};

#endif // CELL_H
