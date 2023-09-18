#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cell.h>

using namespace std;

class Grid
{
public:
    Grid();
    Grid(int nz, int nr, double depth, double radius);

private:
    vector<vector<Cell>> cells;
};

#endif // GRID_H
