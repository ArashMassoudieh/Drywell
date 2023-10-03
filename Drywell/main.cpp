#include <iostream>
#include "grid.h"

using namespace std;

int main()
{
    int nz=10;
    int nr=5;
    Grid G(nz,nr,4,4);
    for (int i=0; i<nz; i++)
    {
        for (int j=0; j<nr; j++)
        {   G.cell(i,0)->Boundary.type = boundaryType::fixedmoisture;
            G.cell(i,0)->Boundary.value = 0.4;
            G.cell(i,nr-1)->Boundary.type = boundaryType::fixedmoisture;
            G.cell(i,nr-1)->Boundary.value = 0.4;
            G.cell(0,j)->Boundary.type = boundaryType::fixedmoisture;
            G.cell(0,j)->Boundary.value = 0.4;
            G.cell(nz-1,j)->Boundary.type = boundaryType::fixedmoisture;
            G.cell(nz-1,j)->Boundary.value = 0.4;
         }
    }

    G.OneStepSolve(0.1);
    G.write_to_vtp("/home/arash/Projects/Drywell_Results/theta.vtp");
    cout<<"done!"<<endl;
}
