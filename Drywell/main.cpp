#include <iostream>
#include "grid.h"

using namespace std;

int main()
{
    Grid G(4,4,1,1);
    for (int i=0; i<4; i++)
    {
        G.cell(i,0)->Boundary.type = boundaryType::fixedmoisture;
        G.cell(i,0)->Boundary.value = 0.4;
        G.cell(i,3)->Boundary.type = boundaryType::fixedmoisture;
        G.cell(i,3)->Boundary.value = 0.4;
        G.cell(0,i)->Boundary.type = boundaryType::fixedmoisture;
        G.cell(0,i)->Boundary.value = 0.4;
        G.cell(3,i)->Boundary.type = boundaryType::fixedmoisture;
        G.cell(3,i)->Boundary.value = 0.4;


    }

    G.OneStepSolve(0.001);

    cout<<"done!"<<endl;
}
