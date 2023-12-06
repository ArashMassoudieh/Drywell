#include <iostream>
#include "grid.h"

using namespace std;

int main()
{
    int nz=20;
    int nr=20;
    double dg=0.7;
    Grid G(nz,nr,4,4);
    for (int i=0; i<nz; i++)
    {
        if (i*4.0/double(nz)<=2)
        {   G.cell(i,0)->Boundary.type = boundaryType::fixedmoisture;
            G.cell(i,0)->Boundary.value = 0.4;
        }
        else
        {   G.cell(i,0)->Boundary.type = boundaryType::symmetry;
            G.cell(i,0)->Boundary.boundary_edge = edge::left;
        }
        G.cell(i,nr-1)->Boundary.type = boundaryType::symmetry;
        G.cell(i,nr-1)->Boundary.boundary_edge = edge::right;
    }

    for (int i=0; i<nz*dg; i++)
        G.cell(i,0)->Boundary.type = boundaryType::fixedpressure;


    for (int j=1; j<nr-1; j++)
    {
        G.cell(0,j)->Boundary.type = boundaryType::symmetry;
        G.cell(0,j)->Boundary.boundary_edge = edge::up;
        G.cell(nz-1,j)->Boundary.type = boundaryType::fixedmoisture;
        G.cell(nz-1,j)->Boundary.value = 0.4;
    }


    G.SetProp("pond_alpha","1000.0");
    G.SetProp("pond_beta","2");
    G.SetProp("alpha","20");
    G.SetProp("inflow","/home/arash/Project_Khiem/Drywell_Results/inflow_zero.txt");
    G.SetProp("well_H","0");
    G.SetProp("well_H_old","0");
    G.SetProp("r_w","0.75");
    G.Solve(0,0.1,10,0.2);
    G.WriteResults("/home/arash/Projects/Drywell_Results/R=20,r_w=0_75/theta_7.vtp");
    G.WaterDepth().make_uniform(0.1).writefile("/home/arash/Projects/Drywell_Results/R=20,r_w=0_75/waterdepth_7.csv");
    cout<<"done!"<<endl;
}
