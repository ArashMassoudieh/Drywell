#include <iostream>
#include "grid.h"


using namespace std;

int main()
{
    cout<<aquiutils::numbertostring(5, 3)<<std::endl;

    int nz=5;
    int nr=5;
    double dg=0.5;
    Grid G(nz,nr,1/dg,1.5);
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
    G.SetProp("alpha","30");
    G.SetProp("n","1.5");
    G.SetProp("inflow","/home/khiem/Projects/Drywell_Result/inflow_zero.txt");
    G.SetProp("well_H","0");
    G.SetProp("well_H_old","0");
    G.SetProp("r_w","0.025");
    G.Solve(0,0.1,10,0.05);
    G.ExtractMoisture(2,2).writefile("/home/khiem/Projects/Drywell_Result/12_07_23/Trial_5_1/moist_5_5.csv");
    G.ExtractMoisture(3,2).writefile("/home/khiem/Projects/Drywell_Result/12_07_23/Trial_5_1/moist_10_5.csv");
    G.ExtractMoisture(2,3).writefile("/home/khiem/Projects/Drywell_Result/12_07_23/Trial_5_1/moist_10_5.csv");
    G.WriteResults("/home/khiem/Projects/Drywell_Result/12_07_23/Trial_5_1/theta5_1.vtp");
    G.WaterDepth().make_uniform(0.1).writefile("/home/khiem/Projects/Drywell_Result/12_07_23/Trial_5_1/waterdepth5_1.csv");
    cout<<"done!"<<endl;
}
