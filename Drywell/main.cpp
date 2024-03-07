#include <iostream>
#include "grid.h"
#include "propertygenerator.h"

using namespace std;

int main()
{

    PropertyGenerator P(20);
    P.correlation_length_scale = 1;
    P.dx = 0.2;
    P.assign_K_gauss();
    P.Normalize_Ksat_normal_scores(0,1);
    P.write("K_sat_normal_score","/home/arash/Projects/Drywell_Results/Test/K_sat_score_1.txt");
    P.Normalize_Ksat_normal_scores(1,0.1);
    P.write("K_sat_normal_score","/home/arash/Projects/Drywell_Results/Test/K_sat_score_2.txt");

    CTimeSeries<double> K_sat_marginal_CDF("/home/arash/Projects/Drywell_Results/K_sat_Marginal_Distribution.csv");
    CTimeSeries<double> alpha_marginal_CDF("/home/arash/Projects/Drywell_Results/K_sat_Marginal_Distribution.csv");
    CTimeSeries<double> n_marginal_CDF("/home/arash/Projects/Drywell_Results/K_sat_Marginal_Distribution.csv");
    P.SetCorr(params::alpha, 0.2);
    P.SetCorr(params::n, 0.2);
    P.SetMarginalDistribution("K_sat",K_sat_marginal_CDF);
    P.SetMarginalDistribution("alpha",alpha_marginal_CDF);
    P.SetMarginalDistribution("n",n_marginal_CDF);
    P.PopulateRealValue("K_sat","K_sat_normal_score");
    P.Populate_Alpha_n_normal_scores(params::alpha);
    P.Populate_Alpha_n_normal_scores(params::n);
    P.PopulateRealValue("alpha","alpha_normal_score");
    P.PopulateRealValue("n","n_normal_score");


    int nz=20;
    int nr=20;


    double dg=0.5;
    Grid G(nz,nr,1/dg,1.5);
    G.AssignProperty(&P);
    for (int i=0; i<nz; i++)
    {
        G.cell(i,0)->Boundary.type = boundaryType::symmetry;
        G.cell(i,0)->Boundary.boundary_edge = edge::left;
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
    G.SetProp("inflow","/home/arash/Projects/Drywell_Results/Test/inflow_zero.txt");
    G.SetProp("well_H","0");
    G.SetProp("well_H_old","0");
    G.SetProp("r_w","0.025");
    G.Solve(0,0.1,10,0.05);
    //G.ExtractMoisture(1,0).writefile("/home/arash/Projects/Drywell_Results/Test/moist_well.csv");
    //G.ExtractMoisture(10,5).writefile("/home/arash/Projects/Drywell_Results/Test/moist_10_5.csv");
    //G.ExtractMoisture(5,10).writefile("/home/arash/Projects/Drywell_Results/Test/moist_5_10.csv");
    G.WriteResults("/home/arash/Projects/Drywell_Results/Test/theta5_1.vtp");
    G.WaterDepth().make_uniform(0.1).writefile("/home/arash/Projects/Drywell_Results/Test/waterdepth5_1.csv");
    G.OutFlow().make_uniform(0.1).writefile("/home/arash/Projects/Drywell_Results/Test/outflow.csv");
    cout<<"done!"<<endl;
}
