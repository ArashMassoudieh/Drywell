#include <iostream>
#include "grid.h"
#include "propertygenerator.h"

using namespace std;

int main()
{

    int nz=10;
    int nr=10;
    PropertyGenerator P(nz);
    P.correlation_length_scale = 1;
    P.dx = 0.2;
    P.assign_K_gauss(); // Assigns normal scores for K_sat
    P.Normalize_Ksat_normal_scores(0,1);
    P.write("K_sat_normal_score","/home/arash/Projects/Drywell_Result/Heterogeneous/Test/K_sat_score_1.txt");
    P.Normalize_Ksat_normal_scores(1,0.1);
    P.write("K_sat_normal_score","/home/arash/Projects/Drywell_Result/Heterogeneous/Test/K_sat_score_2.txt");

    CTimeSeries<double> K_sat_marginal_CDF("/home/arash/Projects/Drywell_Result/Heterogeneous/K_sat_Marginal_Distribution.csv"); // Loading Cummulative Distributions
    CTimeSeries<double> alpha_marginal_CDF("/home/arash/Projects/Drywell_Result/Heterogeneous/alpha_Marginal_Distribution.csv");
    CTimeSeries<double> n_marginal_CDF("/home/arash/Projects/Drywell_Result/Heterogeneous/n_Marginal_Distribution.csv");
    P.SetCorr(params::alpha, 0.378); //Correlation between alpha and K_sat normal scores
    P.SetCorr(params::n, 0.2); //Correlation between n and K_sat normal scores
    P.SetMarginalDistribution("K_sat",K_sat_marginal_CDF); //Assigning CDFs to properties
    P.SetMarginalDistribution("alpha",alpha_marginal_CDF);
    P.SetMarginalDistribution("n",n_marginal_CDF);
    P.PopulateRealValue("K_sat","K_sat_normal_score"); //Assign the actual K_sat
    P.Normalize("K_sat",P.mean("K_sat",true)); //Normalize K_sat by the geometrical mean of K_sat so the geometrical mean is zero
    P.Populate_Alpha_n_normal_scores(params::alpha); //Creates normal scores for alpha
    P.Populate_Alpha_n_normal_scores(params::n); //Creates normal scores for n
    P.PopulateRealValue("alpha","alpha_normal_score"); //Assign the actual values of alpha
    P.Normalize("alpha",0.05);
    P.PopulateRealValue("n","n_normal_score"); //.. n
    P.write("K_sat","/home/arash/Projects/Drywell_Result/Heterogeneous/Test/K_sat.txt");
    P.write("alpha","/home/arash/Projects/Drywell_Result/Heterogeneous/Test/alpha.txt");
    P.write("n","/home/arash/Projects/Drywell_Result/Heterogeneous/Test/n.txt");

    double dg=1;
    Grid G(nz,nr,1/dg,1.5);
    G.AssignProperty(&P); // Assign K_sat, alpha, n based on the Property Generator
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
    G.SetProp("inflow","/home/arash/Projects/Drywell_Result/inflow_zero.txt");
    G.SetProp("well_H","0");
    G.SetProp("well_H_old","0");
    G.SetProp("r_w","0.025");
    G.Solve(0,0.1,10,0.5);
    G.ExtractMoisture(1,0).writefile("/home/arash/Projects/Drywell_Result/Heterogeneous/Test/moist_well.csv");
    //G.ExtractMoisture(9,5).writefile("/home/arash/Projects/Drywell_Results/Test/moist_10_5.csv");
    //G.ExtractMoisture(5,9).writefile("/home/arash/Projects/Drywell_Results/Test/moist_5_10.csv");
    G.WriteResults("/home/arash/Projects/Drywell_Result/Heterogeneous/Test/theta.vtp");
    G.WaterDepth().make_uniform(0.1).writefile("/home/arash/Projects/Drywell_Result/Heterogeneous/Test/waterdepth.csv");
    cout<<"done!"<<endl;
}
