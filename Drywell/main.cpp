#include <iostream>
#include "grid.h"
#include "propertygenerator.h"
#include <QString>
#include <QMap>
#include <QDir>
#include <QDebug>

#ifdef VALGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;

struct ModelParameters{
    QString alpha;
    QString n;
    QString folder;
};


int main()
{
    QMap<QString, ModelParameters> parameter_set;
    ModelParameters case1;
    case1.alpha = "30";
    case1.n = "2";
    case1.folder = "Alpha=30, n=2";
    parameter_set[case1.folder] = case1;

    ModelParameters case2;
    case2.alpha = "10";
    case2.n = "2";
    case2.folder = "Alpha=10, n=2";
    parameter_set[case2.folder] = case2;

    ModelParameters case3;
    case3.alpha = "15";
    case3.n = "2";
    case3.folder = "Alpha=15, n=2";
    parameter_set[case3.folder] = case3;

    enum class _mode {homogeneous, heterogeneous};
    _mode mode = _mode::homogeneous;
    int nz=20;
    int nr=20;

    for (QMap<QString, ModelParameters>::iterator it = parameter_set.begin(); it!= parameter_set.end(); it++)
    {
        string Results_Folder = "/home/arash/Projects/Drywell_Result/Homogeneous/" + it->folder.toStdString();
        //string Results_Folder = "F:/Projects/Drywell_Result";
        QDir dir("/home/arash/Projects/Drywell_Result/Homogeneous/");

        if (!dir.exists(it->folder))
        {
            bool res = dir.mkdir(it->folder);
            qDebug()<<res;
        }
        PropertyGenerator P(nz);
        P.correlation_length_scale = 1;
        P.dx = 0.2;
        P.assign_K_gauss(); // Assigns normal scores for K_sat
        P.Normalize_Ksat_normal_scores(0,1);
        P.write("K_sat_normal_score", Results_Folder + "/K_sat_score_1.txt");
    cout<<"1"<<endl;
        CTimeSeries<double> K_sat_marginal_CDF(Results_Folder + "/K_sat_Marginal_Distribution.csv"); // Loading Cummulative Distributions
        CTimeSeries<double> alpha_marginal_CDF(Results_Folder + "/alpha_Marginal_Distribution.csv");
        CTimeSeries<double> n_marginal_CDF(Results_Folder + "/n_Marginal_Distribution.csv");
        P.SetCorr(params::alpha, 0.378); //Correlation between alpha and K_sat normal scores
        P.SetCorr(params::n, 0.2); //Correlation between n and K_sat normal scores
    cout<<"2"<<endl;
        P.SetMarginalDistribution("K_sat",K_sat_marginal_CDF); //Assigning CDFs to properties
    cout<<"2.1"<<endl;
        P.SetMarginalDistribution("alpha",alpha_marginal_CDF);
    cout<<"2.2"<<endl;
        P.SetMarginalDistribution("n",n_marginal_CDF);
    cout<<"2.3"<<endl;
        P.PopulateRealValue("K_sat","K_sat_normal_score"); //Assign the actual K_sat
    cout<<"2.4"<<endl;
        P.Normalize("K_sat",P.mean("K_sat",true)); //Normalize K_sat by the geometrical mean of K_sat so the geometrical mean is zero
    cout<<"3"<<endl;
        P.Populate_Alpha_n_normal_scores(params::alpha); //Creates normal scores for alpha
        P.Populate_Alpha_n_normal_scores(params::n); //Creates normal scores for n
        P.PopulateRealValue("alpha","alpha_normal_score"); //Assign the actual values of alpha
        P.Normalize("alpha",0.05);
        P.PopulateRealValue("n","n_normal_score"); //.. n
        P.write("K_sat", Results_Folder + "/K_sat.txt");
        P.write("alpha", Results_Folder + "/alpha.txt");
        P.write("n", Results_Folder + "/n.txt");
    cout<<"4"<<endl;
        double dg=0.5;
        Grid G(nz,nr,1/dg,1.5);
        if (mode == _mode::heterogeneous)
            G.AssignProperty(&P); // Assign K_sat, alpha, n based on the Property Generator
    cout<<"5"<<endl;
        for (int i=0; i<nz; i++)
        {
            G.cell(i,0)->Boundary.type = boundaryType::symmetry;
            G.cell(i,0)->Boundary.boundary_edge = edge::left;
            G.cell(i,nr-1)->Boundary.type = boundaryType::symmetry;
            G.cell(i,nr-1)->Boundary.boundary_edge = edge::right;
        }
    cout<<"6"<<endl;
        for (int i=0; i<nz*dg; i++)
            G.cell(i,0)->Boundary.type = boundaryType::fixedpressure;


        for (int j=1; j<nr-1; j++)
        {
            G.cell(0,j)->Boundary.type = boundaryType::symmetry;
            G.cell(0,j)->Boundary.boundary_edge = edge::up;
            G.cell(nz-1,j)->Boundary.type = boundaryType::fixedmoisture;
            G.cell(nz-1,j)->Boundary.value = 0.4;
        }
    cout<<"7"<<endl;

        G.SetProp("pond_alpha","1000.0");
        G.SetProp("pond_beta","2");
    if (mode == _mode::heterogeneous)
    {   G.SetProp("alpha","30");
        G.SetProp("n","1.5");
    }
    cout<<"7.1"<<endl;
        G.SetProp("inflow", "/home/arash/Projects/Drywell_Result/inflow_zero.txt");
    cout<<"7.2"<<endl;
        G.SetProp("well_H","0");
    cout<<"7.3"<<endl;
        G.SetProp("well_H_old","0");
    cout<<"7.4"<<endl;
        G.SetProp("r_w","0.025");
    cout<<"8"<<endl;
        G.Solve(0,0.00001,20,0.1);
    cout<<"9"<<endl;
        G.ExtractMoisture(1,0).writefile(Results_Folder + "/moist_well.csv");
    cout<<"9"<<endl;
        G.WriteResults(Results_Folder + "/theta.vtp");
    cout<<"9.1"<<endl;
        G.WriteResults("alpha", Results_Folder + "/alpha.vtp");
    cout<<"9.2"<<endl;
        G.WriteResults("Ks", Results_Folder + "/Ksat.vtp");
    cout<<"9.3"<<endl;
        G.WriteResults("n", Results_Folder + "/n.vtp");
    cout<<"9.4"<<endl;
        G.WaterDepth().make_uniform(0.1).writefile(Results_Folder + "/waterdepth.csv");
    cout<<"done!"<<endl;
    }
#ifdef VALGRIND
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
}
