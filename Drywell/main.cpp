

#include <iostream>
#include "grid.h"
#include "propertygenerator.h"
#include <QString>
#include <QMap>
#include <QDir>
#include <QDebug>
#include <QTime>

#ifdef VALGRIND
#include <valgrind/callgrind.h>
#endif

//using namespace std;

struct ModelParameters{
    QString alpha;
    QString n;
    QString folder;
    QString rw;
    QString dg;
};


int main()
{
    QMap<QString, ModelParameters> parameter_set;

    /*for (int i=0; i<10; i++)
    {   ModelParameters case1;
        case1.alpha = QString::number((i+1)*5);
        case1.n = "2";
        case1.dg = "0.5";
        case1.rw = "0";
        case1.folder = "Alpha_" + QString::number((i+1)*5);
        parameter_set[case1.folder] = case1;
    }*/

    for (int i=0; i<10; i++)
    {   ModelParameters case1;
        case1.alpha = "30";
        case1.n = "2";
        case1.dg = "0.5";
        case1.rw = "0";
        case1.folder = "Realization_" + QString::number(i+1);
        parameter_set[case1.folder] = case1;
    }



    /*
    ModelParameters case2;
    case2.alpha = "10";
    case2.n = "2";
    case2.dg = "0.5";
    case2.rw = "0";
    case2.folder = "Alpha=10, n=2";
    parameter_set[case2.folder] = case2;

    ModelParameters case3;
    case3.alpha = "15";
    case3.n = "2";
    case3.dg = "0.5";
    case3.rw = "0";
    case3.folder = "Alpha=15, n=2";
    parameter_set[case3.folder] = case3;

    ModelParameters case4;
    case4.alpha = "15";
    case4.n = "2";
    case4.dg = "0.7";
    case4.rw = "0";
    case4.folder = "Alpha=15, n=2, dg=0p7";
    parameter_set[case4.folder] = case4;

    ModelParameters case5;
    case5.alpha = "15";
    case5.n = "2";
    case5.dg = "0.5";
    case5.rw = "0.5";
    case5.folder = "Alpha=15, n=2, rw=0p5";
    parameter_set[case5.folder] = case5;
    */
    enum class _mode {homogeneous, heterogeneous};
    _mode mode = _mode::heterogeneous;
    int nz=20;
    int nr=20;

    for (QMap<QString, ModelParameters>::iterator it = parameter_set.begin(); it!= parameter_set.end(); it++)
    {
        QTime start_time = QTime::currentTime();
        cout<<"Start time = " + start_time.toString().toStdString();
        std::string Results_Folder = "/home/arash/Projects/Drywell_Result/Homogeneous/" + it->folder.toStdString();
        //std::string Results_Folder = "F:/Projects/Drywell_Result";
        QDir dir("/home/arash/Projects/Drywell_Result/Homogeneous");
        std::string SoilDataFolder = "/home/arash/Projects/Drywell_Result/Heterogeneous";
        if (!dir.exists(it->folder))
        {
            bool res = dir.mkdir(it->folder);
            qDebug()<<res;
        }
        PropertyGenerator P(nz);
        P.correlation_length_scale = 1;
        P.dx = 0.2;
        P.assign_K_gauss(); // Assigns normal scores for K_sat
        P.Normalize_Ksat_normal_scores(0,0.5);
        P.write("K_sat_normal_score", Results_Folder + "/K_sat_score_1.txt");
    cout<<"1"<<std::endl;
        CTimeSeries<double> K_sat_marginal_CDF(SoilDataFolder + "/K_sat_Marginal_Distribution.csv"); // Loading Cummulative Distributions
        CTimeSeries<double> alpha_marginal_CDF(SoilDataFolder + "/alpha_Marginal_Distribution.csv");
        CTimeSeries<double> n_marginal_CDF(SoilDataFolder + "/n_Marginal_Distribution.csv");
        P.SetCorr(params::alpha, 0.378); //Correlation between alpha and K_sat normal scores
        P.SetCorr(params::n, 0.2); //Correlation between n and K_sat normal scores
    cout<<"2"<<std::endl;
        P.SetMarginalDistribution("K_sat",K_sat_marginal_CDF); //Assigning CDFs to properties
    cout<<"2.1"<<std::endl;
        P.SetMarginalDistribution("alpha",alpha_marginal_CDF);
    cout<<"2.2"<<std::endl;
        P.SetMarginalDistribution("n",n_marginal_CDF);
    cout<<"2.3"<<std::endl;
        P.PopulateRealValue("K_sat","K_sat_normal_score"); //Assign the actual K_sat
    cout<<"2.4"<<std::endl;
        P.Normalize("K_sat",P.mean("K_sat",true)); //Normalize K_sat by the geometrical mean of K_sat so the geometrical mean is zero
    cout<<"3"<<std::endl;
        P.Populate_Alpha_n_normal_scores(params::alpha); //Creates normal scores for alpha
        P.Populate_Alpha_n_normal_scores(params::n); //Creates normal scores for n
        P.PopulateRealValue("alpha","alpha_normal_score"); //Assign the actual values of alpha
        P.Normalize("alpha",0.05);
        P.PopulateRealValue("n","n_normal_score"); //.. n
        P.write("K_sat", Results_Folder + "/K_sat.txt");
        P.write("alpha", Results_Folder + "/alpha.txt");
        P.write("n", Results_Folder + "/n.txt");
    cout<<"4"<<std::endl;
        double dg=it->dg.toDouble();
        Grid G(nz,nr,1/dg,1.5);
        if (mode == _mode::heterogeneous)
            G.AssignProperty(&P); // Assign K_sat, alpha, n based on the Property Generator
    cout<<"5"<<std::endl;
        for (int i=0; i<nz; i++)
        {
            G.cell(i,0)->Boundary.type = boundaryType::symmetry;
            G.cell(i,0)->Boundary.boundary_edge = edge::left;
            G.cell(i,nr-1)->Boundary.type = boundaryType::symmetry;
            G.cell(i,nr-1)->Boundary.boundary_edge = edge::right;
        }
    cout<<"6"<<std::endl;
        for (int i=0; i<nz*dg; i++)
            G.cell(i,0)->Boundary.type = boundaryType::fixedpressure;


        for (int j=1; j<nr-1; j++)
        {
            G.cell(0,j)->Boundary.type = boundaryType::symmetry;
            G.cell(0,j)->Boundary.boundary_edge = edge::up;
            G.cell(nz-1,j)->Boundary.type = boundaryType::fixedmoisture;
            G.cell(nz-1,j)->Boundary.value = 0.4;
        }
    cout<<"7"<<std::endl;

        G.SetProp("pond_alpha","1000.0");
        G.SetProp("pond_beta","2");
    if (mode == _mode::homogeneous)
    {   G.SetProp("alpha",it->alpha.toStdString());
        G.SetProp("n",it->n.toStdString());
    }

    cout<<"7.1"<<std::endl;
        G.SetProp("inflow", "/home/arash/Projects/Drywell_Result/inflow_zero.txt");
    cout<<"7.2"<<std::endl;
        G.SetProp("well_H","0");
    cout<<"7.3"<<std::endl;
        G.SetProp("well_H_old","0");
    cout<<"7.4"<<std::endl;
        G.SetProp("r_w",it->rw.toStdString());
    cout<<"8"<<std::endl;
        G.Solve(0,0.000001,1,0.005);
    cout<<"9"<<std::endl;
        G.ExtractMoisture(1,0).writefile(Results_Folder + "/moist_well.csv");
    cout<<"9"<<std::endl;
        G.WriteResults(Results_Folder + "/theta.vtp");
    cout<<"9.1"<<std::endl;
        G.WriteResults("alpha", Results_Folder + "/alpha.vtp");
    cout<<"9.2"<<std::endl;
        G.WriteResults("Ks", Results_Folder + "/Ksat.vtp");
    cout<<"9.3"<<std::endl;
        G.WriteResults("n", Results_Folder + "/n.vtp");
    cout<<"9.4"<<endl;
        G.WaterDepth().make_uniform(0.01).writefile(Results_Folder + "/waterdepth.csv");
        G.TotalWaterInSoil().make_uniform(0.01).writefile(Results_Folder + "/totalwaterinsoil.csv");
        G.TotalWaterInWell().make_uniform(0.01).writefile(Results_Folder + "/totalwaterinwell.csv");
    cout<<"done!"<<endl;
    QTime end_time = QTime::currentTime();
    cout<<"End time = " + end_time.toString().toStdString();
    cout<<"Simulation Time = "<<(start_time.secsTo(end_time))/60.0;

    }


#ifdef VALGRIND
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
}
