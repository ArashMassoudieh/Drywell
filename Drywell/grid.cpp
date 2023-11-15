#include "grid.h"
#include "vtk.h"
#include "Utilities.h"

#define pi 3.14159265

Grid::Grid()
{

}

Grid::Grid(int _nz, int _nr, double depth, double radius)
{
    nz = _nz;
    nr = _nr;
    dz = depth/nz;
    dr = radius/nr;
    cells.resize(nz);
    for (unsigned int i=0; i<nz; i++)
        cells[i].resize(nr);
}

CVector_arma Grid::Residual(const CVector_arma &X, const double &dt)
{
    CVector_arma Res(nr*nz+1);
    SetStateVariable(X);
    CMatrix Pressure = H();
    CMatrix Saturation = Se();
    Res[nz*nr] = (2*pi*(r_w+dr/2)*(well_H<0?1:0)+alpha*beta*pow(max(well_H,0.0),beta-1))*(well_H-well_H_old)/dt - inflow.interpol(Solution_State.t);
    for (unsigned int i=0; i<nz; i++)
    {   double lower_el = -dz*(i+1);
        double interface_length = max(min(well_H - lower_el,dz),0.0);
        for (unsigned int j=0; j<nr; j++)
        {
            if (cells[i][j].Boundary.type==boundaryType::none)
            {
                double r = (j+0.5)*dr+r_w;
                Res[j+nr*i] = (cells[i][j].Theta(_time::current)-cells[i][j].Theta(_time::past))/dt;
                Res[j+nr*i] += (r-dr/2)/(r*pow(dr,2))*K(i,j,edge::left)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::left)->H(_time::current));
                Res[j+nr*i] += (r+dr/2)/(r*pow(dr,2))*K(i,j,edge::right)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::right)->H(_time::current));
                Res[j+nr*i] += 1/pow(dz,2)*K(i,j,edge::up)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::up)->H(_time::current));
                Res[j+nr*i] += 1/pow(dz,2)*K(i,j,edge::down)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::down)->H(_time::current));
                Res[j+nr*i] += -1/dz*(K(i,j,edge::up)-K(i,j,edge::down));
            }
            else if (cells[i][j].Boundary.type==boundaryType::fixedmoisture)
            {
                Res[j+nr*i] = cells[i][j].Theta(_time::current) - cells[i][j].Boundary.value;
            }
            else if (cells[i][j].Boundary.type==boundaryType::fixedpressure)
            {
                Res[j+nr*i] = cells[i][j].H(_time::current) - max(well_H+(i+1)*dz,0.0)*(interface_length/dz);
                Res[nz*nr]+=2*interface_length*pi*(r_w+dr/2)*K(i,j,edge::right)*(max(well_H+(i+1)*dz,0.0)-Neighbour(i,j,edge::right)->H(_time::current));
            }
            else if (cells[i][j].Boundary.type==boundaryType::gradient)
            {
                if (cells[i][j].Boundary.boundary_edge == edge::down || cells[i][j].Boundary.boundary_edge == edge::up)
                    Res[j+nr*i] = (cells[i][j].Theta(_time::current)-Neighbour(i,j,cells[i][j].Boundary.boundary_edge, true)->Theta(_time::current))/dz;
                else
                    Res[j+nr*i] = (cells[i][j].Theta(_time::current)-Neighbour(i,j,cells[i][j].Boundary.boundary_edge, true)->Theta(_time::current))/dr;
            }
            else if (cells[i][j].Boundary.type==boundaryType::symmetry)
            {
                double r = (j+0.5)*dr;
                Res[j+nr*i] = (cells[i][j].Theta(_time::current)-cells[i][j].Theta(_time::past))/dt;
                if (cells[i][j].Boundary.boundary_edge!=edge::left && Neighbour(i,j,edge::left))
                {
                    Res[j+nr*i] += (r-dr/2)/(r*pow(dr,2))*K(i,j,edge::left)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::left)->H(_time::current));
                }
                if (cells[i][j].Boundary.boundary_edge!=edge::right && Neighbour(i,j,edge::right))
                {
                    Res[j+nr*i] += (r+dr/2)/(r*pow(dr,2))*K(i,j,edge::right)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::right)->H(_time::current));
                }
                if (cells[i][j].Boundary.boundary_edge!=edge::up && Neighbour(i,j,edge::up))
                {
                    Res[j+nr*i] += 1/pow(dz,2)*K(i,j,edge::up)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::up)->H(_time::current));
                }
                if (cells[i][j].Boundary.boundary_edge!=edge::down && Neighbour(i,j,edge::down))
                {   Res[j+nr*i] += 1/pow(dz,2)*K(i,j,edge::down)*(cells[i][j].H(_time::current)-Neighbour(i,j,edge::down)->H(_time::current));
                    Res[j+nr*i] += -1/dz*(K(i,j,edge::up)-K(i,j,edge::down));
                }

            }

        }
    }
    return Res;
}

void Grid::SetStateVariable(const CVector_arma &X,const _time &t)
{
    for (unsigned int i=0; i<nz; i++)
        for (unsigned int j=0; j<nr; j++)
        {
            if (cells[i][j].Boundary.type == boundaryType::fixedpressure)
            {
                cells[i][j].Boundary.value = min(X[nz*nr]-i*dz,0.0);
            }
            cells[i][j].SetTheta(X[j+nr*i],t);

        }
    if (t==_time::current)
        well_H = X[nz*nr];
    else
        well_H_old = X[nz*nr];

}

CVector_arma Grid::GetStateVariable(const _time &t) const
{
    CVector_arma X(nr*nz+1);
    for (unsigned int i=0; i<nz; i++)
        for (unsigned int j=0; j<nr; j++)
            X[j+nr*i] = cells[i][j].Theta(t);

    X[nr*nz] = well_H;
    return X;
}

double Grid::D(int i,int j,const edge &ej)
{
    return K(i,j,ej)*invC(i,j,ej);
}

double Grid::K(int i,int j,const edge &ej)
{
    Cell* neighbour = Neighbour(i,j,ej);
    if (!neighbour)
        neighbour = &cells[i][j];
    double Se = min(max((max(cells[i][j].Theta(_time::current),neighbour->Theta(_time::current))-getVal(i,j,"theta_r",ej))/(getVal(i,j,"theta_s",ej)-getVal(i,j,"theta_r",ej)),1e-6),0.99999);
    double m = 1.0-1.0/getVal(i,j,"n",ej);
    double K = pow(Se,0.5)*getVal(i,j,"Ks",ej)*pow(1-pow(1-pow(Se,1.0/m),m),2);
    return K;

}

double Grid::invC(int i,int j,const edge &ej)
{
    double C;
    double Se = min(max((max(cells[i][j].Theta(_time::current),Neighbour(i,j,ej)->Theta(_time::current))-getVal(i,j,"theta_r",ej))/(getVal(i,j,"theta_s",ej)-getVal(i,j,"theta_r",ej)),1e-6),0.99999);
    double m = 1.0-1.0/getVal(i,j,"n",ej);
    if (getVal(i,j,"theta",ej)>getVal(i,j,"theta_s",ej))
        C = 1.0/getVal(i,j,"epsilon",ej);
    else
        C = 1.0/getVal(i,j,"alpha",ej)/(getVal(i,j,"theta_s",ej)-getVal(i,j,"theta_r",ej))*pow(pow(Se,-m)-1,-m)*pow(Se,-m-1);
    return C;
}

double Grid::getVal(int i, int j, const string &val, const edge &ej)
{
    Cell* neighbour = Neighbour(i,j,ej);
    if (!neighbour)
        neighbour = &cells[i][j];

    return 0.5*(cells[i][j].getValue(val)+neighbour->getValue(val));

}

Cell* Grid::Neighbour(int i, int j, const edge &ej, bool op)
{
    if (op)
    {   if (ej == edge::down)
        {
            if (i>0)
                return &cells[i-1][j];
            else
                return nullptr;
        }
        else if (ej == edge::up)
        {
            if (i<nz-1)
                return &cells[i+1][j];
            else
                return nullptr;

        }
        else if (ej == edge::left)
        {
            if (j<nr-1)
                return &cells[i][j+1];
            else
                return nullptr;

        }
        else if (ej == edge::right)
        {
            if (j>0)
                return &cells[i][j-1];
            else
                return nullptr;
        }
    }
    else
    {   if (ej == edge::down)
        {
            if (i<nz-1)
                return &cells[i+1][j];
            else
                return nullptr;
        }
        else if (ej == edge::up)
        {
            if (i>0)
                return &cells[i-1][j];
            else
                return nullptr;
        }
        else if (ej == edge::left)
        {
            if (j>0)
                return &cells[i][j-1];
            else
                return nullptr;

        }
        else if (ej == edge::right)
        {
            if (j<nr-1)
                return &cells[i][j+1];
            else
                return nullptr;

        }
    }

}



CMatrix_arma Grid::Jacobian(const CVector_arma &X, const double &dt)
{
    CVector_arma F_base = Residual(X,dt);
    CMatrix_arma M(X.getsize(),X.getsize());
    for (unsigned int i=0; i<X.getsize(); i++)
    {
        CVector_arma X1 = X;
        X1[i]+=1e-6;
        CVector_arma F1 = Residual(X1,dt);
        for (unsigned int j=0; j<X.num; j++)
            M(j,i) = (F1[j]-F_base[j])/1e-6;
    }
    return M;
}

bool Grid::OneStepSolve(const double &dt)
{
    CVector_arma X = GetStateVariable();
    CVector_arma X0 = X;
    CVector_arma Res = Residual(X,dt);
    double err_0 = Res.norm2();
    double err=err_0*10;
    Solution_State.number_of_iterations = 0;
    int count_error_expanding = 0;

    while (err/(err_0+1e-3)>1e-3)
    {
        CMatrix_arma J = Jacobian(X,dt);
        CVector_arma dx = Res/J;
        X-=dx;
        Res = Residual(X,dt);
        double err1=Res.norm2();
        if (err1>err)
        {
            count_error_expanding++;
        }
        err = err1;

        Solution_State.number_of_iterations ++;
        if (Solution_State.number_of_iterations>Solution_State.max_iterations || count_error_expanding>5 || !(err==err))
        {
            //cout<<"Interations: "<<Solution_State.number_of_iterations<<", Error Expanding: "<<count_error_expanding<<", Error: "<<err<< ",Ini Error:"<< err_0<<" dt: "<<dt<< endl;
            SetStateVariable(X0);
            return false;
        }
    }
    return true;
}

bool Grid::Solve(const double &t0, const double &dt0, const double &t_end, const double &write_interval)
{
    Solution_State.t = t0;
    Solution_State.dt = dt0;
    Well_Water_Depth.append(Solution_State.t,well_H);
    while (Solution_State.t + Solution_State.dt<t_end)
    {
        if (!OneStepSolve(Solution_State.dt))
        {
            Solution_State.dt*=Solution_State.dt_scale_factor_fail;
        }
        else
        {
            if (floor(Solution_State.t/write_interval)<floor((Solution_State.t+Solution_State.dt)/write_interval))
            {
                double write_time = floor((Solution_State.t+Solution_State.dt)/write_interval)*write_interval;
                CMatrix snapshot = ((write_time-Solution_State.t)*Theta(_time::past) + (Solution_State.t+Solution_State.dt-write_time)*Theta(_time::current))/Solution_State.dt;
                results.push_back(snapshot);
            }
            Solution_State.t += Solution_State.dt;
            Well_Water_Depth.append(Solution_State.t,well_H);


            if (Solution_State.number_of_iterations>Solution_State.NI_max && Solution_State.dt>1e-5)
            {
                Solution_State.dt*=Solution_State.dt_scale_factor;
            }
            else if (Solution_State.number_of_iterations<Solution_State.NI_min)
            {
                Solution_State.dt/=Solution_State.dt_scale_factor;
            }
        }
    CVector_arma X = GetStateVariable(_time::current);
    SetStateVariable(X,_time::past);

    cout<<Solution_State.t<<endl;
    }
    return true;
}

bool Grid::SetProp(const string &propname, const string &value)
{
    if (propname == "well_H")
    {
        well_H = aquiutils::atof(value);
        return true;
    }
    else if (propname == "well_H_old")
    {
        well_H_old = aquiutils::atof(value);
        return true;
    }
    else if (propname == "r_w")
    {
        r_w = aquiutils::atof(value);
        return true;
    }
    else if (propname == "beta")
    {
        beta = aquiutils::atof(value);
        return true;
    }
    else if (propname == "alpha")
    {
        alpha = aquiutils::atof(value);
        return true;
    }
    else if (propname == "inflow")
    {
        inflow = CTimeSeries<double>(value);
        return true;
    }
    return false;

}

#ifdef use_VTK
void Grid::write_to_vtp(const string &name) const
{
    vtkSmartPointer<vtkPoints> points_3 =
            vtkSmartPointer<vtkPoints>::New();

        double xx, yy, zz;
        vtkSmartPointer<vtkFloatArray> values =
            vtkSmartPointer<vtkFloatArray>::New();

        values->SetNumberOfComponents(1);

        values->SetName("Moisture Content");



        for (unsigned int i = 0; i < cells.size(); i++)
            for (unsigned int j = 0; j < cells[i].size(); j++)
            {
                yy = -dz*(i+0.5);
                xx = dr*(j+0.5)+r_w;
                zz = cells[i][j].Theta(_time::current);


                float t[1] = { float(zz) };
                points_3->InsertNextPoint(xx, yy, zz);
                values->InsertNextTupleValue(t);

            }


        // Add the grid points to a polydata object
        vtkSmartPointer<vtkPolyData> inputPolyData =
            vtkSmartPointer<vtkPolyData>::New();
        inputPolyData->SetPoints(points_3);

        // Triangulate the grid points
        vtkSmartPointer<vtkDelaunay2D> delaunay =
            vtkSmartPointer<vtkDelaunay2D>::New();
    #if VTK_MAJOR_VERSION <= 5
        delaunay->SetInput(inputPolyData);
    #else
        delaunay->SetInputData(inputPolyData);
    #endif
        delaunay->Update();
        vtkPolyData* outputPolyData = delaunay->GetOutput();

        double bounds[6];
        outputPolyData->GetBounds(bounds);

        // Find min and max z
        double minz = bounds[4];
        double maxz = bounds[5];

        // Create the color map
        vtkSmartPointer<vtkLookupTable> colorLookupTable =
            vtkSmartPointer<vtkLookupTable>::New();
        colorLookupTable->SetTableRange(minz, maxz);
        colorLookupTable->Build();

        // Generate the colors for each point based on the color map
        vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
            vtkSmartPointer<vtkUnsignedCharArray>::New();
        colors_2->SetNumberOfComponents(3);
        colors_2->SetName("Colors");


        for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
        {
            double p[3];
            outputPolyData->GetPoint(i, p);

            double dcolor[3];
            colorLookupTable->GetColor(p[2], dcolor);

            unsigned char color[3];
            for (unsigned int j = 0; j < 3; j++)
            {
                color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
            }
            //std::cout << "color: "
            //	<< (int)color[0] << " "
            //	<< (int)color[1] << " "
            //	<< (int)color[2] << std::endl;

            colors_2->InsertNextTupleValue(color);
        }

        outputPolyData->GetPointData()->SetScalars(values);


        //Append the two meshes
        vtkSmartPointer<vtkAppendPolyData> appendFilter =
            vtkSmartPointer<vtkAppendPolyData>::New();
    #if VTK_MAJOR_VERSION <= 5
        appendFilter->AddInputConnection(input1->GetProducerPort());
        appendFilter->AddInputConnection(input2->GetProducerPort());
    #else
        //appendFilter->AddInputData(polydata);
        //appendFilter->AddInputData(polydata_1);
        appendFilter->AddInputData(outputPolyData);
    #endif
        appendFilter->Update();


        // Visualization
        vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
        mapper->SetInputConnection(polydata->GetProducerPort());
    #else
        mapper->SetInputConnection(appendFilter->GetOutputPort());
        //mapper->SetInputData(polydata_1);
    #endif

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetPointSize(5);

        vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(name.c_str());
        writer->SetInputData(mapper->GetInput());
        // This is set so we can see the data in a text editor.
        writer->SetDataModeToAscii();
        writer->Write();

}


void Grid::write_to_vtp(const string &name,const CMatrix &res) const
{
    vtkSmartPointer<vtkPoints> points_3 =
            vtkSmartPointer<vtkPoints>::New();

        double xx, yy, zz;
        vtkSmartPointer<vtkFloatArray> values =
            vtkSmartPointer<vtkFloatArray>::New();

        values->SetNumberOfComponents(1);

        values->SetName("Moisture Content");



        for (unsigned int i = 0; i < cells.size(); i++)
            for (unsigned int j = 0; j < cells[i].size(); j++)
            {
                yy = -dz*(i+0.5);
                xx = dr*(j+0.5);
                zz = res[i][j];


                float t[1] = { float(zz) };
                points_3->InsertNextPoint(xx, yy, zz);
                values->InsertNextTupleValue(t);

            }


        // Add the grid points to a polydata object
        vtkSmartPointer<vtkPolyData> inputPolyData =
            vtkSmartPointer<vtkPolyData>::New();
        inputPolyData->SetPoints(points_3);

        // Triangulate the grid points
        vtkSmartPointer<vtkDelaunay2D> delaunay =
            vtkSmartPointer<vtkDelaunay2D>::New();
    #if VTK_MAJOR_VERSION <= 5
        delaunay->SetInput(inputPolyData);
    #else
        delaunay->SetInputData(inputPolyData);
    #endif
        delaunay->Update();
        vtkPolyData* outputPolyData = delaunay->GetOutput();

        double bounds[6];
        outputPolyData->GetBounds(bounds);

        // Find min and max z
        double minz = bounds[4];
        double maxz = bounds[5];

        // Create the color map
        vtkSmartPointer<vtkLookupTable> colorLookupTable =
            vtkSmartPointer<vtkLookupTable>::New();
        colorLookupTable->SetTableRange(minz, maxz);
        colorLookupTable->Build();

        // Generate the colors for each point based on the color map
        vtkSmartPointer<vtkUnsignedCharArray> colors_2 =
            vtkSmartPointer<vtkUnsignedCharArray>::New();
        colors_2->SetNumberOfComponents(3);
        colors_2->SetName("Colors");


        for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
        {
            double p[3];
            outputPolyData->GetPoint(i, p);

            double dcolor[3];
            colorLookupTable->GetColor(p[2], dcolor);

            unsigned char color[3];
            for (unsigned int j = 0; j < 3; j++)
            {
                color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
            }
            //std::cout << "color: "
            //	<< (int)color[0] << " "
            //	<< (int)color[1] << " "
            //	<< (int)color[2] << std::endl;

            colors_2->InsertNextTupleValue(color);
        }

        outputPolyData->GetPointData()->SetScalars(values);


        //Append the two meshes
        vtkSmartPointer<vtkAppendPolyData> appendFilter =
            vtkSmartPointer<vtkAppendPolyData>::New();
    #if VTK_MAJOR_VERSION <= 5
        appendFilter->AddInputConnection(input1->GetProducerPort());
        appendFilter->AddInputConnection(input2->GetProducerPort());
    #else
        //appendFilter->AddInputData(polydata);
        //appendFilter->AddInputData(polydata_1);
        appendFilter->AddInputData(outputPolyData);
    #endif
        appendFilter->Update();


        // Visualization
        vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
        mapper->SetInputConnection(polydata->GetProducerPort());
    #else
        mapper->SetInputConnection(appendFilter->GetOutputPort());
        //mapper->SetInputData(polydata_1);
    #endif

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetPointSize(5);

        vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(name.c_str());
        writer->SetInputData(mapper->GetInput());
        // This is set so we can see the data in a text editor.
        writer->SetDataModeToAscii();
        writer->Write();

}
#endif

void Grid::WriteResults(const string &filename)
{
    for (unsigned k=0; k<results.size(); k++)
    {
        string name = aquiutils::split(filename,'.')[0]+"_"+aquiutils::numbertostring(k+1)+".vtp";
        write_to_vtp(name,results[k]);
    }
}

CMatrix Grid::H()
{
    CMatrix out(nz,nr);
    for (unsigned int i = 0; i < cells.size(); i++)
        for (unsigned int j = 0; j < cells[i].size(); j++)
        {
            out[i][j] = cells[i][j].H(_time::current);
        }
    return out;
}

CMatrix Grid::Theta(_time t)
{
    CMatrix out(nz,nr);
    for (unsigned int i = 0; i < cells.size(); i++)
        for (unsigned int j = 0; j < cells[i].size(); j++)
        {
            out[i][j] = cells[i][j].Theta(t);
        }
    return out;
}

CMatrix Grid::Se()
{
    CMatrix out(nz,nr);
    for (unsigned int i = 0; i < cells.size(); i++)
        for (unsigned int j = 0; j < cells[i].size(); j++)
        {
            out[i][j] = (cells[i][j].Theta(_time::current)-cells[i][j].getValue("theta_r"))/(cells[i][j].getValue("theta_s")-cells[i][j].getValue("theta_r"));
        }
    return out;
}

