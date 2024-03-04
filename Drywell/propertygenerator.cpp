#include "propertygenerator.h"
#include "gsl/gsl_randist.h"

PropertyGenerator::PropertyGenerator():vector<Propfull>()
{

}


PropertyGenerator::PropertyGenerator(unsigned int n):vector<Propfull>(n)
{

}


correl_mat_vec PropertyGenerator::get_correll_matrix_vec(int i)
{
    correl_mat_vec Correl_Matrix_Vector;
    int num_determined = GetNumberOfPointsDetermined();
    vector<int> determined = Determined();
#ifdef arma
    M.M_22 = CMatrix_arma(num_determined);
    M.V_21 = CVector_arma(num_determined);
    M.V_RHS = CVector_arma(num_determined);
#else
    Correl_Matrix_Vector.M_22 = CMatrix(num_determined);
    Correl_Matrix_Vector.V_21 = CVector(num_determined);
    Correl_Matrix_Vector.V_RHS = CVector(num_determined);
#endif //  arma
    for (int ii = 0; ii < num_determined; ii++)
    {
        Correl_Matrix_Vector.V_21[ii] = exp(-fabs(determined[ii]-i)*dx/correlation_length_scale);
        Correl_Matrix_Vector.V_RHS[ii] = at(determined[ii]).normal_scores.K_sat;
        for (int jj = 0; jj < num_determined; jj++)
        {
#ifdef arma
            Correl_Matrix_Vector.M_22(ii,jj) = exp(-fabs(determined[ii]-jj)*dx/correlation_length_scale);
#else
            Correl_Matrix_Vector.M_22[ii][jj] = exp(-fabs(determined[ii]-jj)*dx/correlation_length_scale);
#endif // arma
        }
    }
    return Correl_Matrix_Vector;
}

void PropertyGenerator::assign_K_gauss(unsigned int i)
{
    correl_mat_vec M = get_correll_matrix_vec(i);
    double mu;
    double sigma;

#ifdef arma
        CMatrix_arma M_inv = inv(M.M_22);
#else
        CMatrix M_inv = Invert(M.M_22);
#endif
        mu = dotproduct(M_inv*M.V_21, M.V_RHS);
        sigma = 1.0 - dotproduct(M_inv*M.V_21, M.V_21);

    double K_gauss = mu + gsl_ran_ugaussian(r)*sigma;
    at(i).k_det = true;
    at(i).normal_scores.K_sat = K_gauss;

}


void PropertyGenerator::assign_K_gauss()
{
    unsigned int n_filled = 0;
    srand(time(NULL));
    int i = gsl_rng_uniform(r)*(size()-1) + 0.5;
    at(i).normal_scores.K_sat = gsl_ran_ugaussian(r);
    at(i).k_det = true;
    n_filled++;
    while (n_filled<size())
    {
        i = gsl_rng_uniform(r)*(size()-1) + 0.5;
        if (!at(i).k_det)
        {
            assign_K_gauss(i);
            n_filled++;
        }
    }

}


unsigned int PropertyGenerator::GetNumberOfPointsDetermined()
{
    int num_determined = 0;
    for (unsigned int i=0; i<size(); i++)
        if (at(i).k_det)
            num_determined++;
    return num_determined;
}

vector<int> PropertyGenerator::Determined()
{
    vector<int> determined;
    for (unsigned int i=0; i<size(); i++)
        if (at(i).k_det)
            determined.push_back(i);
    return determined;
}

