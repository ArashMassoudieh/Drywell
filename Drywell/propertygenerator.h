#ifndef PROPERTYGENERATOR_H
#define PROPERTYGENERATOR_H

#include <vector>
#include <Matrix.h>
#include <Vector.h>
#include <gsl/gsl_rng.h>

using namespace std;

struct correl_mat_vec
{
#ifdef  arma
    CMatrix_arma M_22;
    CVector_arma V_21;
    CVector_arma V_RHS;
#else
    CMatrix M_22;
    CVector V_21;
    CVector V_RHS;
#endif //  arma
};

struct ival
{
    int i;
    double val;
};

struct Prop
{
    double K_sat;
    double alpha;
    double n;
};

struct Propfull
{
    Prop realvalues;
    Prop normal_scores;
    bool k_det = false;
};

class PropertyGenerator: public vector<Propfull>
{
public:
    PropertyGenerator();
    PropertyGenerator(unsigned int n);
    double dx;
    double correlation_length_scale;
    void assign_K_gauss();
private:
    CMatrix K_alpha_n_corr_matrix;
    double K_sat_normal_score_mean;
    double K_sat_normal_score_std;
    correl_mat_vec get_correll_matrix_vec(int i);
    vector<ival> get_top_n(const vector<ival> &vec);
    vector<ival> get_closest_K_dets(unsigned int i);
    void assign_K_gauss(unsigned int i);
    unsigned int GetNumberOfPointsDetermined();
    vector<int> Determined();

    const gsl_rng_type * T;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

};

#endif // PROPERTYGENERATOR_H
