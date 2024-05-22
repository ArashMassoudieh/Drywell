#ifndef PROPERTYGENERATOR_H
#define PROPERTYGENERATOR_H

#include <vector>
#include <Matrix.h>
#include <Vector.h>
#include <gsl/gsl_rng.h>
#include "BTC.h"

using namespace std;

struct correl_mat_vec
{
#ifdef  _arma
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

enum class params {K_sat, alpha, n};

class PropertyGenerator: public vector<Propfull>
{
public:
    PropertyGenerator();
    PropertyGenerator(unsigned int n);
    double dx;
    double correlation_length_scale;
    void assign_K_gauss();
    bool write(const string &quan, const string &filename) const;
    void Normalize_Ksat_normal_scores(const double &mean, const double &std);
    double mean(const string &quan, bool log=false) const;
    double std(const string &quan, bool log=false) const;
    vector<double> vals(const string &quan) const;
    void SetMarginalDistribution(const string &quan, const CTimeSeries<double> series);
    CTimeSeries<double> MarginalDistribution(const string &quan);
    void PopulateRealValue(const string &quan, const string &quanfrom);
    double val(const string &quan, int i) const;
    bool SetVal(const string &quan, int i, const double &value);
    void SetCorr(params, const double &value);
    void Populate_Alpha_n_normal_scores(params p);
    void Normalize(const string &quan, const double &denominator);
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
    map<string,CTimeSeries<double>> marginal_distributions;
    double K_sat_alpha_correlation = 1;
    double K_sat_n_correlation = 1;
    const gsl_rng_type * T;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

};

#endif // PROPERTYGENERATOR_H
