// written by Szymon Winczewski

#include "spherical_harmonics.h"

#include "memory.h"

using namespace std;
using namespace LAMMPS_NS;


double SphericalHarmonics::calculateFactorial(int i)
{
    if ( i == 0 )
        return 1;
    else if ( i == 1 )
        return 1;
    else if ( i == 2 )
        return 2;
    else if ( i == 3 )
        return 6;
    else if ( i == 4 )
        return 24;
    else if ( i == 5 )
        return 120;
    else if ( i == 6 )
        return 720;
    else if ( i == 7 )
        return 5040;
    else if ( i == 8 )
        return 40320;
    else if ( i == 9 )
        return 362880;
    else if ( i == 10 )
        return 3628800;
    else if ( i > 10 )
    {
        double fact = 3628800;
        for (int j = 11; j <= i; j++)
            fact *= j;
        return fact;
    }

    return -1;
}


double SphericalHarmonics::calculateDoubleFactorial(int i)
{
    if ( ( i == 0 ) || ( i == 1 ) )
        return 1.0;
    else
    {
        div_t div_res = div(i, 2);
        if ( div_res.rem == 0 )
            return pow(2, div_res.quot) * calculateFactorial(div_res.quot);
        else
            return calculateFactorial(i) / ( pow(2, div_res.quot) * calculateFactorial(div_res.quot) );
    }
}


double SphericalHarmonics::calculateLegendre(double x, int l, int m)
{
#ifndef SPHERICAL_HARMONICS_FAST_LEGENDRE
    if ( x > 1.0 )
        x = 1.0;
    else if ( x < - 1.0 )
        x = -1.0;
#endif
    if ( l == 0 )
        return 1.0;
    else if ( ( l == 1 ) && ( m == 0 ) )
        return x;
    else
    {
#ifndef SPHERICAL_HARMONICS_FAST_LEGENDRE
        if ( m < 0 )
        {
            m = -m;
            return pow(-1.0, m) * ( vec_fact_[l - m] / vec_fact_[l + m] ) * calculateLegendre(x, l, m);
        }
        else if ( m >= 0 )
        {
#endif
            if ( m == l )
                return pow(-1.0, l) * vec_dfact_[2 * l - 1] * pow((1.0 - x * x), l / 2.0);
            else if ( m == l - 1 )
                return x * (2 * m + 1) * calculateLegendre(x, m, m);
            else
                return (x * (2 * l - 1) * calculateLegendre(x, l - 1, m) - (l + m - 1) * calculateLegendre(x, l - 2, m)) / (l - m);

#ifndef SPHERICAL_HARMONICS_FAST_LEGENDRE
        }
#endif
    }

    return 0.0;
}


void SphericalHarmonics::allocateMemory()
{
    memory->create(vec_fact_, SPHERICAL_HARMONICS_MAX_FACTORIAL, "spherical_harmonics:vec_fact");
    memory->create(vec_dfact_, SPHERICAL_HARMONICS_MAX_DOUBLE_FACTORIAL, "spherical_harmonics:vec_dfact");

    memory->create(vec_Legendre_arg_, Legendre_grid_size_ + 1, "spherical_harmonics:vec_Legendre_arg_");
    memory->create(vec_Legendre_, Legendre_grid_size_ + 1, 24, "spherical_harmonics:vec_Legendre");
    
    memory->create(vec_phi_, phi_grid_size_ + 1, "spherical_harmonics:vec_phi");
    memory->create(vec_trig_mphi_, phi_grid_size_ + 1, 24, "spherical_harmonics:vec_trig_mphi");

    memory->create(ret_vec_, 24, "spherical_harmonics:ret_vec");
}


void SphericalHarmonics::deallocateMemory()
{
    memory->destroy(vec_fact_);
    memory->destroy(vec_dfact_);

    memory->destroy(vec_Legendre_arg_);
    memory->destroy(vec_Legendre_);

    memory->destroy(vec_phi_);
    memory->destroy(vec_trig_mphi_);


    memory->destroy(ret_vec_);
}


void SphericalHarmonics::initFactorials()
{
    for (int i = 0; i < SPHERICAL_HARMONICS_MAX_FACTORIAL; i++)
        vec_fact_[i] = calculateFactorial(i);
    for (int i = 0; i < SPHERICAL_HARMONICS_MAX_DOUBLE_FACTORIAL; i++)
        vec_dfact_[i] = calculateDoubleFactorial(i);
}


void SphericalHarmonics::initLegendre()
{
    int i, j, j_minus_one;

    Legendre_dx_ = 2.0 / Legendre_grid_size_;
    inv_Legendre_dx_ = 1.0 / Legendre_dx_;

    for (i = 0; i < Legendre_grid_size_; i++)
        vec_Legendre_arg_[i] = - 1.0 + i * Legendre_dx_;

    for (i = 0; i < Legendre_grid_size_; i++)
    {
        vec_Legendre_[i][0]   = K_4_0 * calculateLegendre(vec_Legendre_arg_[i], 4, 0);
        vec_Legendre_[i][2]   = K_6_0 * calculateLegendre(vec_Legendre_arg_[i], 6, 0);
        vec_Legendre_[i][4]   = K_4_1 * calculateLegendre(vec_Legendre_arg_[i], 4, 1);
        vec_Legendre_[i][6]   = K_6_1 * calculateLegendre(vec_Legendre_arg_[i], 6, 1);
        vec_Legendre_[i][8]   = K_4_2 * calculateLegendre(vec_Legendre_arg_[i], 4, 2);
        vec_Legendre_[i][10]  = K_6_2 * calculateLegendre(vec_Legendre_arg_[i], 6, 2);
        vec_Legendre_[i][12]  = K_4_3 * calculateLegendre(vec_Legendre_arg_[i], 4, 3);
        vec_Legendre_[i][14]  = K_6_3 * calculateLegendre(vec_Legendre_arg_[i], 6, 3);
        vec_Legendre_[i][16]  = K_4_4 * calculateLegendre(vec_Legendre_arg_[i], 4, 4);
        vec_Legendre_[i][18]  = K_6_4 * calculateLegendre(vec_Legendre_arg_[i], 6, 4);
        vec_Legendre_[i][20]  = K_6_5 * calculateLegendre(vec_Legendre_arg_[i], 6, 5);
        vec_Legendre_[i][22]  = K_6_6 * calculateLegendre(vec_Legendre_arg_[i], 6, 6);
    }

    for (j = 1; j < 24; j += 2)
    {
        j_minus_one = j - 1;
        for (i = 0; i < Legendre_grid_size_ - 1; i++)
            vec_Legendre_[i][j]  = ( vec_Legendre_[i + 1][j_minus_one]  - vec_Legendre_[i][j_minus_one] )  * inv_Legendre_dx_;
    }

    vec_Legendre_[Legendre_grid_size_ - 1][1]  = ( K_4_0 * calculateLegendre(1.0, 4, 0) - vec_Legendre_[Legendre_grid_size_ - 1][0] )  * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][3]  = ( K_6_0 * calculateLegendre(1.0, 6, 0) - vec_Legendre_[Legendre_grid_size_ - 1][2] )  * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][5]  = ( K_4_1 * calculateLegendre(1.0, 4, 1) - vec_Legendre_[Legendre_grid_size_ - 1][4] )  * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][7]  = ( K_6_1 * calculateLegendre(1.0, 6, 1) - vec_Legendre_[Legendre_grid_size_ - 1][6] )  * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][9]  = ( K_4_2 * calculateLegendre(1.0, 4, 2) - vec_Legendre_[Legendre_grid_size_ - 1][8] )  * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][11] = ( K_6_2 * calculateLegendre(1.0, 6, 2) - vec_Legendre_[Legendre_grid_size_ - 1][10] ) * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][13] = ( K_4_3 * calculateLegendre(1.0, 4, 3) - vec_Legendre_[Legendre_grid_size_ - 1][12] ) * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][15] = ( K_6_3 * calculateLegendre(1.0, 6, 3) - vec_Legendre_[Legendre_grid_size_ - 1][14] ) * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][17] = ( K_4_4 * calculateLegendre(1.0, 4, 4) - vec_Legendre_[Legendre_grid_size_ - 1][16] ) * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][19] = ( K_6_4 * calculateLegendre(1.0, 6, 4) - vec_Legendre_[Legendre_grid_size_ - 1][18] ) * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][21] = ( K_6_5 * calculateLegendre(1.0, 6, 5) - vec_Legendre_[Legendre_grid_size_ - 1][20] ) * inv_Legendre_dx_;
    vec_Legendre_[Legendre_grid_size_ - 1][23] = ( K_6_6 * calculateLegendre(1.0, 6, 6) - vec_Legendre_[Legendre_grid_size_ - 1][22] ) * inv_Legendre_dx_;

// problematyczny prawy koniec dziedziny
    vec_Legendre_arg_[Legendre_grid_size_] = vec_Legendre_arg_[Legendre_grid_size_ - 1];
    for (j = 0; j < 24; j++)
        vec_Legendre_[Legendre_grid_size_][j] = vec_Legendre_[Legendre_grid_size_ - 1][j];
}


void SphericalHarmonics::initTrigonometric()
{
    int i, j, j_minus_one;
    double phi_1, phi_2, phi_3, phi_4, phi_5, phi_6;

    dphi_ = TWO_PI / phi_grid_size_;
    inv_dphi_ = 1.0 / dphi_;

    for (i = 0; i < phi_grid_size_; i++)
        vec_phi_[i] = i * dphi_;

    for (i = 0; i < phi_grid_size_; i++)
    {
        phi_1 = 1.0 * vec_phi_[i];
        phi_2 = 2.0 * vec_phi_[i];
        phi_3 = 3.0 * vec_phi_[i];
        phi_4 = 4.0 * vec_phi_[i];
        phi_5 = 5.0 * vec_phi_[i];
        phi_6 = 6.0 * vec_phi_[i];

        vec_trig_mphi_[i][0]  = cos(phi_1);
        vec_trig_mphi_[i][2]  = sin(phi_1);
        vec_trig_mphi_[i][4]  = cos(phi_2);
        vec_trig_mphi_[i][6]  = sin(phi_2);
        vec_trig_mphi_[i][8]  = cos(phi_3);
        vec_trig_mphi_[i][10] = sin(phi_3);
        vec_trig_mphi_[i][12] = cos(phi_4);
        vec_trig_mphi_[i][14] = sin(phi_4);
        vec_trig_mphi_[i][16] = cos(phi_5);
        vec_trig_mphi_[i][18] = sin(phi_5);
        vec_trig_mphi_[i][20] = cos(phi_6);
        vec_trig_mphi_[i][22] = sin(phi_6);
    }

    for (j = 1; j < 24; j += 2)
    {
        j_minus_one = j - 1;
        for (i = 0; i < phi_grid_size_ - 1; i++)
            vec_trig_mphi_[i][j] = ( vec_trig_mphi_[i + 1][j_minus_one] - vec_trig_mphi_[i][j_minus_one] ) * inv_dphi_;
        vec_trig_mphi_[phi_grid_size_ - 1][j] = ( vec_trig_mphi_[0][j_minus_one] - vec_trig_mphi_[phi_grid_size_ - 1][j_minus_one] ) * inv_dphi_;
    }

// problematyczny prawy koniec dziedziny
    vec_phi_[phi_grid_size_] = vec_phi_[phi_grid_size_ - 1];
    for (j = 0; j < 24; j++)
        vec_trig_mphi_[phi_grid_size_][j] = vec_trig_mphi_[phi_grid_size_ - 1][j];
}


SphericalHarmonics::SphericalHarmonics(int n_knots_Legendre, int n_knots_trigonometric, Memory* memory)
 : memory(memory)
{
    /*
    if ( n_knots_Legendre < 2 )
        raiseError(ERR_SPHERICAL_HARMONICS, 1, "SphericalHarmonics", "incorrect number of knots for ALPs");

    if ( n_knots_trigonometric < 2 )
        raiseError(ERR_SPHERICAL_HARMONICS, 2, "SphericalHarmonics", "incorrect number of knots for TFs");
    */

    Legendre_grid_size_ = n_knots_Legendre;
    phi_grid_size_      = n_knots_trigonometric;

    allocateMemory();

    initFactorials();
    initLegendre();
    initTrigonometric();
}


SphericalHarmonics::~SphericalHarmonics()
{
    deallocateMemory();
}


cmpx SphericalHarmonics::calculateSphHarm(double cos_theta, double phi, int l, int m)
{
    cmpx Y_lm;
    double mult;
    double mphi;

    mult = sqrt( ( 2 * l + 1 ) * INV_OF_FOUR_PI * vec_fact_[l - m] / vec_fact_[l + m] ) * calculateLegendre(cos_theta, l, m);
    mphi = m * phi;
    Y_lm.re = mult * cos(mphi);
    Y_lm.im = mult * sin(mphi);

    return Y_lm;
}


void SphericalHarmonics::getSphHarmVecExact(double cos_theta, double phi, cmpx *vec_Ylm)
{
    double sphPlm;
    double mphi;
    double sin_mphi, cos_mphi;

// l = 4, m = 0
    vec_Ylm[0].re = K_4_0 * calculateLegendre(cos_theta, 4, 0);
    vec_Ylm[0].im = 0.0;
// l = 6, m = 0
    vec_Ylm[5].re = K_6_0 * calculateLegendre(cos_theta, 6, 0);
    vec_Ylm[5].im = 0.0;

// m = 1
    sin_mphi = sin(phi);
    cos_mphi = cos(phi);
// l = 4, m = 1
    sphPlm = K_4_1 * calculateLegendre(cos_theta, 4, 1);
    vec_Ylm[1].re = sphPlm * cos_mphi;
    vec_Ylm[1].im = sphPlm * sin_mphi;
// l = 6, m = 1
    sphPlm = K_6_1 * calculateLegendre(cos_theta, 6, 1);
    vec_Ylm[6].re = sphPlm * cos_mphi;
    vec_Ylm[6].im = sphPlm * sin_mphi;

// m = 2
    mphi = 2 * phi;
    sin_mphi = sin(mphi);
    cos_mphi = cos(mphi);
// l = 4, m = 2
    sphPlm = K_4_2 * calculateLegendre(cos_theta, 4, 2);
    vec_Ylm[2].re = sphPlm * cos_mphi;
    vec_Ylm[2].im = sphPlm * sin_mphi;
// l = 6, m = 2
    sphPlm = K_6_2 * calculateLegendre(cos_theta, 6, 2);
    vec_Ylm[7].re = sphPlm * cos_mphi;
    vec_Ylm[7].im = sphPlm * sin_mphi;

// m = 3
    mphi = 3 * phi;
    sin_mphi = sin(mphi);
    cos_mphi = cos(mphi);
// l = 4, m = 3
    sphPlm = K_4_3 * calculateLegendre(cos_theta, 4, 3);
    vec_Ylm[3].re = sphPlm * cos_mphi;
    vec_Ylm[3].im = sphPlm * sin_mphi;
// l = 6, m = 3
    sphPlm = K_6_3 * calculateLegendre(cos_theta, 6, 3);
    vec_Ylm[8].re = sphPlm * cos_mphi;
    vec_Ylm[8].im = sphPlm * sin_mphi;

// m = 4
    mphi = 4 * phi;
    sin_mphi = sin(mphi);
    cos_mphi = cos(mphi);
// l = 4, m = 4
    sphPlm = K_4_4 * calculateLegendre(cos_theta, 4, 4);
    vec_Ylm[4].re = sphPlm * cos_mphi;
    vec_Ylm[4].im = sphPlm * sin_mphi;
// l = 6, m = 4
    sphPlm = K_6_4 * calculateLegendre(cos_theta, 6, 4);
    vec_Ylm[9].re = sphPlm * cos_mphi;
    vec_Ylm[9].im = sphPlm * sin_mphi;

// m = 5
    mphi = 5 * phi;
// l = 6, m = 5
    sphPlm = K_6_5 * calculateLegendre(cos_theta, 6, 5);
    vec_Ylm[10].re = sphPlm * cos(mphi);
    vec_Ylm[10].im = sphPlm * sin(mphi);

// m = 6
    mphi = 6 * phi;
// l = 6, m = 6
    sphPlm = K_6_6 * calculateLegendre(cos_theta, 6, 6);
    vec_Ylm[11].re = sphPlm * cos(mphi);
    vec_Ylm[11].im = sphPlm * sin(mphi);
}


void SphericalHarmonics::getSphHarmVecFSI(double cos_theta, double phi, cmpx *vec_Ylm)
{
    int interval;
    double delta;

// pobieramy interpolowane wartosci stowarzyszonych wielomianow Legendre'a oraz
// funkcji trygonometrycznych sin i cos

// interpolujemy stowarzyszone wielomiany Legendre'a
    interval = int( floor ( (cos_theta + 1.0) * inv_Legendre_dx_ ) );
    delta = cos_theta - vec_Legendre_arg_[interval];

    ret_vec_[0] = vec_Legendre_[interval][0] + vec_Legendre_[interval][1] * delta;
    ret_vec_[1] = vec_Legendre_[interval][2] + vec_Legendre_[interval][3] * delta;

    ret_vec_[2] = vec_Legendre_[interval][4] + vec_Legendre_[interval][5] * delta;
    ret_vec_[3] = vec_Legendre_[interval][6] + vec_Legendre_[interval][7] * delta;

    ret_vec_[6] = vec_Legendre_[interval][8] + vec_Legendre_[interval][9] * delta;
    ret_vec_[7] = vec_Legendre_[interval][10] + vec_Legendre_[interval][11] * delta;

    ret_vec_[10] = vec_Legendre_[interval][12] + vec_Legendre_[interval][13] * delta;
    ret_vec_[11] = vec_Legendre_[interval][14] + vec_Legendre_[interval][15] * delta;

    ret_vec_[14] = vec_Legendre_[interval][16] + vec_Legendre_[interval][17] * delta;
    ret_vec_[15] = vec_Legendre_[interval][18] + vec_Legendre_[interval][19] * delta;

    ret_vec_[18] = vec_Legendre_[interval][20] + vec_Legendre_[interval][21] * delta;

    ret_vec_[21] = vec_Legendre_[interval][22] + vec_Legendre_[interval][23] * delta;

// interpolujemy funkcje trygonometryczne
    interval = int( floor ( phi * inv_dphi_ ) );
    delta = phi - vec_phi_[interval];

    ret_vec_[4] = vec_trig_mphi_[interval][0] + vec_trig_mphi_[interval][1] * delta;
    ret_vec_[5] = vec_trig_mphi_[interval][2] + vec_trig_mphi_[interval][3] * delta;

    ret_vec_[8] = vec_trig_mphi_[interval][4] + vec_trig_mphi_[interval][5] * delta;
    ret_vec_[9] = vec_trig_mphi_[interval][6] + vec_trig_mphi_[interval][7] * delta;

    ret_vec_[12] = vec_trig_mphi_[interval][8] + vec_trig_mphi_[interval][9] * delta;
    ret_vec_[13] = vec_trig_mphi_[interval][10] + vec_trig_mphi_[interval][11] * delta;

    ret_vec_[16] = vec_trig_mphi_[interval][12] + vec_trig_mphi_[interval][13] * delta;
    ret_vec_[17] = vec_trig_mphi_[interval][14] + vec_trig_mphi_[interval][15] * delta;

    ret_vec_[19] = vec_trig_mphi_[interval][16] + vec_trig_mphi_[interval][17] * delta;
    ret_vec_[20] = vec_trig_mphi_[interval][18] + vec_trig_mphi_[interval][19] * delta;

    ret_vec_[22] = vec_trig_mphi_[interval][20] + vec_trig_mphi_[interval][21] * delta;
    ret_vec_[23] = vec_trig_mphi_[interval][22] + vec_trig_mphi_[interval][23] * delta;

// obliczamy harmoniki sferyczne

// l = 4 m = 0 i l = 6 m = 0
    vec_Ylm[0].re = ret_vec_[0];
    vec_Ylm[0].im = 0.0;
    vec_Ylm[5].re = ret_vec_[1];
    vec_Ylm[5].im = 0.0;

// l = 4 m = 1 i l = 6 m = 1
    vec_Ylm[1].re = ret_vec_[2] * ret_vec_[4];
    vec_Ylm[6].re = ret_vec_[3] * ret_vec_[4];
    vec_Ylm[1].im = ret_vec_[2] * ret_vec_[5];
    vec_Ylm[6].im = ret_vec_[3] * ret_vec_[5];

// l = 4 m = 2 i l = 6 m = 2
    vec_Ylm[2].re = ret_vec_[6] * ret_vec_[8];
    vec_Ylm[7].re = ret_vec_[7] * ret_vec_[8];
    vec_Ylm[2].im = ret_vec_[6] * ret_vec_[9];
    vec_Ylm[7].im = ret_vec_[7] * ret_vec_[9];

// l = 4 m = 3 i l = 6 m = 3
    vec_Ylm[3].re = ret_vec_[10] * ret_vec_[12];
    vec_Ylm[8].re = ret_vec_[11] * ret_vec_[12];
    vec_Ylm[3].im = ret_vec_[10] * ret_vec_[13];
    vec_Ylm[8].im = ret_vec_[11] * ret_vec_[13];

// l = 4 m = 4 i l = 6 m = 4
    vec_Ylm[4].re = ret_vec_[14] * ret_vec_[16];
    vec_Ylm[9].re = ret_vec_[15] * ret_vec_[16];
    vec_Ylm[4].im = ret_vec_[14] * ret_vec_[17];
    vec_Ylm[9].im = ret_vec_[15] * ret_vec_[17];

// l = 6 m = 5
    vec_Ylm[10].re = ret_vec_[18] * ret_vec_[19];
    vec_Ylm[10].im = ret_vec_[18] * ret_vec_[20];

// l = 6 m = 6
    vec_Ylm[11].re = ret_vec_[21] * ret_vec_[22];
    vec_Ylm[11].im = ret_vec_[21] * ret_vec_[23];
}


#ifdef SPHERICAL_HARMONICS_HAVE_GSL
void SphericalHarmonics::getSphHarmVecGSL(double cos_theta, double phi, cmpx *vec_Ylm)
{
    double sphPlm;
    double mphi;
    double sin_mphi, cos_mphi;

// m = 0
// l = 4, m = 0
    sphPlm = gsl_sf_legendre_sphPlm(4, 0, cos_theta);
    vec_Ylm[0].re = sphPlm;
    vec_Ylm[0].im = 0.0;
// l = 6, m = 0
    sphPlm = gsl_sf_legendre_sphPlm(6, 0, cos_theta);
    vec_Ylm[5].re = sphPlm;
    vec_Ylm[5].im = 0.0;

// m = 1
    sin_mphi = sin(phi);
    cos_mphi = cos(phi);
// l = 4, m = 1
    sphPlm = gsl_sf_legendre_sphPlm(4, 1, cos_theta);
    vec_Ylm[1].re = sphPlm * cos_mphi;
    vec_Ylm[1].im = sphPlm * sin_mphi;
// l = 6, m = 1
    sphPlm = gsl_sf_legendre_sphPlm(6, 1, cos_theta);
    vec_Ylm[6].re = sphPlm * cos_mphi;
    vec_Ylm[6].im = sphPlm * sin_mphi;

// m = 2
    mphi = 2 * phi;
    sin_mphi = sin(mphi);
    cos_mphi = cos(mphi);
// l = 4, m = 2
    sphPlm = gsl_sf_legendre_sphPlm(4, 2, cos_theta);
    vec_Ylm[2].re = sphPlm * cos_mphi;
    vec_Ylm[2].im = sphPlm * sin_mphi;
// l = 6, m = 2
    sphPlm = gsl_sf_legendre_sphPlm(6, 2, cos_theta);
    vec_Ylm[7].re = sphPlm * cos_mphi;
    vec_Ylm[7].im = sphPlm * sin_mphi;

// m = 3
    mphi = 3 * phi;
    sin_mphi = sin(mphi);
    cos_mphi = cos(mphi);
// l = 4, m = 3
    sphPlm = gsl_sf_legendre_sphPlm(4, 3, cos_theta);
    vec_Ylm[3].re = sphPlm * cos_mphi;
    vec_Ylm[3].im = sphPlm * sin_mphi;
// l = 6, m = 3
    sphPlm = gsl_sf_legendre_sphPlm(6, 3, cos_theta);
    vec_Ylm[8].re = sphPlm * cos_mphi;
    vec_Ylm[8].im = sphPlm * sin_mphi;

// m = 4
    mphi = 4 * phi;
    sin_mphi = sin(mphi);
    cos_mphi = cos(mphi);
// l = 4, m = 4
    sphPlm = gsl_sf_legendre_sphPlm(4, 4, cos_theta);
    vec_Ylm[4].re = sphPlm * cos_mphi;
    vec_Ylm[4].im = sphPlm * sin_mphi;
// l = 6, m = 4
    sphPlm = gsl_sf_legendre_sphPlm(6, 4, cos_theta);
    vec_Ylm[9].re = sphPlm * cos_mphi;
    vec_Ylm[9].im = sphPlm * sin_mphi;

// m = 5
    mphi = 5 * phi;
// l = 6, m = 5
    sphPlm = gsl_sf_legendre_sphPlm(6, 5, cos_theta);
    vec_Ylm[10].re = sphPlm * cos(mphi);
    vec_Ylm[10].im = sphPlm * sin(mphi);

// m = 6
    mphi = 6 * phi;
// l = 6, m = 6
    sphPlm = gsl_sf_legendre_sphPlm(6, 6, cos_theta);
    vec_Ylm[11].re = sphPlm * cos(mphi);
    vec_Ylm[11].im = sphPlm * sin(mphi);
}
#endif
