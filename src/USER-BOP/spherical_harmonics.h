// written by Szymon Winczewski

#ifndef spherical_harmonics_h
#define spherical_harmonics_h

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "complex_numbers.h"
#include "constants.h"

#ifdef SPHERICAL_HARMONICS_HAVE_GSL
#include <gsl/gsl_sf_legendre.h>
#endif


// UWAGA:
// moze byc wymagana biblioteka GSL (Gnu Scientific Library)
// konfiguracja Code::Blocks, w "Project"->"Build options"
// 1) w "Linker settings" w "Other linker options" dodać "-lgsl" oraz "-lgslcblas"
// w przypadku gdy biblioteka jest niedostepna nalezy w defines.h
// zakomentowac linie:
// #define SPHERICAL_HARMONICS_HAVE_GSL

// klasa sluzy do szybkiego obliczania harmonik sferycznych,
// dla zadanych wartosci katow cos_theta oraz phi zwraca
// wektor harmonik:
//   0   Y_40(cos_theta, phi)
//   1   Y_41(cos_theta, phi)
//   2   Y_42(cos_theta, phi)
//   3   Y_43(cos_theta, phi)
//   4   Y_44(cos_theta, phi)
//   5   Y_60(cos_theta, phi)
//   6   Y_61(cos_theta, phi)
//   7   Y_62(cos_theta, phi)
//   8   Y_63(cos_theta, phi)
//   9   Y_64(cos_theta, phi)
//  10   Y_65(cos_theta, phi)
//  11   Y_66(cos_theta, phi)

// wspolczynniki normalizacyjne harmonik sferycznych
// (zakodowane na twardo i wlaczone w stowarzyszone wielomiany Legendre'a)
// postac wspolczynnika:
//   K_l_m = sqrt( ( 2 * l + 1 ) / ( 4 * PI ) * (l - m)! / (l + m)! );

#define K_4_0 0.846284375321634474431675698724575340747833251953125
#define K_4_1 0.18923493915151201605340247624553740024566650390625
#define K_4_2 0.0446031029038192750046931678298278711736202239990234375
#define K_4_3 0.01192068067522240176758785423771769274026155471801757812
#define K_4_4 0.004214597070904596912144235432151617715135216712951660156
#define K_6_0 1.0171072362820547940742699211114086210727691650390625
#define K_6_1 0.1569430538290060017647675749685731716454029083251953125
#define K_6_2 0.02481487565210345469512986937843379564583301544189453125
#define K_6_3 0.004135812608683909694096136888674664078280329704284667969
#define K_6_4 0.0007550926197968212474909144305001973407343029975891113281
#define K_6_5 0.0001609862874555168616927391944670944212703034281730651855
#define K_6_6 4.647273819914056675658281525542747658619191497564315796e-05

// tablice interpolacyjne stowarzyszonych wielomianow Legendre'a
// vec_Legendre_[interval_x][]
//  0 - P40    1 - dP40
//  2 - P60    3 - dP60
//  4 - P41    5 - dP41
//  6 - P61    7 - dP61
//  8 - P42    9 - dP42
// 10 - P62   11 - dP62
// 12 - P43   13 - dP43
// 14 - P63   15 - dP63
// 16 - P44   17 - dP44
// 18 - P64   19 - dP64
// 20 - P65   21 - dP65
// 22 - P66   23 - dP66

// tablice interpolacyjne funkcji trygonometrycznych
// vec_trig_mphi_[interval_phi][]
//  0 - cos(phi)      1 - dcos(phi)
//  2 - sin(phi)      3 - dsin(phi)
//  4 - cos(2 phi)    5 - dcos(2 phi)
//  6 - sin(2 phi)    7 - dsin(2 phi)
//  8 - cos(3 phi)    9 - dcos(3 phi)
// 10 - sin(3 phi)   11 - dsin(3 phi)
// 12 - cos(4 phi)   13 - dcos(4 phi)
// 14 - sin(4 phi)   15 - dsin(4 phi)
// 16 - cos(5 phi)   17 - dcos(5 phi)
// 18 - sin(5 phi)   19 - dsin(5 phi)
// 20 - cos(6 phi)   21 - dcos(6 phi)
// 22 - sin(6 phi)   23 - dsin(6 phi)

// ret_vec_[]
// 0  - P40(x)    1 - P60(x)
// 2  - P41(x)    3 - P61(x)    4 - cos(1 phi)    5 - sin(1 phi)
// 6  - P42(x)    7 - P62(x)    8 - cos(2 phi)    9 - sin(2 phi)
// 10 - P43(x)   11 - P63(x)   12 - cos(3 phi)   13 - sin(3 phi)
// 14 - P44(x)   15 - P64(x)   16 - cos(4 phi)   17 - sin(4 phi)
//               18 - P65(x)   19 - cos(5 phi)   20 - sin(5 phi)
//               21 - P66(x)   22 - cos(6 phi)   23 - sin(6 phi)


// lookupy silni wypelniane do ponizszej wartosci
const int SPHERICAL_HARMONICS_MAX_FACTORIAL = 20;

// lookupy podwojnej silni wypelniane do ponizszej wartosci
const int SPHERICAL_HARMONICS_MAX_DOUBLE_FACTORIAL = 40;

namespace LAMMPS_NS
{
  class Memory;
}

class SphericalHarmonics
{
protected:
    double *vec_fact_;
    double *vec_dfact_;

    int Legendre_grid_size_;
    double Legendre_dx_;
    double inv_Legendre_dx_;
    double *vec_Legendre_arg_;
    double **vec_Legendre_;

    int phi_grid_size_;
    double dphi_;
    double inv_dphi_;
    double *vec_phi_;
    double **vec_trig_mphi_;

    double *ret_vec_;

    LAMMPS_NS::Memory* memory;

    double calculateFactorial(int i);
    double calculateDoubleFactorial(int i);
    double calculateLegendre(double x, int l, int m);

    void allocateMemory();
    void deallocateMemory();

    void initFactorials();
    void initLegendre();
    void initTrigonometric();


public:
    SphericalHarmonics(int n_knots_Legendre, int n_knots_trigonometric, LAMMPS_NS::Memory* memory);
    ~SphericalHarmonics();

    cmpx calculateSphHarm(double cos_theta, double phi, int l, int m);

    void getSphHarmVecExact(double cos_theta, double phi, cmpx *vec_Ylm);
    void getSphHarmVecFSI(double cos_theta, double phi, cmpx *vec_Ylm);

    #ifdef SPHERICAL_HARMONICS_HAVE_GSL
    void getSphHarmVecGSL(double cos_theta, double phi, cmpx *vec_Ylm);
    #endif
};

#endif
