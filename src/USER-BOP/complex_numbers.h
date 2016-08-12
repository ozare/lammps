// written by Szymon Winczewski

#ifndef complex_numbers_h
#define complex_numbers_h

#include <cmath>


typedef struct
{
    double re, im;
} cmpx;
// liczba zespolona


inline cmpx CmpxSet(double re, double im);
inline cmpx CmpxAdd(cmpx op1, cmpx op2);
inline cmpx CmpxSub(cmpx op1, cmpx op2);
inline cmpx CmpxMult(cmpx op1, cmpx op2);
inline cmpx CmpxMult(double op1, cmpx op2);
inline double CmpxNorm(cmpx op);
inline double CmpxAbs(cmpx op);
inline cmpx CmpxConjugate(cmpx op);


inline cmpx CmpxSet(double re, double im)
{
  cmpx res;
  res.re = re;
  res.im = im;
  return res;
}


inline cmpx CmpxAdd(cmpx op1, cmpx op2)
{
  cmpx res;
  res.re = op1.re + op2.re;
  res.im = op1.im + op2.im;
  return res;
}


inline cmpx CmpxSub(cmpx op1, cmpx op2)
{
  cmpx res;
  res.re = op1.re - op2.re;
  res.im = op1.im - op2.im;
  return res;
}


inline cmpx CmpxMult(cmpx op1, cmpx op2)
{
  cmpx res;
  res.re = op1.re * op2.re - op1.im * op2.im;
  res.im = op1.re * op2.im + op2.re * op1.im;
  return res;
}


inline cmpx CmpxMult(double op1, cmpx op2)
{
  cmpx res;
  res.re = op1 * op2.re;
  res.im = op1 * op2.im;
  return res;
}


inline double CmpxNorm(cmpx op)
{
  return ( op.re * op.re + op.im * op.im );
}


inline double CmpxAbs(cmpx op)
{
  return sqrt( op.re * op.re + op.im * op.im );
}


inline cmpx CmpxConjugate(cmpx op)
{
  cmpx res;
  res.re =   op.re;
  res.im = - op.im;
  return res;
}


inline double Cmpx_Re_Z1_Z2_Z3Star(cmpx z1, cmpx z2, cmpx z3)
{
    return z1.re * z2.re * z3.re + z1.re * z2.im * z3.im + z2.re * z1.im * z3.im - z3.re * z1.im * z2.im;
}


inline double Cmpx_Re_Z1_Z1_Z2Star(cmpx z1, cmpx z2)
{
    return z1.re * z1.re * z2.re + 2.0 * z1.re * z1.im * z2.im - z2.re * z1.im * z1.im;
}

#endif
