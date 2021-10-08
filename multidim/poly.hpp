#pragma once
#include <complex>
#include <array>

using Complex = std::complex<double>;

using CCoeffs = std::array<Complex, 4>;
using RCoeffs = std::array<double, 10>;

double norm_sqr(Complex z);

class Polynomial {
    CCoeffs m_coeffs;
    RCoeffs m_real, m_imag;
    RCoeffs m_rdx, m_rdy;
    RCoeffs m_idx, m_idy;

    void recalc(); 
public:
    Polynomial(const CCoeffs &coeffs);

    double eval(Complex z) const;

    double eval_der_x(Complex z) const;
    double eval_der_y(Complex z) const;

    Complex divide(Complex z0);
};
