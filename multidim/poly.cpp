#include "poly.hpp"
#include <cstdio>

double norm_sqr(Complex z) {
    return z.real() * z.real() + z.imag() * z.imag();
}

RCoeffs real(const CCoeffs &a) {
    return {
        a[0].real(), //x^3
        -3 * a[0].imag(), //x^2y
        -3 * a[0].real(), //xy^2
        a[0].imag(), //y^3
        a[1].real(), //x^2
        -2 * a[1].imag(), //xy
        -a[1].real(), //y^2
        a[2].real(), //x
        -a[2].imag(), //y
        a[3].real() //1
    };
}

RCoeffs imag(const CCoeffs &a) {
    return {
        a[0].imag(), //x^3
        3 * a[0].real(), //x^2y
        -3 * a[0].imag(), //xy^2
        -a[0].real(), //y^3
        a[1].imag(), //x^2
        2 * a[1].real(), //xy
        -a[1].imag(), //y^2
        a[2].imag(), //x
        a[2].real(), //y
        a[3].imag() //1
    };
}

RCoeffs der_x(const RCoeffs &a) {
    return {
        0, //x^3
        0, //x^2y
        0, //xy^2
        0, //y^3
        3 * a[0], //x^2
        2 * a[1], //xy
        a[2], //y^2
        2 * a[4], //x
        a[5], //y
        a[7], //1
    };
}

RCoeffs der_y(const RCoeffs &a) {
    return {
        0, //x^3
        0, //x^2y
        0, //xy^2
        0, //y^3
        a[1], //x^2
        2 * a[2], //xy
        3 * a[3], //y^2
        a[5], //x
        2 * a[6], //y
        a[8], //1
    };
}

double polyval2d(const RCoeffs &a, Complex z) {
    double x = z.real(), y = z.imag();
    return a[0] * x * x * x
         + a[1] * x * x * y
         + a[2] * x * y * y
         + a[3] * y * y * y
         + a[4] * x * x
         + a[5] * x * y
         + a[6] * y * y
         + a[7] * x
         + a[8] * y
         + a[9];
}

Polynomial::Polynomial(const CCoeffs &coeffs) 
    : m_coeffs(coeffs) 
{
    recalc();
}

void Polynomial::recalc() {
    m_real = real(m_coeffs); m_imag = imag(m_coeffs);
    m_rdx = der_x(m_real); m_rdy = der_y(m_real);
    m_idx = der_x(m_imag); m_idy = der_y(m_imag);
}


double Polynomial::eval(Complex z) const {
    return norm_sqr(
        m_coeffs[0] * z * z * z
      + m_coeffs[1] * z * z
      + m_coeffs[2] * z
      + m_coeffs[3]
    );
}

double Polynomial::eval_der_x(Complex z) const {
    return 2 * (polyval2d(m_real, z) * polyval2d(m_rdx, z)
        + polyval2d(m_imag, z) * polyval2d(m_idx, z));
}

double Polynomial::eval_der_y(Complex z) const {
    return 2 * (polyval2d(m_real, z) * polyval2d(m_rdy, z)
        + polyval2d(m_imag, z) * polyval2d(m_idy, z));
}

Complex Polynomial::divide(Complex z0) {
    CCoeffs q;
    Complex k = 0;
    for (int i = 0; i <= 2; ++i) {
        k = m_coeffs[i] + k * z0;
        q[i + 1] = k;
    }

    Complex remainder = m_coeffs[3] + k * z0;
    m_coeffs = q;
    recalc();
    return remainder;
}

