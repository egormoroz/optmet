#pragma once

#include <array>
#include <complex>
using Complex = std::complex<double>;

double norm_sqr(Complex z) {
    return z.real() * z.real() + z.imag() * z.imag();
}

struct Polynomial {
    Complex coeffs[4];

    double eval(Complex z) const;
    double eval_der_x(Complex z) const;
    double eval_der_y(Complex z) const;

    //divide by z-z0
    std::pair<Polynomial, Complex> divide(Complex z0) const;
};

double Polynomial::eval(Complex z) const {
    Complex result = 0;
    for (int i = 0; i <= 3; ++i) {
        result += coeffs[i] * pow(z, 3 - i);
    }

    return norm_sqr(result);
}

double Polynomial::eval_der_x(Complex z) const {
    const double h = 1e-15;
    const Complex t(h, 0);
    return (eval(z + h) - eval(z - h)) / (2*h);
}

double Polynomial::eval_der_y(Complex z) const {
    const double h = 1e-15;
    const Complex t(0, h);
    return (eval(z + t) - eval(z - t)) / (2*h);
}

std::pair<Polynomial, Complex> Polynomial::divide(Complex z0) const {
    Polynomial result;
    Complex k = 0;
    for (int i = 0; i <= 2; ++i) {
        k = coeffs[i] + k * z0;
        result.coeffs[i + 1] = k;
    }
    Complex remainder = coeffs[3] + k * z0;

    return {result, remainder};
}
