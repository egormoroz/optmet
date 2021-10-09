#pragma once
#include <complex>
#include <array>

using Complex = std::complex<double>;
using Coeffs = std::array<Complex, 4>;

double norm_sqr(Complex z) {
    return z.real() * z.real() + z.imag() * z.imag();
}

class Polynomial {
    Coeffs m_coeffs;
public:
    Polynomial(const std::array<Complex, 3> &r) {
        m_coeffs[0] = 1;
        m_coeffs[1] = -r[0] - r[1] - r[2];
        m_coeffs[2] = r[0]*r[1] + r[0]*r[2] + r[1]*r[2];
        m_coeffs[3] = -r[0] * r[1] * r[2];
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
        const double h = 1e-7;
        return (eval(z + h) - eval(z - h)) / (2 * h);
    }

    double Polynomial::eval_der_y(Complex z) const {
        const double h = 1e-7;
        const Complex t(0, h);
        return (eval(z + t) - eval(z - t)) / (2*h);
    }

    Complex Polynomial::get_last_root() const {
        return -m_coeffs[3] / m_coeffs[2];
    }

    Complex Polynomial::divide(Complex z0) {
        Coeffs q;
        Complex k = 0;
        for (int i = 0; i <= 2; ++i) {
            k = m_coeffs[i] + k * z0;
            q[i + 1] = k;
        }

        Complex remainder = m_coeffs[3] + k * z0;
        m_coeffs = q;
        return remainder;
    }
};
