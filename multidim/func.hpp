#pragma once

#include <complex>
#include <cmath>

using Complex = std::complex<double>;


template<size_t N>
Complex polyval(const Complex (&p)[N], Complex z) {
    Complex val = 0;
    for (int i = 0; i < N; ++i) {
        val = val * z + p[i];
    }
    return val;
}

//Скалярное произведение (в типе Complex хранятся векторы!)
double muls(Complex a, Complex b) {
    return a.real() * b.real() + a.imag() * b.imag();
}

//Квадрат нормы а
double norm_sqr(Complex a) {
    return muls(a, a);
}

//Норма а
double norm(Complex a) {
    return sqrt(norm_sqr(a));
}

/*
* Описание целевой функции с счётчиком вызовов:
* Квадрат модуля полинома с корнями z1, z2, z3
* */
class Function {
    Complex m_poly[4], m_der[3];
    int m_counter;

    void recalc_der() {
        for (int i = 0; i < 3; ++i)
            m_der[i] = double(3 - i) * m_poly[i];
    }

public:
    Function()
        : Function({4, 5}, {4, 9}, {9, 4}) {}

    Function(Complex z1, Complex z2, Complex z3) {
        m_counter = 0;
        m_poly[0] = 1.0;
        m_poly[1] = -(z1 + z2 + z3);
        m_poly[2] = z1 * z2 + z1 * z3 + z2 * z3;
        m_poly[3] = -z1 * z2 * z3;
        recalc_der();
    }

    //Вычислить значение в точке z
    double eval(Complex z) {
        m_counter++;
        return dbg_eval(z);
    }

    //Вычислить градиент в точке z и вернуь его в виде (df/dx, i*df/dy)
    Complex eval_grad(Complex z) {
        m_counter += 2;
        return dbg_eval_grad(z);
    }

    //Вычислить значение функции, не увеличивая счётчик (используется для отладки)
    double dbg_eval(Complex z) const {
        return norm_sqr(polyval(m_poly, z));
    }

    //Вычислить градиент функции, не увеличивая счётчик
    Complex dbg_eval_grad(Complex z) const {
        return 2.0 * std::conj(polyval(m_der, z)) * polyval(m_poly, z);
    }

    //Сбросить кол-во вызовов
    void reset() { m_counter = 0; }
    //Получить количество вызовов
    int num_calls() const { return m_counter; }

    //Разделить полином на корень z
    //(делит изначальный комплексный полином третьей степени, а не квадрат модуля)
    Complex divide(Complex z) {
        Complex k = 0;
        Complex new_coeffs[4];
        new_coeffs[0] = 0;
        for (int i = 0; i < 3; ++i) {
            k = k * z + m_poly[i];
            new_coeffs[i + 1] = k;
        }
        Complex remainder = k * z + m_poly[3];

        for (int i = 0; i < 4; ++i)
            m_poly[i] = new_coeffs[i];
        recalc_der();
        return remainder;
    }

    //Получить последний корень, когда полином выглядит так: az+b=0
    Complex get_last_root() const {
        return -m_poly[3] / m_poly[2];
    }
};

