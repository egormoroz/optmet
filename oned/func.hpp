#pragma once
#include <cmath>

const double A = 4;
const double B = 5;
const double C = 9;

class Function {
    int m_counter = 0;
public:
    Function() = default;

    double eval(double x) {
        m_counter++;
        return exp(B * x) + exp(-A * x) + A * x * x + C * x + B - x * x * x;
    }

    double eval_der(double x) {
        m_counter += 1;
        return B * exp(B * x) - A * exp(-A * x) + 2 * A * x + C - 3 * x * x;
    }

    double eval_der2(double x) {
        m_counter += 1;
        return B * B * exp(B * x) + A * A * exp(-A * x) + 2 * A - 6 * x;
    }

    int num_of_calls() const { return m_counter; }
};

