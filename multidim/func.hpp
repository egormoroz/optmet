#pragma once
#include "poly.hpp"

class Function {
public:
    Function() : m_poly{{ 1, {-17, -18}, {-13, 214}, {485, -388} } } 
    {
    }

    double eval(Complex z) {
        m_counter++;
        return m_poly.eval(z);
    }

    //not really a comlex number; just a nice way to represent a vector
    Complex eval_grad(Complex z) {
        m_counter += 2;
        return {m_poly.eval_der_x(z), m_poly.eval_der_y(z)};
    }

    void reset() {
        m_counter = 0;
    }

    Complex exclude_root(Complex z) {
        return m_poly.divide(z);
    }

    Complex get_last_root() const {
        return m_poly.get_last_root();
    }

    int num_calls() const { return m_counter; }

    double dbg_eval(Complex z) const { return m_poly.eval(z); }
    Complex dbg_eval_grad(Complex z) const {
        return {m_poly.eval_der_x(z), m_poly.eval_der_y(z)};
    }

private:
    Polynomial m_poly;
    int m_counter = 0;
};
