#include <cstdio>
#include <iostream>
#include "poly.hpp"

using namespace std;


class Function {
public:
    Function() : m_poly{ { 1, {-17, -18}, {-13, 214}, {485, -388} } } 
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

    Complex exclude_root(Complex z) {
        return m_poly.divide(z);
    }

    int num_calls() const { return m_counter; }

private:
    Polynomial m_poly;
    int m_counter = 0;
};

bool try_advance(Function &f, Complex &p, double &z, Complex step) {
    auto pp = p + step;
    double zz = f.eval(pp);
    if (zz < z) {
        do {
            p = pp; z = zz;
            pp = p + step; zz = f.eval(pp);
        } while (zz < z);
        return true;
    }
    return false;
}

Complex coord_descent(Function &f, Complex p, Complex step, double eps) {
    Complex dx(step.real(), 0.), dy(0., step.imag()), pp;
    double z = f.eval(p), zz;

    while (true) {
        pp = p; zz = z;
        int flag = 0;
        double old_dx = dx.real(), old_dy = dy.imag();
        if (!try_advance(f, p, z, dx) && !try_advance(f, p, z, -dx)) {
            dx /= 2;
            flag++;
        }
        dx /= 1.1;
        if (!try_advance(f, p, z, dy) && !try_advance(f, p, z, -dy)) {
            dy /= 2;
            flag++;
        }
        dy /= 1.1;

        if (old_dx < eps && old_dy < eps)
            break;
    } 

    return p;
}

Complex grad_descent(Function &f, Complex z0, double alpha, double delta) {
    const double EPS = 0.1;
    const double L = 0.5;

    const bool OVERJUMP_GUARD = false;

    Complex grad = f.eval_grad(z0);
    double val = f.eval(z0);
    while (norm_sqr(grad) >= delta) {
        Complex z = z0 - alpha*grad;
        double new_val = f.eval(z);
        if (new_val - val > -alpha * EPS * norm_sqr(grad)) {
            alpha *= L;
            continue;
        }

        Complex new_grad = f.eval_grad(z);
        if (OVERJUMP_GUARD) {
            //если градиент ВНЕЗАПНО развернулся
            if (grad.real() * new_grad.real() + grad.imag() * new_grad.imag() < 0) {
                alpha *= L;
                continue;
            }
        }
        z0 = z;
        grad = new_grad;
        val = new_val;
    }

    return z0;
}

void test_coord_descent() {
    Function f;
    for (int i = 0; i < 3; ++i) {
        auto z0 = coord_descent(f, {0, 0}, {1, 1}, 1e-4);
        auto r = f.exclude_root(z0);
        printf("f(%.5f + %.5fi) = %12.5e + %12.5ei; ", 
            z0.real(), z0.imag(), r.real(), r.imag());
        printf("%d calls\n", f.num_calls());
    }
}

void test_grad_descent() {
    Function f;
    for (int i = 0; i < 3; ++i) {
        auto z0 = grad_descent(f, {0, 0}, 1, 1e-4);
        auto r = f.exclude_root(z0);
        printf("f(%.5f + %.5fi) = %12.5e + %12.5ei; ", 
            z0.real(), z0.imag(), r.real(), r.imag());
        printf("%d calls\n", f.num_calls());
    }
}

int main() {
    test_coord_descent();
    test_grad_descent();
}

