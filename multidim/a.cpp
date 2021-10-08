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

    Complex exclude_root(Complex z) {
        auto result = m_poly.divide(z);
        m_poly = result.first;
        return result.second;
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

int main() {
    Function f;
    for (int i = 0; i < 3; ++i) {
        auto z0 = coord_descent(f, {1.0, 1.0}, {0.33, 0.33}, 1e-4);
        auto r = f.exclude_root(z0);
        printf("f(%.5f + %.5fi) = %12.5e + %12.5ei; ", 
            z0.real(), z0.imag(), r.real(), r.imag());
        printf("%d calls\n", f.num_calls());
    }
}
