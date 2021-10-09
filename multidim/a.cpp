#include <cstdio>
#include <iostream>
#include "poly.hpp"
#include "oned.hpp"

#define uwu(s) #s, s

using namespace std;

double muls(Complex x, Complex y) {
    return x.real() * y.real() + x.imag() * y.imag();
}

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

void log_iter(int k, Complex z0, double r, int calls) {
 //   printf("%4d) f(%.5f + %.5fi) = %12.5e;\t%d calls\n",
 //       k, z0.real(), z0.imag(), r, calls);
}

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

    int k = 0;

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

        log_iter(k++, p, z, f.num_calls());

        if (old_dx < eps && old_dy < eps)
            break;
    } 

    return p;
}

Complex grad_descent_div(Function &f, Complex z0, double alpha, double delta) {
    const double EPS = 0.1;
    const double L = 0.5;

    int k = 0;

    Complex grad = f.eval_grad(z0);
    double val = f.eval(z0);
    while (norm_sqr(grad) >= delta) {
        Complex z = z0 - alpha*grad;
        double new_val = f.eval(z);
        if (new_val - val > -alpha * EPS * norm_sqr(grad)) {
            alpha *= L;
            continue;
        }

        z0 = z;
        grad = f.eval_grad(z);
        val = new_val;

        log_iter(k++, z0, val, f.num_calls());
    }

    return z0;
}

Complex grad_descent_const(Function &f, Complex z, double alpha, double delta) {
    Complex grad = f.eval_grad(z);
    double val = f.eval(z);
    int counter = 0;
    int k = 0;
    while (norm_sqr(grad) >= delta && counter < 10) {
        z = z - alpha*grad;
        auto new_val = f.eval(z);
        if (new_val > val) {
            counter++;
        }
        val = new_val;
        grad = f.eval_grad(z);
        log_iter(k++, z, val, f.num_calls());
    }

    return z;
}

Complex grad_descent_seq(Function &f, Complex z, int k, double delta) {
    Complex grad = f.eval_grad(z);
    double val = f.eval(z);
    int j = 0;
    while (norm_sqr(grad) >= delta) {
        z = z - (1.0/k++) * grad;
        val = f.eval(z);
        grad = f.eval_grad(z);
        log_iter(j++, z, val, f.num_calls());
    }

    return z;
}

Complex grad_descent_fastest(Function &f, Complex z, double delta) {
    const double EPS = 1e-5;
    Complex grad = f.eval_grad(z);
    double val = f.eval(z);
    int k = 0;
    while (norm_sqr(grad) >= delta) {
        double alpha = golden_ratio([&](double alpha) {
            return f.eval(z - alpha * grad);
        }, 0, 1, EPS);
        z = z - alpha * grad;
        val = f.eval(z);
        grad = f.eval_grad(z);
        log_iter(k++, z, val, f.num_calls());
    }

    return z;
}

template<typename Method, typename... Args>
void test_method(const char *name, Method m, Complex z0, Args&& ...args) {
    Function f;
    printf("%s:\n", name);
    for (int i = 0; i < 3; ++i) {
        Complex r;
        if (i < 2) {
          z0 = m(f, z0, args...);
          r = f.exclude_root(z0);
        } else {
            z0 = f.get_last_root();
            r = 0;
        }
        printf("f(%.5f + %.5fi) = %12.5e + %12.5ei; ", 
            z0.real(), z0.imag(), r.real(), r.imag());
        printf("%d calls\n", f.num_calls());
    }
    printf("\n");
}


int main() {
    using namespace std::complex_literals;
    auto start = 10.0 + 10i;
    test_method(uwu(coord_descent), 0, start, 1e-3);
    test_method(uwu(grad_descent_div), start, 1, 1e-3);
    test_method(uwu(grad_descent_const), start, 1e-4, 1e-3);
    test_method(uwu(grad_descent_seq), start, 10000, 1e-3);
    test_method(uwu(grad_descent_fastest), start, 1e-3);
}

