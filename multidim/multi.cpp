#include <cstdio>
#include <array>
#include "func.hpp"
#include "oned.hpp"

using namespace std;
using namespace std::complex_literals;

struct Result {
    Complex z;
    int calls, iterations;
};

void log_gd(int k, int calls, Complex z, Complex ag, double val) {
    printf("%d, %d, (%f; %f), (%f; %f), %f\n", ++k, calls,
        z.real(), z.imag(), ag.real(), ag.imag(), val);
}

void print_gd_results(const Result r[3]) {
    static Function f;
    for (int j = 0; j < 3; ++j) {
        auto& i = r[j];
        printf("%d, %d, (%.3f; %.3f), 0, %f\n", i.iterations, i.calls, 
            i.z.real(), i.z.imag(), f.eval(i.z));
    }
}


Complex estm_min_loc(double x0, double xx, double y0, double yy, double step) {
    Function f;
    Complex z(x0, y0);
    double min_val = f.eval(z);

    for (double x = x0; x < xx; x += step) {
        for (double y = y0; y < yy; y += step) {
            double val = f.eval({x, y});
            if (val < min_val) {
                z = {x, y};
                min_val = val;
            }
        }
    }

    return z;
}

bool try_advance(Function &f, Complex &z, double &val, Complex step, int k) {
    auto zz = z + step;
    double new_val = f.eval(zz);
    printf("%d, %d, %f, %f, (%f; %f), %f\n", k, f.num_calls(), 
        step.real(), step.imag(), zz.real(), zz.imag(), new_val);
    if (new_val >= val)
        return false;
    do {
        z = zz; val = new_val;
        zz = z + step; new_val = f.eval(zz);
        printf("%d, %d, %f, %f, (%f; %f), %f\n", k, f.num_calls(), 
            step.real(), step.imag(), zz.real(), zz.imag(), new_val);
    } while (new_val < val);
    return true;
}

int coord_descent(Function &f, Complex &z, Complex &step, int r) {
    Complex dx(step.real(), 0.), dy(0., step.imag()), pp;
    double val = f.eval(z);

    double K = 1.0 / (2.0 * pow(10, r)),
        old_dx, old_dy;

    int k = 1;

    do {
        old_dx = abs(dx.real()); old_dy = abs(dy.imag());

        if (old_dx > K * abs(z.real())) {
            if (!try_advance(f, z, val, dx, k)) {
                dx = -dx;
                try_advance(f, z, val, dx, k);
            }
            dx /= 2;
        }

        if (old_dy > K * abs(z.imag())) {
            if (!try_advance(f, z, val, dy, k)) {
                dy = -dy;
                try_advance(f, z, val, dy, k);
            }
            dy /= 2;
        }

        ++k;
    } while (old_dx > K * abs(z.real()) || old_dy > K * abs(z.imag()));

    step = dx + dy;

    return k;
}

int grad_descent_div(Function &f, Complex &z0, double alpha, double delta) {
    const double EPS = 0.1;
    const double L = 0.5;

    int k = 0;

    Complex grad = f.eval_grad(z0);
    double val = f.eval(z0);
    while (norm_sqr(grad) >= delta) {
        Complex z = z0 - alpha*grad;
        double new_val = f.eval(z);
        if (new_val - val > -alpha * EPS * norm_sqr(grad)) {
            printf("*");
            log_gd(k, f.num_calls(), z, alpha * grad, new_val);
            alpha *= L;
            continue;
        }
        log_gd(++k, f.num_calls(), z, alpha * grad, new_val);

        z0 = z;
        grad = f.eval_grad(z);
        val = new_val;
    }

    return k;
}

template<bool CAREFUL = false>
int grad_descent_const(Function &f, Complex &z, double alpha, double delta) {
    Complex grad = f.eval_grad(z);
    double val;
    if (CAREFUL) val = f.eval(z);
    int k = 0;
    while (norm_sqr(grad) >= delta) {
        auto ag = alpha * grad;
        log_gd(++k, f.num_calls(), z, ag, f.dbg_eval(z - ag));

        if (CAREFUL) {
            double new_val = f.eval(z - alpha*grad);
            if (new_val > val)
                break;
            val = new_val;
        }
        z = z - alpha*grad;
        grad = f.eval_grad(z);
    }

    return k;
}

int grad_descent_seq(Function &f, Complex &z, int &n, double delta) {
    Complex grad = f.eval_grad(z);
    double val = f.eval(z);
    int k = 0;
    while (norm_sqr(grad) >= delta) {
        grad /= norm(grad);
        double alpha = 1.0/n++;
        auto ag = alpha * grad;
        log_gd(++k, f.num_calls(), z, ag, f.dbg_eval(z - ag));

        z = z - alpha * grad;
        val = f.eval(z);
        grad = f.eval_grad(z);
    }

    return k;
}

int grad_descent_fastest(Function &f, Complex &z, double delta) {
    const double EPS = 1e-5;
    Complex grad = f.eval_grad(z);
    int k = 1;
    while (norm_sqr(grad) >= delta) {
        double alpha = xyz([&](double alpha) {
            double val = f.eval(z - alpha * grad);
            log_gd(k, f.num_calls(), z, alpha * grad, val);
            return val;
        }, 0.1, EPS);

        z = z - alpha * grad;
        grad = f.eval_grad(z);
        ++k;
    }

    return k;
}


void test_coord_descent(Complex z0, Complex stp, int R) {
    Result r[3];
    freopen("coord_descent.csv", "w", stdout);
    printf("Iteration, Calls, StepX, StepY, Z, f(z)\n");

    Complex step;
    int its, lits, sits;
    Function orig, f;

    //first root
    step = stp;
    its = coord_descent(f, z0, step, R);
    r[0] = { z0, f.num_calls(), its };
    f.divide(z0);

    //second root
    f.reset(); orig.reset();
    step = stp;
    lits = coord_descent(f, z0, step, R - 1);
    sits = coord_descent(orig, z0, step, R);
    r[1] = { z0, f.num_calls() + orig.num_calls(), lits + sits };
    f.divide(z0);

    //third root
    orig.reset();
    double K = 1.0 / (2.0 * pow(10, R));
    z0 = f.get_last_root();
    step = Complex(abs(z0.real()), abs(z0.imag())) * K;
    its = coord_descent(orig, z0, step, R);
    r[2] = { z0, orig.num_calls(), its };

    for (auto& i : r) {
        printf("%d, %d, 0, 0, (%.3f; %.3f), %f\n", i.iterations, i.calls, 
            i.z.real(), i.z.imag(), orig.dbg_eval(i.z));
    }
}

void test_gd_div(Complex z0, double alpha, double delta) {
    freopen("gd_div.csv", "w", stdout);
    printf("Iteration, Calls, Z, alpha * Gradient, f(z - alpha * gradient)\n");

    int its, lits, sits;
    Function orig, f;
    Result r[3];

    //first root
    its = grad_descent_div(f, z0, alpha, delta);
    r[0] = { z0, f.num_calls(), its };
    f.divide(z0);

    //second root
    f.reset(); orig.reset();
    lits = grad_descent_div(f, z0, alpha, delta);
    sits = grad_descent_div(orig, z0, alpha, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.divide(z0);

    //third root
    orig.reset();
    z0 = f.get_last_root();
    its = grad_descent_div(orig, z0, alpha, delta);
    r[2] = {z0, orig.num_calls(), its};

    print_gd_results(r);
}

void test_gd_const(Complex z0, double alpha, double delta) {
    freopen("gd_const.csv", "w", stdout);
    printf("Iteration, Calls, Z, alpha * Gradient, f(z - alpha * gradient)\n");

    int its, lits, sits;
    Function orig, f;
    Result r[3];

    //first root
    its = grad_descent_const(f, z0, alpha, delta);
    r[0] = { z0, f.num_calls(), its };
    f.divide(z0);

    //second root
    f.reset(); orig.reset();
    lits = grad_descent_const(f, z0, alpha, delta);
    sits = grad_descent_const<true>(orig, z0, alpha, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.divide(z0);

    //third root
    orig.reset();
    z0 = f.get_last_root();
    its = grad_descent_const<true>(orig, z0, alpha, delta);
    r[2] = {z0, orig.num_calls(), its};

    print_gd_results(r);
}

void test_gd_seq(Complex z0, int n, double delta) {
    freopen("gd_seq.csv", "w", stdout);
    printf("Iteration, Calls, Z, alpha * Gradient, f(z - alpha * gradient)\n");

    int its, lits, sits;
    Function orig, f;
    Result r[3];

    //first root
    int m = n;
    its = grad_descent_seq(f, z0, m, delta);
    r[0] = { z0, f.num_calls(), its };
    f.divide(z0);

    //second root
    f.reset(); orig.reset();
    m = n;
    lits = grad_descent_seq(f, z0, m, delta);
    sits = grad_descent_seq(orig, z0, m, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.divide(z0);

    //third root
    orig.reset();
    z0 = f.get_last_root();
    its = grad_descent_seq(orig, z0, n, delta);
    r[2] = {z0, orig.num_calls(), its};

    print_gd_results(r);
}

void test_gd_fastest(Complex z0, double delta) {
    freopen("gd_fastest.csv", "w", stdout);
    printf("Iteration, Calls, Z, alpha * Gradient, f(z - alpha * gradient)\n");

    int its, lits, sits;
    Function orig, f;
    Result r[3];

    //first root
    its = grad_descent_fastest(f, z0, delta);
    r[0] = { z0, f.num_calls(), its };
    f.divide(z0);

    //second root
    f.reset(); orig.reset();
    lits = grad_descent_fastest(f, z0, delta * 10.0);
    sits = grad_descent_fastest(orig, z0, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.divide(z0);

    //third root
    orig.reset();
    z0 = f.get_last_root();
    its = grad_descent_fastest(orig, z0, delta);
    r[2] = {z0, orig.num_calls(), its};

    print_gd_results(r);
}

int main() {
    Complex start = estm_min_loc(-10, 10, -10, 10, 4);
    printf("%f, %f\n", start.real(), start.imag());
    const int R = 4;
    test_coord_descent(start, 0.9 + 0.9i, R);
    test_gd_div(start, 0.01, 0.01);
    test_gd_const(start, 1e-3, 0.01);
    test_gd_seq(start, 1, 0.01);
    test_gd_fastest(start, 0.01);
}
