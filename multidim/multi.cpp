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
        printf("%d, %d, (%f; %f), 0, %f\n", i.iterations, i.calls, 
            i.z.real(), i.z.imag(), f.eval(i.z));
    }
}


double muls(Complex x, Complex y) {
    return x.real() * y.real() + x.imag() * y.imag();
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

int coord_descent(Function &f, Complex &p, Complex &step, int r) {
    Complex dx(step.real(), 0.), dy(0., step.imag()), pp;
    double z = f.eval(p), zz;

    double K = 1.0 / (2.0 * pow(10, r)),
        old_dx, old_dy;

    int k = 0;
    double D = 2.0;

    do {
        pp = p; zz = z;
        old_dx = dx.real(); old_dy = dy.imag();

        int flag = 0;
        if (!try_advance(f, p, z, dx) && !try_advance(f, p, z, -dx)) {
            dx /= D;
            flag++;
        }

        if (!try_advance(f, p, z, dy) && !try_advance(f, p, z, -dy)) {
            dy /= D;
            flag++;
        }

        D = flag == 2 ? D * 2.0 : 2.0;

        printf("%d, %d, %f, %f, (%f; %f), %f\n", ++k, f.num_calls(), 
            dx.real(), dy.imag(), p.real(), p.imag(), f.dbg_eval(p));

        dx /= 1.1;
        dy /= 1.1;
    } while (old_dx > K * abs(p.real()) || old_dy > K * abs(p.imag()));

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
            alpha *= L;
            printf("*");
            continue;
        }
        log_gd(++k, f.num_calls(), z0, alpha * grad, new_val);

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
    const double EPS = 1e-5; //if bigger then method diverges
    Complex grad = f.eval_grad(z);
    int k = 0;
    while (norm_sqr(grad) >= delta) {
        double alpha = golden_ratio([&](double alpha) {
            return f.eval(z - alpha * grad);
        }, 0, 0.05, EPS);
        auto ag = alpha * grad;
        log_gd(++k, f.num_calls(), z, ag, f.dbg_eval(z - ag));

        z = z - alpha * grad;
        grad = f.eval_grad(z);
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
    f.exclude_root(z0);

    //second root
    f.reset(); orig.reset();
    step = stp;
    lits = coord_descent(f, z0, step, R - 1);
    sits = coord_descent(orig, z0, step, R);
    r[1] = { z0, f.num_calls() + orig.num_calls(), lits + sits };
    f.exclude_root(z0);

    //third root
    orig.reset();
    double K = 1.0 / (2.0 * pow(10, R));
    z0 = f.get_last_root();
    step = Complex(abs(z0.real()), abs(z0.imag())) * K;
    its = coord_descent(orig, z0, step, R);
    r[2] = { z0, orig.num_calls(), its };

    for (auto& i : r) {
        printf("%d, %d, 0, 0, (%f; %f), %f\n", i.iterations, i.calls, 
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
    f.exclude_root(z0);

    //second root
    f.reset(); orig.reset();
    lits = grad_descent_div(f, z0, alpha, delta);
    sits = grad_descent_div(orig, z0, alpha, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.exclude_root(z0);

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
    f.exclude_root(z0);

    //second root
    f.reset(); orig.reset();
    lits = grad_descent_const(f, z0, alpha, delta);
    sits = grad_descent_const<true>(orig, z0, alpha, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.exclude_root(z0);

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
    f.exclude_root(z0);

    //second root
    f.reset(); orig.reset();
    m = n;
    lits = grad_descent_seq(f, z0, m, delta);
    sits = grad_descent_seq(orig, z0, m, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.exclude_root(z0);

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
    f.exclude_root(z0);

    //second root
    f.reset(); orig.reset();
    lits = grad_descent_fastest(f, z0, delta * 10.0);
    sits = grad_descent_fastest(orig, z0, delta);
    r[1] = {z0, f.num_calls() + orig.num_calls(), lits + sits};
    f.exclude_root(z0);

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
    test_gd_seq(start, 1000, 0.01);
    test_gd_fastest(start, 0.01);
}
