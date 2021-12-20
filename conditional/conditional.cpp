#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include "func.hpp"
#include <algorithm>
using namespace std;

double g(Complex z) {
    const double R = 6;
    double x = z.real(), y = z.imag();
    return x * x + y * y - R * R;
}

Complex g_grad(Complex z) {
    return z * 2.;
}

double my_cos(Complex a, Complex b) {
    return muls(a, b) / (norm(a) * norm(b));
}

void exterior_grad(Function &f, Complex &z, double step, double min_step, double eps) {
    double r = 1.0;
    int k = 1;
    auto phi = [&](Complex z) {
        double t = max(0.0, g(z));
        double val = f.eval(z), penalty = r * t * t;

        auto fg = f.dbg_eval_grad(z), gg = g_grad(z);
        double cos_val = my_cos(-fg, gg);
        double angle = acos(cos_val) * 180 / M_PI;
        printf("%d; %d; (%f, %f); %f; %f; %f; (%f, %f); (%f, %f); %f; %f\n", 
            k, f.num_calls(), z.real(), z.imag(), val, g(z), val + penalty,
            fg.real(), fg.imag(), gg.real(), gg.imag(), cos_val, angle);

        return val + penalty;
    };

    while(true) {
        auto fg = f.eval_grad(z), gg = g_grad(z);
        auto grad = fg + 2 * max(0.0, g(z)) * gg;
        grad /= norm(grad);
        double st = step, alpha = 0, val = phi(z);
        while (st > min_step) {
            double beta = alpha + st, new_val = phi(z - beta * grad);
            if (new_val > val) {
                st /= 2;
                continue;
            }

            alpha = beta;
            val = new_val;
        }
        z -= alpha * grad;

        gg = g_grad(z);
        fg = f.dbg_eval_grad(z);
        double cos_val = my_cos(-fg, gg);
        double angle = acos(cos_val) * 180 / M_PI;
        printf("*%d; %d; (%f, %f); %f; %f; %f; (%f, %f); (%f, %f); %f; %f\n", 
            k, f.num_calls(), z.real(), z.imag(), f.dbg_eval(z), g(z), val,
            fg.real(), fg.imag(), gg.real(), gg.imag(), cos_val, angle);

        if (g(z) < eps && my_cos(-f.eval_grad(z), gg) > 0.9994)
            break;

        ++k;
        r *= 2;
        z -= gg / norm(gg) * step;
    }
}

void interior_grad(Function &f, Complex &z, double step, double min_step) {
    double t = 1.0;
    int k = 1;
    auto psi = [&](Complex z) {
        double val = f.eval(z), penalty = -t / g(z);

        Complex fg = f.dbg_eval_grad(z);
        Complex gg = g_grad(z);
        double cos_val = my_cos(-fg, gg);
        double angle = acos(cos_val) * 180 / M_PI;
        printf("%d; %d; (%f, %f); %f; %f; %f; (%f, %f); (%f, %f); %f; %f\n", 
            k, f.num_calls(), z.real(), z.imag(), f.dbg_eval(z), g(z), val,
            fg.real(), fg.imag(), gg.real(), gg.imag(), my_cos(fg, gg), angle);

        return val + penalty;
    };

    while(true) {
        auto fg = f.eval_grad(z), gg = g_grad(z);
        auto grad = fg + t / (g(z)*g(z));
        grad /= norm(grad);

        double st = step, alpha = 0, val = psi(z);
        while (st > min_step) {
            double beta = alpha + st;
            auto zz = z - beta * grad;
            if (g(zz) > 0) {
                st /= 2;
                continue;
            }
            double new_val = psi(zz);
            if (new_val > val) {
                st /= 2;
                continue;
            }

            alpha = beta;
            val = new_val;
        }

        z -= alpha * grad;

        gg = g_grad(z);
        fg = f.dbg_eval_grad(z);
        double cos_val = my_cos(-fg, gg);
        double angle = acos(cos_val) * 180 / M_PI;
        printf("*%d; %d; (%f, %f); %f; %f; %f; (%f, %f); (%f, %f); %f; %f\n", 
            k, f.num_calls(), z.real(), z.imag(), f.dbg_eval(z), g(z), val,
            fg.real(), fg.imag(), gg.real(), gg.imag(), cos_val, angle);

        if (my_cos(-f.eval_grad(z), gg) > 0.9994)
            break;

        ++k;
        t /= 2;
        z -= gg / norm(gg) * step;
    }
}

void surf_border(Function& f, Complex& z, double eps) {
    const double EPS = 1e-5;
    const double L = 0.5;
    const double APPROACH_STEP = 1.;
    int k = 1;

    Complex fg = f.eval_grad(z);
    while (abs(g(z)) > EPS) {
        fg /= norm(fg);
        if (g(z - APPROACH_STEP * fg) > EPS) {
            z -= g(z) / muls(g_grad(z), fg) * fg;
        } else {
            z -= APPROACH_STEP * fg;
        }
        fg = f.eval_grad(z);
        printf("%d; %d; (%f, %f); %f; %f; (%f, %f); -; -; -; -\n", 
            k, f.num_calls(), z.real(), z.imag(), f.dbg_eval(z), g(z),
            fg.real(), fg.imag());
        ++k;
    }

    Complex gg = g_grad(z), prev_tan = 0.;
    double cos = my_cos(-fg, gg);
    const double TARGET_COS = 0.9994;
    double alpha = 0.1;
    do {
        Complex tan = -fg - gg * muls(-fg, gg) / norm_sqr(gg);
        if (muls(tan, prev_tan) < 0)
            alpha /= 2;
        prev_tan = tan;

        tan /= norm(tan);
        z += tan * alpha;

        if (g(z) > EPS)
            z -= g(z) / norm_sqr(gg) * gg;

        fg = f.eval_grad(z);
        gg = g_grad(z);
        cos = my_cos(-fg, gg);
        ++k;

        printf("%d; %d; (%f, %f); %f; %f; (%f, %f); (%f, %f); (%f, %f); %f; %f\n", 
            k, f.num_calls(), z.real(), z.imag(), f.dbg_eval(z), g(z),
            fg.real(), fg.imag(), gg.real(), gg.imag(), tan.real(), tan.imag(), cos,
            acos(cos) * 180 / M_PI);
    } while (cos < TARGET_COS);
}


int main() {
    Complex z(0.);
    Function f({4, 5}, {4, 9}, {9, 4});

    /* freopen("exterior.csv", "w", stdout); */
    /* printf("Iteration; Calls; (x, y); f(x, y); g(x, y);" */
    /*         "phi(x, y); grad_f; grad_g; cos; angle\n"); */
    /* exterior_grad(f, z, 0.1, 1e-3, 1e-3); */

    /* f.reset(); */
    /* z = 0.; */
    /* freopen("interior.csv", "w", stdout); */
    /* printf("Iteration; Calls; (x, y); f(x, y); g(x, y);" */
    /*        "psi(x, y); grad_f; grad_g; cos; angle\n"); */
    /* interior_grad(f, z, 0.1, 1e-3); */

    /* f.reset(); */
    /* z = 0.; */
    /* freopen("surf_border.csv", "w", stdout); */
    printf("Iteration; Calls; (x, y); f(x, y); g(x, y); grad_f; grad_g;"
        "alpha * tan; cos(-grad_f, grad_g); angle\n");
    surf_border(f, z, 0.);
}

