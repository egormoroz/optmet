#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include "display.hpp"
#include "func.hpp"

using namespace std;

vector<it::Passive> pits;
vector<it::Localizer> bits;
vector<it::Localizer> gits;
vector<it::Localizer> fits;
vector<it::Tangent> tanits; //o_O
vector<it::Newton> nits;
vector<it::Newton> sits;

double passive(double a, double b, double step);
double dichotomy(double a, double b, int r);
double golden_ratio(double a, double b, int r);
double fibonacci(double a, double b, double eps);

double tangent(double a, double b, int r);
double newton(double a, double b, double x0, double eps);
double secant(double a, double b, double x1, double x2, double eps);

int main() {
    const double LEFT = -1.0;
    const double RIGHT = 1.0;
    const double EPS = 1e-5;

    const int DIGITS = 5;

    passive(LEFT, RIGHT, EPS);
    dichotomy(LEFT, RIGHT, DIGITS);
    golden_ratio(LEFT, RIGHT, DIGITS);
    fibonacci(LEFT, RIGHT, EPS);

    tangent(LEFT, RIGHT, DIGITS);
    newton(LEFT, RIGHT, RIGHT, EPS);
    secant(LEFT, RIGHT, RIGHT, RIGHT - EPS, EPS);

    print_iteration_table(pits, bits, gits, fits, tanits, nits, sits);
    printf("\n\n");
    auto results = extract_results(pits, bits, gits, fits, tanits, nits, sits);
    print_results(results.data(), results.size());
}

double passive(double a, double b, double step) {
    Function f;

    double min_x = a, min_y = f.eval(a);
    pits.push_back({ min_x, f.num_of_calls() });

    
    for (a += step; a <= b; a += step) {
        double y = f.eval(a);
        if (y < min_y) {
            min_x = a;
            min_y = y;
        } 
        pits.push_back({ min_x, f.num_of_calls() });
    }

    return min_x;
}

double dichotomy(double a, double b, int r) {
    Function f;
    double middle;
    double c, d, cy, dy;

    double K = 1.0 / (2.0 * pow(10, r));

    do {
        middle = (a + b) / 2;
        double delta = (b - a) / 20;

        c = middle - delta, d = middle + delta,
            cy = f.eval(c), dy = f.eval(d);

        bits.push_back({ a, b, c, d, f.num_of_calls() });

        if (cy < dy)
            b = d;
        else
            a = c;
    } while ((b - a) / 2 > K * abs(middle));

    bits.push_back({ a, b, c, d, f.num_of_calls() }); //just for debug info, not an actual iteration

    return (a + b) / 2;
}

double golden_ratio(double a, double b, int digits) {
    Function f;

    double K = 1.0 / (2.0 * pow(10, digits));
    double P = (3 - sqrt(5)) / 2,
           Q = (sqrt(5) - 1) / 2;

    double c = P * (b - a) + a, cy = f.eval(c),
           d = Q * (b - a) + a, dy = f.eval(d);
    gits.push_back({ a, b, c, d, f.num_of_calls() });
    if (cy < dy) {
        b = d;
        d = c; dy = cy;
        c = P * (b - a) + a; cy = f.eval(c);
    } else {
        a = c;
        c = d;cy = dy;
        d = Q * (b - a) + a; dy = f.eval(d);
    }

    while (true) {
        gits.push_back({ a, b, c, d, f.num_of_calls() });
        if (cy < dy) {
            b = d;
            d = c; dy = cy;
            c = P * (b - a) + a; 
            if (b - a < K * abs(a + b))
                break;
            cy = f.eval(c);
        } else {
            a = c;
            c = d;cy = dy;
            d = Q * (b - a) + a; 
            if (b - a < K * abs(a + b))
                break;
            dy = f.eval(d);
        }
    }

    gits.push_back({ a, b, c, d, f.num_of_calls() }); //fake iteration just for debug purposes

    return (a + b) / 2;
}

double fibonacci(double a, double b, double eps) {
    Function f;
    vector<double> fib = { 1.0, 1.0 };
    int64_t p = 1, q = 1;
    while (fib.back() < 11.0 / 20.0 * (b-a) / eps) {
        size_t n = fib.size();
        int64_t t = p + q;
        p = q; q = t;
        fib.push_back(q);
    }

    size_t n = fib.size() - 1;
    double c = a + (b - a) * fib[n - 2] / fib[n],
           cy = f.eval(c),
           d = a + (b - a) * fib[n - 1] / fib[n],
           dy = f.eval(d);
    fits.push_back({ a, b, c, d, f.num_of_calls() });
    n--;
    if (cy < dy) {
        b = d;
        d = c; dy = cy;
        c = a + (b - a) * fib[n - 2] / fib[n];
        cy = f.eval(c);
    } else {
        a = c;
        c = d; cy = dy;
        d = a + (b - a) * fib[n - 1] / fib[n],
        dy = f.eval(d);
    }

    while (true) {
        fits.push_back({ a, b, c, d, f.num_of_calls() });
        n--;
        if (cy < dy) {
            b = d;
            d = c; dy = cy;
            c = a + (b - a) * fib[n - 2] / fib[n];
            if (n <= 2)
                break;
            cy = f.eval(c);
        } else {
            a = c;
            c = d; cy = dy;
            d = a + (b - a) * fib[n - 1] / fib[n];
            if (n <= 2)
                break;
            dy = f.eval(d);
        }
    }
    double middle = c,
           delta = 0.1 * (b - a) / fib.back();
    c = middle + delta;
    d = middle - delta;

    if (f.eval(c) < f.eval(d)) {
        b = d;
    } else {
        a = c;
    }

    fits.push_back({ a, b, c, d, f.num_of_calls() });

    return (a + b) / 2;
}

double tangent(double a, double b, int r) {
    Function f;
    double K = 1.0 / (2.0 * pow(10, r));

    double ay = f.eval(a), by = f.eval(b),
           ady = f.eval_der(a), bdy = f.eval_der(b),
           c, cdy;
    while (true) {
        c = (by - b*bdy - ay + a*ady) / (ady - bdy);
        cdy = f.eval_der(c);
        tanits.push_back({ a, b, c, f.num_of_calls() });
        if (cdy > 0) {
            b = c;
            if (b - a < K * abs(a + b))
                break;
            by = f.eval(c);
            bdy = cdy;
        } else if (cdy < 0) {
            a = c;
            if (b - a < K * abs(a + b))
                break;
            ay = f.eval(c);
            ady = cdy;
        } else { //got lucky
            break;
        }
    }

    tanits.push_back({ a, b, c, f.num_of_calls() });

    return c;
}

double newton(double a, double b, double x, double eps) {
    Function f;

    double d = f.eval_der(x); 
    while (abs(d) > eps) {
        double dd = f.eval_der2(x);
        x = x - d/dd;
        d = f.eval_der(x);
        nits.push_back({ x, f.num_of_calls() });
    }

    return x;
}

double secant(double a, double b, double x1, double x2, double eps) {
    Function f;

    double d1 = f.eval_der(x1), d2 = f.eval_der(x2);
    while (abs(d2) > eps) {
        double dd = (d2 - d1) / (x2 - x1);
        x1 = x2; d1 = d2;
        x2 = x2 - d2/dd;
        d2 = f.eval_der(x2);
        sits.push_back({ x2, f.num_of_calls() });
    }

    return x2;
}

