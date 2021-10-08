#pragma once

template<typename F>
double golden_ratio(F f, double a, double b, double eps) {
    double P = (3 - sqrt(5)) / 2,
           Q = (sqrt(5) - 1) / 2;

    double c = P * (b - a) + a, cy = f(c),
           d = Q * (b - a) + a, dy = f(d);
    if (cy < dy) {
        b = d;
        d = c; dy = cy;
        c = P * (b - a) + a; cy = f(c);
    } else {
        a = c;
        c = d;cy = dy;
        d = Q * (b - a) + a; dy = f(d);
    }

    while (true) {
        if (cy < dy) {
            b = d;
            d = c; dy = cy;
            c = P * (b - a) + a; 
            if (b - a < 2 * eps)
                break;
            cy = f(c);
        } else {
            a = c;
            c = d;cy = dy;
            d = Q * (b - a) + a; 
            if (b - a < 2 * eps)
                break;
            dy = f(d);
        }
    }

    return (a + b) / 2;
}

