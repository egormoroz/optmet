#pragma once
/*
* Одномерная минимизация для МНГС
* Представляет из себя пассиный поиск, 
* делающий убывающие по длине шаги
* */
template<typename F>
double xyz(F f, double step, double eps) {
    double x = 0., y = f(x);
    const double L = 0.5;

    while (step > eps) {
        double xx = x + step, yy = f(xx);
        if (yy > y) {
            step *= L;
            continue;
        }
        x = xx;
        y = yy;
    }

    return x;
}

