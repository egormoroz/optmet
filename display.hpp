#pragma once

#include <cstdio>
#include <vector>
#include <algorithm>
#include <array>
#include "func.hpp"

namespace it {
    struct Passive {
        double x;
        int calls;
    };

    struct Localizer {
        double a, b;
        double c, d;
        int calls;
    };

    struct Tangent {
        double a, b;
        double c;
        int calls;
    };

    struct Newton {
        double x;
        int calls;
    };
}

void print_iteration_table(
        const std::vector<it::Passive> &ps, 
        const std::vector<it::Localizer> &bs,
        const std::vector<it::Localizer> &gs,
        const std::vector<it::Localizer> &fs,
        const std::vector<it::Tangent> &ts,
        const std::vector<it::Newton> &ns,
        const std::vector<it::Newton> &ss)
{
    Function f;
    printf("passive search\n");
    printf("     # |  Calls |        (x, f(x))    \n");
    printf("--------------------------------------\n");
    int step = (ps.size() + 99) / 100;
    for (int k = 0; k < ps.size(); k += step) {
        auto &i = ps[k];
        printf("%6d | %6d | (%9f, %9f)\n",
               k+1, i.calls, i.x, f.eval(i.x));
    }
    printf("------------------------------------------------------------------|");
    printf("----------------------------------------------------------|");
    printf("---------------------------------------------------------\n");
    printf("                         dichotomy                                |");
    printf("                        golden ratio                      |");
    printf("                          fibonacci                               \n");
    printf("     # |  Calls |         [a; b]         |         (x, f(x))      | ");
    printf(" Calls |         [a; b]         |         (x, f(x))      |");
    printf(" Calls |         [a; b]         |         (x, f(x))      \n");
    printf("------------------------------------------------------------------|");
    printf("----------------------------------------------------------|");
    printf("---------------------------------------------------------\n");
    size_t n = std::max({bs.size(), gs.size(), fs.size()});
    for (int k = 0; k < n; ++k) {
        printf("%6d | ", k + 1);
        if (k < bs.size()) {
            auto &i = bs[k];
            double x = (i.a + i.b) / 2;
            printf("%6d | [%9f; %9f] | (%9f, %9f) | ",
                   i.calls, i.a, i.b, x, f.eval(x));
        } else {
            printf("       |                        |                        | ");
        }
        if (k < gs.size()) {
            auto &i = gs[k];
            double x = (i.a + i.b) / 2;
            printf("%6d | [%9f; %9f] | (%9f, %9f) |",
                    i.calls, i.a, i.b, x, f.eval(x));
        } else {
            printf("       |                        |                        |");
        }
        if (k < fs.size()) {
            auto &i = fs[k];
            double x = (i.a + i.b) / 2;
            printf("%6d | [%9f; %9f] | (%9f, %9f)",
                    i.calls, i.a, i.b, x, f.eval(x));
        }
        printf("\n");
    }
    
    printf("\n");
    printf("------------------------------------------------------------------|");
    printf("---------------------------------|");
    printf("--------------------------------\n");
    printf("                          tangent                                 |");
    printf("             newton              |");
    printf("             secant             \n");
    printf("     # |  Calls |         [a; b]         |         (c, f(c))      | ");
    printf(" Calls |         (x, f(x))      |");
    printf(" Calls |         (x, f(x))      \n");
    printf("------------------------------------------------------------------|");
    printf("---------------------------------|");
    printf("--------------------------------\n");
    n = std::max({ts.size(), ns.size(), ss.size()});
    for (int k = 0; k < n; ++k) {
        printf("%6d | ", k + 1);
        if (k < ts.size()) {
            auto &i = ts[k];
            printf("%6d | [%9f; %9f] | (%9f, %9f) | ",
                   i.calls, i.a, i.b, i.c, f.eval(i.c));
        } else {
            printf("       |                        |                        | ");
        }
        if (k < ns.size()) {
            auto &i = ns[k];
            printf("%6d | (%9f, %9f) |",
                    i.calls, i.x, f.eval(i.x));
        } else {
            printf("       |                        |");
        }
        if (k < ss.size()) {
            auto &i = ss[k];
            printf("%6d | (%9f, %9f)",
                    i.calls, i.x, f.eval(i.x));
        }
        printf("\n");
    }
 
}

struct Result {
    const char *name;
    size_t its;
    int calls;
    double x;
};

std::array<Result, 7> extract_results (
        const std::vector<it::Passive> &ps, 
        const std::vector<it::Localizer> &bs,
        const std::vector<it::Localizer> &gs,
        const std::vector<it::Localizer> &fs,
        const std::vector<it::Tangent> &ts,
        const std::vector<it::Newton> &ns,
        const std::vector<it::Newton> &ss)

{
    auto x = [](auto &i) { return (i.a + i.b) / 2; };
    return {
        Result{ "passive", ps.size(), ps.back().calls, ps.back().x },
        { "dichotomy", bs.size(), bs.back().calls, x(bs.back()) },
        { "golden ratio", gs.size(), gs.back().calls, x(gs.back()) },
        { "fibonacci", fs.size(), fs.back().calls, x(fs.back()) },
        { "tangent", ts.size(), ts.back().calls, ts.back().c },
        { "newton", ns.size(), ns.back().calls, ns.back().x },
        { "secant", ss.size(), ss.back().calls, ss.back().x }
    };
}

void print_results(const Result *rs, int n) {
    Function f;
    printf("----------------------------------------------------------------\n");
    printf("     name     | iterations |   calls   |          (x, y)      \n");
    printf("----------------------------------------------------------------\n");
    for (int i = 0; i < n; ++i) {
        auto &r = rs[i];
        printf(" %12s | %10d | %9d | (%.3e, %.3e)\n",
               r.name, int(r.its), r.calls, r.x, f.eval(r.x));
    }
}

