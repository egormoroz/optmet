#pragma once

#include <cstdio>
#include <vector>
#include <algorithm>
#include <array>
#include "func.hpp"

/*
* Информация об итерации соответсвующего метода
* (a; b) - локализирующий интервал
* x - приблизительная точка минимума
* calls - количество вызовов функций f(x), f'(x), f''(x) (т.е. трудоёмкость)
*/
namespace it {
    //Пассивный поиск
    struct Passive {
        double x;
        int calls;
    };

    //Дихотомия, ЗС, фибоначчи
    struct Localizer {
        double a, b;
        double c, d;
        int calls;
    };

    //Касательные
    struct Tangent {
        double a, b;
        double c;
        int calls;
    };

    //Ньютон-Рафсон, секущие
    struct Newton {
        double x;
        int calls;
    };
}

/*
* Принятые корткие обозначения массивов:
* ps - passive search (пассивный поиск)
* bs - binary search (дихотомия)
* gs - golden search (ЗС)
* fs - fibonacci search (фибоначчи)
* ts - tangent search (касательные)
* ns - newton search (ньютон)
* ss - secant search (секущие)
* 
* Они везде одни и те же
*/

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

    /*
    * Как понять, что выводится в printf'ах:
    * 1. Посмотреть на printf выше с названиями столбцов
    * 2. Посмотреть на printf ниже (он находится в цикле):
    * в нём значения, разделённые прямой палочкой - |, идут в том же порядке, что и названия столбцов
    */

    printf("passive search\n");
    printf("     # |  Calls |        (x, f(x))    \n"); //названия столбцов, всего 4 числа
    printf("--------------------------------------\n");
    int step = (ps.size() + 99) / 100;
    for (int k = 0; k < ps.size(); k += step) {
        auto &i = ps[k];
        printf("%6d | %6d | (%9f, %9f)\n",
               k+1, i.calls, i.x, f.eval(i.x)); //номер итерации, кол-во вызовов, (точка, значение функции)
    }

    printf("------------------------------------------------------------------|");
    printf("----------------------------------------------------------|");
    printf("---------------------------------------------------------\n");
    printf("                         dichotomy                                |");
    printf("                        golden ratio                      |");
    printf("                          fibonacci                               \n");
    printf("     # |  Calls |         [a; b]         |         (x, f(x))      | "); //первые 4 стоблца (дихотомия)
    printf(" Calls |         [a; b]         |         (x, f(x))      |"); //след. 3 столбца (ЗС)
    printf(" Calls |         [a; b]         |         (x, f(x))      \n");//последние 3 столбца (фибоначчи)
    printf("------------------------------------------------------------------|");
    printf("----------------------------------------------------------|");
    printf("---------------------------------------------------------\n");
    size_t n = std::max({bs.size(), gs.size(), fs.size()});
    for (int k = 0; k < n; ++k) {
        printf("%6d | ", k + 1); //номер итерации
        //итерации дихотомия
        if (k < bs.size()) {
            auto &i = bs[k];
            double x = (i.a + i.b) / 2;
            printf("%6d | [%9f; %9f] | (%9f, %9f) | ",
                   i.calls, i.a, i.b, x, f.eval(x));
        } else {
            printf("       |                        |                        | ");
        }
        //итерации золотого сечения
        if (k < gs.size()) {
            auto &i = gs[k];
            double x = (i.a + i.b) / 2;
            printf("%6d | [%9f; %9f] | (%9f, %9f) |",
                    i.calls, i.a, i.b, x, f.eval(x));
        } else {
            printf("       |                        |                        |");
        }
        //итерации фибоначчи
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
    printf("     # |  Calls |         [a; b]         |         (c, f(c))      | ");//касательные
    printf(" Calls |         (x, f(x))      |"); //ньютон-рафсон
    printf(" Calls |         (x, f(x))      \n"); //секущие
    printf("------------------------------------------------------------------|");
    printf("---------------------------------|");
    printf("--------------------------------\n");
    n = std::max({ts.size(), ns.size(), ss.size()});
    for (int k = 0; k < n; ++k) {
        printf("%6d | ", k + 1);
        //касательные
        if (k < ts.size()) {
            auto &i = ts[k];
            printf("%6d | [%9f; %9f] | (%9f, %9f) | ",
                   i.calls, i.a, i.b, i.c, f.eval(i.c));
        } else {
            printf("       |                        |                        | ");
        }
        //ньютон-рафсон
        if (k < ns.size()) {
            auto &i = ns[k];
            printf("%6d | (%9f, %9f) |",
                    i.calls, i.x, f.eval(i.x));
        } else {
            printf("       |                        |");
        }
        //секущие
        if (k < ss.size()) {
            auto &i = ss[k];
            printf("%6d | (%9f, %9f)",
                    i.calls, i.x, f.eval(i.x));
        }
        printf("\n");
    }
 
}

struct Result {
    const char *name; //название метода
    size_t its; //кол-во итераций
    int calls; //трудоёмкость
    double x; //результат (т.е. примерный минимум)
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
    /*
    * 1. Из всех массивов итераций изъяли информацию с последней итерации
    * (последняя итерация всегда содержит самую последнюю информацию: для этого были созданы итерации болванки
    * в которых ничего не вычислялось (calls не изменилось))
    * 2. Изъятую инфу положили в структуру Result с соотв. именем
    * 3. Вернули все Result'ы в качестве массива
    */
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
        /*
        * %.3e - числа в формате X.YZWeN, где X,Y,Z,W - цифры, 
        * N - степень 10 (3 символа: + или - и число из двух цифр).
        * Так же printf автоматически округляет число по правилам математики.
        * Пример:
        * printf("%.3e", 314.159);
        * напечатает 3.142e+02
        */ 
        printf(" %12s | %10d | %9d | (%.3e, %.3e)\n",
               r.name, int(r.its), r.calls, r.x, f.eval(r.x));
    }
}

