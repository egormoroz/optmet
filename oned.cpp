#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class Function {
    int m_counter = 0;
public:
    Function() = default;

    double eval(double x) {
        m_counter++;
        return eval_without_inc(x);
    }

    double eval_without_inc(double x) const {
        const double A = 4;
        const double B = 5;
        const double C = 9;

        //can be improved
        return exp(B * x) + exp(-A * x) + A * x * x + C * x + B - x * x * x;
    }

    void reset() { m_counter = 0; }

    int num_of_calls() const { return m_counter; }

} f;

struct Result {
    const char* name = 0;
    double x = 0.0;
    int num_calls = 0;
    int num_iterations = 0;
};


Result search_passive(double a, double b, double step);
Result bin_search(double a, double b, double delta, double tol);

void print_results(const Result* r, int n, int p, int q);

int main() {
    const double LEFT = -1.0;
    const double RIGHT = 1.0;
    const double TOLERANCE = 1e-5; //four digits + one extra in case of rounding error
    const double DELTA = 1e-6;

    ios::sync_with_stdio(false);

    Result r[2];

    r[0] = search_passive(LEFT, RIGHT, TOLERANCE);
    r[1] = bin_search(LEFT, RIGHT, DELTA, TOLERANCE);

    printf("\n\n");
    print_results(r, size(r), 2, 4);
}

void print_iter_info(int i, double a, double b, double x, double y, bool chk = false) {
    printf("%10d) bounds: [%8.5f; %8.5f], point: (%8.5f, %8.5f)%s\n",
        i, a, b, x, y, chk ? "\tx" : "");
}

Result search_passive(double a, double b, double step) {
    f.reset();
    Result r;
    r.name = "passive search";
    printf("%s\n", r.name);

    double min_x = a, min_y = f.eval(a);
    print_iter_info(1, a, b, min_x, min_y);
    r.num_iterations++;
    
    //stop after having scanned the interval
    for (a += step; a <= b; a += step, r.num_iterations++) {
        double y = f.eval(a);
        bool should_print = false;
        bool chk = false;
        if (y < min_y) {
            min_x = a;
            min_y = y;
            chk = true;
            should_print = r.num_iterations % 1'00 == 0;
        } 
        should_print |= r.num_iterations % 1'000 == 0;

        if (should_print)
            print_iter_info(r.num_iterations + 1, a, b, min_x, min_y, chk);
    }

    r.x = min_x;
    r.num_calls = f.num_of_calls();
    return r;
}

Result bin_search(double a, double b, double delta, double tol) {
    f.reset();
    Result r;
    r.name = "binary search";
    printf("%s\n", r.name);

    double middle;
    double c, d, cy, dy;

    do {
        middle = (a + b) / 2;
        c = middle - delta / 2, d = middle + delta / 2,
            cy = f.eval(c), dy = f.eval(d);
        print_iter_info(r.num_iterations + 1, a, b, middle, f.eval_without_inc(middle));
        
        if (cy < dy)
            b = d;
        else
            a = c;

        r.num_iterations++;
    } while ((b - a) / 2 > tol);
    //we are done when half of the remaining interval is less than tolerance

    r.x = middle;
    r.num_calls = f.num_of_calls();
    return r;
}

void print_results(const Result* r, int n, int p, int q) {
    const char header[] = "  #  |       Name       | Iterations |   Calls    |     Result    ";
    char sep[size(header)];
    memset(sep, '-', size(header));
    sep[size(header) - 1] = 0;
    printf("%s\n", header);
    int total = p + q + 2;
    for (int i = 0; i < n; ++i) {
        double y = f.eval_without_inc(r[i].x);
        printf("%s\n", sep);
        printf("  %2d | %16s | %10d | %10d | (%*.*f, %*.*f)\n",
            i + 1, r[i].name, r[i].num_iterations, r[i].num_calls,
            total, q, r[i].x, total, q, y);
    }
}

