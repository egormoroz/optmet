#include <cstdio>
#include "vec.hpp"

using V2 = Vec2<double>;

class Function {
public:
    double eval(const V2 &p) {
        m_counter++;
        return 9 * p.x * p.x + p.y * p.y;
    }

private:
    int m_counter = 0;
};

bool try_advance(Function &f, V2 &p, double &z, V2 step) {
    V2 pp = p + step;
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

V2 coord_descent(V2 p, V2 step, double eps) {
    Function f;
    V2 dx{step.x, 0.}, dy{0., step.y}, pp;
    double z = f.eval(p), zz;
    double K = 2.0;

    while (true) {
        pp = p; zz = z;
        int flag = 0;
        if (!try_advance(f, p, z, dx) && !try_advance(f, p, z, -dx)) {
            dx /= K;
            flag++;
        }
        dx /= 1.1;
        if (!try_advance(f, p, z, dy) && !try_advance(f, p, z, -dy)) {
            dy /= K;
            flag++;
        }
        dy /= 1.1;

        if ((dx + dy).norm() < eps)
            break;

        if (flag == 2) {
            K *= 2; //If the step is too big, we decrease it exponentially
        } else if (zz - z < eps) {
            break;
        } else {
            K = 2.0;
        }
    } 

    return p;
}

int main() {
    coord_descent({1.0, 1.0}, {0.33, 0.33}, 1e-4);
}
