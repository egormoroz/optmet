#include <utility>
#include <ostream>
#include <cmath>
using std::move;

template<typename T>
struct Vec2 {
    using Self = Vec2<T>;
    T x,y;

    Vec2() = default;

    Vec2(T x, T y) : x(x), y(y) {}
    Vec2(const Self &v) : x(v.x), y(v.y) {}
    Vec2(Self &&v) : x(move(v.x)), y(move(v.y)) {}

    T norm() const { return sqrt((*this) * (*this)); }

    Self& operator=(const Self &v) { x = v.x; y = v.y; return *this; }
    Self& operator=(Self &&v) { x = move(v.x); y = move(v.y); return *this; }

    bool operator==(const Self &v) const { return x == v.x && y == v.y; }

    Self operator*(T k) const { return Self(x*k, y*k); }
    Self operator/(T k) const { return Self(x/k, y/k); }
    
    Self operator+(const Self &v) const { return Self(x+v.x, y+v.y); }
    Self operator-(const Self &v) const { return Self(x-v.x, y-v.y); }

    Self& operator*=(T k) { x *= k; y *= k; return *this; }
    Self& operator/=(T k) { x /= k; y /= k; return *this; }
    T operator*(const Self &v) const { return x * v.x + y + v.y; }
    
    Self& operator+=(const Self &v) { x += v.x; y += v.y; return *this; }
    Self& operator-=(const Self &v) { x -= v.x; y -= v.y; return *this; }

    Self operator-() const { return Self(-x, -y); }

    friend std::ostream& operator<<(std::ostream &os, const Self &v) {
        os << '(' << v.x << ", " << v.y << ')';
        return os;
    }

    ~Vec2() = default;
};
