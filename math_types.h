#pragma once
// math_types.h — Vec3 and Mat3, nothing fancy.
// No Eigen, no GLM. Just enough for rigid body mechanics on flat ground.

#include <cmath>
#include <iostream>
#include <algorithm>

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    Vec3 cross(const Vec3& v) const {
        return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
    }

    double length_sq() const { return x * x + y * y + z * z; }
    double length() const { return std::sqrt(length_sq()); }

    Vec3 normalized() const {
        double len = length();
        if (len < 1e-12) return {0, 0, 0};
        return *this / len;
    }
};

inline Vec3 operator*(double s, const Vec3& v) { return v * s; }

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

// Row-major 3x3. m[row][col].
struct Mat3 {
    double m[3][3] = {};

    Mat3() = default;

    Mat3(const Vec3& r0, const Vec3& r1, const Vec3& r2) {
        m[0][0] = r0.x; m[0][1] = r0.y; m[0][2] = r0.z;
        m[1][0] = r1.x; m[1][1] = r1.y; m[1][2] = r1.z;
        m[2][0] = r2.x; m[2][1] = r2.y; m[2][2] = r2.z;
    }

    static Mat3 identity() {
        Mat3 I;
        I.m[0][0] = I.m[1][1] = I.m[2][2] = 1.0;
        return I;
    }

    static Mat3 diagonal(double a, double b, double c) {
        Mat3 D;
        D.m[0][0] = a; D.m[1][1] = b; D.m[2][2] = c;
        return D;
    }

    Vec3 operator*(const Vec3& v) const {
        return {
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
        };
    }

    Mat3 operator*(const Mat3& B) const {
        Mat3 C;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    C.m[i][j] += m[i][k] * B.m[k][j];
        return C;
    }

    Mat3 transposed() const {
        Mat3 T;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                T.m[i][j] = m[j][i];
        return T;
    }

    // General 3x3 inverse via cofactors. Overkill for diagonal inertia tensors
    // but we'll want it once products of inertia come into play.
    Mat3 inverse() const {
        double det =
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

        if (std::abs(det) < 1e-15) return Mat3::identity(); // singular fallback

        double id = 1.0 / det;
        Mat3 inv;
        inv.m[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * id;
        inv.m[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * id;
        inv.m[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id;
        inv.m[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * id;
        inv.m[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id;
        inv.m[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * id;
        inv.m[2][0] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * id;
        inv.m[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * id;
        inv.m[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id;
        return inv;
    }

    // Yaw-only rotation (about Z). SIMPLIFICATION: no pitch/roll — flat ground.
    static Mat3 from_yaw(double yaw_rad) {
        double c = std::cos(yaw_rad);
        double s = std::sin(yaw_rad);
        Mat3 R;
        R.m[0][0] =  c; R.m[0][1] = -s;
        R.m[1][0] =  s; R.m[1][1] =  c;
        R.m[2][2] = 1;
        return R;
    }
};