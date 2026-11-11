// camera.hpp: Camera model, intrinsics, 3D math and RadTan distortion inversion
#pragma once

#include <cmath>
#include <cstddef>
#include <algorithm>

namespace evimo2 {

// Lekki wektor 3D (double)
struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const {
        double inv = (std::abs(s) > 1e-15) ? (1.0 / s) : 0.0;
        return {x * inv, y * inv, z * inv};
    }

    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    double norm() const { return std::sqrt(x * x + y * y + z * z); }
};

// Lekka macierz 3x3 (row-major: m[row * 3 + col])
struct Mat3 {
    double m[9] = {1.0, 0.0, 0.0,
                   0.0, 1.0, 0.0,
                   0.0, 0.0, 1.0};

    Mat3() = default;
    Mat3(std::initializer_list<double> list) {
        size_t idx = 0;
        for (double v : list) {
            if (idx < 9) m[idx++] = v;
        }
    }

    static Mat3 identity() {
        return Mat3{1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0};
    }

    static Mat3 zeros() {
        return Mat3{0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0};
    }

    double operator()(int r, int c) const { return m[r * 3 + c]; }
    double& operator()(int r, int c) { return m[r * 3 + c]; }

    Mat3 operator+(const Mat3& o) const {
        Mat3 res;
        for (int i = 0; i < 9; ++i) res.m[i] = m[i] + o.m[i];
        return res;
    }

    Mat3 operator*(double s) const {
        Mat3 res;
        for (int i = 0; i < 9; ++i) res.m[i] = m[i] * s;
        return res;
    }

    Vec3 operator*(const Vec3& v) const {
        return {
            m[0] * v.x + m[1] * v.y + m[2] * v.z,
            m[3] * v.x + m[4] * v.y + m[5] * v.z,
            m[6] * v.x + m[7] * v.y + m[8] * v.z
        };
    }

    Mat3 operator*(const Mat3& o) const {
        Mat3 res = Mat3::zeros();
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                res.m[r * 3 + c] = m[r * 3 + 0] * o.m[0 * 3 + c] +
                                   m[r * 3 + 1] * o.m[1 * 3 + c] +
                                   m[r * 3 + 2] * o.m[2 * 3 + c];
            }
        }
        return res;
    }

    Mat3 transpose() const {
        return Mat3{
            m[0], m[3], m[6],
            m[1], m[4], m[7],
            m[2], m[5], m[8]
        };
    }
};

// Parametry wewnętrzne kamery (model pinhole + dystorsja RadTan / Brown-Conrady)
struct CameraIntrinsics {
    double fx = 0.0, fy = 0.0, cx = 0.0, cy = 0.0;
    double k1 = 0.0, k2 = 0.0, k3 = 0.0, p1 = 0.0, p2 = 0.0;
};

// Wzór Rodriguesa: wektor obrotu axis-angle rv -> macierz obrotu SO(3)
inline Mat3 rodrigues(const Vec3& rv) {
    double th = rv.norm();
    double a, b;
    if (th < 1e-6) {
        a = 1.0 - th * th / 6.0;
        b = 0.5 - th * th / 24.0;
    } else {
        a = std::sin(th) / th;
        b = (1.0 - std::cos(th)) / (th * th);
    }

    // Macierz skośno-symetryczna K:
    // [  0,  -z,   y ]
    // [  z,   0,  -x ]
    // [ -y,   x,   0 ]
    Mat3 K = Mat3::zeros();
    K.m[1] = -rv.z; K.m[2] = rv.y;
    K.m[3] = rv.z;  K.m[5] = -rv.x;
    K.m[6] = -rv.y; K.m[7] = rv.x;

    Mat3 K2 = K * K;
    return Mat3::identity() + K * a + K2 * b;
}

// Inwersja modelu radtan OpenCV (Brown-Conrady) metodą punktu stałego
// Zwraca znormalizowane współrzędne promienia (out_x, out_y) na płaszczyźnie rzutowej z=1
inline void undistort_point(double px, double py, const CameraIntrinsics& intr,
                            double& out_x, double& out_y, int iters = 8) {
    double xd = (px - intr.cx) / intr.fx;
    double yd = (py - intr.cy) / intr.fy;
    double x = xd, y = yd;
    for (int iter = 0; iter < iters; ++iter) {
        double r2 = x * x + y * y;
        double rad = 1.0 + intr.k1 * r2 + intr.k2 * r2 * r2 + intr.k3 * r2 * r2 * r2;
        x = (xd - (2.0 * intr.p1 * x * y + intr.p2 * (r2 + 2.0 * x * x))) / rad;
        y = (yd - (intr.p1 * (r2 + 2.0 * y * y) + 2.0 * intr.p2 * x * y)) / rad;
    }
    out_x = x;
    out_y = y;
}

} // namespace evimo2
