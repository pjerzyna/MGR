// rotcomp.hpp: Motion compensation (SO(3), Rodrigues integration, warp, IWE accumulation)
#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

#include "camera.hpp"
#include "calibration.hpp"
#include <cstdint>

namespace evimo2 {

constexpr int DEFAULT_W = 640;
constexpr int DEFAULT_H = 480;

struct WarpedEvent {
    double t = 0.0;          // Czas zdarzenia [s]
    float px = 0.0f;         // Oryginalna współrzędna x [px]
    float py = 0.0f;         // Oryginalna współrzędna y [px]
    Vec3 bearing;            // Promień bez dystorsji: (bx, by, 1.0)
    uint8_t p = 0;           // Polaryzacja (0 lub 1)
};

struct EventWindow {
    double tc = 0.0;                     // Środek okna czasowego [s]
    double wlen = 0.0;                   // Szerokość okna [s]
    std::vector<WarpedEvent> events;     // Zdarzenia w oknie
    std::vector<double> grid;            // Siatka całkowania [s]
    int kref = 0;                        // Indeks chwili odniesienia w siatce
    CameraIntrinsics intr;               // Parametry kamery
};

// Konfiguracja adaptacyjnego okienkowania hybrydowego (Bounded Adaptive Windowing)
// Opisana w:
//   * Gallego et al. (IEEE TPAMI 2022) - 'Event-based Vision: A Survey',
//   * Gehrig et al. (IJCV 2020) - 'EKLT: Asynchronous Photometric Feature Tracking using Events and Frames',
//   * Rosinol et al. (IEEE RA-L 2018) - 'Ultimate SLAM? Combining Events, Images, and IMU'.
struct AdaptiveWindowConfig {
    double dt_min = 0.010;   // Minimalny czas okna [s] (np. 10 ms -> max 100 FPS)
    double dt_max = 0.040;   // Maksymalny czas okna / timeout [s] (np. 40 ms -> min 25 FPS, ochrona IMU)
    size_t n_target = 60000; // Docelowa liczba zdarzeń w oknie (stała złożoność O(N))
    size_t n_min = 5000;     // Minimalna liczba zdarzeń do uznania okna za poprawne
};

struct WindowSlice {
    size_t i0 = 0;                   // Początkowy indeks zdarzenia
    size_t i1 = 0;                   // Końcowy indeks zdarzenia (exclusive)
    double t_start = 0.0;            // Czas początku okna [s]
    double t_end = 0.0;              // Czas końca okna [s]
    double tc = 0.0;                 // Środek czasowy okna [s]
    double wlen = 0.0;               // Długość okna wlen = t_end - t_start [s]
    size_t count = 0;                // Liczba zdarzeń w oknie
    bool triggered_by_count = false; // true = N >= n_target; false = dt >= dt_max (timeout)
};

// Asynchroniczny generator okien adaptacyjnych sterowany podwójnym wyzwalaczem (Dual-Trigger)
template <typename TimeIterator>
inline std::vector<WindowSlice> generate_adaptive_windows(
    TimeIterator t_begin, TimeIterator t_end,
    double t_start_limit, double t_end_limit,
    const AdaptiveWindowConfig& cfg, int max_windows = -1)
{
    std::vector<WindowSlice> slices;
    size_t total = std::distance(t_begin, t_end);
    if (total == 0) return slices;

    auto it_start = std::lower_bound(t_begin, t_end, t_start_limit);
    size_t i0 = std::distance(t_begin, it_start);

    while (i0 < total) {
        double t0 = *(t_begin + i0);
        if (t0 >= t_end_limit) break;

        size_t i1 = i0;
        bool triggered_by_count = false;
        bool triggered = false;

        while (i1 < total) {
            double ti = *(t_begin + i1);
            if (ti > t_end_limit) break;

            double dt = ti - t0;
            size_t n = i1 - i0;

            // Warunek A: Osiągnięto N_target (po upływie co najmniej dt_min)
            if (n >= cfg.n_target && dt >= cfg.dt_min) {
                triggered_by_count = true;
                triggered = true;
                break;
            }
            // Warunek B: Osiągnięto timeout dt_max
            if (dt >= cfg.dt_max) {
                if (n >= cfg.n_min) {
                    triggered = true;
                }
                break;
            }
            ++i1;
        }

        size_t n = i1 - i0;
        if (triggered && n >= cfg.n_min) {
            WindowSlice s;
            s.i0 = i0;
            s.i1 = i1;
            s.t_start = t0;
            s.t_end = *(t_begin + (i1 < total ? i1 : total - 1));
            s.wlen = std::max(s.t_end - s.t_start, 1e-4);
            s.tc = 0.5 * (s.t_start + s.t_end);
            s.count = n;
            s.triggered_by_count = triggered_by_count;
            slices.push_back(s);

            if (max_windows > 0 && static_cast<int>(slices.size()) >= max_windows) {
                break;
            }
            i0 = i1; // Asynchroniczne okna ciągłe rozłączne
        } else {
            i0 += std::max<size_t>(1000, cfg.n_min / 2);
        }
    }
    return slices;
}

// Tradycyjny generator okien o stałej długości czasowej
template <typename TimeIterator>
inline std::vector<WindowSlice> generate_fixed_windows(
    TimeIterator t_begin, TimeIterator t_end,
    double t_start_limit, double t_end_limit,
    double wlen, double step, int max_windows = -1, size_t n_min = 5000)
{
    std::vector<WindowSlice> slices;
    size_t total = std::distance(t_begin, t_end);
    if (total == 0) return slices;

    for (double tc = t_start_limit; tc <= t_end_limit; tc += step) {
        double t_lo = tc - wlen * 0.5;
        double t_hi = tc + wlen * 0.5;

        auto it0 = std::lower_bound(t_begin, t_end, t_lo);
        auto it1 = std::lower_bound(t_begin, t_end, t_hi);

        size_t idx0 = std::distance(t_begin, it0);
        size_t idx1 = std::distance(t_begin, it1);
        size_t n = idx1 > idx0 ? (idx1 - idx0) : 0;

        if (n >= n_min) {
            WindowSlice s;
            s.i0 = idx0;
            s.i1 = idx1;
            s.tc = tc;
            s.wlen = wlen;
            s.t_start = t_lo;
            s.t_end = t_hi;
            s.count = n;
            s.triggered_by_count = false;
            slices.push_back(s);
            if (max_windows > 0 && static_cast<int>(slices.size()) >= max_windows) {
                break;
            }
        }
    }
    return slices;
}

// Liniowa interpolacja wektora 3D w zadanym punkcie czasowym t_query
inline Vec3 interp_vec3(double t_query, const double* ts, const Vec3* vals, size_t n) {
    if (n == 0) return Vec3{0, 0, 0};
    if (t_query <= ts[0]) return vals[0];
    if (t_query >= ts[n - 1]) return vals[n - 1];

    const double* it = std::upper_bound(ts, ts + n, t_query);
    size_t idx1 = it - ts;
    size_t idx0 = idx1 - 1;

    double t0 = ts[idx0];
    double t1 = ts[idx1];
    double dt = std::max(t1 - t0, 1e-12);
    double a = std::clamp((t_query - t0) / dt, 0.0, 1.0);
    return vals[idx0] * (1.0 - a) + vals[idx1] * a;
}

// Przygotowanie okna zdarzeń: wycięcie z bufora czasowego, odwrócenie dystorsji i budowa siatki całkowania
template <typename TimeIterator, typename XIterator, typename YIterator, typename PolIterator>
inline EventWindow prep_window(double tc, double wlen,
                               TimeIterator t_begin, TimeIterator t_end,
                               XIterator x_begin, YIterator y_begin, PolIterator p_begin,
                               const CameraIntrinsics& intr,
                               const UndistortLUT* lut = nullptr) {
    EventWindow win;
    win.tc = tc;
    win.wlen = wlen;
    win.intr = intr;

    double t_start = tc - wlen / 2.0;
    double t_stop  = tc + wlen / 2.0;

    auto it0 = std::lower_bound(t_begin, t_end, t_start);
    auto it1 = std::lower_bound(t_begin, t_end, t_stop);

    size_t idx0 = std::distance(t_begin, it0);
    size_t idx1 = std::distance(t_begin, it1);
    size_t count = idx1 > idx0 ? idx1 - idx0 : 0;

    win.events.resize(count);
    if (lut && !lut->empty()) {
        for (size_t i = 0; i < count; ++i) {
            size_t src_idx = idx0 + i;
            double t = *(t_begin + src_idx);
            float x = static_cast<float>(*(x_begin + src_idx));
            float y = static_cast<float>(*(y_begin + src_idx));
            uint8_t p = static_cast<uint8_t>(*(p_begin + src_idx));

            win.events[i].t = t;
            win.events[i].px = x;
            win.events[i].py = y;
            win.events[i].bearing = lut->lookup(static_cast<int>(x), static_cast<int>(y));
            win.events[i].p = p;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            size_t src_idx = idx0 + i;
            double t = *(t_begin + src_idx);
            float x = static_cast<float>(*(x_begin + src_idx));
            float y = static_cast<float>(*(y_begin + src_idx));
            uint8_t p = static_cast<uint8_t>(*(p_begin + src_idx));

            double bx = 0.0, by = 0.0;
            undistort_point(x, y, intr, bx, by, 8);

            win.events[i].t = t;
            win.events[i].px = x;
            win.events[i].py = y;
            win.events[i].bearing = Vec3{bx, by, 1.0};
            win.events[i].p = p;
        }
    }

    int nst = std::max(2, static_cast<int>(std::round(wlen / 0.001)));
    win.grid.resize(nst + 1);
    double dt_step = wlen / static_cast<double>(nst);
    for (int k = 0; k <= nst; ++k) {
        win.grid[k] = t_start + k * dt_step;
    }
    win.kref = nst / 2;

    return win;
}

// Przygotowanie EventWindow na podstawie wcześniej wyznaczonego WindowSlice
template <typename TimeIterator, typename XIterator, typename YIterator, typename PolIterator>
inline EventWindow prep_window_slice(const WindowSlice& slice,
                                     TimeIterator t_begin,
                                     XIterator x_begin, YIterator y_begin, PolIterator p_begin,
                                     const CameraIntrinsics& intr,
                                     const UndistortLUT* lut = nullptr) {
    EventWindow win;
    win.tc = slice.tc;
    win.wlen = slice.wlen;
    win.intr = intr;

    size_t count = slice.count;
    size_t idx0 = slice.i0;

    win.events.resize(count);
    if (lut && !lut->empty()) {
        for (size_t i = 0; i < count; ++i) {
            size_t src_idx = idx0 + i;
            double t = *(t_begin + src_idx);
            float x = static_cast<float>(*(x_begin + src_idx));
            float y = static_cast<float>(*(y_begin + src_idx));
            uint8_t p = static_cast<uint8_t>(*(p_begin + src_idx));

            win.events[i].t = t;
            win.events[i].px = x;
            win.events[i].py = y;
            win.events[i].bearing = lut->lookup(static_cast<int>(x), static_cast<int>(y));
            win.events[i].p = p;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            size_t src_idx = idx0 + i;
            double t = *(t_begin + src_idx);
            float x = static_cast<float>(*(x_begin + src_idx));
            float y = static_cast<float>(*(y_begin + src_idx));
            uint8_t p = static_cast<uint8_t>(*(p_begin + src_idx));

            double bx = 0.0, by = 0.0;
            undistort_point(x, y, intr, bx, by, 8);

            win.events[i].t = t;
            win.events[i].px = x;
            win.events[i].py = y;
            win.events[i].bearing = Vec3{bx, by, 1.0};
            win.events[i].p = p;
        }
    }

    int nst = std::max(2, static_cast<int>(std::round(slice.wlen / 0.001)));
    win.grid.resize(nst + 1);
    double dt_step = slice.wlen / static_cast<double>(nst);
    for (int k = 0; k <= nst; ++k) {
        win.grid[k] = slice.t_start + k * dt_step;
    }
    win.kref = nst / 2;

    return win;
}

// Wyznacza macierze obrotu M_i = R(t_ref)^T * R(t_i) dla każdego zdarzenia w oknie
inline std::vector<Mat3> rot_matrices(const EventWindow& win,
                                      const double* cal_ts, const Vec3* cal_w, size_t n_imu,
                                      double d) {
    size_t nst = win.grid.size() - 1;
    std::vector<Vec3> om(nst);
    std::vector<Mat3> steps(nst);

    for (size_t k = 0; k < nst; ++k) {
        double t_mid = 0.5 * (win.grid[k] + win.grid[k + 1]) + d;
        om[k] = interp_vec3(t_mid, cal_ts, cal_w, n_imu);
        double dt = win.grid[k + 1] - win.grid[k];
        steps[k] = rodrigues(om[k] * dt);
    }

    std::vector<Mat3> Rg(nst + 1);
    Rg[0] = Mat3::identity();
    for (size_t k = 0; k < nst; ++k) {
        Rg[k + 1] = Rg[k] * steps[k];
    }

    Mat3 R_ref_inv = Rg[win.kref].transpose();

    size_t n_ev = win.events.size();
    std::vector<Mat3> M(n_ev);

    for (size_t i = 0; i < n_ev; ++i) {
        double te = win.events[i].t;
        auto it = std::upper_bound(win.grid.begin(), win.grid.end(), te);
        int ki = static_cast<int>(std::distance(win.grid.begin(), it)) - 1;
        ki = std::clamp(ki, 0, static_cast<int>(nst) - 1);

        double sub_dt = te - win.grid[ki];
        Mat3 E = rodrigues(om[ki] * sub_dt);
        Mat3 Ri = Rg[ki] * E;
        M[i] = R_ref_inv * Ri;
    }

    return M;
}

// Rzutowanie (warping) promieni zdarzeń do chwili środkowej okna
// Jeśli d jest NaN, zwracane są współrzędne bez kompensacji rotacji (wyprostowane promienie)
inline void warp(const EventWindow& win,
                 const double* cal_ts, const Vec3* cal_w, size_t n_imu,
                 double d,
                 float* out_xs, float* out_ys) {
    size_t n = win.events.size();
    double fx = win.intr.fx;
    double fy = win.intr.fy;
    double cx = win.intr.cx;
    double cy = win.intr.cy;

    if (std::isnan(d)) {
        for (size_t i = 0; i < n; ++i) {
            const Vec3& b = win.events[i].bearing;
            out_xs[i] = static_cast<float>(fx * b.x + cx);
            out_ys[i] = static_cast<float>(fy * b.y + cy);
        }
        return;
    }

    std::vector<Mat3> M = rot_matrices(win, cal_ts, cal_w, n_imu, d);
    for (size_t i = 0; i < n; ++i) {
        Vec3 b2 = M[i] * win.events[i].bearing;
        double z_clamped = (std::abs(b2.z) > 1e-12) ? b2.z : (b2.z >= 0.0 ? 1e-12 : -1e-12);
        double inv_z = 1.0 / z_clamped;
        out_xs[i] = static_cast<float>(fx * b2.x * inv_z + cx);
        out_ys[i] = static_cast<float>(fy * b2.y * inv_z + cy);
    }
}

// Akumulacja biliniowa zdarzeń do obrazu (Image of Warped Events - IWE)
inline void iwe(const float* xs, const float* ys, size_t n,
                float* img, int W = DEFAULT_W, int H = DEFAULT_H) {
    std::fill(img, img + W * H, 0.0f);
    for (size_t i = 0; i < n; ++i) {
        float x = xs[i];
        float y = ys[i];
        int x0 = static_cast<int>(std::floor(x));
        int y0 = static_cast<int>(std::floor(y));

        if (x0 < 0 || x0 >= W - 1 || y0 < 0 || y0 >= H - 1) {
            continue;
        }

        float ax = x - static_cast<float>(x0);
        float ay = y - static_cast<float>(y0);

        float w00 = (1.0f - ax) * (1.0f - ay);
        float w01 = ax * (1.0f - ay);
        float w10 = (1.0f - ax) * ay;
        float w11 = ax * ay;

        size_t idx = static_cast<size_t>(y0) * W + x0;
        img[idx]         += w00;
        img[idx + 1]     += w01;
        img[idx + W]     += w10;
        img[idx + W + 1] += w11;
    }
}

} // namespace evimo2
