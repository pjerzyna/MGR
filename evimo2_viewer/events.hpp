// Event representation + loader for EVIMO2 dataset_events_{xy,t,p}.npy
#pragma once
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "npy.hpp"

struct Event {
    int64_t t_us;   // microseconds since the first event in the file
    uint16_t x, y;
    uint8_t p;      // 1 = ON (brightness increase), 0 = OFF
};

struct LoadOptions {
    double t_begin_s = 0.0;                                   // relative to first event in file
    double t_end_s = std::numeric_limits<double>::infinity(); // relative to first event in file
    bool invert_polarity = false;                             // needed for Samsung DVS Gen3
};

struct FileInfo {
    size_t total_events = 0;
    double t0_abs_s = 0.0;       // absolute timestamp of first event
    double duration_s = 0.0;
    std::string xy_desc, t_desc, p_desc;
};

inline std::vector<Event> load_events(const std::string& dir, const LoadOptions& o,
                                      FileInfo* info = nullptr) {
    namespace fs = std::filesystem;
    NpyArray xy((fs::path(dir) / "dataset_events_xy.npy").string());
    NpyArray t((fs::path(dir) / "dataset_events_t.npy").string());
    NpyArray p((fs::path(dir) / "dataset_events_p.npy").string());

    const size_t n = t.size();
    if (n == 0) throw std::runtime_error("empty event file");
    if (xy.size() != 2 * n || p.size() != n)
        throw std::runtime_error("size mismatch between xy/t/p arrays");

    auto time_at = [&](size_t i) {
        double v = 0;
        t.visit([&](auto* ptr) { v = static_cast<double>(ptr[i]); });
        return v;
    };
    const double t0 = time_at(0);
    auto lower = [&](double target) {   // first index with time >= target (timestamps are sorted)
        size_t lo = 0, hi = n;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (time_at(mid) < target) lo = mid + 1; else hi = mid;
        }
        return lo;
    };

    if (info) {
        info->total_events = n; info->t0_abs_s = t0; info->duration_s = time_at(n - 1) - t0;
        info->xy_desc = xy.describe(); info->t_desc = t.describe(); info->p_desc = p.describe();
    }

    const size_t i0 = lower(t0 + o.t_begin_s);
    const size_t i1 = std::isinf(o.t_end_s) ? n : lower(t0 + o.t_end_s);
    const size_t m = i1 > i0 ? i1 - i0 : 0;
    std::vector<Event> ev(m);

    t.visit([&](auto* ptr) {
        for (size_t k = 0; k < m; ++k)
            ev[k].t_us = std::llround((static_cast<double>(ptr[i0 + k]) - t0) * 1e6);
    });
    xy.visit([&](auto* ptr) {
        for (size_t k = 0; k < m; ++k) {
            ev[k].x = static_cast<uint16_t>(ptr[2 * (i0 + k)]);      // negative -> huge -> caught as out-of-range
            ev[k].y = static_cast<uint16_t>(ptr[2 * (i0 + k) + 1]);
        }
    });
    p.visit([&](auto* ptr) {
        for (size_t k = 0; k < m; ++k) {
            const uint8_t v = ptr[i0 + k] != 0;
            ev[k].p = o.invert_polarity ? !v : v;
        }
    });
    return ev;
}
