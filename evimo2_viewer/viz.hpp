// Simple event representations + image dumping (PPM/PGM, no external libs).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "events.hpp"

// Signed polarity count image: +1 per ON event, -1 per OFF event.
inline std::vector<int16_t> accumulate_signed(const Event* ev, size_t n, int w, int h) {
    std::vector<int16_t> img(static_cast<size_t>(w) * h, 0);
    for (size_t i = 0; i < n; ++i) {
        if (ev[i].x >= w || ev[i].y >= h) continue;
        img[static_cast<size_t>(ev[i].y) * w + ev[i].x] += ev[i].p ? 1 : -1;
    }
    return img;
}

// Exponential time surface in [0,1]: exp(-(t_ref - t_last)/tau), 0 where no event.
inline std::vector<float> time_surface(const Event* ev, size_t n, int64_t t_ref_us, double tau_us,
                                       int w, int h) {
    const int64_t none = std::numeric_limits<int64_t>::min();
    std::vector<int64_t> last(static_cast<size_t>(w) * h, none);
    for (size_t i = 0; i < n; ++i) {
        if (ev[i].x >= w || ev[i].y >= h) continue;
        last[static_cast<size_t>(ev[i].y) * w + ev[i].x] = ev[i].t_us;
    }
    std::vector<float> ts(last.size(), 0.f);
    for (size_t i = 0; i < last.size(); ++i)
        if (last[i] != none) ts[i] = static_cast<float>(std::exp(-double(t_ref_us - last[i]) / tau_us));
    return ts;
}

// ON = red, OFF = blue, no events = white.
inline void write_ppm_polarity(const std::string& path, const std::vector<int16_t>& img, int w, int h) {
    std::ofstream f(path, std::ios::binary);
    f << "P6\n" << w << " " << h << "\n255\n";
    for (int16_t v : img) {
        unsigned char px[3] = {255, 255, 255};
        if (v > 0) { px[1] = px[2] = static_cast<unsigned char>(std::max(0, 255 - 80 * v)); }
        else if (v < 0) { px[0] = px[1] = static_cast<unsigned char>(std::max(0, 255 + 80 * v)); }
        f.write(reinterpret_cast<char*>(px), 3);
    }
}

inline void write_pgm(const std::string& path, const std::vector<float>& img, int w, int h) {
    std::ofstream f(path, std::ios::binary);
    f << "P5\n" << w << " " << h << "\n255\n";
    for (float v : img) {
        unsigned char b = static_cast<unsigned char>(std::clamp(v, 0.f, 1.f) * 255.f);
        f.write(reinterpret_cast<char*>(&b), 1);
    }
}
