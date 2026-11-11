// viz.hpp: Event representations, IMO cluster rendering, PPM image export, and side-by-side composite
#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "events.hpp"
#include "detector.hpp"

namespace evimo2 {

// Podstawowa struktura koloru RGB 24-bit
struct ColorRGB {
    uint8_t r = 0, g = 0, b = 0;
};

// Paleta wyrazistych barw dla głównych klastrów IMO
inline ColorRGB get_cluster_color(int id) {
    static const ColorRGB palette[] = {
        {255, 59, 48},    // 1: Koralowy czerwony
        {52, 199, 89},    // 2: Jaskrawy zielony
        {0, 199, 230},    // 3: Neonowy cyan
        {255, 149, 0},    // 4: Bursztynowy pomarańcz
        {175, 82, 222},   // 5: Fioletowy / Purpura
        {255, 214, 10},   // 6: Złoty żółty
        {255, 45, 85},    // 7: Różowy
        {48, 209, 88},    // 8: Szmaragdowy
        {100, 210, 255},  // 9: Lodowy błękit
        {218, 165, 32}    // 10: Ciepły brąz
    };
    if (id <= 0) return {255, 59, 48};
    return palette[(id - 1) % 10];
}

// 5x7 ASCII bitmap font (ASCII 32 ' ' do 95 '_')
inline const uint8_t FONT_5X7[64][5] = {
    {0x00, 0x00, 0x00, 0x00, 0x00}, /* 32: " " */
    {0x00, 0x00, 0x5f, 0x00, 0x00}, /* 33: "!" */
    {0x00, 0x07, 0x00, 0x07, 0x00}, /* 34: "\"" */
    {0x14, 0x7f, 0x14, 0x7f, 0x14}, /* 35: "#" */
    {0x24, 0x2a, 0x7f, 0x2a, 0x12}, /* 36: "$" */
    {0x23, 0x13, 0x08, 0x64, 0x62}, /* 37: "%" */
    {0x36, 0x49, 0x55, 0x22, 0x50}, /* 38: "&" */
    {0x00, 0x05, 0x03, 0x00, 0x00}, /* 39: "'" */
    {0x00, 0x1c, 0x22, 0x41, 0x00}, /* 40: "(" */
    {0x00, 0x41, 0x22, 0x1c, 0x00}, /* 41: ")" */
    {0x14, 0x08, 0x3e, 0x08, 0x14}, /* 42: "*" */
    {0x08, 0x08, 0x3e, 0x08, 0x08}, /* 43: "+" */
    {0x00, 0x50, 0x30, 0x00, 0x00}, /* 44: "," */
    {0x08, 0x08, 0x08, 0x08, 0x08}, /* 45: "-" */
    {0x00, 0x60, 0x60, 0x00, 0x00}, /* 46: "." */
    {0x20, 0x10, 0x08, 0x04, 0x02}, /* 47: "/" */
    {0x3e, 0x51, 0x49, 0x45, 0x3e}, /* 48: "0" */
    {0x00, 0x42, 0x7f, 0x40, 0x00}, /* 49: "1" */
    {0x42, 0x61, 0x51, 0x49, 0x46}, /* 50: "2" */
    {0x21, 0x41, 0x45, 0x4b, 0x31}, /* 51: "3" */
    {0x18, 0x14, 0x12, 0x7f, 0x10}, /* 52: "4" */
    {0x27, 0x45, 0x45, 0x45, 0x39}, /* 53: "5" */
    {0x3c, 0x4a, 0x49, 0x49, 0x30}, /* 54: "6" */
    {0x01, 0x71, 0x09, 0x05, 0x03}, /* 55: "7" */
    {0x36, 0x49, 0x49, 0x49, 0x36}, /* 56: "8" */
    {0x06, 0x49, 0x49, 0x29, 0x1e}, /* 57: "9" */
    {0x00, 0x36, 0x36, 0x00, 0x00}, /* 58: ":" */
    {0x00, 0x56, 0x36, 0x00, 0x00}, /* 59: ";" */
    {0x08, 0x14, 0x22, 0x41, 0x00}, /* 60: "<" */
    {0x14, 0x14, 0x14, 0x14, 0x14}, /* 61: "=" */
    {0x00, 0x41, 0x22, 0x14, 0x08}, /* 62: ">" */
    {0x02, 0x01, 0x51, 0x09, 0x06}, /* 63: "?" */
    {0x32, 0x49, 0x79, 0x41, 0x3e}, /* 64: "@" */
    {0x7c, 0x12, 0x11, 0x12, 0x7c}, /* 65: "A" */
    {0x7f, 0x49, 0x49, 0x49, 0x36}, /* 66: "B" */
    {0x3e, 0x41, 0x41, 0x41, 0x22}, /* 67: "C" */
    {0x7f, 0x41, 0x41, 0x22, 0x1c}, /* 68: "D" */
    {0x7f, 0x49, 0x49, 0x49, 0x41}, /* 69: "E" */
    {0x7f, 0x09, 0x09, 0x09, 0x01}, /* 70: "F" */
    {0x3e, 0x41, 0x49, 0x49, 0x7a}, /* 71: "G" */
    {0x7f, 0x08, 0x08, 0x08, 0x7f}, /* 72: "H" */
    {0x00, 0x41, 0x7f, 0x41, 0x00}, /* 73: "I" */
    {0x20, 0x40, 0x41, 0x3f, 0x01}, /* 74: "J" */
    {0x7f, 0x08, 0x14, 0x22, 0x41}, /* 75: "K" */
    {0x7f, 0x40, 0x40, 0x40, 0x40}, /* 76: "L" */
    {0x7f, 0x02, 0x0c, 0x02, 0x7f}, /* 77: "M" */
    {0x7f, 0x04, 0x08, 0x10, 0x7f}, /* 78: "N" */
    {0x3e, 0x41, 0x41, 0x41, 0x3e}, /* 79: "O" */
    {0x7f, 0x09, 0x09, 0x09, 0x06}, /* 80: "P" */
    {0x3e, 0x41, 0x51, 0x21, 0x5e}, /* 81: "Q" */
    {0x7f, 0x09, 0x19, 0x29, 0x46}, /* 82: "R" */
    {0x46, 0x49, 0x49, 0x49, 0x31}, /* 83: "S" */
    {0x01, 0x01, 0x7f, 0x01, 0x01}, /* 84: "T" */
    {0x3f, 0x40, 0x40, 0x40, 0x3f}, /* 85: "U" */
    {0x1f, 0x20, 0x40, 0x20, 0x1f}, /* 86: "V" */
    {0x3f, 0x40, 0x38, 0x40, 0x3f}, /* 87: "W" */
    {0x63, 0x14, 0x08, 0x14, 0x63}, /* 88: "X" */
    {0x07, 0x08, 0x70, 0x08, 0x07}, /* 89: "Y" */
    {0x61, 0x51, 0x49, 0x45, 0x43}, /* 90: "Z" */
    {0x00, 0x7f, 0x41, 0x41, 0x00}, /* 91: "[" */
    {0x02, 0x04, 0x08, 0x10, 0x20}, /* 92: "\\" */
    {0x00, 0x41, 0x41, 0x7f, 0x00}, /* 93: "]" */
    {0x04, 0x02, 0x01, 0x02, 0x04}, /* 94: "^" */
    {0x40, 0x40, 0x40, 0x40, 0x40}, /* 95: "_" */
};

// Bufor obrazu RGB 24-bit z funkcjami rysowania
class ImageRGB {
public:
    int w = 0;
    int h = 0;
    std::vector<uint8_t> data;

    ImageRGB() = default;
    ImageRGB(int width, int height, ColorRGB bg = {0, 0, 0})
        : w(width), h(height), data(static_cast<size_t>(width) * height * 3) {
        fill(bg);
    }

    void fill(ColorRGB c) {
        for (size_t i = 0; i < static_cast<size_t>(w) * h; ++i) {
            data[3 * i + 0] = c.r;
            data[3 * i + 1] = c.g;
            data[3 * i + 2] = c.b;
        }
    }

    inline void set_pixel(int x, int y, ColorRGB c) {
        if (x < 0 || x >= w || y < 0 || y >= h) return;
        size_t idx = (static_cast<size_t>(y) * w + x) * 3;
        data[idx + 0] = c.r;
        data[idx + 1] = c.g;
        data[idx + 2] = c.b;
    }

    inline ColorRGB get_pixel(int x, int y) const {
        if (x < 0 || x >= w || y < 0 || y >= h) return {0, 0, 0};
        size_t idx = (static_cast<size_t>(y) * w + x) * 3;
        return {data[idx + 0], data[idx + 1], data[idx + 2]};
    }

    void fill_rect(int x0, int y0, int x1, int y1, ColorRGB c) {
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);
        x0 = std::clamp(x0, 0, w - 1);
        x1 = std::clamp(x1, 0, w - 1);
        y0 = std::clamp(y0, 0, h - 1);
        y1 = std::clamp(y1, 0, h - 1);

        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                set_pixel(x, y, c);
            }
        }
    }

    void draw_rect(int x0, int y0, int x1, int y1, ColorRGB c, int thickness = 1) {
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);

        for (int t = 0; t < thickness; ++t) {
            int cx0 = std::clamp(x0 - t, 0, w - 1);
            int cx1 = std::clamp(x1 + t, 0, w - 1);
            int cy0 = std::clamp(y0 - t, 0, h - 1);
            int cy1 = std::clamp(y1 + t, 0, h - 1);

            for (int x = cx0; x <= cx1; ++x) {
                set_pixel(x, cy0, c);
                set_pixel(x, cy1, c);
            }
            for (int y = cy0; y <= cy1; ++y) {
                set_pixel(cx0, y, c);
                set_pixel(cx1, y, c);
            }
        }
    }

    void draw_text(int x0, int y0, const std::string& text, ColorRGB color, int scale = 1) {
        int cur_x = x0;
        for (char c : text) {
            char uc = (c >= 'a' && c <= 'z') ? static_cast<char>(c - 'a' + 'A') : c;
            int idx = static_cast<int>(uc) - 32;
            if (idx < 0 || idx >= 64) idx = 0;
            const uint8_t* cols = FONT_5X7[idx];

            for (int col = 0; col < 5; ++col) {
                uint8_t bits = cols[col];
                for (int row = 0; row < 7; ++row) {
                    if (bits & (1 << row)) {
                        for (int sy = 0; sy < scale; ++sy) {
                            for (int sx = 0; sx < scale; ++sx) {
                                set_pixel(cur_x + col * scale + sx, y0 + row * scale + sy, color);
                            }
                        }
                    }
                }
            }
            cur_x += (5 + 1) * scale;
        }
    }

    bool write_ppm(const std::string& path) const {
        std::ofstream f(path, std::ios::binary);
        if (!f.is_open()) return false;
        f << "P6\n" << w << " " << h << "\n255\n";
        f.write(reinterpret_cast<const char*>(data.data()), data.size());
        return true;
    }
};

struct VizOptions {
    std::string theme = "dark";        // "dark", "light"
    std::string color_mode = "unified"; // "unified" (jednolity wyrazisty kolor IMO), "clusters" (różne barwy per klaster)
    std::string boxes = "major";       // "none" (czysta maska), "major" (tylko duże obiekty >= min_box_events), "all" (wszystkie)
    int min_box_events = 2500;         // Min. liczba zdarzeń w klastrze, by narysować ramkę w trybie "major"
    int max_boxes = 5;                 // Maks. liczba rysowanych ramek w trybie "major"
};

// Renderowanie czytelnej nakładki zdarzeń IMO
inline ImageRGB render_imo_frame(
    const float* xs, const float* ys,
    const uint8_t* pol,
    const std::vector<bool>& is_imo,
    const int* event_cluster,
    size_t n,
    const std::vector<ClusterInfo>& clusters,
    const std::string& title,
    const std::string& info_line,
    const VizOptions& opt,
    int W = 640, int H = 480
) {
    bool dark = (opt.theme == "dark");
    ColorRGB bg_col    = dark ? ColorRGB{14, 16, 20}   : ColorRGB{246, 247, 249};
    ColorRGB bar_bg    = dark ? ColorRGB{24, 27, 34}   : ColorRGB{228, 232, 238};
    ColorRGB bar_line  = dark ? ColorRGB{42, 47, 56}   : ColorRGB{200, 205, 215};
    ColorRGB text_col  = dark ? ColorRGB{240, 245, 255} : ColorRGB{25, 30, 40};
    ColorRGB sub_col   = dark ? ColorRGB{150, 160, 175} : ColorRGB{90, 100, 115};

    // Główny wyrazisty kolor dla wykrytego ruchu (IMO)
    ColorRGB imo_unified = dark ? ColorRGB{255, 60, 60} : ColorRGB{220, 25, 30};

    ImageRGB img(W, H, bg_col);

    // 1. Górny pasek statusu / metryk (24 px)
    img.fill_rect(0, 0, W - 1, 23, bar_bg);
    img.draw_rect(0, 23, W - 1, 23, bar_line, 1);

    img.draw_text(8, 6, title, text_col, 1);
    if (!info_line.empty()) {
        int info_w = static_cast<int>(info_line.size()) * 6;
        img.draw_text(std::max(8, W - info_w - 10), 6, info_line, sub_col, 1);
    }

    // 2. Rysowanie zdarzeń tła (ruch własny kamery / ego-motion)
    for (size_t i = 0; i < n; ++i) {
        if (!is_imo.empty() && is_imo[i]) continue;

        int px = static_cast<int>(std::round(xs[i]));
        int py = static_cast<int>(std::round(ys[i]));
        if (px < 0 || px >= W || py < 24 || py >= H) continue;

        ColorRGB c;
        if (dark) {
            if (pol && pol[i] == 1)      c = {60, 85, 120};  // ON: subtelny błękit
            else if (pol && pol[i] == 0) c = {100, 75, 70};  // OFF: subtelny łupek
            else                         c = {75, 80, 90};   // neutralny szary
        } else {
            if (pol && pol[i] == 1)      c = {140, 160, 190};
            else if (pol && pol[i] == 0) c = {190, 155, 145};
            else                         c = {175, 180, 185};
        }
        img.set_pixel(px, py, c);
    }

    // 3. Rysowanie zdarzeń IMO (niezależny ruch obiektów)
    for (size_t i = 0; i < n; ++i) {
        if (is_imo.empty() || !is_imo[i]) continue;

        int px = static_cast<int>(std::round(xs[i]));
        int py = static_cast<int>(std::round(ys[i]));
        if (px < 0 || px >= W || py < 24 || py >= H) continue;

        ColorRGB c = imo_unified;
        if (opt.color_mode == "clusters") {
            int cid = (event_cluster && event_cluster[i] > 0) ? event_cluster[i] : 1;
            c = get_cluster_color(cid);
        }
        img.set_pixel(px, py, c);
    }

    // 4. Rysowanie Bounding Boxów (opcjonalne i inteligentnie filtrowane)
    if (opt.boxes != "none" && !clusters.empty()) {
        std::vector<ClusterInfo> sorted_cl = clusters;
        std::sort(sorted_cl.begin(), sorted_cl.end(), [](const ClusterInfo& a, const ClusterInfo& b) {
            return a.event_count > b.event_count;
        });

        int drawn = 0;
        for (const auto& cl : sorted_cl) {
            if (opt.boxes == "major") {
                if (cl.event_count < opt.min_box_events) break;
                if (drawn >= opt.max_boxes) break;
            }

            ColorRGB cl_col = (opt.color_mode == "clusters") ? get_cluster_color(cl.id) : ColorRGB{255, 220, 40}; // złotożółty
            int bx0 = std::clamp(cl.x_min - 2, 0, W - 1);
            int bx1 = std::clamp(cl.x_max + 2, 0, W - 1);
            int by0 = std::clamp(cl.y_min - 2, 24, H - 1);
            int by1 = std::clamp(cl.y_max + 2, 24, H - 1);

            // Subtelna ramka 1 px (lub 2 px)
            img.draw_rect(bx0, by0, bx1, by1, cl_col, 2);

            // Zwięzły znacznik tekstowy
            std::ostringstream tag_oss;
            tag_oss << "#" << cl.id << " (" << cl.event_count << " EV)";
            std::string tag = tag_oss.str();

            int tag_w = static_cast<int>(tag.size()) * 6 + 4;
            int tag_h = 10;
            int tag_x = bx0;
            int tag_y = (by0 >= 24 + tag_h) ? (by0 - tag_h) : (by0 + 2);
            if (tag_x + tag_w >= W) tag_x = W - 1 - tag_w;

            img.fill_rect(tag_x, tag_y, tag_x + tag_w - 1, tag_y + tag_h - 1, dark ? ColorRGB{20, 22, 26} : ColorRGB{250, 250, 250});
            img.draw_rect(tag_x, tag_y, tag_x + tag_w - 1, tag_y + tag_h - 1, cl_col, 1);
            img.draw_text(tag_x + 2, tag_y + 1, tag, cl_col, 1);

            ++drawn;
        }
    }

    return img;
}

// Składa dwa obrazy o identycznej wysokości obok siebie w jeden panel (np. Sensor | Warped SO(3))
inline ImageRGB make_side_by_side(const ImageRGB& left, const ImageRGB& right, ColorRGB sep_color = {55, 60, 75}) {
    int H = std::max(left.h, right.h);
    int W = left.w + right.w;
    ImageRGB res(W, H);

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < left.w; ++x) {
            res.set_pixel(x, y, left.get_pixel(x, y));
        }
        for (int x = 0; x < right.w; ++x) {
            res.set_pixel(left.w + x, y, right.get_pixel(x, y));
        }
    }

    // Pionowa linia rozdzielająca (2 px)
    res.draw_rect(left.w - 1, 0, left.w, H - 1, sep_color, 1);
    return res;
}

// Funkcje zachowane dla kompatybilności z aplikacją viewer
inline std::vector<int16_t> accumulate_signed(const Event* ev, size_t n, int w, int h) {
    std::vector<int16_t> img(static_cast<size_t>(w) * h, 0);
    for (size_t i = 0; i < n; ++i) {
        if (ev[i].x >= w || ev[i].y >= h) continue;
        img[static_cast<size_t>(ev[i].y) * w + ev[i].x] += ev[i].p ? 1 : -1;
    }
    return img;
}

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

} // namespace evimo2
