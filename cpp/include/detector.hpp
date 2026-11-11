// detector.hpp: IMO detection pipeline (residual map, box blur via SAT, CCL, morphological ops)
#pragma once

#include <cmath>
#include <vector>
#include <deque>
#include <algorithm>
#include <numeric>

#include "rotcomp.hpp"

namespace evimo2 {

constexpr int CELL_SIZE = 4;

struct DetectParams {
    int blur = 3;             // Promień filtru pudełkowego (2 przebiegi)
    float tau = 0.5f;         // Próg mapy niespójności D
    float support = 1.0f;     // Współczynnik progu gęstości S (support * mean(S))
    int min_cells = 8;        // Minimalna liczba komórek 4x4 w klastrze
    bool fill = false;        // Wypełnianie dziur wewnątrz klastrów
    int close = 0;            // Domknięcie morfologiczne siatki komórek (promień komórek)
};

struct ClusterInfo {
    int id = 0;
    int x_min = 100000;
    int y_min = 100000;
    int x_max = -1;
    int y_max = -1;
    int event_count = 0;
    int cell_count = 0;
};

struct DetectResult {
    std::vector<bool> is_imo;           // Etykieta IMO per zdarzenie (N)
    std::vector<int> event_cluster;     // ID klastra (0 = tło, 1..K = klastry) per zdarzenie (N)
    int n_clusters = 0;                 // Liczba zachowanych klastrów
    std::vector<ClusterInfo> clusters;  // Metadane klastrów wraz z Bounding Boxami
};

// Szybki filtr pudełkowy 2D przez tablicę sum skumulowanych (SAT / Integral Image)
// z replikacją krawędzi (odpowiednik mode="edge" w NumPy)
inline void box_blur(const float* img, float* out, int W, int H, int r) {
    int k = 2 * r + 1;
    int W_pad = W + 2 * r;
    int H_pad = H + 2 * r;

    // 1. Tablica sum integralnych SAT o rozmiarze (H_pad + 1) x (W_pad + 1)
    std::vector<double> sat((H_pad + 1) * (W_pad + 1), 0.0);
    size_t sat_stride = W_pad + 1;

    for (int y = 0; y < H_pad; ++y) {
        int orig_y = std::clamp(y - r, 0, H - 1);
        double row_sum = 0.0;
        for (int x = 0; x < W_pad; ++x) {
            int orig_x = std::clamp(x - r, 0, W - 1);
            row_sum += img[orig_y * W + orig_x];
            sat[(y + 1) * sat_stride + (x + 1)] = sat[y * sat_stride + (x + 1)] + row_sum;
        }
    }

    // 2. Odczyt sum w oknie k x k i normalizacja przez k^2
    double inv_area = 1.0 / (static_cast<double>(k) * k);
    for (int y = 0; y < H; ++y) {
        int y1 = y;
        int y2 = y + k;
        for (int x = 0; x < W; ++x) {
            int x1 = x;
            int x2 = x + k;
            double sum = sat[y2 * sat_stride + x2]
                       - sat[y1 * sat_stride + x2]
                       - sat[y2 * sat_stride + x1]
                       + sat[y1 * sat_stride + x1];
            out[y * W + x] = static_cast<float>(sum * inv_area);
        }
    }
}

// Generuje mapę niespójności D = |E - L| / (E + L + eps) pomiędzy wczesną i późną połówką okna
inline bool residual_map(const float* xs, const float* ys, const double* te,
                         size_t n, double tc, int r,
                         std::vector<float>& out_D, std::vector<float>& out_S,
                         int W = DEFAULT_W, int H = DEFAULT_H, float eps = 1e-12f) {
    out_D.assign(W * H, 0.0f);
    out_S.assign(W * H, 0.0f);

    std::vector<float> early_x, early_y;
    std::vector<float> late_x, late_y;
    early_x.reserve(n); early_y.reserve(n);
    late_x.reserve(n); late_y.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        if (te[i] < tc) {
            early_x.push_back(xs[i]);
            early_y.push_back(ys[i]);
        } else {
            late_x.push_back(xs[i]);
            late_y.push_back(ys[i]);
        }
    }

    if (early_x.empty() || late_x.empty()) {
        return false;
    }

    std::vector<float> E(W * H);
    std::vector<float> L(W * H);
    iwe(early_x.data(), early_y.data(), early_x.size(), E.data(), W, H);
    iwe(late_x.data(), late_y.data(), late_x.size(), L.data(), W, H);

    double sum_E = 0.0;
    double sum_L = 0.0;
    for (int i = 0; i < W * H; ++i) {
        sum_E += E[i];
        sum_L += L[i];
    }

    if (sum_E < 1.0 || sum_L < 1.0) {
        return false;
    }

    float inv_E = static_cast<float>(1.0 / sum_E);
    float inv_L = static_cast<float>(1.0 / sum_L);
    for (int i = 0; i < W * H; ++i) {
        E[i] *= inv_E;
        L[i] *= inv_L;
    }

    // Dwa przebiegi filtru pudełkowego na obu połówkach okna
    std::vector<float> tmp(W * H);
    box_blur(E.data(), tmp.data(), W, H, r);
    box_blur(tmp.data(), E.data(), W, H, r);

    box_blur(L.data(), tmp.data(), W, H, r);
    box_blur(tmp.data(), L.data(), W, H, r);

    for (int i = 0; i < W * H; ++i) {
        float s = E[i] + L[i];
        float d = std::abs(E[i] - L[i]) / (s + eps);
        out_S[i] = s;
        out_D[i] = d;
    }

    return true;
}

// Szybka dylatacja logiczna maski binarnej (8-sąsiedztwo, odpowiednik scipy binary_dilation rank-2)
inline std::vector<bool> dilate_mask(const std::vector<bool>& m, int W, int H, int r) {
    if (r <= 0) return m;
    std::vector<bool> cur = m;
    std::vector<bool> nxt(W * H);

    for (int step = 0; step < r; ++step) {
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                bool val = cur[y * W + x];
                for (int dy = -1; dy <= 1 && !val; ++dy) {
                    for (int dx = -1; dx <= 1 && !val; ++dx) {
                        int yy = y + dy;
                        int xx = x + dx;
                        if (yy >= 0 && yy < H && xx >= 0 && xx < W) {
                            if (cur[yy * W + xx]) {
                                val = true;
                                break;
                            }
                        }
                    }
                }
                nxt[y * W + x] = val;
            }
        }
        cur = nxt;
    }
    return cur;
}

// Domknięcie morfologiczne (dylatacja, a następnie erozja)
inline std::vector<bool> close_mask(const std::vector<bool>& m, int W, int H, int r) {
    if (r <= 0) return m;
    // Dylatacja
    std::vector<bool> d = dilate_mask(m, W, H, r);
    // Erozja: ~dilate(~d, r)
    std::vector<bool> inv_d(W * H);
    for (int i = 0; i < W * H; ++i) inv_d[i] = !d[i];
    std::vector<bool> d_inv = dilate_mask(inv_d, W, H, r);
    std::vector<bool> res(W * H);
    for (int i = 0; i < W * H; ++i) res[i] = !d_inv[i];
    return res;
}

// Składowe spójne (8-sąsiedztwo) maski logicznej (indeks 0 = tło)
inline int label_cc(const std::vector<bool>& g, int W, int H,
                    std::vector<int>& lab, std::vector<int>& areas) {
    lab.assign(W * H, 0);
    areas.clear();
    areas.push_back(0); // Indeks 0 = tło

    int n_comp = 0;
    std::deque<std::pair<int, int>> q;

    for (int y0 = 0; y0 < H; ++y0) {
        for (int x0 = 0; x0 < W; ++x0) {
            int idx0 = y0 * W + x0;
            if (!g[idx0] || lab[idx0] != 0) {
                continue;
            }

            ++n_comp;
            lab[idx0] = n_comp;
            int count = 0;
            q.emplace_back(y0, x0);

            while (!q.empty()) {
                auto [y, x] = q.front();
                q.pop_front();
                ++count;

                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int yy = y + dy;
                        int xx = x + dx;
                        if (yy >= 0 && yy < H && xx >= 0 && xx < W) {
                            int n_idx = yy * W + xx;
                            if (g[n_idx] && lab[n_idx] == 0) {
                                lab[n_idx] = n_comp;
                                q.emplace_back(yy, xx);
                            }
                        }
                    }
                }
            }
            areas.push_back(count);
        }
    }
    return n_comp;
}

// Wypełnianie dziur w masce binarnej (flood-fill 4-sąsiedztwo od krawędzi siatki)
inline std::vector<bool> fill_holes(const std::vector<bool>& g, int W, int H) {
    std::vector<bool> reach(W * H, false);
    std::deque<std::pair<int, int>> q;

    for (int y = 0; y < H; ++y) {
        for (int x : {0, W - 1}) {
            int idx = y * W + x;
            if (!g[idx] && !reach[idx]) {
                reach[idx] = true;
                q.emplace_back(y, x);
            }
        }
    }

    for (int x = 0; x < W; ++x) {
        for (int y : {0, H - 1}) {
            int idx = y * W + x;
            if (!g[idx] && !reach[idx]) {
                reach[idx] = true;
                q.emplace_back(y, x);
            }
        }
    }

    constexpr int dys[4] = {1, -1, 0, 0};
    constexpr int dxs[4] = {0, 0, 1, -1};

    while (!q.empty()) {
        auto [y, x] = q.front();
        q.pop_front();

        for (int i = 0; i < 4; ++i) {
            int yy = y + dys[i];
            int xx = x + dxs[i];
            if (yy >= 0 && yy < H && xx >= 0 && xx < W) {
                int n_idx = yy * W + xx;
                if (!g[n_idx] && !reach[n_idx]) {
                    reach[n_idx] = true;
                    q.emplace_back(yy, xx);
                }
            }
        }
    }

    std::vector<bool> out(W * H);
    for (int i = 0; i < W * H; ++i) {
        out[i] = g[i] || !reach[i];
    }
    return out;
}

// Główna funkcja detekcji IMO v0
inline DetectResult detect_imo(const float* xs, const float* ys,
                               const double* te, size_t n, double tc,
                               const DetectParams& params,
                               int W = DEFAULT_W, int H = DEFAULT_H) {
    DetectResult result;
    result.is_imo.assign(n, false);
    result.event_cluster.assign(n, 0);
    result.n_clusters = 0;

    std::vector<float> D, S;
    if (!residual_map(xs, ys, te, n, tc, params.blur, D, S, W, H)) {
        return result;
    }

    double s_mean = 0.0;
    for (int i = 0; i < W * H; ++i) s_mean += S[i];
    s_mean /= static_cast<double>(W * H);

    float s_thresh = static_cast<float>(params.support * s_mean);

    // Redukcja na siatkę komórek 4x4 px (downsampling)
    int Wg = W / CELL_SIZE; // 160
    int Hg = H / CELL_SIZE; // 120
    std::vector<bool> g(Wg * Hg, false);

    for (int cy = 0; cy < Hg; ++cy) {
        for (int cx = 0; cx < Wg; ++cx) {
            int count = 0;
            for (int dy = 0; dy < CELL_SIZE; ++dy) {
                int py = cy * CELL_SIZE + dy;
                for (int dx = 0; dx < CELL_SIZE; ++dx) {
                    int px = cx * CELL_SIZE + dx;
                    int pidx = py * W + px;
                    if (D[pidx] >= params.tau && S[pidx] >= s_thresh) {
                        ++count;
                    }
                }
            }
            // mean >= 0.3 na bloku 16 pikseli oznacza count >= 5
            if (count >= 5) {
                g[cy * Wg + cx] = true;
            }
        }
    }

    if (params.close > 0) {
        g = close_mask(g, Wg, Hg, params.close);
    }

    std::vector<int> lab, areas;
    label_cc(g, Wg, Hg, lab, areas);

    std::vector<bool> keep_label(areas.size(), false);
    std::vector<int> remap_id(areas.size(), 0);
    int kept_count = 0;
    for (size_t k = 1; k < areas.size(); ++k) {
        if (areas[k] >= params.min_cells) {
            keep_label[k] = true;
            ++kept_count;
            remap_id[k] = kept_count;
        }
    }

    if (kept_count == 0) {
        return result;
    }

    result.n_clusters = kept_count;
    result.clusters.resize(kept_count);
    for (int k = 0; k < kept_count; ++k) {
        result.clusters[k].id = k + 1;
    }

    std::vector<int> cell_cluster(Wg * Hg, 0);
    for (int cy = 0; cy < Hg; ++cy) {
        for (int cx = 0; cx < Wg; ++cx) {
            int idx = cy * Wg + cx;
            int l = lab[idx];
            if (l > 0 && keep_label[l]) {
                int cid = remap_id[l];
                cell_cluster[idx] = cid;
                auto& cl = result.clusters[cid - 1];
                cl.cell_count++;
                int px0 = cx * CELL_SIZE;
                int py0 = cy * CELL_SIZE;
                int px1 = std::min(W - 1, (cx + 1) * CELL_SIZE - 1);
                int py1 = std::min(H - 1, (cy + 1) * CELL_SIZE - 1);
                cl.x_min = std::min(cl.x_min, px0);
                cl.y_min = std::min(cl.y_min, py0);
                cl.x_max = std::max(cl.x_max, px1);
                cl.y_max = std::max(cl.y_max, py1);
            }
        }
    }

    // Wypełnianie dziur wewnątrz klastrów (jeśli włączone) - per klaster, aby przypisanie było spójne
    if (params.fill) {
        for (int cid = 1; cid <= kept_count; ++cid) {
            std::vector<bool> cl_mask(Wg * Hg, false);
            for (int i = 0; i < Wg * Hg; ++i) {
                if (cell_cluster[i] == cid) cl_mask[i] = true;
            }
            std::vector<bool> filled = fill_holes(cl_mask, Wg, Hg);
            auto& cl = result.clusters[cid - 1];
            for (int cy = 0; cy < Hg; ++cy) {
                for (int cx = 0; cx < Wg; ++cx) {
                    int idx = cy * Wg + cx;
                    if (filled[idx] && cell_cluster[idx] == 0) {
                        cell_cluster[idx] = cid;
                        cl.cell_count++;
                        int px0 = cx * CELL_SIZE;
                        int py0 = cy * CELL_SIZE;
                        int px1 = std::min(W - 1, (cx + 1) * CELL_SIZE - 1);
                        int py1 = std::min(H - 1, (cy + 1) * CELL_SIZE - 1);
                        cl.x_min = std::min(cl.x_min, px0);
                        cl.y_min = std::min(cl.y_min, py0);
                        cl.x_max = std::max(cl.x_max, px1);
                        cl.y_max = std::max(cl.y_max, py1);
                    }
                }
            }
        }
    }

    // Przypisanie etykiety IMO do pojedynczych zdarzeń
    for (size_t i = 0; i < n; ++i) {
        int cx = static_cast<int>(std::floor(xs[i] / static_cast<float>(CELL_SIZE)));
        int cy = static_cast<int>(std::floor(ys[i] / static_cast<float>(CELL_SIZE)));
        if (cx >= 0 && cx < Wg && cy >= 0 && cy < Hg) {
            int idx = cy * Wg + cx;
            int cid = cell_cluster[idx];
            if (cid > 0) {
                result.is_imo[i] = true;
                result.event_cluster[i] = cid;
                result.clusters[cid - 1].event_count++;
            }
        }
    }

    return result;
}

} // namespace evimo2
