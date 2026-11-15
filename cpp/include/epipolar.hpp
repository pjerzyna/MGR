// epipolar.hpp: Geometria epipolarna i weryfikacja Focus of Expansion (FOE) dla detektora IMO v2
#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

#include "camera.hpp"
#include "detector.hpp"

namespace evimo2 {

// Lekki wektor 2D dla płaszczyzny matrycy
struct Vec2 {
    double x = 0.0, y = 0.0;

    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_) {}

    Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(double s) const { return {x * s, y * s}; }
    Vec2 operator/(double s) const {
        double inv = (std::abs(s) > 1e-15) ? (1.0 / s) : 0.0;
        return {x * inv, y * inv};
    }

    double dot(const Vec2& o) const { return x * o.x + y * o.y; }
    double norm() const { return std::hypot(x, y); }
};

// Parametry filtru epipolarnego FOE
struct EpipolarParams {
    double min_speed = 0.05;         // Min. prędkość translacji kamery [m/s] do aktywacji filtru
    double angle_thresh_deg = 35.0;  // Próg kątowy w stopniach dla stożka paralaksy epipolarnej
    int max_shift = 8;               // Maksymalny promień przeszukiwania przesunięcia [px]
    int margin = 4;                  // Margines rozszerzenia bounding boxa klastra [px]
    bool enabled = true;             // Czy filtr epipolarny jest włączony
};

// Wynik klasyfikacji epipolarnej per klaster
struct ClusterEpipolarInfo {
    int cluster_id = 0;              // ID klastra (1..K)
    int dx = 0;                      // Przesunięcie w poziomie (wczesny -> późny) [px]
    int dy = 0;                      // Przesunięcie w pionie (wczesny -> późny) [px]
    double shift_len = 0.0;          // Długość wektora przesunięcia [px]
    float corr = 0.0f;               // Maksymalna znormalizowana korelacja wzajemna (NCC)
    double angle_deg = 0.0;          // Kąt odchylenia od lokalnego wektora epipolarnego [deg]
    bool is_imo = false;             // True = potwierdzony obiekt IMO, False = odrzucona paralaksa tła
};

// Wyznacza teoretyczny kierunek wektora pozornego ruchu tła na matrycy w punkcie (u, v)
// na podstawie prędkości liniowej kamery v_cam = (vx, vy, vz) w układzie kamery
inline std::pair<Vec2, double> epipolar_direction(
    double u, double v,
    const Vec3& v_cam,
    const CameraIntrinsics& intr
) {
    double x_norm = (u - intr.cx) / intr.fx;
    double y_norm = (v - intr.cy) / intr.fy;

    double u_dot = intr.fx * (-v_cam.x + x_norm * v_cam.z);
    double v_dot = intr.fy * (-v_cam.y + y_norm * v_cam.z);

    double norm = std::hypot(u_dot, v_dot);
    if (norm < 1e-7) {
        return {{0.0, 0.0}, 0.0};
    }
    return {{u_dot / norm, v_dot / norm}, norm};
}

// Wyznacza optymalne dyskretne przesunięcie (dx, dy) między wczesnym obrazem Es a późnym Ls
// wewnątrz wycinka klastra metodą maksymalizacji znormalizowanej korelacji wzajemnej (NCC)
inline std::pair<std::pair<int, int>, float> estimate_cluster_shift(
    const float* Es,
    const float* Ls,
    int x0, int x1, int y0, int y1,
    int W = DEFAULT_W, int H = DEFAULT_H,
    int max_shift = 8
) {
    x0 = std::clamp(x0, 0, W);
    x1 = std::clamp(x1, 0, W);
    y0 = std::clamp(y0, 0, H);
    y1 = std::clamp(y1, 0, H);

    if (x1 <= x0 || y1 <= y0) {
        return {{0, 0}, 0.0f};
    }

    // Obliczenie energii wycinka E
    float norm_E_sq = 0.0f;
    int bw = x1 - x0;
    int bh = y1 - y0;
    if (bw <= 0 || bh <= 0) return {{0, 0}, 0.0f};

    // Subsampling dla bardzo dużych klastrów (>64 px) dla zachowania ultra-wysokiej wydajności
    int step = 1;
    if (bw > 64 || bh > 64) step = 2;

    for (int y = y0; y < y1; y += step) {
        const float* __restrict__ pE = Es + y * W;
        for (int x = x0; x < x1; x += step) {
            float v = pE[x];
            norm_E_sq += v * v;
        }
    }
    float norm_E = std::sqrt(norm_E_sq);
    if (norm_E < 1e-4f) {
        return {{0, 0}, 0.0f};
    }

    float best_corr = -1.0f;
    std::pair<int, int> best_shift = {0, 0};

    for (int dy = -max_shift; dy <= max_shift; ++dy) {
        int y0_l = y0 + dy;
        int y1_l = y1 + dy;
        if (y0_l < 0 || y1_l > H) continue;

        for (int dx = -max_shift; dx <= max_shift; ++dx) {
            int x0_l = x0 + dx;
            int x1_l = x1 + dx;
            if (x0_l < 0 || x1_l > W) continue;

            float dot = 0.0f;
            float norm_L_sq = 0.0f;

            for (int y = y0; y < y1; y += step) {
                const float* __restrict__ pE = Es + y * W;
                const float* __restrict__ pL = Ls + (y + dy) * W + dx;
                for (int x = x0; x < x1; x += step) {
                    float ve = pE[x];
                    float vl = pL[x];
                    dot += ve * vl;
                    norm_L_sq += vl * vl;
                }
            }

            float norm_L = std::sqrt(norm_L_sq);
            if (norm_L < 1e-4f) continue;

            float corr = dot / (norm_E * norm_L);
            if (corr > best_corr) {
                best_corr = corr;
                best_shift = {dx, dy};
            }
        }
    }

    return {best_shift, best_corr};
}

// Główny filtr epipolarny FOE dla klastrów
// Weryfikuje każdy klaster z DetectResult pod kątem zgodności z ruchem paralaksy tła
inline std::vector<ClusterEpipolarInfo> filter_clusters_epipolar(
    DetectResult& result,
    const float* Es,
    const float* Ls,
    const Vec3& v_cam,
    const CameraIntrinsics& intr,
    const EpipolarParams& params = EpipolarParams(),
    int W = DEFAULT_W, int H = DEFAULT_H
) {
    std::vector<ClusterEpipolarInfo> cluster_info;
    if (result.n_clusters == 0) {
        return cluster_info;
    }

    double speed = v_cam.norm();
    cluster_info.resize(result.clusters.size());

    // Jeśli filtr jest wyłączony lub kamera porusza się zbyt wolno (czysty obrót / spoczynek):
    // akceptujemy wszystkie decyzje z detektora v1
    if (!params.enabled || speed < params.min_speed) {
        for (size_t k = 0; k < result.clusters.size(); ++k) {
            cluster_info[k].cluster_id = result.clusters[k].id;
            cluster_info[k].is_imo = true;
        }
        return cluster_info;
    }

    for (size_t k = 0; k < result.clusters.size(); ++k) {
        const auto& cl = result.clusters[k];
        auto& ci = cluster_info[k];
        ci.cluster_id = cl.id;

        int x0 = std::max(0, cl.x_min - params.margin);
        int x1 = std::min(W, cl.x_max + params.margin + 1);
        int y0 = std::max(0, cl.y_min - params.margin);
        int y1 = std::min(H, cl.y_max + params.margin + 1);

        auto [shift, corr] = estimate_cluster_shift(Es, Ls, x0, x1, y0, y1, W, H, params.max_shift);
        ci.dx = shift.first;
        ci.dy = shift.second;
        ci.shift_len = std::hypot(static_cast<double>(ci.dx), static_cast<double>(ci.dy));
        ci.corr = corr;

        double xm = 0.5 * (x0 + x1);
        double ym = 0.5 * (y0 + y1);
        auto [u_epi, epi_norm] = epipolar_direction(xm, ym, v_cam, intr);

        if (epi_norm > 1e-4 && ci.shift_len > 0) {
            Vec2 u_s{static_cast<double>(ci.dx) / ci.shift_len, static_cast<double>(ci.dy) / ci.shift_len};
            double cos_th = std::clamp(u_s.dot(u_epi), -1.0, 1.0);
            ci.angle_deg = std::acos(cos_th) * 180.0 / M_PI;

            // Warunek IMO:
            // Ruch odchyla się od linii epipolarnej o co najmniej angle_thresh_deg LUB
            // jest skierowany przeciwnie do wektora epipolarnego (cos_th < 0)
            ci.is_imo = (ci.angle_deg >= params.angle_thresh_deg);
        } else {
            // Brak wyrazistego przesunięcia lub osobliwość FOE - zachowujemy decyzję v1
            ci.angle_deg = 0.0;
            ci.is_imo = true;
        }
    }

    // Aktualizacja etykiet zdarzeń i listy klastrów w DetectResult:
    // Klastry odrzucone jako paralaksa tła zostają usunięte, a zdarzenia oznaczone jako tło
    std::vector<int> remap(result.clusters.size() + 1, 0);
    std::vector<ClusterInfo> kept_clusters;
    int new_id = 1;
    for (size_t k = 0; k < cluster_info.size(); ++k) {
        if (cluster_info[k].is_imo) {
            remap[cluster_info[k].cluster_id] = new_id;
            ClusterInfo cl = result.clusters[k];
            cl.id = new_id;
            kept_clusters.push_back(cl);
            new_id++;
        }
    }
    result.clusters = std::move(kept_clusters);
    result.n_clusters = static_cast<int>(result.clusters.size());

    size_t n_events = result.is_imo.size();
    for (size_t i = 0; i < n_events; ++i) {
        if (result.is_imo[i]) {
            int cid = result.event_cluster[i];
            if (cid > 0 && cid < static_cast<int>(remap.size()) && remap[cid] > 0) {
                result.event_cluster[i] = remap[cid];
            } else {
                result.event_cluster[i] = 0;
                result.is_imo[i] = false;
            }
        } else {
            result.event_cluster[i] = 0;
        }
    }

    return cluster_info;
}

} // namespace evimo2
