// imo_detect: Standalone C++ IMO detection and performance benchmark on EVIMO2 sequences
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "camera.hpp"
#include "calibration.hpp"
#include "rotcomp.hpp"
#include "detector.hpp"
#include "viz.hpp"
#include "npy.hpp"

using namespace evimo2;

struct AppArgs {
    std::string seq_dir;
    double window_ms = 40.0;
    double step_ms = 40.0;
    int max_windows = 40;
    int blur = 3;
    float tau = 0.5f;
    float support = 1.0f;
    int min_cells = 8;
    double lag_ms = -15.7;
    std::string out_dir = "";

    // Wizualizacja
    std::string viz_dir = "";
    std::string viz_mode = "both";   // "both", "sensor", "warped"
    std::string viz_theme = "dark";  // "dark", "light"
    int viz_stride = 1;
    bool viz_only_imo = false;
    bool png = true;
};

static bool parse_args(int argc, char** argv, AppArgs& a) {
    if (argc < 2) return false;
    a.seq_dir = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        auto need = [&](int k) { return i + k < argc; };
        if (s == "--window-ms" && need(1)) a.window_ms = std::atof(argv[++i]);
        else if (s == "--step-ms" && need(1)) a.step_ms = std::atof(argv[++i]);
        else if (s == "--max-windows" && need(1)) a.max_windows = std::atoi(argv[++i]);
        else if (s == "--blur" && need(1)) a.blur = std::atoi(argv[++i]);
        else if (s == "--tau" && need(1)) a.tau = static_cast<float>(std::atof(argv[++i]));
        else if (s == "--support" && need(1)) a.support = static_cast<float>(std::atof(argv[++i]));
        else if (s == "--min-cells" && need(1)) a.min_cells = std::atoi(argv[++i]);
        else if (s == "--lag-ms" && need(1)) a.lag_ms = std::atof(argv[++i]);
        else if (s == "--out" && need(1)) a.out_dir = argv[++i];
        else if (s == "--viz" && need(1)) a.viz_dir = argv[++i];
        else if (s == "--viz-mode" && need(1)) a.viz_mode = argv[++i];
        else if (s == "--viz-theme" && need(1)) a.viz_theme = argv[++i];
        else if (s == "--viz-stride" && need(1)) a.viz_stride = std::atoi(argv[++i]);
        else if (s == "--viz-only-imo") a.viz_only_imo = true;
        else if (s == "--no-png") a.png = false;
        else return false;
    }
    return true;
}

// Prosty parser JSON dla pliku dataset_intrinsics.json
static bool load_intrinsics(const std::string& path, CameraIntrinsics& intr) {
    std::ifstream f(path);
    if (!f.is_open()) return false;
    std::string line;
    auto extract_val = [](const std::string& l) -> double {
        size_t colon = l.find(':');
        if (colon == std::string::npos) return 0.0;
        size_t comma = l.find(',', colon);
        std::string s = l.substr(colon + 1, comma - colon - 1);
        return std::stod(s);
    };
    while (std::getline(f, line)) {
        if (line.find("\"fx\"") != std::string::npos) intr.fx = extract_val(line);
        else if (line.find("\"fy\"") != std::string::npos) intr.fy = extract_val(line);
        else if (line.find("\"cx\"") != std::string::npos) intr.cx = extract_val(line);
        else if (line.find("\"cy\"") != std::string::npos) intr.cy = extract_val(line);
        else if (line.find("\"k1\"") != std::string::npos) intr.k1 = extract_val(line);
        else if (line.find("\"k2\"") != std::string::npos) intr.k2 = extract_val(line);
        else if (line.find("\"k3\"") != std::string::npos) intr.k3 = extract_val(line);
        else if (line.find("\"p1\"") != std::string::npos) intr.p1 = extract_val(line);
        else if (line.find("\"p2\"") != std::string::npos) intr.p2 = extract_val(line);
    }
    return true;
}

int main(int argc, char** argv) {
    AppArgs a;
    if (!parse_args(argc, argv, a)) {
        std::cerr << "Użycie: " << argv[0] << " <katalog_sekwencji> [opcje]\n"
                  << "  --window-ms N   Długość okna w ms (domyślnie 40.0)\n"
                  << "  --step-ms N     Krok przesuwania okna w ms (domyślnie 40.0)\n"
                  << "  --max-windows N Maksymalna liczba okien (domyślnie 40)\n"
                  << "  --blur N        Promień filtru pudełkowego (domyślnie 3)\n"
                  << "  --tau F         Próg niespójności D (domyślnie 0.5)\n"
                  << "  --support F     Mnożnik progu wsparcia S (domyślnie 1.0)\n"
                  << "  --min-cells N   Min. liczba komórek 4x4 w klastrze (domyślnie 8)\n"
                  << "  --lag-ms F      Opóźnienie IMU w ms (domyślnie -15.7)\n"
                  << "  --out DIR       Katalog wyjściowy wyników numerycznych\n"
                  << "  --viz DIR       Katalog zapisu wizualizacji (włącza rendering klatek)\n"
                  << "  --viz-mode M    Tryb: both (1280x480, domyślny), sensor (640x480), warped (640x480)\n"
                  << "  --viz-theme T   Motyw: dark (neonowy, domyślny), light (papierowy)\n"
                  << "  --viz-stride N  Zapisuj co N-te okno czasowe (domyślnie 1)\n"
                  << "  --viz-only-imo  Zapisuj tylko okna, w których wykryto klastry IMO\n"
                  << "  --no-png        Pomiń automatyczną konwersję ImageMagick PPM do PNG\n";
        return 1;
    }

    namespace fs = std::filesystem;
    fs::path seq_p(a.seq_dir);

    // 1. Intrinsics
    CameraIntrinsics intr;
    fs::path intr_p = seq_p / "dataset_intrinsics.json";
    if (!load_intrinsics(intr_p.string(), intr)) {
        std::cerr << "Błąd: Nie znaleziono " << intr_p << ". Uruchom python/scripts/prepare_cpp_data.py\n";
        return 2;
    }

    // 2. Wczytanie zdarzeń przez mmap
    NpyArray npy_t((seq_p / "dataset_events_t.npy").string());
    NpyArray npy_xy((seq_p / "dataset_events_xy.npy").string());
    NpyArray npy_p((seq_p / "dataset_events_p.npy").string());

    size_t total_events = npy_t.size();
    const float* ev_t_f = (npy_t.itemsize() == 4) ? static_cast<const float*>(npy_t.data()) : nullptr;
    const double* ev_t_d = (npy_t.itemsize() == 8) ? static_cast<const double*>(npy_t.data()) : nullptr;

    auto get_time = [&](size_t idx) -> double {
        return ev_t_f ? static_cast<double>(ev_t_f[idx]) : ev_t_d[idx];
    };

    const uint16_t* ev_xy = static_cast<const uint16_t*>(npy_xy.data());
    const uint8_t* ev_p = static_cast<const uint8_t*>(npy_p.data());

    // 3. Wczytanie IMU
    fs::path imu_p = seq_p / "dataset_imu.npy";
    if (!fs::exists(imu_p)) {
        std::cerr << "Błąd: Brak pliku " << imu_p << ". Uruchom python/scripts/prepare_cpp_data.py\n";
        return 3;
    }
    NpyArray npy_imu(imu_p.string());
    size_t n_imu = npy_imu.shape()[0];
    const double* imu_raw = static_cast<const double*>(npy_imu.data());

    std::vector<double> imu_ts(n_imu);
    std::vector<Vec3> cal_w(n_imu);
    Mat3 R = default_R_imu2cam();
    Vec3 bias = default_gyro_bias();

    for (size_t i = 0; i < n_imu; ++i) {
        imu_ts[i] = imu_raw[i * 7 + 0];
        Vec3 raw_gyro{imu_raw[i * 7 + 1], imu_raw[i * 7 + 2], imu_raw[i * 7 + 3]};
        cal_w[i] = calibrated_gyro(raw_gyro, R, bias);
    }

    double wlen = a.window_ms * 1e-3;
    double step = a.step_ms * 1e-3;
    double lag_s = a.lag_ms * 1e-3;

    double t_first_ev = get_time(0);
    double t_last_ev  = get_time(total_events - 1);
    double t_start = std::max(t_first_ev, imu_ts[0]) + wlen / 2.0 + 0.06;
    double t_end   = std::min(t_last_ev, imu_ts[n_imu - 1]) - wlen / 2.0 - 0.06;

    std::vector<double> window_centers;
    for (double tc = t_start; tc <= t_end; tc += step) {
        window_centers.push_back(tc);
        if (static_cast<int>(window_centers.size()) >= a.max_windows) break;
    }

    std::cout << "================================================================================\n";
    std::cout << "  EVIMO2 IMO DETECTOR v0 (C++ Engine) | Sekwencja: " << seq_p.filename().string() << "\n";
    std::cout << "  Zdarzeń w pliku : " << total_events << " (" << std::fixed << std::setprecision(2)
              << (t_last_ev - t_first_ev) << " s)\n";
    std::cout << "  Próbek IMU      : " << n_imu << " | Opóźnienie d = " << a.lag_ms << " ms\n";
    std::cout << "  Konfiguracja    : okno " << a.window_ms << " ms, krok " << a.step_ms
              << " ms, max okien: " << window_centers.size() << "\n";
    std::cout << "  Parametry       : blur " << a.blur << " (2x), tau " << a.tau
              << ", support " << a.support << ", min_cells " << a.min_cells << "\n";
    std::cout << "================================================================================\n";
    std::cout << std::setw(6)  << "Okno"
              << std::setw(12) << "t_c [s]"
              << std::setw(12) << "Zdarzeń"
              << std::setw(10) << "Klastry"
              << std::setw(12) << "Zdarz. IMO"
              << std::setw(12) << "Undistort"
              << std::setw(10) << "Warp"
              << std::setw(10) << "Detect"
              << std::setw(12) << "Suma [ms]" << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";

    DetectParams det_params;
    det_params.blur = a.blur;
    det_params.tau = a.tau;
    det_params.support = a.support;
    det_params.min_cells = a.min_cells;

    double total_time_ms = 0.0;
    double total_undist_ms = 0.0;
    double total_warp_ms = 0.0;
    double total_det_ms = 0.0;
    size_t processed_events = 0;
    int processed_windows = 0;

    // Bufory wielokrotnego użytku per okno
    std::vector<float> xs, ys;
    std::vector<double> te;

    for (size_t win_idx = 0; win_idx < window_centers.size(); ++win_idx) {
        double tc = window_centers[win_idx];

        // 1. Undistort & prep_window
        auto t0 = std::chrono::high_resolution_clock::now();

        double t_lo = tc - wlen / 2.0;
        double t_hi = tc + wlen / 2.0;

        size_t i0 = 0, i1 = 0;
        if (ev_t_f) {
            float lo_f = static_cast<float>(t_lo);
            float hi_f = static_cast<float>(t_hi);
            i0 = std::lower_bound(ev_t_f, ev_t_f + total_events, lo_f) - ev_t_f;
            i1 = std::lower_bound(ev_t_f, ev_t_f + total_events, hi_f) - ev_t_f;
        } else {
            i0 = std::lower_bound(ev_t_d, ev_t_d + total_events, t_lo) - ev_t_d;
            i1 = std::lower_bound(ev_t_d, ev_t_d + total_events, t_hi) - ev_t_d;
        }

        size_t n = (i1 > i0) ? (i1 - i0) : 0;
        if (n < 5000) continue; // Pomiń zbyt małe okna

        EventWindow win;
        win.tc = tc;
        win.wlen = wlen;
        win.intr = intr;
        win.events.resize(n);

        for (size_t i = 0; i < n; ++i) {
            size_t src = i0 + i;
            double t = get_time(src);
            float x = static_cast<float>(ev_xy[2 * src]);
            float y = static_cast<float>(ev_xy[2 * src + 1]);
            double bx = 0.0, by = 0.0;
            undistort_point(x, y, intr, bx, by, 8);

            win.events[i].t = t;
            win.events[i].px = x;
            win.events[i].py = y;
            win.events[i].bearing = Vec3{bx, by, 1.0};
            win.events[i].p = ev_p[src];
        }

        int nst = std::max(2, static_cast<int>(std::round(wlen / 0.001)));
        win.grid.resize(nst + 1);
        double dt_step = wlen / static_cast<double>(nst);
        for (int k = 0; k <= nst; ++k) win.grid[k] = t_lo + k * dt_step;
        win.kref = nst / 2;

        auto t1 = std::chrono::high_resolution_clock::now();

        // 2. Warp
        xs.resize(n);
        ys.resize(n);
        te.resize(n);
        for (size_t i = 0; i < n; ++i) te[i] = win.events[i].t;

        warp(win, imu_ts.data(), cal_w.data(), n_imu, lag_s, xs.data(), ys.data());

        auto t2 = std::chrono::high_resolution_clock::now();

        // 3. Detect IMO
        DetectResult res = detect_imo(xs.data(), ys.data(), te.data(), n, tc, det_params);

        auto t3 = std::chrono::high_resolution_clock::now();

        double d_undist = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double d_warp   = std::chrono::duration<double, std::milli>(t2 - t1).count();
        double d_det    = std::chrono::duration<double, std::milli>(t3 - t2).count();
        double d_sum    = d_undist + d_warp + d_det;

        total_undist_ms += d_undist;
        total_warp_ms   += d_warp;
        total_det_ms    += d_det;
        total_time_ms   += d_sum;
        processed_events += n;
        ++processed_windows;

        size_t imo_count = 0;
        for (bool b : res.is_imo) imo_count += b;

        std::cout << std::setw(6)  << processed_windows
                  << std::setw(12) << std::fixed << std::setprecision(4) << tc
                  << std::setw(12) << n
                  << std::setw(10) << res.n_clusters
                  << std::setw(12) << imo_count
                  << std::setw(11) << std::setprecision(2) << d_undist << " ms"
                  << std::setw(9)  << std::setprecision(2) << d_warp << " ms"
                  << std::setw(9)  << std::setprecision(2) << d_det << " ms"
                  << std::setw(11) << std::setprecision(2) << d_sum << " ms" << "\n";

        // 4. Wizualizacja klastrów IMO (jeśli podano --viz)
        if (!a.viz_dir.empty()) {
            bool should_save = true;
            if (a.viz_only_imo && res.n_clusters == 0) should_save = false;
            if ((processed_windows - 1) % a.viz_stride != 0) should_save = false;

            if (should_save) {
                // Wyznaczenie bounding boxów klastrów we współrzędnych sensora
                std::vector<ClusterInfo> sensor_clusters = res.clusters;
                for (auto& cl : sensor_clusters) {
                    cl.x_min = DEFAULT_W; cl.y_min = DEFAULT_H;
                    cl.x_max = -1; cl.y_max = -1;
                }
                for (size_t i = 0; i < n; ++i) {
                    int cid = res.event_cluster[i];
                    if (cid > 0 && cid <= static_cast<int>(sensor_clusters.size())) {
                        int px = static_cast<int>(std::round(win.events[i].px));
                        int py = static_cast<int>(std::round(win.events[i].py));
                        auto& cl = sensor_clusters[cid - 1];
                        cl.x_min = std::min(cl.x_min, px);
                        cl.y_min = std::min(cl.y_min, py);
                        cl.x_max = std::max(cl.x_max, px);
                        cl.y_max = std::max(cl.y_max, py);
                    }
                }

                std::vector<float> orig_px(n), orig_py(n);
                std::vector<uint8_t> pol(n);
                for (size_t i = 0; i < n; ++i) {
                    orig_px[i] = win.events[i].px;
                    orig_py[i] = win.events[i].py;
                    pol[i] = win.events[i].p;
                }

                std::ostringstream ss_info;
                ss_info << "WIN #" << processed_windows
                        << " | T: " << std::fixed << std::setprecision(3) << tc << "S"
                        << " | IMO: " << imo_count << " EV (" << std::setprecision(1) << (100.0 * imo_count / n) << "%)"
                        << " | KLASTRY: " << res.n_clusters;

                ImageRGB out_img;
                if (a.viz_mode == "sensor") {
                    out_img = render_imo_frame(orig_px.data(), orig_py.data(), pol.data(),
                                               res.is_imo, res.event_cluster.data(), n,
                                               sensor_clusters, "WIDOK SENSORA (RAW EVENTS + IMO)",
                                               ss_info.str(), a.viz_theme);
                } else if (a.viz_mode == "warped") {
                    out_img = render_imo_frame(xs.data(), ys.data(), pol.data(),
                                               res.is_imo, res.event_cluster.data(), n,
                                               res.clusters, "KOMPENSACJA ROTACJI (SO3 WARPED)",
                                               ss_info.str(), a.viz_theme);
                } else { // "both"
                    ImageRGB left = render_imo_frame(orig_px.data(), orig_py.data(), pol.data(),
                                                     res.is_imo, res.event_cluster.data(), n,
                                                     sensor_clusters, "WIDOK SENSORA (RAW EVENTS)",
                                                     ss_info.str(), a.viz_theme);
                    ImageRGB right = render_imo_frame(xs.data(), ys.data(), pol.data(),
                                                      res.is_imo, res.event_cluster.data(), n,
                                                      res.clusters, "KOMPENSACJA ROTACJI SO(3)",
                                                      ss_info.str(), a.viz_theme);
                    out_img = make_side_by_side(left, right);
                }

                fs::create_directories(a.viz_dir);
                std::ostringstream fn_ppm, fn_png;
                fn_ppm << a.viz_dir << "/frame_" << std::setfill('0') << std::setw(4) << processed_windows
                       << "_t" << std::fixed << std::setprecision(4) << tc << "_imo.ppm";
                fn_png << a.viz_dir << "/frame_" << std::setfill('0') << std::setw(4) << processed_windows
                       << "_t" << std::fixed << std::setprecision(4) << tc << "_imo.png";

                out_img.write_ppm(fn_ppm.str());
                if (a.png) {
                    std::string conv_cmd = "convert " + fn_ppm.str() + " " + fn_png.str() + " 2>/dev/null";
                    int rc = std::system(conv_cmd.c_str());
                    (void)rc;
                }
            }
        }
    }

    std::cout << "================================================================================\n";
    std::cout << "  PODSUMOWANIE BENCHMARKU WYDAJNOŚCI C++:\n";
    std::cout << "  Przetworzonych okien     : " << processed_windows << "\n";
    std::cout << "  Przetworzonych zdarzeń   : " << processed_events << "\n";
    std::cout << "  Łączny czas obliczeń     : " << std::fixed << std::setprecision(2) << total_time_ms << " ms\n";
    std::cout << "  Średni czas per okno     : " << (total_time_ms / processed_windows) << " ms "
              << " (undistort: " << (total_undist_ms / processed_windows) << " ms, "
              << "warp: " << (total_warp_ms / processed_windows) << " ms, "
              << "detect: " << (total_det_ms / processed_windows) << " ms)\n";

    double throughput_mev = (processed_events / (total_time_ms * 1e-3)) * 1e-6;
    std::cout << "  Średnia przepustowość     : " << std::setprecision(2) << throughput_mev << " Mev/s\n";
    double rt_factor = (processed_windows * a.step_ms) / total_time_ms;
    std::cout << "  Współczynnik czasu rzecz.: " << std::setprecision(2) << rt_factor
              << "x czasu rzeczywistego (przetwarza " << rt_factor << "s danych w 1s)\n";
    if (!a.viz_dir.empty()) {
        std::cout << "  Wizualizacje zapisano do : " << a.viz_dir << "/ (tryb: " << a.viz_mode
                  << ", motyw: " << a.viz_theme << (a.png ? ", format: PPM + PNG" : ", format: PPM") << ")\n";
    }
    std::cout << "================================================================================\n";

    return 0;
}
