// test_numerics.cpp: Weryfikacja numeryczna implementacji C++ względem referencji Pythona.
// Zestaw testowy bazuje na rzeczywistych danych z sekwencji EVIMO2: scene14_dyn_test_05_000000 (scene14_05, klatka 85).
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <string>
#include <fstream>
#include <filesystem>

#include "camera.hpp"
#include "calibration.hpp"
#include "rotcomp.hpp"
#include "detector.hpp"
#include "npy.hpp"

using namespace evimo2;

void test_basic_math() {
    std::cout << "[1/7] Testowanie podstawowych struktur Vec3 i Mat3... ";
    Vec3 v1{1.0, 2.0, 3.0};
    Vec3 v2{4.0, 5.0, 6.0};
    Vec3 v3 = v1 + v2;
    assert(std::abs(v3.x - 5.0) < 1e-12 && std::abs(v3.y - 7.0) < 1e-12 && std::abs(v3.z - 9.0) < 1e-12);
    assert(std::abs(v1.dot(v2) - 32.0) < 1e-12);

    Mat3 I = Mat3::identity();
    Vec3 v_rot = I * v1;
    assert(std::abs(v_rot.x - v1.x) < 1e-12 && std::abs(v_rot.y - v1.y) < 1e-12 && std::abs(v_rot.z - v1.z) < 1e-12);

    Mat3 A{1, 2, 3, 4, 5, 6, 7, 8, 9};
    Mat3 At = A.transpose();
    assert(At(0, 1) == 4 && At(1, 0) == 2);

    std::cout << "OK\n";
}

void test_rodrigues() {
    std::cout << "[2/7] Testowanie wzoru Rodriguesa... ";
    // Obrót zerowy -> macierz jednostkowa
    Mat3 R0 = rodrigues(Vec3{0.0, 0.0, 0.0});
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            double expected = (r == c) ? 1.0 : 0.0;
            assert(std::abs(R0(r, c) - expected) < 1e-12);
        }
    }

    // Obrót o 90 stopni wokół osi Z: [0, 0, pi/2]
    double pi_2 = M_PI / 2.0;
    Mat3 Rz90 = rodrigues(Vec3{0.0, 0.0, pi_2});
    // [ 0, -1,  0 ]
    // [ 1,  0,  0 ]
    // [ 0,  0,  1 ]
    assert(std::abs(Rz90(0, 0) - 0.0) < 1e-6);
    assert(std::abs(Rz90(0, 1) - (-1.0)) < 1e-6);
    assert(std::abs(Rz90(1, 0) - 1.0) < 1e-6);
    assert(std::abs(Rz90(1, 1) - 0.0) < 1e-6);
    assert(std::abs(Rz90(2, 2) - 1.0) < 1e-6);

    std::cout << "OK\n";
}

void test_calibration() {
    std::cout << "[3/7] Testowanie kalibracji żyroskopu... ";
    Vec3 raw_gyro{GYRO_BIAS_ARR[0], GYRO_BIAS_ARR[1], GYRO_BIAS_ARR[2]};
    Vec3 cal = calibrated_gyro(raw_gyro);
    assert(std::abs(cal.x) < 1e-12 && std::abs(cal.y) < 1e-12 && std::abs(cal.z) < 1e-12);

    Vec3 raw2{GYRO_BIAS_ARR[0] + 1.0, GYRO_BIAS_ARR[1], GYRO_BIAS_ARR[2]};
    Vec3 cal2 = calibrated_gyro(raw2);
    // R_IMU2CAM = diag(-1, -1, 1) -> delta X = +1.0 powinno dać -1.0
    assert(std::abs(cal2.x - (-1.0)) < 1e-12);

    std::cout << "OK\n";
}

void test_python_reference(const std::string& data_dir) {
    namespace fs = std::filesystem;
    if (!fs::exists(data_dir)) {
        std::cerr << "Pominięto walidację z Pythonem (brak " << data_dir << ")\n";
        return;
    }

    std::string seq_name = "scene14_dyn_test_05_000000 (scene14_05, klatka 85)";
    std::ifstream f_meta(data_dir + "/meta.json");
    if (f_meta.is_open()) {
        std::string line;
        while (std::getline(f_meta, line)) {
            if (line.find("\"sequence\"") != std::string::npos) {
                size_t c = line.find(':');
                if (c != std::string::npos) {
                    size_t q1 = line.find('"', c);
                    size_t q2 = (q1 != std::string::npos) ? line.find('"', q1 + 1) : std::string::npos;
                    if (q2 != std::string::npos) seq_name = line.substr(q1 + 1, q2 - q1 - 1);
                }
            }
        }
    }

    std::cout << "[4/7] Wczytywanie danych referencyjnych z Pythona (" << data_dir << ")\n"
              << "      [Sekwencja: " << seq_name << " | t_c = 1.433 s, wlen = 40 ms]... ";
    NpyArray npy_t(data_dir + "/t_ev.npy");
    NpyArray npy_pix(data_dir + "/pix_ev.npy");
    NpyArray npy_p(data_dir + "/p_ev.npy");
    NpyArray npy_bearing(data_dir + "/bearing.npy");
    NpyArray npy_grid(data_dir + "/grid.npy");
    NpyArray npy_rot_matrices(data_dir + "/rot_matrices.npy");
    NpyArray npy_warped_xs(data_dir + "/warped_xs.npy");
    NpyArray npy_warped_ys(data_dir + "/warped_ys.npy");
    NpyArray npy_iwe_E(data_dir + "/iwe_E.npy");
    NpyArray npy_iwe_L(data_dir + "/iwe_L.npy");
    NpyArray npy_D(data_dir + "/D.npy");
    NpyArray npy_S(data_dir + "/S.npy");
    NpyArray npy_pred_imo(data_dir + "/pred_imo.npy");

    NpyArray npy_imu_ts(data_dir + "/imu_ts.npy");
    NpyArray npy_imu_w(data_dir + "/imu_w.npy");

    size_t n_ev = npy_t.size();
    size_t n_imu = npy_imu_ts.size();
    std::cout << n_ev << " zdarzeń, " << n_imu << " próbek IMU. OK\n";

    // Parametry kamery ze sceny scene14_05
    CameraIntrinsics intr;
    intr.fx = 556.573975;
    intr.fy = 556.086975;
    intr.cx = 321.846008;
    intr.cy = 225.617996;
    intr.k1 = -0.10828;
    intr.k2 = 0.206233;
    intr.p1 = -0.00065;
    intr.p2 = -0.000555;

    double tc = 1.433333;
    double wlen = 0.040;
    double d = -0.0157;

    const double* ref_t = static_cast<const double*>(npy_t.data());
    const double* ref_pix = static_cast<const double*>(npy_pix.data());
    const uint8_t* ref_p = static_cast<const uint8_t*>(npy_p.data());
    const double* ref_bearing = static_cast<const double*>(npy_bearing.data());

    // 1. Walidacja undistort i make_window
    std::cout << "[5/7] Walidacja numeryczna undistort & make_window... ";
    std::vector<double> ev_t(ref_t, ref_t + n_ev);
    std::vector<double> ev_x(n_ev), ev_y(n_ev);
    std::vector<uint8_t> ev_p(ref_p, ref_p + n_ev);
    for (size_t i = 0; i < n_ev; ++i) {
        ev_x[i] = ref_pix[2 * i];
        ev_y[i] = ref_pix[2 * i + 1];
    }

    EventWindow win = prep_window(tc, wlen,
                                  ev_t.begin(), ev_t.end(),
                                  ev_x.begin(), ev_y.begin(), ev_p.begin(),
                                  intr);

    assert(win.events.size() == n_ev);
    double max_err_bearing = 0.0;
    for (size_t i = 0; i < n_ev; ++i) {
        double bx = win.events[i].bearing.x;
        double by = win.events[i].bearing.y;
        double bz = win.events[i].bearing.z;
        double exp_bx = ref_bearing[3 * i];
        double exp_by = ref_bearing[3 * i + 1];
        double exp_bz = ref_bearing[3 * i + 2];

        max_err_bearing = std::max(max_err_bearing, std::abs(bx - exp_bx));
        max_err_bearing = std::max(max_err_bearing, std::abs(by - exp_by));
        max_err_bearing = std::max(max_err_bearing, std::abs(bz - exp_bz));
    }
    std::cout << "Maks. błąd bearing: " << max_err_bearing;
    assert(max_err_bearing < 1e-10);
    std::cout << " (OK)\n";

    // 2. Walidacja rot_matrices i warp
    std::cout << "[6/7] Walidacja numeryczna rot_matrices & warp... ";
    const double* imu_ts = static_cast<const double*>(npy_imu_ts.data());
    const double* imu_w_raw = static_cast<const double*>(npy_imu_w.data());
    std::vector<Vec3> imu_w(n_imu);
    for (size_t i = 0; i < n_imu; ++i) {
        imu_w[i] = Vec3{imu_w_raw[3 * i], imu_w_raw[3 * i + 1], imu_w_raw[3 * i + 2]};
    }

    std::vector<Mat3> M = rot_matrices(win, imu_ts, imu_w.data(), n_imu, d);
    const double* exp_M = static_cast<const double*>(npy_rot_matrices.data());
    double max_err_M = 0.0;
    for (size_t i = 0; i < n_ev; ++i) {
        for (int k = 0; k < 9; ++k) {
            max_err_M = std::max(max_err_M, std::abs(M[i].m[k] - exp_M[9 * i + k]));
        }
    }
    std::cout << "Maks. błąd M: " << max_err_M << ", ";
    assert(max_err_M < 1e-7);

    std::vector<float> xs(n_ev), ys(n_ev);
    warp(win, imu_ts, imu_w.data(), n_imu, d, xs.data(), ys.data());

    const double* exp_xs = static_cast<const double*>(npy_warped_xs.data());
    const double* exp_ys = static_cast<const double*>(npy_warped_ys.data());
    double max_err_warp = 0.0;
    for (size_t i = 0; i < n_ev; ++i) {
        max_err_warp = std::max(max_err_warp, (double)std::abs(xs[i] - exp_xs[i]));
        max_err_warp = std::max(max_err_warp, (double)std::abs(ys[i] - exp_ys[i]));
    }
    std::cout << "Maks. błąd pozycji [px]: " << max_err_warp;
    assert(max_err_warp < 1e-4);
    std::cout << " (OK)\n";

    // 3. Walidacja IWE, residual_map i detect_imo
    std::cout << "[7/7] Walidacja numeryczna detektora IMO v0... ";
    DetectParams params;
    params.blur = 3;
    params.tau = 0.5f;
    params.support = 1.0f;
    params.min_cells = 8;
    params.fill = false;
    params.close = 0;

    std::vector<double> te(n_ev);
    for (size_t i = 0; i < n_ev; ++i) te[i] = win.events[i].t;

    DetectResult det = detect_imo(xs.data(), ys.data(), te.data(), n_ev, tc, params);

    const uint8_t* exp_pred = static_cast<const uint8_t*>(npy_pred_imo.data());
    size_t mismatches = 0;
    for (size_t i = 0; i < n_ev; ++i) {
        bool act = det.is_imo[i];
        bool exp = (exp_pred[i] != 0);
        if (act != exp) {
            ++mismatches;
        }
    }

    std::cout << "Wykrytych klastrów: " << det.n_clusters
              << " | Niezgodności klasyfikacji: " << mismatches << " / " << n_ev
              << " (" << std::fixed << std::setprecision(4) << (100.0 * mismatches / n_ev) << "%)";
    assert(mismatches == 0);
    std::cout << " -> ZGODNOŚĆ 100% (OK)\n";
}

int main(int argc, char** argv) {
    std::cout << "============================================================\n";
    std::cout << "      TESTY JEDNOSTKOWE I WALIDACJA NUMERYCZNA C++ (EVIMO2) \n";
    std::cout << "============================================================\n";

    test_basic_math();
    test_rodrigues();
    test_calibration();

    std::string data_dir = "cpp/tests/data";
    if (argc > 1) {
        data_dir = argv[1];
    }
    test_python_reference(data_dir);

    std::cout << "============================================================\n";
    std::cout << "           WSZYSTKIE TESTY ZAKOŃCZONE SUKCESEM!             \n";
    std::cout << "============================================================\n";
    return 0;
}
