// evimo2_viewer: load an EVIMO2 sequence, print stream statistics, dump event frames.
//
// Usage: evimo2_viewer <seq_dir> [--cam W H] [--t0 SEC] [--dt MS] [--frames N]
//                      [--invert-pol] [--tau MS] [--out DIR]
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>

#include "events.hpp"
#include "viz.hpp"

struct Args {
    std::string dir, out = "out";
    int w = 640, h = 480;
    double t0 = 0.0, dt_ms = 10.0, tau_ms = 30.0;
    int frames = 10;
    bool invert = false;
};

static bool parse(int argc, char** argv, Args& a) {
    if (argc < 2) return false;
    a.dir = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        auto need = [&](int k) { return i + k < argc; };
        if (s == "--cam" && need(2)) { a.w = std::atoi(argv[++i]); a.h = std::atoi(argv[++i]); }
        else if (s == "--t0" && need(1)) a.t0 = std::atof(argv[++i]);
        else if (s == "--dt" && need(1)) a.dt_ms = std::atof(argv[++i]);
        else if (s == "--tau" && need(1)) a.tau_ms = std::atof(argv[++i]);
        else if (s == "--frames" && need(1)) a.frames = std::atoi(argv[++i]);
        else if (s == "--out" && need(1)) a.out = argv[++i];
        else if (s == "--invert-pol") a.invert = true;
        else return false;
    }
    return true;
}

int main(int argc, char** argv) {
    Args a;
    if (!parse(argc, argv, a)) {
        std::fprintf(stderr,
            "Usage: %s <seq_dir> [--cam W H] [--t0 SEC] [--dt MS] [--frames N] "
            "[--invert-pol] [--tau MS] [--out DIR]\n", argv[0]);
        return 1;
    }
    try {
        LoadOptions lo;
        lo.t_begin_s = a.t0;
        lo.t_end_s = a.t0 + a.frames * a.dt_ms * 1e-3;
        lo.invert_polarity = a.invert;

        FileInfo info;
        const auto tic = std::chrono::steady_clock::now();
        std::vector<Event> ev = load_events(a.dir, lo, &info);
        const double load_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - tic).count();

        std::printf("File : %zu events, %.3f s total\n", info.total_events, info.duration_s);
        std::printf("dtype: xy=%s t=%s p=%s\n", info.xy_desc.c_str(), info.t_desc.c_str(), info.p_desc.c_str());
        std::printf("Window [%.3f, %.3f] s: loaded %zu events in %.3f s (%zu MB)\n",
                    lo.t_begin_s, lo.t_end_s, ev.size(), load_s, ev.size() * sizeof(Event) >> 20);
        if (ev.empty()) { std::fprintf(stderr, "No events in window.\n"); return 2; }

        // Sanity statistics
        size_t on = 0, oob = 0, nonmono = 0;
        const int64_t t_first = ev.front().t_us;
        std::vector<uint32_t> per_ms(static_cast<size_t>((ev.back().t_us - t_first) / 1000) + 1, 0);
        for (size_t i = 0; i < ev.size(); ++i) {
            on += ev[i].p;
            oob += (ev[i].x >= a.w || ev[i].y >= a.h);
            if (i && ev[i].t_us < ev[i - 1].t_us) ++nonmono;
            ++per_ms[static_cast<size_t>((ev[i].t_us - t_first) / 1000)];
        }
        const double span_s = std::max(1e-6, (ev.back().t_us - t_first) * 1e-6);
        std::printf("ON ratio: %.3f | out-of-range: %zu | non-monotonic: %zu\n",
                    double(on) / ev.size(), oob, nonmono);
        std::printf("Mean rate: %.3f Mev/s | peak 1ms-bin rate: %.3f Mev/s\n",
                    ev.size() / span_s * 1e-6, *std::max_element(per_ms.begin(), per_ms.end()) * 1e-3);

        // Dump frames
        std::filesystem::create_directories(a.out);
        const int64_t win_start = static_cast<int64_t>(std::llround(a.t0 * 1e6));
        for (int k = 0; k < a.frames; ++k) {
            const int64_t ta = win_start + static_cast<int64_t>(std::llround(k * a.dt_ms * 1e3));
            const int64_t tb = win_start + static_cast<int64_t>(std::llround((k + 1) * a.dt_ms * 1e3));
            auto cmp = [](const Event& e, int64_t t) { return e.t_us < t; };
            const Event* b = &*std::lower_bound(ev.begin(), ev.end(), ta, cmp);
            const Event* e = ev.data() + (std::lower_bound(ev.begin(), ev.end(), tb, cmp) - ev.begin());
            if (e < b) e = b;
            const size_t n = static_cast<size_t>(e - b);
            char name[256];
            std::snprintf(name, sizeof name, "%s/frame_%03d_pol.ppm", a.out.c_str(), k);
            write_ppm_polarity(name, accumulate_signed(b, n, a.w, a.h), a.w, a.h);
            std::snprintf(name, sizeof name, "%s/frame_%03d_ts.pgm", a.out.c_str(), k);
            write_pgm(name, time_surface(b, n, tb, a.tau_ms * 1e3, a.w, a.h), a.w, a.h);
        }
        std::printf("Wrote %d frame pairs to %s/\n", a.frames, a.out.c_str());
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "Error: %s\n", ex.what());
        return 3;
    }
    return 0;
}
