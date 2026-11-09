// Minimal, dependency-free reader for NumPy .npy files (memory-mapped, POSIX).
// Supports: format v1/v2/v3, little-endian, C-order, dtypes f4 f8 i1 i2 i4 i8 u1 u2 u4 u8 b1.
#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

class NpyArray {
public:
    explicit NpyArray(const std::string& path) {
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0) throw std::runtime_error("Cannot open " + path);
        struct stat st {};
        if (::fstat(fd_, &st) != 0) fail(path, "fstat failed");
        map_len_ = static_cast<size_t>(st.st_size);
        if (map_len_ < 12) fail(path, "file too small");
        map_ = ::mmap(nullptr, map_len_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (map_ == MAP_FAILED) { map_ = nullptr; fail(path, "mmap failed"); }
        parse_header(path);
    }
    ~NpyArray() { release(); }
    NpyArray(const NpyArray&) = delete;
    NpyArray& operator=(const NpyArray&) = delete;

    const std::vector<size_t>& shape() const { return shape_; }
    size_t size() const { return count_; }          // total number of elements
    char kind() const { return kind_; }
    size_t itemsize() const { return itemsize_; }
    const void* data() const { return static_cast<const unsigned char*>(map_) + data_offset_; }

    // Calls f(const T*) with T matching the stored dtype. Lets callers write one generic lambda.
    template <class F>
    void visit(F&& f) const {
        const void* p = data();
        switch (kind_) {
            case 'f':
                if (itemsize_ == 4) return f(static_cast<const float*>(p));
                if (itemsize_ == 8) return f(static_cast<const double*>(p));
                break;
            case 'i':
                if (itemsize_ == 1) return f(static_cast<const int8_t*>(p));
                if (itemsize_ == 2) return f(static_cast<const int16_t*>(p));
                if (itemsize_ == 4) return f(static_cast<const int32_t*>(p));
                if (itemsize_ == 8) return f(static_cast<const int64_t*>(p));
                break;
            case 'u':
                if (itemsize_ == 1) return f(static_cast<const uint8_t*>(p));
                if (itemsize_ == 2) return f(static_cast<const uint16_t*>(p));
                if (itemsize_ == 4) return f(static_cast<const uint32_t*>(p));
                if (itemsize_ == 8) return f(static_cast<const uint64_t*>(p));
                break;
            case 'b':
                return f(static_cast<const uint8_t*>(p));
            default: break;
        }
        throw std::runtime_error(std::string("Unsupported dtype: ") + descr_);
    }

    std::string describe() const {
        std::string s = descr_ + " shape=(";
        for (size_t i = 0; i < shape_.size(); ++i) s += (i ? "," : "") + std::to_string(shape_[i]);
        return s + (shape_.size() == 1 ? ",)" : ")");
    }

private:
    int fd_ = -1;
    void* map_ = nullptr;
    size_t map_len_ = 0, data_offset_ = 0, count_ = 0, itemsize_ = 0;
    char kind_ = 0;
    std::string descr_;
    std::vector<size_t> shape_;

    void release() {
        if (map_) { ::munmap(map_, map_len_); map_ = nullptr; }
        if (fd_ >= 0) { ::close(fd_); fd_ = -1; }
    }
    [[noreturn]] void fail(const std::string& path, const std::string& msg) {
        release();
        throw std::runtime_error(path + ": " + msg);
    }

    static std::string value_after(const std::string& h, const std::string& key) {
        auto p = h.find("'" + key + "'");
        if (p == std::string::npos) throw std::runtime_error("npy header: missing key " + key);
        p = h.find(':', p) + 1;
        while (p < h.size() && h[p] == ' ') ++p;
        return h.substr(p);
    }

    void parse_header(const std::string& path) {
        const unsigned char* b = static_cast<const unsigned char*>(map_);
        if (std::memcmp(b, "\x93NUMPY", 6) != 0) fail(path, "not a .npy file");
        const int major = b[6];
        size_t hlen, hstart;
        if (major == 1) { hlen = b[8] | (b[9] << 8); hstart = 10; }
        else { hlen = b[8] | (b[9] << 8) | (b[10] << 16) | (static_cast<size_t>(b[11]) << 24); hstart = 12; }
        if (hstart + hlen > map_len_) fail(path, "truncated header");
        const std::string h(reinterpret_cast<const char*>(b + hstart), hlen);
        data_offset_ = hstart + hlen;

        std::string d = value_after(h, "descr");
        descr_ = d.substr(1, d.find(d[0], 1) - 1);               // e.g. "<f8"
        if (descr_.size() < 3) fail(path, "bad descr " + descr_);
        if (descr_[0] == '>') fail(path, "big-endian not supported");
        kind_ = descr_[1];
        itemsize_ = static_cast<size_t>(std::stoi(descr_.substr(2)));

        std::string s = value_after(h, "shape");
        s = s.substr(1, s.find(')') - 1);
        count_ = 1;
        size_t pos = 0;
        while (pos < s.size()) {
            while (pos < s.size() && (s[pos] == ' ' || s[pos] == ',')) ++pos;
            if (pos >= s.size()) break;
            size_t used = 0;
            size_t v = std::stoull(s.substr(pos), &used);
            shape_.push_back(v); count_ *= v; pos += used;
        }
        if (value_after(h, "fortran_order").rfind("True", 0) == 0) {
            size_t big = 0; for (auto v : shape_) big += (v > 1);
            if (big > 1) fail(path, "fortran_order arrays not supported");
        }
        if (data_offset_ + count_ * itemsize_ > map_len_) fail(path, "file shorter than header claims");
    }
};
