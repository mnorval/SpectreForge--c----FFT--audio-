#pragma once
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>
namespace spectre {
using cx = std::complex<double>;
constexpr double kPi = 3.14159265358979323846;
inline bool is_pow2(std::size_t n) { return n && ((n & (n - 1)) == 0); }
inline void hann(std::vector<double>& x) {
  const std::size_t n = x.size();
  for (std::size_t i = 0; i < n; ++i)
    x[i] *= 0.5 * (1.0 - std::cos(2.0 * kPi * double(i) / (n - 1)));
}
inline void fft(std::vector<cx>& a) {
  const std::size_t n = a.size();
  if (!is_pow2(n)) return;
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (std::size_t len = 2; len <= n; len <<= 1) {
    const double ang = -2.0 * kPi / double(len);
    const cx wlen(std::cos(ang), std::sin(ang));
    for (std::size_t i = 0; i < n; i += len) {
      cx w(1.0, 0.0);
      for (std::size_t j = 0; j < len / 2; ++j) {
        cx u = a[i + j];
        cx v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
}
inline std::vector<double> magnitude_db(const std::vector<cx>& a) {
  std::vector<double> db(a.size() / 2);
  for (std::size_t k = 0; k < db.size(); ++k) {
    const double mag = std::abs(a[k]) / double(a.size());
    db[k] = 20.0 * std::log10(mag + 1e-12);
  }
  return db;
}
inline double goertzel(const std::vector<double>& x, double freq_hz, double fs) {
  const std::size_t n = x.size();
  const double k = std::round(n * freq_hz / fs);
  const double omega = 2.0 * kPi * k / double(n);
  const double coeff = 2.0 * std::cos(omega);
  double s0 = 0, s1 = 0, s2 = 0;
  for (double v : x) { s0 = v + coeff * s1 - s2; s2 = s1; s1 = s0; }
  const double re = s1 - s2 * std::cos(omega);
  const double im = s2 * std::sin(omega);
  return std::hypot(re, im) / double(n);
}
struct Biquad {
  double b0{}, b1{}, b2{}, a1{}, a2{}, z1{}, z2{};
  double tick(double x) {
    const double y = b0 * x + z1;
    z1 = b1 * x - a1 * y + z2;
    z2 = b2 * x - a2 * y;
    return y;
  }
  static Biquad lowpass(double fc, double fs, double q) {
    const double w0 = 2.0 * kPi * fc / fs;
    const double c = std::cos(w0), s = std::sin(w0);
    const double alpha = s / (2.0 * q);
    const double a0 = 1.0 + alpha;
    Biquad f;
    f.b0 = ((1.0 - c) / 2.0) / a0;
    f.b1 = (1.0 - c) / a0;
    f.b2 = f.b0;
    f.a1 = (-2.0 * c) / a0;
    f.a2 = (1.0 - alpha) / a0;
    return f;
  }
};
}
