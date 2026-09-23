#include "spectre.hpp"
#include <cstdio>
#include <vector>
int main() {
  constexpr int N = 256;
  constexpr double fs = 8000.0;
  std::vector<double> x(N);
  for (int n = 0; n < N; ++n) {
    double t = n / fs;
    x[n] = 0.7 * std::sin(2 * spectre::kPi * 440 * t) +
           0.3 * std::sin(2 * spectre::kPi * 1200 * t);
  }
  auto lp = spectre::Biquad::lowpass(800, fs, 0.707);
  auto y = x;
  for (double& s : y) s = lp.tick(s);
  auto win = y; spectre::hann(win);
  std::vector<spectre::cx> spec(N);
  for (int i = 0; i < N; ++i) spec[i] = {win[i], 0};
  spectre::fft(spec);
  auto db = spectre::magnitude_db(spec);
  std::printf("SpectreForge N=%d fs=%.0f\n", N, fs);
  std::printf("Goertzel LPF: 440=%.4f 1200=%.4f\n",
    spectre::goertzel(y,440,fs), spectre::goertzel(y,1200,fs));
  for (size_t k = 1; k < db.size(); ++k)
    if (db[k] > -25)
      std::printf("  bin %zu  %.1f Hz  %.1f dB\n", k, k*fs/N, db[k]);
}
