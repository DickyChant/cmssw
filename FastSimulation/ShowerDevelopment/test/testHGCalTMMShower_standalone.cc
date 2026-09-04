// Standalone cross-validation of the C++ TMM sampler vs the python
// population: draws N decks at x = ln(50), prints per-block means/stds.
#include "FastSimulation/ShowerDevelopment/interface/HGCalTMMShower.h"
#include <cstdio>
#include <random>
#include <vector>

int main(int argc, char** argv) {
  HGCalTMMShower model(argv[1]);
  std::mt19937_64 eng(42);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::normal_distribution<double> G(0.0, 1.0);
  auto flat = [&]() { return U(eng); };
  auto gauss = [&]() { return G(eng); };
  const double x = std::log(50.0);
  const int N = 20000;
  int K;
  std::vector<double> th;
  // accumulate stats for the dominant K
  std::vector<double> s1, s2;
  int nAcc = 0, kDom = -1;
  std::vector<int> kCount(20, 0);
  for (int i = 0; i < N; ++i) {
    model.drawDeck(x, flat, gauss, K, th);
    if (K < 20) kCount[K]++;
    if (kDom < 0) kDom = K;
    if (K != kDom) continue;
    if (s1.empty()) { s1.assign(th.size(), 0); s2.assign(th.size(), 0); }
    for (size_t d = 0; d < th.size(); ++d) {
      s1[d] += th[d];
      s2[d] += th[d] * th[d];
    }
    ++nAcc;
  }
  printf("K counts:");
  for (int k = 0; k < 20; ++k)
    if (kCount[k]) printf(" K=%d:%d", k, kCount[k]);
  printf("\nstats for K=%d (n=%d):\n", kDom, nAcc);
  for (size_t d = 0; d < s1.size(); ++d) {
    double m = s1[d] / nAcc;
    double sd = std::sqrt(std::max(s2[d] / nAcc - m * m, 0.0));
    printf("dim %2zu: mean %+.4f std %.4f\n", d, m, sd);
  }
  // density sanity: integrate over a coarse grid at layer 8
  double sum = 0.0;
  model.drawDeck(x, flat, gauss, K, th);
  for (double cx = -15; cx <= 15; cx += 0.25)
    for (double cy = -15; cy <= 15; cy += 0.25)
      sum += model.density(th, K, 8, cx, cy) * 0.0625;
  printf("density layer-8 integral (one deck): %.4f  eSpot(8)=%.2e "
         "ring(3)=%.3f\n", sum, model.eSpot(x, 8), model.ringCal(x, 3));
  return 0;
}
