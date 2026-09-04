// TMM splatting shower model -- standalone sampler + renderer core.
// See interface/HGCalTMMShower.h and doc/TMM_SPLAT_INTEGRATION.md.
// Pure standard C++; parameters from the flat file written by
// test/export_tmm_params.py.

#include "FastSimulation/ShowerDevelopment/interface/HGCalTMMShower.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {
  double phiCdf(double t) { return 0.5 * (1.0 + std::erf(t * M_SQRT1_2)); }

  double tLogPdf(double d2, double sig2, double nu) {
    if (nu > 1e5)
      return -d2 / (2.0 * sig2) - std::log(2.0 * M_PI * sig2);
    return std::lgamma(0.5 * (nu + 2.0)) - std::lgamma(0.5 * nu) - std::log(M_PI * nu * sig2) -
           0.5 * (nu + 2.0) * std::log1p(d2 / (nu * sig2));
  }

  std::vector<double> readRow(std::istream& in, size_t n) {
    std::vector<double> v(n);
    for (size_t i = 0; i < n; ++i)
      if (!(in >> v[i]))
        throw std::runtime_error("HGCalTMMShower: truncated param file");
    return v;
  }

  void expectTag(std::istream& in, const std::string& tag) {
    std::string s;
    in >> s;
    if (s != tag)
      throw std::runtime_error("HGCalTMMShower: expected '" + tag + "', got '" + s + "'");
  }

  // interpolation helpers on an ascending grid
  void bracket(const std::vector<double>& g, double x, int& j, double& t) {
    j = int(std::upper_bound(g.begin(), g.end(), x) - g.begin()) - 1;
    j = std::max(0, std::min<int>(j, int(g.size()) - 2));
    t = (x - g[j]) / (g[j + 1] - g[j]);
    t = std::max(0.0, std::min(1.0, t));
  }

  // in-place Cholesky of a symmetric PSD matrix (row-major, dim D);
  // small ridge for numerical safety
  void cholesky(std::vector<double>& A, int D) {
    for (int i = 0; i < D; ++i)
      A[i * D + i] += 1e-9;
    for (int i = 0; i < D; ++i) {
      for (int j = 0; j <= i; ++j) {
        double s = A[i * D + j];
        for (int k = 0; k < j; ++k)
          s -= A[i * D + k] * A[j * D + k];
        if (i == j) {
          A[i * D + i] = std::sqrt(std::max(s, 1e-12));
        } else {
          A[i * D + j] = s / A[j * D + j];
        }
      }
      for (int j = i + 1; j < D; ++j)
        A[i * D + j] = 0.0;
    }
  }
}  // namespace

HGCalTMMShower::HGCalTMMShower(const std::string& paramFile) { loadParams(paramFile); }

void HGCalTMMShower::loadParams(const std::string& paramFile) {
  std::ifstream in(paramFile);
  if (!in)
    throw std::runtime_error("HGCalTMMShower: cannot open " + paramFile);
  std::string tok;
  while (in >> tok) {
    if (tok[0] == '#') {
      std::string line;
      std::getline(in, line);
      continue;
    }
    if (tok == "ZLAYERS") {
      int nl;
      in >> nl;
      for (int l = 0; l < nl; ++l) {
        double a, b;
        in >> a >> b;
        if (l < kNLay) {
          zLo_[l] = a;
          zHi_[l] = b;
        }
      }
    } else if (tok == "PK") {
      int nx, km;
      in >> nx >> km;
      pKx_.resize(nx);
      pK_.assign(size_t(nx) * km, 0.0);
      kMax_ = km;
      for (int i = 0; i < nx; ++i) {
        in >> pKx_[i];
        for (int k = 0; k < km; ++k)
          in >> pK_[size_t(i) * km + k];
      }
    } else if (tok == "POP") {
      PopulationK p;
      int D, NB, NQ;
      in >> p.K >> D >> NB >> NQ;
      expectTag(in, "MEAN");
      p.meanCoef = readRow(in, size_t(3) * D);
      expectTag(in, "XC");
      p.xCenters = readRow(in, NB);
      expectTag(in, "SIG");
      p.sigma = readRow(in, size_t(NB) * D);
      expectTag(in, "COR");
      p.corr = readRow(in, size_t(NB) * D * D);
      expectTag(in, "QEW");
      p.qEw = readRow(in, size_t(NB) * NQ);
      expectTag(in, "MREW");
      p.mrEw = readRow(in, NB);
      p.D = D;
      p.NQ = NQ;
      pops_.push_back(std::move(p));
    } else if (tok == "CALIB") {
      int nx, nl, nr;
      in >> nx >> nl >> nr;
      nRing_ = nr;
      expectTag(in, "XGRID");
      calib_.xGrid = readRow(in, nx);
      expectTag(in, "ESPOT");
      calib_.eSpot = readRow(in, size_t(nx) * nl);
      expectTag(in, "RING");
      calib_.ringCal = readRow(in, size_t(nx) * nr);
      expectTag(in, "OFFX");
      calib_.layOffX = readRow(in, size_t(nx) * nl);
      expectTag(in, "OFFY");
      calib_.layOffY = readRow(in, size_t(nx) * nl);
    } else {
      throw std::runtime_error("HGCalTMMShower: unknown tag " + tok);
    }
  }
  if (pops_.empty())
    throw std::runtime_error("HGCalTMMShower: no populations loaded");
}

void HGCalTMMShower::drawDeck(double x, UniformRng flat, GaussRng gauss, int& K, std::vector<double>& theta) const {
  // 1. K ~ P(K | x)
  int j;
  double t;
  bracket(pKx_, x, j, t);
  std::vector<double> pk(kMax_);
  double tot = 0.0;
  for (int k = 0; k < kMax_; ++k) {
    pk[k] = (1 - t) * pK_[size_t(j) * kMax_ + k] + t * pK_[size_t(j + 1) * kMax_ + k];
    tot += pk[k];
  }
  double u = flat() * tot, acc = 0.0;
  K = pops_.front().K;
  for (int k = 0; k < kMax_; ++k) {
    acc += pk[k];
    if (u <= acc) {
      K = k + 1;
      break;
    }
  }
  const PopulationK* pop = nullptr;
  for (const auto& p : pops_)
    if (p.K == K)
      pop = &p;
  if (!pop) {
    pop = &pops_.front();
    K = pop->K;
  }
  const int D = pop->D;
  // 2. interpolated moments + Cholesky draw
  int jb;
  double tb;
  bracket(pop->xCenters, x, jb, tb);
  std::vector<double> sig(D), C(size_t(D) * D);
  for (int d = 0; d < D; ++d)
    sig[d] = (1 - tb) * pop->sigma[size_t(jb) * D + d] + tb * pop->sigma[size_t(jb + 1) * D + d];
  for (int a = 0; a < D; ++a)
    for (int b = 0; b < D; ++b)
      C[size_t(a) * D + b] =
          (1 - tb) * pop->corr[(size_t(jb) * D + a) * D + b] + tb * pop->corr[(size_t(jb + 1) * D + a) * D + b];
  cholesky(C, D);
  std::vector<double> z(D), zc(D);
  for (int d = 0; d < D; ++d)
    z[d] = gauss();
  for (int a = 0; a < D; ++a) {
    double s = 0.0;
    for (int b = 0; b <= a; ++b)
      s += C[size_t(a) * D + b] * z[b];
    zc[a] = s;
  }
  theta.assign(D, 0.0);
  for (int d = 0; d < D; ++d) {
    double mean = pop->meanCoef[d] + pop->meanCoef[size_t(D) + d] * x + pop->meanCoef[size_t(2) * D + d] * x * x;
    theta[d] = mean + zc[d] * sig[d];
  }
  // 3. last dim: empirical quantile marginal for ln Ew
  double uq = phiCdf(zc[D - 1]);
  uq = std::max(1e-4, std::min(1.0 - 1e-4, uq));
  const int NQ = pop->NQ;
  double pos = 0.001 + uq * 0.998;  // grid spans [0.001, 0.999]
  double fi = (pos - 0.001) / 0.998 * (NQ - 1);
  int qi = std::max(0, std::min(NQ - 2, int(fi)));
  double qt = fi - qi;
  double q0 = (1 - tb) * pop->qEw[size_t(jb) * NQ + qi] + tb * pop->qEw[size_t(jb + 1) * NQ + qi];
  double q1 = (1 - tb) * pop->qEw[size_t(jb) * NQ + qi + 1] + tb * pop->qEw[size_t(jb + 1) * NQ + qi + 1];
  double mr = (1 - tb) * pop->mrEw[jb] + tb * pop->mrEw[jb + 1];
  double meanEw =
      pop->meanCoef[D - 1] + pop->meanCoef[size_t(D) + D - 1] * x + pop->meanCoef[size_t(2) * D + D - 1] * x * x;
  theta[D - 1] = meanEw + mr + (1 - qt) * q0 + qt * q1;
}

double HGCalTMMShower::density(const std::vector<double>& theta, int K, int layer, double cx, double cy) const {
  if (K <= 0 || theta.empty() || layer < 0 || layer >= kNLay)
    return 0.0;
  // softmax over the ln pi block for exact normalization
  double pmax = -1e300;
  for (int k = 0; k < K; ++k)
    pmax = std::max(pmax, theta[k]);
  double psum = 0.0;
  for (int k = 0; k < K; ++k)
    psum += std::exp(theta[k] - pmax);
  double g = 0.0;
  for (int k = 0; k < K; ++k) {
    const double pi = std::exp(theta[k] - pmax) / psum;
    const double mux = theta[K + k];
    const double muy = theta[2 * K + k];
    const double muz = theta[3 * K + k];
    const double sr = std::exp(theta[4 * K + k]);
    const double sz = std::exp(theta[5 * K + k]);
    const double nu = std::exp(theta[6 * K + k]);
    const double zm = phiCdf((zHi_[layer] - muz) / sz) - phiCdf((zLo_[layer] - muz) / sz);
    const double d2 = (cx - mux) * (cx - mux) + (cy - muy) * (cy - muy);
    g += pi * zm * std::exp(tLogPdf(d2, sr * sr, nu));
  }
  return g;
}

double HGCalTMMShower::eSpot(double x, int layer) const {
  if (calib_.xGrid.empty())
    return 2e-4;
  int j;
  double t;
  bracket(calib_.xGrid, x, j, t);
  const int nl = kNLay;
  return (1 - t) * calib_.eSpot[size_t(j) * nl + layer] + t * calib_.eSpot[size_t(j + 1) * nl + layer];
}

double HGCalTMMShower::ringCal(double x, int ring) const {
  if (calib_.ringCal.empty() || ring < 0 || ring >= nRing_)
    return 1.0;
  int j;
  double t;
  bracket(calib_.xGrid, x, j, t);
  return (1 - t) * calib_.ringCal[size_t(j) * nRing_ + ring] + t * calib_.ringCal[size_t(j + 1) * nRing_ + ring];
}

// ---------------------------------------------------------------------------
// Spot generation: the deck sampled into HGCalShowerSpot quanta.

namespace {
  // Marsaglia-Tsang Gamma(alpha, 1) sampler (alpha >= 0.5 here)
  template <class Flat, class Gauss>
  double gammaDraw(double alpha, Flat& flat, Gauss& gauss) {
    if (alpha < 1.0) {
      const double u = std::max(flat(), 1e-12);
      return gammaDraw(alpha + 1.0, flat, gauss) * std::pow(u, 1.0 / alpha);
    }
    const double d = alpha - 1.0 / 3.0;
    const double c = 1.0 / std::sqrt(9.0 * d);
    while (true) {
      double x = gauss();
      double v = 1.0 + c * x;
      if (v <= 0.0)
        continue;
      v = v * v * v;
      const double u = std::max(flat(), 1e-12);
      if (std::log(u) < 0.5 * x * x + d - d * v + d * std::log(v))
        return d * v;
    }
  }
}  // namespace

void HGCalTMMShower::computeSpots(double e0,
                                  const double entry[3],
                                  const double dir[3],
                                  const double* zLayers,
                                  const double* layerWeights,
                                  UniformRng flat,
                                  GaussRng gauss,
                                  std::vector<HGCalShowerSpot>& spots) const {
  const double x = std::log(std::max(e0, 1e-3));
  int K = 0;
  std::vector<double> th;
  drawDeck(x, flat, gauss, K, th);
  if (K <= 0)
    return;
  const double Ew = std::exp(std::min(th[7 * K], 20.0));
  // per-component mixture weights and per-layer interval masses
  std::vector<double> pi(K);
  double pmax = -1e300;
  for (int k = 0; k < K; ++k)
    pmax = std::max(pmax, th[k]);
  double psum = 0.0;
  for (int k = 0; k < K; ++k) {
    pi[k] = std::exp(th[k] - pmax);
    psum += pi[k];
  }
  for (int k = 0; k < K; ++k)
    pi[k] /= psum;
  std::vector<double> mass(size_t(kNLay) * K), mL(kNLay, 0.0);
  for (int l = 0; l < kNLay; ++l)
    for (int k = 0; k < K; ++k) {
      const double muz = th[3 * K + k];
      const double sz = std::exp(th[5 * K + k]);
      const double zm = phiCdf((zHi_[l] - muz) / sz) - phiCdf((zLo_[l] - muz) / sz);
      mass[size_t(l) * K + k] = pi[k] * zm;
      mL[l] += pi[k] * zm;
    }
  double wm = 0.0;
  for (int l = 0; l < kNLay; ++l)
    wm += layerWeights[l] * mL[l];
  if (wm <= 0.0)
    return;
  const double cnorm = Ew / wm;  // raw energy per unit layer mass
  // transverse basis orthogonal to dir
  double e1[3], e2[3];
  const double az = std::abs(dir[2]);
  if (az < 0.99) {
    e1[0] = -dir[1];
    e1[1] = dir[0];
    e1[2] = 0.0;
  } else {
    e1[0] = 1.0;
    e1[1] = 0.0;
    e1[2] = 0.0;
  }
  double n1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (double& v : e1)
    v /= n1;
  e2[0] = dir[1] * e1[2] - dir[2] * e1[1];
  e2[1] = dir[2] * e1[0] - dir[0] * e1[2];
  e2[2] = dir[0] * e1[1] - dir[1] * e1[0];
  const double zsign = (dir[2] < 0.0) ? -1.0 : 1.0;
  double totRaw = 0.0, totW = 0.0;
  const size_t i0 = spots.size();
  for (int l = 0; l < kNLay; ++l) {
    const double eL = cnorm * mL[l];
    const double es = eSpot(x, l);
    if (eL <= 0.0 || es <= 0.0)
      continue;
    const double lam = eL / es;
    long n;
    if (lam > 50.0) {
      n = std::lround(lam + std::sqrt(lam) * gauss());
    } else {  // Knuth
      const double Lp = std::exp(-lam);
      double p = 1.0;
      n = -1;
      do {
        ++n;
        p *= flat();
      } while (p > Lp);
    }
    if (n <= 0)
      continue;
    const double zg = zsign * zLayers[l];
    const double s = (std::abs(dir[2]) > 1e-6) ? (zg - entry[2]) / dir[2] : 0.0;
    const double px = entry[0] + dir[0] * s;
    const double py = entry[1] + dir[1] * s;
    for (long i = 0; i < n; ++i) {
      // component ~ mass(l, k)
      double u = flat() * mL[l], acc = 0.0;
      int k = K - 1;
      for (int kk = 0; kk < K; ++kk) {
        acc += mass[size_t(l) * K + kk];
        if (u <= acc) {
          k = kk;
          break;
        }
      }
      const double sr = std::exp(th[4 * K + k]);
      const double nu = std::exp(th[6 * K + k]);
      double scale = 1.0;
      if (nu < 1e5)
        scale = std::sqrt(nu / (2.0 * gammaDraw(0.5 * nu, flat, gauss)));
      const double uu = th[K + k] + sr * scale * gauss();
      const double vv = th[2 * K + k] + sr * scale * gauss();
      HGCalShowerSpot sp;
      sp.layer = l + 1;
      sp.x = px + uu * e1[0] + vv * e2[0];
      sp.y = py + uu * e1[1] + vv * e2[1];
      sp.z = zg;
      sp.energy = es;
      const double dx = sp.x - entry[0], dy = sp.y - entry[1], dz = sp.z - entry[2];
      sp.t = std::sqrt(dx * dx + dy * dy + dz * dz) / 29.9792458;
      spots.push_back(sp);
      totRaw += es;
      totW += layerWeights[l] * es;
    }
  }
  // exact weighted-energy closure (the python renderer's final rescale)
  if (totW > 0.0) {
    const double f = Ew / totW;
    for (size_t i = i0; i < spots.size(); ++i)
      spots[i].energy *= f;
  }
}
