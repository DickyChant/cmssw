#include "FastSimulation/ShowerDevelopment/interface/HGCalHadronModel.h"
#include "FastSimulation/CaloGeometryTools/interface/HGCalFastGeometry.h"
#include "FastSimulation/CalorimeterProperties/interface/HGCalProperties.h"
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <Math/SpecFuncMathCore.h>

#include <algorithm>
#include <cmath>

namespace {
  constexpr double kCLight = 29.9792458;  // cm/ns

  inline double gammaCdf(double a, double x) { return (x <= 0.) ? 0. : ROOT::Math::inc_gamma(a, x); }
}  // namespace

HGCalHadronModel::HGCalHadronModel(const edm::ParameterSet& ps,
                                   const HGCalFastGeometry* cee,
                                   const HGCalFastGeometry* ceh,
                                   const HGCalProperties* props)
    : cee_(cee), ceh_(ceh), props_(props) {
  const edm::ParameterSet lon = ps.getParameter<edm::ParameterSet>("Longitudinal");
  const std::vector<double> muS = lon.getParameter<std::vector<double>>("muSlope");
  const std::vector<double> muC = lon.getParameter<std::vector<double>>("muConst");
  const std::vector<double> sgS = lon.getParameter<std::vector<double>>("sigmaSlope");
  const std::vector<double> sgC = lon.getParameter<std::vector<double>>("sigmaConst");
  for (int i = 0; i < kNPar; ++i) {
    muSlope_[i] = muS.at(i);
    muConst_[i] = muC.at(i);
    sigmaSlope_[i] = sgS.at(i);
    sigmaConst_[i] = sgC.at(i);
  }
  const std::vector<double> ch = lon.getParameter<std::vector<double>>("cholCorr");
  for (int i = 0; i < kNPar * (kNPar + 1) / 2; ++i)
    chol_[i] = ch.at(i);
  cehDepthPerLayer_ = lon.getParameter<double>("cehDepthPerLayer");

  eh0_ = ps.getParameter<double>("ehSlope");
  eh1_ = ps.getParameter<double>("ehConst");
  ecrit_ = ps.getParameter<double>("criticalEnergy");

  const edm::ParameterSet tr = ps.getParameter<edm::ParameterSet>("Transverse");
  r68Slope_ = tr.getParameter<double>("r68Slope");
  r68Const_ = tr.getParameter<double>("r68Const");
  coreOverR68_ = tr.getParameter<double>("coreOverR68");
  tailOverCore_ = tr.getParameter<double>("tailOverCore");
  coreFraction_ = tr.getParameter<double>("coreFraction");

  spotEnergyGeV_ = ps.getParameter<double>("spotEnergy");
  maxSpots_ = ps.getParameter<unsigned>("maxSpots");
}

double HGCalHadronModel::ehRatio(double y) const {
  // Measured 0.623 / 0.557 / 0.548 at 5 / 50 / 500 GeV: falls with energy as the
  // neutral-pion fraction grows, so a single constant will not do.
  return std::clamp(eh0_ * std::log(y) + eh1_, 0.2, 1.0);
}

void HGCalHadronModel::compute(double e0,
                               const double entry[3],
                               const double dir[3],
                               const RandomEngineAndDistribution* random,
                               std::vector<HGCalShowerSpot>& spots) const {
  spots.clear();
  if (e0 <= 0. || cee_ == nullptr || props_ == nullptr)
    return;

  const double y = e0 / ecrit_;
  if (y <= 1.)
    return;
  const double lny = std::log(y);

  // ---- correlated draw of the 5 shape parameters -------------------------
  // theta = mu(ln y) + D(sigma(ln y)) . L . z. The correlations matter: aE-TE
  // at +0.81 (early first interaction -> both compartments shift together) and
  // TH-fH at +0.52 (late showers put more into CE-H).
  double z[kNPar];
  for (int i = 0; i < kNPar; ++i)
    z[i] = random->gaussShoot(0., 1.);
  double theta[kNPar];
  int k = 0;
  for (int i = 0; i < kNPar; ++i) {
    double corr = 0.;
    for (int j = 0; j <= i; ++j)
      corr += chol_[k++] * z[j];
    const double sig = std::clamp(sigmaSlope_[i] * lny + sigmaConst_[i], 0.05, 5.0);
    theta[i] = (muSlope_[i] * lny + muConst_[i]) + sig * corr;
  }
  const double aE = std::clamp(std::exp(theta[0]), 1.05, 200.);
  const double TE = std::clamp(std::exp(theta[1]), 1.0, 60.);
  const double aH = std::clamp(std::exp(theta[2]), 0.2, 200.);
  const double TH = std::clamp(std::exp(theta[3]), 1.0, 80.);
  const double fH = 1. / (1. + std::exp(-theta[4]));
  // MEAN-depth convention (moment estimator): beta = alpha / T.
  const double bE = aE / TE;
  const double bH = aH / TH;

  // Visible energy: only a fraction of the incident energy becomes silicon
  // signal, and less than for an electromagnetic shower of the same energy.
  const double eVisible = e0 * ehRatio(y);

  // ---- layer fractions: two compartment Gammas ---------------------------
  // CE-E on its per-layer X0 table; CE-H on a uniform local depth axis, both
  // exactly as in the parameter derivation. Each compartment is normalized on
  // its own support, then weighted by (1-fH) / fH: fH IS the CE-H energy
  // share, not a tail integral.
  std::vector<double> eLayer(kNTotal + 1, 0.);
  double normE = 0., normH = 0.;
  for (int L = 1; L <= kNCEE; ++L) {
    const double f = props_->frontX0(L);
    const double b = f + props_->layerX0(L);
    eLayer[L] = std::max(gammaCdf(aE, bE * b) - gammaCdf(aE, bE * f), 0.);
    normE += eLayer[L];
  }
  for (int L = kNCEE + 1; L <= kNTotal; ++L) {
    const double f = (L - kNCEE - 1) * cehDepthPerLayer_;
    const double b = f + cehDepthPerLayer_;
    eLayer[L] = std::max(gammaCdf(aH, bH * b) - gammaCdf(aH, bH * f), 0.);
    normH += eLayer[L];
  }
  if (normE <= 0. && normH <= 0.)
    return;
  for (int L = 1; L <= kNCEE; ++L)
    eLayer[L] *= (normE > 0.) ? eVisible * (1. - fH) / normE : 0.;
  for (int L = kNCEE + 1; L <= kNTotal; ++L)
    eLayer[L] *= (normH > 0.) ? eVisible * fH / normH : 0.;

  // Basis perpendicular to the shower axis (as in the electromagnetic model).
  double u1[3], u2[3];
  {
    const double ax = std::abs(dir[0]), ay = std::abs(dir[1]), az = std::abs(dir[2]);
    double ref[3] = {0., 0., 0.};
    if (ax <= ay && ax <= az)
      ref[0] = 1.;
    else if (ay <= az)
      ref[1] = 1.;
    else
      ref[2] = 1.;
    u1[0] = ref[1] * dir[2] - ref[2] * dir[1];
    u1[1] = ref[2] * dir[0] - ref[0] * dir[2];
    u1[2] = ref[0] * dir[1] - ref[1] * dir[0];
    const double n1 = std::sqrt(u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2]);
    for (int i = 0; i < 3; ++i)
      u1[i] /= (n1 > 0. ? n1 : 1.);
    u2[0] = dir[1] * u1[2] - dir[2] * u1[1];
    u2[1] = dir[2] * u1[0] - dir[0] * u1[2];
    u2[2] = dir[0] * u1[1] - dir[1] * u1[0];
  }

  const int zside = (dir[2] >= 0.) ? 1 : -1;

  // Transverse scale: tau relative to the energy-weighted global mean depth,
  // expressed in layers as before. Convert the two compartment means back to a
  // global mean layer for the tau axis.
  const double meanLayerE = TE / std::max(props_->totalX0() / kNCEE, 1e-3);
  const double meanLayerH = kNCEE + TH / cehDepthPerLayer_;
  const double meanLayer = std::clamp((1. - fH) * meanLayerE + fH * meanLayerH, 1., double(kNTotal));

  for (int L = 1; L <= kNTotal; ++L) {
    if (eLayer[L] <= 0.)
      continue;

    // Layer position from the geometry: CE-E layers 1..26, CE-H 27..47 which the
    // CE-H geometry numbers from 1 again.
    const HGCalFastGeometry* geom = (L <= kNCEE) ? cee_ : ceh_;
    const int localLayer = (L <= kNCEE) ? L : (L - kNCEE);
    if (geom == nullptr)
      continue;
    const double zLayer = geom->layerZ(localLayer);
    if (!(zLayer > 0.) || std::abs(dir[2]) < 1e-6)
      continue;

    const double cz = zside * zLayer;
    const double s = (cz - entry[2]) / dir[2];
    if (s < 0.)
      continue;
    const double cx = entry[0] + dir[0] * s;
    const double cy = entry[1] + dir[1] * s;

    unsigned n = static_cast<unsigned>(std::ceil(eLayer[L] / std::max(spotEnergyGeV_, 1e-9)));
    n = std::max(1u, std::min(n, maxSpots_));
    const double eSpot = eLayer[L] / n;

    const double tau = static_cast<double>(L) / meanLayer;
    const double rc = std::max(r68Slope_ * tau + r68Const_, 0.1) * coreOverR68_;
    const double rt = tailOverCore_ * rc;

    for (unsigned is = 0; is < n; ++is) {
      const double R = (random->flatShoot() < coreFraction_) ? rc : rt;
      const double u = std::min(std::max(random->flatShoot(), 1e-9), 1. - 1e-9);
      const double r = R * std::sqrt(u / (1. - u));
      const double phi = 2. * M_PI * random->flatShoot();
      const double qx = r * (std::cos(phi) * u1[0] + std::sin(phi) * u2[0]);
      const double qy = r * (std::cos(phi) * u1[1] + std::sin(phi) * u2[1]);
      const double qz = r * (std::cos(phi) * u1[2] + std::sin(phi) * u2[2]);
      const double lambda = -qz / dir[2];  // slide back onto the layer plane

      HGCalShowerSpot sp;
      sp.layer = L;  // GLOBAL layer; the caller splits CE-E from CE-H
      sp.x = cx + qx + lambda * dir[0];
      sp.y = cy + qy + lambda * dir[1];
      sp.z = cz;
      sp.energy = eSpot;
      sp.t = std::sqrt(sp.x * sp.x + sp.y * sp.y + sp.z * sp.z) / kCLight;
      spots.push_back(sp);
    }
  }
}
