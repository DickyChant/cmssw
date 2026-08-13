#include "FastSimulation/ShowerDevelopment/interface/HGCalHadronModel.h"
#include "FastSimulation/CaloGeometryTools/interface/HGCalFastGeometry.h"
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <Math/SpecFuncMathCore.h>

#include <algorithm>
#include <cmath>

namespace {
  constexpr double kCLight = 29.9792458;  // cm/ns

  inline double gammaCdf(double a, double x) {
    return (x <= 0.) ? 0. : ROOT::Math::inc_gamma(a, x);
  }
}  // namespace

HGCalHadronModel::HGCalHadronModel(const edm::ParameterSet& ps,
                                   const HGCalFastGeometry* cee,
                                   const HGCalFastGeometry* ceh)
    : cee_(cee), ceh_(ceh) {
  const edm::ParameterSet lon = ps.getParameter<edm::ParameterSet>("Longitudinal");
  a0_ = lon.getParameter<double>("alphaSlope");
  a1_ = lon.getParameter<double>("alphaConst");
  t0_ = lon.getParameter<double>("tSlope");
  t1_ = lon.getParameter<double>("tConst");
  sigmaLnAlpha_ = lon.getParameter<double>("sigmaLnAlpha");
  sigmaLnT_ = lon.getParameter<double>("sigmaLnT");

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

double HGCalHadronModel::meanAlpha(double y) const { return a0_ * std::log(y) + a1_; }
double HGCalHadronModel::meanT(double y) const { return t0_ * std::log(y) + t1_; }

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
  if (e0 <= 0. || cee_ == nullptr)
    return;

  const double y = e0 / ecrit_;
  if (y <= 1.)
    return;

  // Fluctuate the shape. Hadronic showers fluctuate more than electromagnetic
  // ones, which is carried by the configured widths.
  const double alpha = std::max(std::exp(std::log(std::max(meanAlpha(y), 1.05)) +
                                         sigmaLnAlpha_ * random->gaussShoot(0., 1.)),
                                1.05);
  const double T = std::max(std::exp(std::log(std::max(meanT(y), 1.0)) +
                                     sigmaLnT_ * random->gaussShoot(0., 1.)),
                            1.0);
  const double beta = (alpha - 1.) / T;

  // Visible energy: only a fraction of the incident energy becomes silicon
  // signal, and less than for an electromagnetic shower of the same energy.
  const double eVisible = e0 * ehRatio(y);

  // Layer fractions from the Gamma integrated over each layer, in layer units.
  std::vector<double> eLayer(kNTotal + 1, 0.);
  double norm = 0.;
  for (int L = 1; L <= kNTotal; ++L) {
    const double f = gammaCdf(alpha, beta * L) - gammaCdf(alpha, beta * (L - 1));
    eLayer[L] = std::max(f, 0.);
    norm += eLayer[L];
  }
  if (norm <= 0.)
    return;
  for (int L = 1; L <= kNTotal; ++L)
    eLayer[L] *= eVisible / norm;

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

    const double tau = static_cast<double>(L) / std::max(T, 1.0);
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
