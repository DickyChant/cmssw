#include "FastSimulation/ShowerDevelopment/interface/HGCalGFlashModel.h"
#include "FastSimulation/CalorimeterProperties/interface/HGCalProperties.h"
#include "FastSimulation/CaloGeometryTools/interface/HGCalFastGeometry.h"
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <Math/SpecFuncMathCore.h>

#include <algorithm>
#include <cmath>

namespace {
  constexpr double kCLight = 29.9792458;  // cm/ns
  constexpr double kConvX0 = 9.0 / 7.0;   // mean photon conversion depth, X0

  /// Regularised lower incomplete gamma P(a, x): the Gamma CDF.
  inline double gammaCdf(double a, double x) {
    if (x <= 0.)
      return 0.;
    return ROOT::Math::inc_gamma(a, x);
  }
}  // namespace

HGCalGFlashModel::HGCalGFlashModel(const edm::ParameterSet& ps,
                                   const HGCalProperties* props,
                                   const HGCalFastGeometry* geom)
    : props_(props), geom_(geom) {
  const edm::ParameterSet lon = ps.getParameter<edm::ParameterSet>("Longitudinal");
  a0_ = lon.getParameter<double>("alphaSlope");
  a1_ = lon.getParameter<double>("alphaConst");
  t0_ = lon.getParameter<double>("tSlope");
  t1_ = lon.getParameter<double>("tConst");
  sigmaLnAlpha_ = lon.getParameter<double>("sigmaLnAlpha");
  sigmaLnT_ = lon.getParameter<double>("sigmaLnT");
  rhoLnAlphaT_ = lon.getParameter<double>("rhoLnAlphaT");

  const edm::ParameterSet cw = ps.getParameter<edm::ParameterSet>("CassetteSplit");
  wa0_ = cw.getParameter<double>("aSlope");
  wa1_ = cw.getParameter<double>("aConst");
  wb0_ = cw.getParameter<double>("bSlope");
  wb1_ = cw.getParameter<double>("bConst");
  applyW_ = cw.getParameter<bool>("apply");

  const edm::ParameterSet tr = ps.getParameter<edm::ParameterSet>("Transverse");
  r68Slope_ = tr.getParameter<double>("r68Slope");
  r68Const_ = tr.getParameter<double>("r68Const");
  coreOverR68_ = tr.getParameter<double>("coreOverR68");
  tailOverCore_ = tr.getParameter<double>("tailOverCore");
  coreFrac0_ = tr.getParameter<double>("coreFractionSlope");
  coreFrac1_ = tr.getParameter<double>("coreFractionConst");
  applyTransverse_ = tr.getParameter<bool>("apply");

  applyConversion_ = lon.getParameter<bool>("applyPhotonConversion");

  spotEnergyGeV_ = ps.getParameter<double>("spotEnergy");
  maxSpots_ = ps.getParameter<unsigned>("maxSpots");
}

double HGCalGFlashModel::meanAlpha(double y) const { return a0_ * std::log(y) + a1_; }
double HGCalGFlashModel::meanT(double y) const { return t0_ * std::log(y) + t1_; }

double HGCalGFlashModel::wSplit(double x0, double y) const {
  const double lny = std::log(y);
  const double a = wa0_ * lny + wa1_;
  const double b = wb0_ * lny + wb1_;
  return a * (1. - std::exp(b * x0));
}

double HGCalGFlashModel::coreRadius(double tau) const {
  const double r68 = std::max(r68Slope_ * tau + r68Const_, 0.05);
  return coreOverR68_ * r68;
}

double HGCalGFlashModel::tailRadius(double tau) const { return tailOverCore_ * coreRadius(tau); }

double HGCalGFlashModel::coreFraction(double tau) const {
  return std::clamp(coreFrac0_ * tau + coreFrac1_, 0.05, 1.0);
}

void HGCalGFlashModel::compute(double e0,
                               const double entry[3],
                               const double dir[3],
                               bool isGamma,
                               const RandomEngineAndDistribution* random,
                               std::vector<HGCalShowerSpot>& spots) const {
  spots.clear();
  if (!props_ || e0 <= 0.)
    return;

  const double ecrit = props_->criticalEnergyGeV();
  const double y = e0 / ecrit;
  if (y <= 1.)
    return;

  // ---- 1. longitudinal: fluctuate (ln alpha, ln T) with their correlation ----
  const double mLnAlpha = std::log(std::max(meanAlpha(y), 1.05));
  const double mLnT = std::log(std::max(meanT(y), 0.2));

  const double z1 = random->gaussShoot(0., 1.);
  const double z2 = random->gaussShoot(0., 1.);
  const double lnAlpha = mLnAlpha + sigmaLnAlpha_ * z1;
  const double lnT =
      mLnT + sigmaLnT_ * (rhoLnAlphaT_ * z1 + std::sqrt(std::max(1. - rhoLnAlphaT_ * rhoLnAlphaT_, 0.)) * z2);

  const double alpha = std::max(std::exp(lnAlpha), 1.05);
  const double T = std::max(std::exp(lnT), 0.2);
  const double beta = (alpha - 1.) / T;

  // Photons start after a conversion; this is sampled explicitly rather than
  // absorbed into the shape, because fitting it as a free shift is degenerate.
  double x0Offset = 0.;
  if (isGamma && applyConversion_) {
    const double u = std::max(random->flatShoot(), 1e-9);
    x0Offset = -std::log(u) * kConvX0;
  }

  // Path-length scaling: the CE-E layers are planes perpendicular to z, so a
  // particle at angle theta traverses more X0 per layer.
  const double cth = std::abs(dir[2]);
  const double pathScale = (cth > 1e-6) ? 1. / cth : 1.;

  // ---- 2. energy per layer: integral of the Gamma over each layer ----
  const unsigned nl = props_->nLayers();
  std::vector<double> eLayer(nl + 1, 0.);
  double norm = 0.;
  for (unsigned L = 1; L <= nl; ++L) {
    // frontX0/layerX0 are normal-incidence depths, so they scale with pathScale.
    // x0Offset is already a distance along the trajectory and must NOT be scaled
    // again -- doing so delayed the conversion by ~3.7% at eta=2, 10.5% at 1.5.
    const double f = props_->frontX0(L) * pathScale - x0Offset;
    const double b = (props_->frontX0(L) + props_->layerX0(L)) * pathScale - x0Offset;
    if (b <= 0.)
      continue;
    const double lo = gammaCdf(alpha, beta * std::max(f, 0.));
    const double hi = gammaCdf(alpha, beta * b);
    eLayer[L] = std::max(hi - lo, 0.);
    norm += eLayer[L];
  }
  if (norm <= 0.)
    return;
  for (unsigned L = 1; L <= nl; ++L)
    eLayer[L] *= e0 / norm;

  // ---- 3. cassette split ----
  // Cassettes are the layer pairs (2,3), (4,5), ... ; layer 1 (half-lead) and the
  // last layer are unpaired. w redistributes each pair without changing its sum.
  if (applyW_) {
    for (unsigned L = 2; L + 1 <= nl; L += 2) {
      const double pair = eLayer[L] + eLayer[L + 1];
      if (pair <= 0.)
        continue;
      const double x0mid = props_->frontX0(L) + 0.5 * (props_->layerX0(L) + props_->layerX0(L + 1));
      const double w = std::clamp(wSplit(x0mid, y), 0.05, 0.95);
      eLayer[L] = w * pair;
      eLayer[L + 1] = (1. - w) * pair;
    }
  }

  // ---- 4. spots ----
  // Orthonormal basis perpendicular to the shower axis, so radii are sampled in
  // the plane perpendicular to the particle and then land where the trajectory
  // crosses each layer (an ellipse in the layer plane, as in FullSim).
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

  const double tmax = std::max(T, 0.2);
  const int zside = (dir[2] >= 0.) ? 1 : -1;

  for (unsigned L = 1; L <= nl; ++L) {
    if (eLayer[L] <= 0.)
      continue;

    // number of spots for this layer, bounded so a very high energy shower
    // cannot blow up the event size
    unsigned n = static_cast<unsigned>(std::ceil(eLayer[L] / std::max(spotEnergyGeV_, 1e-9)));
    n = std::max(1u, std::min(n, maxSpots_));
    const double eSpot = eLayer[L] / n;

    // depth of this layer along the trajectory, and the local radial scale
    const double x0mid = (props_->frontX0(L) + 0.5 * props_->layerX0(L)) * pathScale;
    const double tau = x0mid / tmax;
    const double rc = coreRadius(tau);
    const double rt = tailRadius(tau);
    const double pcore = coreFraction(tau);

    // Where the trajectory crosses this layer. The physical z spacing is NOT
    // X0 * radLenIncm: the configured stack is 25.6 X0 = 19.7 cm while the real
    // CE-E envelope is ~40.9 cm, so using the X0 length compresses the shower
    // toward the front and lands spots in the wrong cells (~5.8 cm error at the
    // back at eta=2). Take the true plane position from the geometry; the X0
    // table still sets the shower age and the energy sharing.
    double cx, cy, cz;
    if (geom_ != nullptr && geom_->layerZ(static_cast<int>(L)) > 0. && std::abs(dir[2]) > 1e-6) {
      cz = zside * geom_->layerZ(static_cast<int>(L));
      const double sAlongCm = (cz - entry[2]) / dir[2];
      cx = entry[0] + dir[0] * sAlongCm;
      cy = entry[1] + dir[1] * sAlongCm;
    } else {
      const double sAlongCm = (props_->frontX0(L) + 0.5 * props_->layerX0(L)) * props_->radLenIncm() * pathScale;
      cx = entry[0] + dir[0] * sAlongCm;
      cy = entry[1] + dir[1] * sAlongCm;
      cz = entry[2] + dir[2] * sAlongCm;
    }

    for (unsigned is = 0; is < n; ++is) {
      double r = 0.;
      if (applyTransverse_) {
        const double R = (random->flatShoot() < pcore) ? rc : rt;
        const double u = std::min(std::max(random->flatShoot(), 1e-9), 1. - 1e-9);
        r = R * std::sqrt(u / (1. - u));
      }
      const double phi = 2. * M_PI * random->flatShoot();
      const double cphi = std::cos(phi), sphi = std::sin(phi);

      // The offset q is perpendicular to the shower axis, so it has a z
      // component. Adding it directly would take the spot off the layer plane,
      // narrow the in-plane footprint by |dir_z|^2 (7% at eta=2, 18% at 1.5) and
      // could even flip a far-tail spot to the other endcap. Slide back along the
      // axis until z is on the plane again -- that is the intended ellipse.
      const double qx = r * (cphi * u1[0] + sphi * u2[0]);
      const double qy = r * (cphi * u1[1] + sphi * u2[1]);
      const double qz = r * (cphi * u1[2] + sphi * u2[2]);
      const double lambda = (std::abs(dir[2]) > 1e-6) ? (-qz / dir[2]) : 0.;

      HGCalShowerSpot sp;
      sp.layer = static_cast<int>(L);
      sp.x = cx + qx + lambda * dir[0];
      sp.y = cy + qy + lambda * dir[1];
      sp.z = cz;
      sp.energy = eSpot;
      sp.t = std::sqrt(sp.x * sp.x + sp.y * sp.y + sp.z * sp.z) / kCLight;
      spots.push_back(sp);
    }
  }
}
