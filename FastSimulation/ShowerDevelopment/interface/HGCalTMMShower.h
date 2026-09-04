#ifndef FastSimulation_ShowerDevelopment_HGCalTMMShower_H
#define FastSimulation_ShowerDevelopment_HGCalTMMShower_H

// TMM splatting: table-driven learned-mixture shower model for HGCAL.
// A shower = K anisotropic 3D components with per-component Student-t
// transverse tails, drawn from an energy-conditional Gaussian copula
// (all constants from FastSimulation/ShowerDevelopment/data/*.params).
// Pure standard C++; RNG is injected so the CMSSW random service
// adapter stays a thin wrapper. See doc/TMM_SPLAT_INTEGRATION.md.

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "FastSimulation/CaloHitMakers/interface/HGCalShowerSpot.h"

class HGCalTMMShower {
public:
  struct PopulationK {          // exact-K copula tables
    int K = 0;                  // components; deck dim D = 7*K + 1
    int D = 0;                  // 7K+1
    int NQ = 0;                 // quantile points of the Ew marginal
    std::vector<double> meanCoef;   // (3, D) quadratic-in-x means
    std::vector<double> xCenters;   // (NB) bin centres in x = ln E
    std::vector<double> sigma;      // (NB, D)
    std::vector<double> corr;       // (NB, D, D)
    std::vector<double> qEw;        // (NB, NQ) empirical ln Ew marginal
    std::vector<double> mrEw;       // (NB)
  };
  struct Calib {                // per-x-grid calibration tables
    std::vector<double> xGrid;
    std::vector<double> eSpot;      // (NX, NLAY)
    std::vector<double> ringCal;    // (NX, NRING)
    std::vector<double> layOffX;    // (NX, NLAY)
    std::vector<double> layOffY;    // (NX, NLAY)
  };

  using UniformRng = std::function<double()>;   // U(0,1)
  using GaussRng = std::function<double()>;     // N(0,1)

  explicit HGCalTMMShower(const std::string& paramFile);

  // Draw one shower's deck at ln E = x: fills K, theta (7K+1).
  void drawDeck(double x, UniformRng flat, GaussRng gauss, int& K,
                std::vector<double>& theta) const;

  // Expected energy density of the drawn deck at (cell x, y, layer),
  // in the incidence frame. Caller applies quanta + transplant.
  double density(const std::vector<double>& theta, int K, int layer,
                 double cx, double cy) const;

  // Full spot generation: deck -> HGCalShowerSpot quanta, transplanted
  // to the track frame (entry, dir). zLayers/layerWeights: 47 values
  // from the geometry and the reverse calibration.
  void computeSpots(double e0, const double entry[3], const double dir[3],
                    const double* zLayers, const double* layerWeights,
                    UniformRng flat, GaussRng gauss,
                    std::vector<HGCalShowerSpot>& spots) const;

  // Per-layer quantum and calibration lookups at ln E = x.
  double eSpot(double x, int layer) const;
  double ringCal(double x, int ring) const;

  static constexpr int kNLay = 47;

private:
  std::vector<PopulationK> pops_;
  std::vector<double> pK_;          // P(K | x) table, (NXK, KMAX)
  std::vector<double> pKx_;         // x grid for pK_
  int kMax_ = 0;
  int nRing_ = 0;
  Calib calib_;
  std::array<double, kNLay> zLo_{};  // layer depth intervals [X0]
  std::array<double, kNLay> zHi_{};

  void loadParams(const std::string& paramFile);
};

#endif
