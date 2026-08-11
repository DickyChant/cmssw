// Validation and benchmark for HGCalFastGeometry.
//
// Mode 1 (self-consistency, default): for every valid DetId of the loaded
// geometry, take HGCalGeometry::getPosition(id) and push it back through
// HGCalFastGeometry::toDetId(). The recovered DetId must equal the original.
// This is the same contract Geometry/HGCalGeometry/test/HGCalGeomLocatorTester
// checks for getClosestCell(), so the two are directly comparable -- including
// the fact that getClosestCell() itself is not perfectly invertible (see the
// known-mismatch lists in Geometry/HGCalGeometry/data/).
//
// Mode 2 (benchmark): time toDetId() against HGCalGeometry::getClosestCell()
// on the same set of points.

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"

#include "FastSimulation/CaloGeometryTools/interface/HGCalFastGeometry.h"

#include <chrono>
#include <map>
#include <string>
#include <vector>

class HGCalFastGeometryTester : public edm::one::EDAnalyzer<> {
public:
  explicit HGCalFastGeometryTester(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  const std::string nameSense_;
  const bool benchmark_;
  const int maxCells_;
  const edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> dddToken_;
  const edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> geomToken_;
};

HGCalFastGeometryTester::HGCalFastGeometryTester(const edm::ParameterSet& ps)
    : nameSense_(ps.getParameter<std::string>("nameSense")),
      benchmark_(ps.getParameter<bool>("benchmark")),
      maxCells_(ps.getParameter<int>("maxCells")),
      dddToken_(esConsumes<HGCalDDDConstants, IdealGeometryRecord>(edm::ESInputTag{"", nameSense_})),
      geomToken_(esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", nameSense_})) {}

void HGCalFastGeometryTester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("nameSense", "HGCalEESensitive");
  desc.add<bool>("benchmark", true);
  desc.add<int>("maxCells", -1);  // -1 = all
  descriptions.add("hgcalFastGeometryTester", desc);
}

void HGCalFastGeometryTester::analyze(const edm::Event&, const edm::EventSetup& iSetup) {
  const HGCalDDDConstants& ddd = iSetup.getData(dddToken_);
  const HGCalGeometry& geom = iSetup.getData(geomToken_);

  const DetId::Detector det = (nameSense_ == "HGCalEESensitive")          ? DetId::HGCalEE
                              : (nameSense_ == "HGCalHESiliconSensitive") ? DetId::HGCalHSi
                                                                          : DetId::HGCalHSc;

  auto t0 = std::chrono::steady_clock::now();
  HGCalFastGeometry fast(ddd, det);
  auto t1 = std::chrono::steady_clock::now();
  const double buildMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

  edm::LogVerbatim("HGCalFastGeom") << "[" << nameSense_ << "] built in " << buildMs << " ms; layers "
                                    << fast.firstLayer() << ".." << fast.lastLayer() << " ("
                                    << fast.nLayers() << "), cached wafers " << fast.nCachedWafers();

  const std::vector<DetId>& ids = geom.getValidDetIds();
  edm::LogVerbatim("HGCalFastGeom") << "valid DetIds: " << ids.size();

  unsigned nTested = 0, nNull = 0, nMismatch = 0;
  std::map<int, unsigned> mismatchByLayer;
  std::vector<GlobalPoint> points;
  points.reserve(ids.size());

  for (const auto& id : ids) {
    if (maxCells_ > 0 && static_cast<int>(nTested) >= maxCells_)
      break;
    const GlobalPoint gp = geom.getPosition(id);
    points.push_back(gp);

    HGCSiliconDetId sid(id);
    const int layer = sid.layer();
    const int zside = sid.zside();

    const DetId back = fast.toDetId(layer, gp.x(), gp.y(), zside);
    ++nTested;
    if (back == DetId()) {
      ++nNull;
    } else if (back.rawId() != id.rawId()) {
      ++nMismatch;
      ++mismatchByLayer[layer];
      if (nMismatch <= 10) {
        HGCSiliconDetId b(back);
        edm::LogVerbatim("HGCalFastGeom")
            << "  mismatch at (" << gp.x() << ", " << gp.y() << ", " << gp.z() << ")"
            << " expected L" << layer << " w(" << sid.waferU() << "," << sid.waferV() << ")"
            << " c(" << sid.cellU() << "," << sid.cellV() << ") t" << sid.type() << " | got L" << b.layer()
            << " w(" << b.waferU() << "," << b.waferV() << ") c(" << b.cellU() << "," << b.cellV() << ") t"
            << b.type();
      }
    }
  }

  const double frac = nTested ? 100. * (nNull + nMismatch) / nTested : 0.;
  edm::LogVerbatim("HGCalFastGeom") << "RESULT [" << nameSense_ << "] tested " << nTested << "  null " << nNull
                                    << "  mismatched " << nMismatch << "  -> " << frac << "% disagree";
  for (const auto& m : mismatchByLayer)
    edm::LogVerbatim("HGCalFastGeom") << "   layer " << m.first << ": " << m.second;

  if (benchmark_ && !points.empty()) {
    // toDetId
    auto b0 = std::chrono::steady_clock::now();
    unsigned long long sink = 0;
    for (const auto& gp : points) {
      const int layer = ddd.getLayer(std::abs(gp.z()), true);
      sink += fast.toDetId(layer, gp.x(), gp.y(), (gp.z() > 0 ? 1 : -1)).rawId();
    }
    auto b1 = std::chrono::steady_clock::now();

    // getClosestCell, on a subset -- it is orders of magnitude slower
    const std::size_t nSlow = std::min<std::size_t>(points.size(), 2000);
    auto c0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < nSlow; ++i)
      sink += geom.getClosestCell(points[i]).rawId();
    auto c1 = std::chrono::steady_clock::now();

    const double fastUs = std::chrono::duration<double, std::micro>(b1 - b0).count() / points.size();
    const double slowUs = std::chrono::duration<double, std::micro>(c1 - c0).count() / nSlow;
    edm::LogVerbatim("HGCalFastGeom") << "BENCH [" << nameSense_ << "] toDetId " << fastUs << " us/call ("
                                      << points.size() << " calls); getClosestCell " << slowUs << " us/call ("
                                      << nSlow << " calls); speedup x" << (slowUs / fastUs) << " [sink " << sink
                                      << "]";
  }
}

DEFINE_FWK_MODULE(HGCalFastGeometryTester);
