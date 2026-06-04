// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      CentralityTableProducer
//
/**\class CentralityTableProducer CentralityTableProducer.cc PhysicsTools/NanoAOD/plugins/CentralityTableProducer.cc

 Event-level heavy-ion event-activity table for NanoAOD.

 Tabulates the contents of a reco::Centrality object (the same information
 stored by the HiForest HiEvtAnalyzer: HF tower/hit sums, ECAL sums,
 pixel/track multiplicities and ZDC sums) into a singleton nanoaod::FlatTable,
 optionally together with the integer centrality bin produced by
 RecoHI/HiCentralityAlgos/CentralityBinProducer.

 For UPC running the centrality *bin* is usually not defined: leave `srcBin`
 with an empty label and the `hiBin` column is filled with -1 while the HF/ZDC
 event-activity columns are still written. This is the recommended UPC setup.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano)
//

#include <memory>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class CentralityTableProducer : public edm::stream::EDProducer<> {
public:
  explicit CentralityTableProducer(const edm::ParameterSet&);
  ~CentralityTableProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const std::string name_;
  const std::string doc_;
  const int precision_;
  const edm::EDGetTokenT<reco::Centrality> centralityToken_;
  bool hasBin_;
  edm::EDGetTokenT<int> binToken_;
};

CentralityTableProducer::CentralityTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      precision_(iConfig.getParameter<int>("precision")),
      centralityToken_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("src"))) {
  const auto binTag = iConfig.getParameter<edm::InputTag>("srcBin");
  hasBin_ = !binTag.label().empty();
  if (hasBin_)
    binToken_ = consumes<int>(binTag);
  produces<nanoaod::FlatTable>();
}

void CentralityTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  auto out = std::make_unique<nanoaod::FlatTable>(1, name_, /*singleton=*/true, /*extension=*/false);
  out->setDoc(doc_);

  int hiBin = -1;
  if (hasBin_) {
    edm::Handle<int> binHandle;
    iEvent.getByToken(binToken_, binHandle);
    if (binHandle.isValid())
      hiBin = *binHandle;
  }
  out->addColumnValue<int>("hiBin", hiBin, "Centrality bin (0-199, 0.5% bins of HFtowers); -1 if not computed (e.g. UPC)");

  edm::Handle<reco::Centrality> cenHandle;
  iEvent.getByToken(centralityToken_, cenHandle);
  const bool ok = cenHandle.isValid();
  const int p = precision_;

  // member-function pointers: the reco::Centrality method is only evaluated when
  // the handle is valid, so we never read an uninitialized object.
  typedef double (reco::Centrality::*CenFun)() const;
  auto addF = [&](const char* n, CenFun m, const char* d) {
    out->addColumnValue<float>(n, ok ? static_cast<float>((cenHandle.product()->*m)()) : -999.f, d, p);
  };
  auto addI = [&](const char* n, CenFun m, const char* d) {
    out->addColumnValue<int>(n, ok ? static_cast<int>((cenHandle.product()->*m)()) : -1, d);
  };

  // HF (the key UPC / event-activity observables)
  addF("hiHF", &reco::Centrality::EtHFtowerSum, "Sum of Et in HF towers");
  addF("hiHFplus", &reco::Centrality::EtHFtowerSumPlus, "Sum of Et in HF+ towers");
  addF("hiHFminus", &reco::Centrality::EtHFtowerSumMinus, "Sum of Et in HF- towers");
  addF("hiHFECut", &reco::Centrality::EtHFtowerSumECut, "Sum of Et in HF towers (energy cut)");
  addF("hiHFhit", &reco::Centrality::EtHFhitSum, "Sum of E in HF rechits");
  addF("hiHFhitPlus", &reco::Centrality::EtHFhitSumPlus, "Sum of E in HF+ rechits");
  addF("hiHFhitMinus", &reco::Centrality::EtHFhitSumMinus, "Sum of E in HF- rechits");
  // ECAL
  addF("hiEB", &reco::Centrality::EtEBSum, "Sum of Et in EB");
  addF("hiEE", &reco::Centrality::EtEESum, "Sum of Et in EE");
  addF("hiET", &reco::Centrality::EtMidRapiditySum, "Sum of Et at mid-rapidity");
  // ZDC
  addF("hiZDC", &reco::Centrality::zdcSum, "Sum of E in ZDC");
  addF("hiZDCplus", &reco::Centrality::zdcSumPlus, "Sum of E in ZDC+");
  addF("hiZDCminus", &reco::Centrality::zdcSumMinus, "Sum of E in ZDC-");
  // Multiplicities
  addI("hiNpix", &reco::Centrality::multiplicityPixel, "Number of pixel hits");
  addI("hiNpixelTracks", &reco::Centrality::NpixelTracks, "Number of pixel tracks");
  addI("hiNtracks", &reco::Centrality::Ntracks, "Number of tracks");

  iEvent.put(std::move(out));
}

void CentralityTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hiCentrality"))
      ->setComment("reco::Centrality collection (event activity / HF / ZDC sums)");
  desc.add<edm::InputTag>("srcBin", edm::InputTag("centralityBin", "HFtowers"))
      ->setComment("int centrality bin; leave the label empty to skip (e.g. UPC)");
  desc.add<std::string>("name", "Cent")->setComment("name of the flat table / branch prefix");
  desc.add<std::string>("doc", "Heavy-ion event-activity (centrality) information");
  desc.add<int>("precision", 10)->setComment("mantissa bits for float columns");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(CentralityTableProducer);
