// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      HIRhoMapTableProducer
//
/**\class HIRhoMapTableProducer HIRhoMapTableProducer.cc PhysicsTools/NanoAOD/plugins/HIRhoMapTableProducer.cc

 Per-eta underlying-event density (rho) table for NanoAOD.

 Ports the HiForest hiFJRhoAnalyzer content: one row per eta strip of the
 hiFJRhoProducer map (etaMin, etaMax, rho, rhom). This is the event-activity
 input for the per-eta rho correction of PF isolation and jets, as done in PbPb.
 The inputs are the std::vector<double> products mapEtaEdges (n+1 edges),
 mapToRho and mapToRhoM (n values) of hiFJRhoProducer.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano)
//

#include <algorithm>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

class HIRhoMapTableProducer : public edm::stream::EDProducer<> {
public:
  explicit HIRhoMapTableProducer(const edm::ParameterSet&);
  ~HIRhoMapTableProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const std::string name_;
  const std::string doc_;
  const int precision_;
  const edm::EDGetTokenT<std::vector<double>> etaMapToken_;
  const edm::EDGetTokenT<std::vector<double>> rhoToken_;
  const edm::EDGetTokenT<std::vector<double>> rhomToken_;
};

HIRhoMapTableProducer::HIRhoMapTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      precision_(iConfig.getParameter<int>("precision")),
      etaMapToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("etaMap"))),
      rhoToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("rho"))),
      rhomToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("rhom"))) {
  produces<nanoaod::FlatTable>();
}

void HIRhoMapTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  std::vector<float> etaMin, etaMax, rho, rhom;

  edm::Handle<std::vector<double>> etaMap, rhoMap, rhomMap;
  iEvent.getByToken(etaMapToken_, etaMap);
  iEvent.getByToken(rhoToken_, rhoMap);
  iEvent.getByToken(rhomToken_, rhomMap);
  if (etaMap.isValid() && rhoMap.isValid() && etaMap->size() > 1) {
    const unsigned int n = std::min(etaMap->size() - 1, rhoMap->size());
    for (unsigned int i = 0; i < n; ++i) {
      etaMin.push_back((*etaMap)[i]);
      etaMax.push_back((*etaMap)[i + 1]);
      rho.push_back((*rhoMap)[i]);
      rhom.push_back(rhomMap.isValid() && i < rhomMap->size() ? (*rhomMap)[i] : -1.f);
    }
  }

  auto table = std::make_unique<nanoaod::FlatTable>(etaMin.size(), name_, false, false);
  table->setDoc(doc_);
  table->addColumn<float>("etaMin", etaMin, "lower edge of the eta strip", precision_);
  table->addColumn<float>("etaMax", etaMax, "upper edge of the eta strip", precision_);
  table->addColumn<float>("rho", rho, "rho (median pt/area) in the eta strip", precision_);
  table->addColumn<float>("rhom", rhom, "rho_m in the eta strip (-1 if not available)", precision_);
  iEvent.put(std::move(table));
}

void HIRhoMapTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("etaMap", edm::InputTag("hiFJRhoProducer", "mapEtaEdges"));
  desc.add<edm::InputTag>("rho", edm::InputTag("hiFJRhoProducer", "mapToRho"));
  desc.add<edm::InputTag>("rhom", edm::InputTag("hiFJRhoProducer", "mapToRhoM"));
  desc.add<std::string>("name", "HIRho");
  desc.add<std::string>("doc", "per-eta rho map (hiFJRhoProducer)");
  desc.add<int>("precision", 10);
  descriptions.add("hiRhoMapTable", desc);
}

DEFINE_FWK_MODULE(HIRhoMapTableProducer);
