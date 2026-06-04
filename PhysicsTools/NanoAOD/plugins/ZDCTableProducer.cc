// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      ZDCTableProducer
//
/**\class ZDCTableProducer ZDCTableProducer.cc PhysicsTools/NanoAOD/plugins/ZDCTableProducer.cc

 Zero-Degree-Calorimeter (ZDC) rechit table for NanoAOD.

 Ports the HiForest ZDC ntuple (HeavyIonsAnalysis/ZDCAnalysis/ZDCRecHitAnalyzerHC)
 to NanoAOD FlatTables. Produces:
   - a variable-length table ("ZDC") with one row per ZDC rechit
     (zside, section, channel, depth, energy, time, ...),
   - a singleton table ("ZDCsum") with the per-side EM+HAD energy sums,
     which are the key observables for UPC neutron (0n/Xn) classification.

 The input is a ZDCRecHitCollection, typically produced on the fly by the
 central RecoLocalCalo/HcalRecProducers ZdcHitReconstructor_Run3 from the
 hcalDigis:ZDC digis that are kept in the heavy-ion MiniAOD.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano)
//

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class ZDCTableProducer : public edm::stream::EDProducer<> {
public:
  explicit ZDCTableProducer(const edm::ParameterSet&);
  ~ZDCTableProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const std::string name_;
  const std::string doc_;
  const int precision_;
  const edm::EDGetTokenT<ZDCRecHitCollection> hitsToken_;
};

ZDCTableProducer::ZDCTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      precision_(iConfig.getParameter<int>("precision")),
      hitsToken_(consumes<ZDCRecHitCollection>(iConfig.getParameter<edm::InputTag>("src"))) {
  produces<nanoaod::FlatTable>();          // per-rechit table  -> name_
  produces<nanoaod::FlatTable>("sums");    // per-side sums      -> name_ + "sum"
}

void ZDCTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  std::vector<int> zside, section, channel, depth;
  std::vector<float> energy, time, tdcTime, cwTime, eSOIp1, rSOIp1, lowGainE;
  float sumPlus = 0.f, sumMinus = 0.f;

  edm::Handle<ZDCRecHitCollection> hits;
  iEvent.getByToken(hitsToken_, hits);
  if (hits.isValid()) {
    for (const auto& hit : *hits) {
      const HcalZDCDetId id = hit.id();
      const int sec = static_cast<int>(id.section());
      zside.push_back(id.zside());
      section.push_back(sec);
      channel.push_back(id.channel());
      depth.push_back(id.depth());
      energy.push_back(hit.energy());
      time.push_back(hit.time());
      tdcTime.push_back(hit.TDCtime());
      cwTime.push_back(hit.chargeWeightedTime());
      eSOIp1.push_back(hit.energySOIp1());
      rSOIp1.push_back(hit.ratioSOIp1());
      lowGainE.push_back(hit.lowGainEnergy());
      // per-side energy sum over the calorimetric sections (EM + HAD), as in the forest
      if (id.section() == HcalZDCDetId::EM || id.section() == HcalZDCDetId::HAD) {
        if (id.zside() > 0)
          sumPlus += hit.energy();
        else
          sumMinus += hit.energy();
      }
    }
  }

  auto tab = std::make_unique<nanoaod::FlatTable>(zside.size(), name_, /*singleton=*/false, /*extension=*/false);
  tab->setDoc(doc_);
  tab->addColumn<int>("zside", zside, "ZDC side (+1 = plus / -1 = minus)");
  tab->addColumn<int>("section", section, "ZDC section (1=EM, 2=HAD, 3=LUM, 4=RPD)");
  tab->addColumn<int>("channel", channel, "ZDC channel within the section");
  tab->addColumn<int>("depth", depth, "ZDC depth");
  tab->addColumn<float>("energy", energy, "ZDC rechit energy", precision_);
  tab->addColumn<float>("time", time, "ZDC rechit time", precision_);
  tab->addColumn<float>("TDCtime", tdcTime, "ZDC TDC time", precision_);
  tab->addColumn<float>("chargeWeightedTime", cwTime, "ZDC charge-weighted time", precision_);
  tab->addColumn<float>("energySOIp1", eSOIp1, "ZDC energy of (slice-of-interest + 1)", precision_);
  tab->addColumn<float>("ratioSOIp1", rSOIp1, "ZDC ratio E(SOI)/E(SOI+1)", precision_);
  tab->addColumn<float>("lowGainEnergy", lowGainE, "ZDC low-gain energy", precision_);
  iEvent.put(std::move(tab));

  auto sums = std::make_unique<nanoaod::FlatTable>(1, name_ + "sum", /*singleton=*/true, /*extension=*/false);
  sums->setDoc("ZDC per-side EM+HAD energy sums");
  sums->addColumnValue<float>("Plus", sumPlus, "Sum of ZDC+ EM+HAD rechit energy", precision_);
  sums->addColumnValue<float>("Minus", sumMinus, "Sum of ZDC- EM+HAD rechit energy", precision_);
  iEvent.put(std::move(sums), "sums");
}

void ZDCTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("zdcrecoRun3"))->setComment("ZDCRecHitCollection input");
  desc.add<std::string>("name", "ZDC")->setComment("name of the per-rechit flat table");
  desc.add<std::string>("doc", "ZDC rechits");
  desc.add<int>("precision", -1)->setComment("mantissa bits for float columns (-1 = full float)");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(ZDCTableProducer);
