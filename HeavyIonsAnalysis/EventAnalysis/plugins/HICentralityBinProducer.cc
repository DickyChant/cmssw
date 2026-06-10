// -*- C++ -*-
//
// Package:    HICentralityBinProducer
// Class:      HICentralityBinProducer
//
/**\class HICentralityBinProducer HICentralityBinProducer.cc RecoHI/HICentralityBinProducer/src/HICentralityBinProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yetkin Yilmaz
//         Created:  Thu Aug 12 05:34:11 EDT 2010
//
//

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

//
// class declaration
//

class HICentralityBinProducer : public edm::stream::EDProducer<> {
  enum VariableType {
    HFtowers = 0,
    HFtowersPlus = 1,
    HFtowersMinus = 2,
    HFtowersTrunc = 3,
    HFtowersPlusTrunc = 4,
    HFtowersMinusTrunc = 5,
    HFhits = 6,
    PixelHits = 7,
    PixelTracks = 8,
    Tracks = 9,
    EB = 10,
    EE = 11,
    ZDChitsPlus = 12,
    ZDChitsMinus = 13,
    PFhf = 14,
    PFhfPlus = 15,
    PFhfMinus = 16,
    Missing = 17
  };

public:
  explicit HICentralityBinProducer(const edm::ParameterSet&);
  ~HICentralityBinProducer() override {};

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------

  const edm::EDGetTokenT<reco::Centrality> token_;
  const std::vector<double> table_;

  const std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;
  VariableType varType_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HICentralityBinProducer::HICentralityBinProducer(const edm::ParameterSet& iConfig) :
  token_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("Centrality"))),
  table_(iConfig.getParameter<std::vector<double>>("table")),
  centralityVariable_(iConfig.getParameter<std::string>("centralityVariable")),
  varType_(Missing) {

  if (centralityVariable_ == "HFtowers")
    varType_ = HFtowers;
  if (centralityVariable_ == "HFtowersPlus")
    varType_ = HFtowersPlus;
  if (centralityVariable_ == "HFtowersMinus")
    varType_ = HFtowersMinus;
  if (centralityVariable_ == "HFtowersTrunc")
    varType_ = HFtowersTrunc;
  if (centralityVariable_ == "HFtowersPlusTrunc")
    varType_ = HFtowersPlusTrunc;
  if (centralityVariable_ == "HFtowersMinusTrunc")
    varType_ = HFtowersMinusTrunc;
  if (centralityVariable_ == "HFhits")
    varType_ = HFhits;
  if (centralityVariable_ == "PixelHits")
    varType_ = PixelHits;
  if (centralityVariable_ == "PixelTracks")
    varType_ = PixelTracks;
  if (centralityVariable_ == "Tracks")
    varType_ = Tracks;
  if (centralityVariable_ == "EB")
    varType_ = EB;
  if (centralityVariable_ == "EE")
    varType_ = EE;
  if (centralityVariable_ == "ZDChitsPlus")
    varType_ = ZDChitsPlus;
  if (centralityVariable_ == "ZDChitsMinus")
    varType_ = ZDChitsMinus;
  if (centralityVariable_ == "PFhf")
    varType_ = PFhf;
  if (centralityVariable_ == "PFhfPlus")
    varType_ = PFhfPlus;
  if (centralityVariable_ == "PFhfMinus")
    varType_ = PFhfMinus;
  if (varType_ == Missing) {
    std::string errorMessage =
        "Requested Centrality variable does not exist : " + centralityVariable_ + "\n" + "Supported variables are: \n" +
        "HFtowers HFtowersPlus HFtowersMinus HFtowersTrunc HFtowersPlusTrunc HFtowersMinusTrunc "
        "HFhits PixelHits PixelTracks Tracks EB EE ZDChitsPlus ZDChitsMinus PFhf PFhfPlus PFhfMinus" +
        "\n";
    throw cms::Exception("Configuration", errorMessage);
  }

  if (iConfig.exists("nonDefaultGlauberModel")) {
    centralityMC_ = iConfig.getParameter<std::string>("nonDefaultGlauberModel");
  }
  centralityLabel_ = centralityVariable_ + centralityMC_;

  produces<int>(centralityVariable_);
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HICentralityBinProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& input = iEvent.get(token_);

  double value = 0;
  switch (varType_) {
    case HFtowers:
      value = input.EtHFtowerSum();
      break;
    case HFtowersPlus:
      value = input.EtHFtowerSumPlus();
      break;
    case HFtowersMinus:
      value = input.EtHFtowerSumMinus();
      break;
    case HFhits:
      value = input.EtHFhitSum();
      break;
    case HFtowersTrunc:
      value = input.EtHFtruncated();
      break;
    case HFtowersPlusTrunc:
      value = input.EtHFtruncatedPlus();
      break;
    case HFtowersMinusTrunc:
      value = input.EtHFtruncatedMinus();
      break;
    case PixelHits:
      value = input.multiplicityPixel();
      break;
    case PixelTracks:
      value = input.NpixelTracks();
      break;
    case Tracks:
      value = input.Ntracks();
      break;
    case EB:
      value = input.EtEBSum();
      break;
    case EE:
      value = input.EtEESum();
      break;
    case ZDChitsPlus:
      value = input.zdcSumPlus();
      break;
    case ZDChitsMinus:
      value = input.zdcSumMinus();
      break;
    case PFhf:
      value = input.EtPFhfSum();
      break;
    case PFhfPlus:
      value = input.EtPFhfSumPlus();
      break;
    case PFhfMinus:
      value = input.EtPFhfSumMinus();
      break;
    default:
      throw cms::Exception("HICentralityBinProducer", "Centrality variable not recognized.");
  }

  int bin(-1);
  assert(table_.size() == 201);
  for(size_t i = 0; i < 200; i++)
    if(value >= table_[199 - i]) {
      bin = i;
      break;
    }

  iEvent.put(std::make_unique<int>(bin), centralityVariable_);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HICentralityBinProducer);
