#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace pat {

  class PATCountFilter : public edm::global::EDFilter<> {
  public:
    explicit PATCountFilter(const edm::ParameterSet& iConfig);
    ~PATCountFilter() override;

  private:
    bool filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  private:
    const edm::EDGetTokenT<edm::View<Electron> > electronToken_;
    const edm::EDGetTokenT<edm::View<Muon> > muonToken_;
    const edm::EDGetTokenT<edm::View<Tau> > tauToken_;
    const edm::EDGetTokenT<edm::View<Photon> > photonToken_;
    const edm::EDGetTokenT<edm::View<reco::Jet> > jetToken_;
    const edm::EDGetTokenT<reco::JetTagCollection> jetTagToken_;
    const double jetMinPt_;
    const bool countElectrons_;
    const bool countMuons_;
    const bool countTaus_;
    const unsigned int minNumber_;
    const unsigned int maxNumber_;
    const unsigned int minJets_;
  };

}  // namespace pat

using namespace pat;

PATCountFilter::PATCountFilter(const edm::ParameterSet& iConfig)
    : electronToken_(mayConsume<edm::View<Electron> >(iConfig.getParameter<edm::InputTag>("electronSource"))),
      muonToken_(mayConsume<edm::View<Muon> >(iConfig.getParameter<edm::InputTag>("muonSource"))),
      tauToken_(mayConsume<edm::View<Tau> >(iConfig.getParameter<edm::InputTag>("tauSource"))),
      photonToken_(mayConsume<edm::View<Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("photonSource", {}))),
      jetToken_(mayConsume<edm::View<reco::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetSource", {}))),
      jetTagToken_(mayConsume<reco::JetTagCollection>(iConfig.getUntrackedParameter<edm::InputTag>("jetTagSource", {}))),
      jetMinPt_(iConfig.getUntrackedParameter<double>("jetMinPt", 0.)),
      countElectrons_(iConfig.getParameter<bool>("countElectrons")),
      countMuons_(iConfig.getParameter<bool>("countMuons")),
      countTaus_(iConfig.getParameter<bool>("countTaus")),
      minNumber_(iConfig.getParameter<unsigned int>("minNumber")),
      maxNumber_(iConfig.getParameter<unsigned int>("maxNumber")),
      minJets_(iConfig.getUntrackedParameter<unsigned int>("minJets", 0)) {}

PATCountFilter::~PATCountFilter() {}

bool PATCountFilter::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<edm::View<Electron> > electrons;
  if (countElectrons_)
    iEvent.getByToken(electronToken_, electrons);
  edm::Handle<edm::View<Muon> > muons;
  if (countMuons_)
    iEvent.getByToken(muonToken_, muons);
  edm::Handle<edm::View<Tau> > taus;
  if (countTaus_)
    iEvent.getByToken(tauToken_, taus);
  unsigned int nrLeptons = 0;
  nrLeptons += (countElectrons_ ? electrons->size() : 0);
  nrLeptons += (countMuons_ ? muons->size() : 0);
  nrLeptons += (countTaus_ ? taus->size() : 0);
  bool passLepton = nrLeptons >= minNumber_ && nrLeptons <= maxNumber_;
  const auto& photons = iEvent.getHandle(photonToken_);
  bool passPhoton = photons.isValid() ? photons->size() >= minNumber_ : false;
  bool passJets(false);
  if (minJets_ > 0) {
    size_t nJets(0);
    const auto& jetTag = iEvent.getHandle(jetTagToken_);
    for (const auto& jet : iEvent.get(jetToken_)) {
      auto jetPt = jet.pt();
      if (jetTag.isValid()) {
        double ptCorr(-999.), maxDR2(0.09);
        for (const auto& t : *jetTag) {
          const auto& dR2 = reco::deltaR2(jet, *(t.first));
          if (dR2 < maxDR2) {
            maxDR2 = dR2;
            ptCorr = t.second;
          }
        }
        if (std::isnan(ptCorr))
          ptCorr = -999.;
        jetPt *= ptCorr;
      }
      if (jetPt >= jetMinPt_ && std::abs(jet.eta()) <= 2.1)
        nJets += 1;
    }
    passJets = nJets >= minJets_;
  }
  return (passLepton || passPhoton || passJets);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATCountFilter);
