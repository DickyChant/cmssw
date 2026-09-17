#include <cmath>
#include <random>
#include <vector>

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// Synthetic events for benchmarks: a few collimated hard sprays on top of a
// uniform soft background, with a tunable multiplicity.  Deterministic per
// event (seeded from the event number).
class FlashJetRandomCandidateProducer : public edm::global::EDProducer<> {
public:
  explicit FlashJetRandomCandidateProducer(edm::ParameterSet const& config)
      : putToken_{produces()},
        nSoft_{config.getParameter<int>("nSoft")},
        nJets_{config.getParameter<int>("nJets")},
        nPerJet_{config.getParameter<int>("nPerJet")},
        jetPt_{config.getParameter<double>("jetPt")},
        jetWidth_{config.getParameter<double>("jetWidth")},
        etaMax_{config.getParameter<double>("etaMax")},
        seed_{config.getParameter<unsigned>("seed")} {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    std::mt19937_64 rng(seed_ ^ (event.id().event() * 0x9E3779B97F4A7C15ULL));
    std::uniform_real_distribution<double> eta(-etaMax_, etaMax_), phi(-M_PI, M_PI), unit(0., 1.);
    std::exponential_distribution<double> softPt(1. / 0.7);
    std::normal_distribution<double> spread(0., jetWidth_);
    std::vector<reco::LeafCandidate> out;
    out.reserve(nSoft_ + nJets_ * nPerJet_);
    auto add = [&out](double pt, double eta, double phi) {
      const double px = pt * std::cos(phi), py = pt * std::sin(phi), pz = pt * std::sinh(eta);
      out.emplace_back(0, reco::Particle::LorentzVector(px, py, pz, std::sqrt(px * px + py * py + pz * pz)));
    };
    for (int k = 0; k < nSoft_; ++k)
      add(0.1 + softPt(rng), eta(rng), phi(rng));
    for (int j = 0; j < nJets_; ++j) {
      const double eta0 = eta(rng) * 0.5, phi0 = phi(rng);
      for (int k = 0; k < nPerJet_; ++k) {
        // steeply falling momentum fractions, a few hard cores
        const double pt = jetPt_ * std::pow(unit(rng), 4.) * 4. / nPerJet_ + 0.2;
        add(pt, eta0 + spread(rng), phi0 + spread(rng));
      }
    }
    event.emplace(putToken_, std::move(out));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<int>("nSoft", 1500)->setComment("uniform soft particles (pileup-like)");
    desc.add<int>("nJets", 6)->setComment("hard sprays");
    desc.add<int>("nPerJet", 80);
    desc.add<double>("jetPt", 400.);
    desc.add<double>("jetWidth", 0.15);
    desc.add<double>("etaMax", 4.7);
    desc.add<unsigned>("seed", 1234);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDPutTokenT<std::vector<reco::LeafCandidate>> putToken_;
  const int nSoft_, nJets_, nPerJet_;
  const double jetPt_, jetWidth_, etaMax_;
  const unsigned seed_;
};

DEFINE_FWK_MODULE(FlashJetRandomCandidateProducer);
