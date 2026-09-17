#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/FlashJet/interface/FlashJetHostCollection.h"
#include "DataFormats/FlashJet/interface/alpaka/FlashJetDeviceCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "FlashJetAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Reclusters the constituents of every jet of `src` (above jetPtMin) in a
  // single batched kernel call: one entry per jet, entry.source = index of the
  // jet in `src`, particle candIdx = daughter index within the jet.  With
  // softDrop.enable the hardest reclustered jet of each entry is groomed on
  // the device as well (see FlashJetSoftDropProducer for the host side).
  class FlashJetReclusterProducer : public global::EDProducer<> {
  public:
    FlashJetReclusterProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          srcToken_{consumes(config.getParameter<edm::InputTag>("src"))},
          putToken_{produces()},
          jetPtMin_{config.getParameter<double>("jetPtMin")},
          algo_{config} {}

    void produce(edm::StreamID, device::Event& event, device::EventSetup const&) const override {
      auto const& jets = event.get(srcToken_);
      std::vector<int32_t> selected;
      int32_t nParticles = 0;
      int32_t maxEntry = 0;
      for (size_t j = 0; j < jets.size(); ++j) {
        auto const& jet = jets[j];
        if (jet.pt() < jetPtMin_ || jet.numberOfDaughters() == 0)
          continue;
        selected.push_back(j);
        nParticles += jet.numberOfDaughters();
        maxEntry = std::max(maxEntry, static_cast<int32_t>(jet.numberOfDaughters()));
      }
      const int32_t nEntries = selected.size();

      ::flashjet::FlashJetHostCollection host{event.queue(), nParticles, nEntries};
      auto particles = host.view().particles();
      auto entries = host.view().entries();
      int32_t row = 0;
      for (int32_t b = 0; b < nEntries; ++b) {
        auto const& jet = jets[selected[b]];
        const int32_t n = jet.numberOfDaughters();
        entries.offset()[b] = row;
        entries.size()[b] = n;
        entries.source()[b] = selected[b];
        entries.nJets()[b] = 0;
        for (int32_t k = 0; k < n; ++k, ++row) {
          auto const* c = jet.daughter(k);
          particles.px()[row] = c->px();
          particles.py()[row] = c->py();
          particles.pz()[row] = c->pz();
          particles.e()[row] = c->energy();
          particles.candIdx()[row] = k;
        }
      }

      flashjet::FlashJetDeviceCollection device{event.queue(), nParticles, nEntries};
      alpaka::memcpy(event.queue(), device.buffer(), host.const_buffer());
      algo_.cluster(event.queue(), device, maxEntry);
      event.emplace(putToken_, std::move(device));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src", edm::InputTag("ak8PFJetsPuppi"));
      desc.add<double>("jetPtMin", 0.)->setComment("only recluster jets with at least this pt");
      // C/A with FastJet's max_allowable_R: every jet reclusters into one, as
      // in fastjet::contrib::Recluster
      FlashJetAlgo::fillPSetDescription(desc, "CambridgeAachen", 1000.);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    const edm::EDGetTokenT<edm::View<reco::Jet>> srcToken_;
    const device::EDPutToken<flashjet::FlashJetDeviceCollection> putToken_;
    const double jetPtMin_;
    const FlashJetAlgo algo_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(FlashJetReclusterProducer);
