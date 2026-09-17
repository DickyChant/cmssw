#ifndef RecoJets_FlashJet_interface_FlashJetValueMaps_h
#define RecoJets_FlashJet_interface_FlashJetValueMaps_h

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ProducesCollector.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

namespace flashjet {

  // Soft-drop observables stored as ValueMap<float> keyed by the source jets;
  // jets that were not groomed get -1.
  class SoftDropValueMaps {
  public:
    static constexpr std::array<const char*, 5> names = {{"mass", "pt", "zg", "rg", "nDropped"}};
    enum { kMass, kPt, kZg, kRg, kNDropped };

    explicit SoftDropValueMaps(edm::ProducesCollector producer) {
      for (size_t k = 0; k < names.size(); ++k)
        tokens_[k] = producer.produces<edm::ValueMap<float>>(names[k]);
    }

    static std::array<std::vector<float>, 5> make(size_t nJets) {
      std::array<std::vector<float>, 5> values;
      for (auto& v : values)
        v.assign(nJets, -1.f);
      return values;
    }

    void put(edm::Event& event,
             edm::Handle<edm::View<reco::Jet>> const& jets,
             std::array<std::vector<float>, 5> const& values) const {
      for (size_t k = 0; k < names.size(); ++k) {
        edm::ValueMap<float> map;
        edm::ValueMap<float>::Filler filler(map);
        filler.insert(jets, values[k].begin(), values[k].end());
        filler.fill();
        event.emplace(tokens_[k], std::move(map));
      }
    }

  private:
    std::array<edm::EDPutTokenT<edm::ValueMap<float>>, 5> tokens_;
  };

}  // namespace flashjet

#endif  // RecoJets_FlashJet_interface_FlashJetValueMaps_h
