#ifndef RecoJets_FlashJet_interface_FlashJetInputs_h
#define RecoJets_FlashJet_interface_FlashJetInputs_h

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace flashjet {

  // generalized-kt exponent p for the FastjetJetProducer algorithm names
  inline double exponentOf(std::string const& algorithm) {
    if (algorithm == "AntiKt")
      return -1.;
    if (algorithm == "Kt")
      return 1.;
    if (algorithm == "CambridgeAachen")
      return 0.;
    throw cms::Exception("Configuration")
        << "FlashJet: unsupported jetAlgorithm '" << algorithm << "' (use AntiKt, Kt or CambridgeAachen)";
  }

  // indices of the candidates that enter the clustering, in View order
  inline std::vector<int32_t> selectInputs(edm::View<reco::Candidate> const& cands, double ptMin) {
    std::vector<int32_t> selected;
    selected.reserve(cands.size());
    for (size_t k = 0; k < cands.size(); ++k) {
      auto const& c = cands[k];
      if (!std::isfinite(c.px()) || !std::isfinite(c.py()) || !std::isfinite(c.pz()) || !std::isfinite(c.energy()))
        continue;
      if (c.pt() < ptMin)
        continue;
      selected.push_back(static_cast<int32_t>(k));
    }
    return selected;
  }

}  // namespace flashjet

#endif  // RecoJets_FlashJet_interface_FlashJetInputs_h
