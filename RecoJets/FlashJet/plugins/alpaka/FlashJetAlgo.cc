#include <cmath>
#include <string>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoJets/FlashJet/interface/FlashJetInputs.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "FlashJetAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  FlashJetAlgo::FlashJetAlgo(edm::ParameterSet const& config)
      : R_(config.getParameter<double>("rParam")),
        p_(::flashjet::exponentOf(config.getParameter<std::string>("jetAlgorithm"))),
        entriesPerThread_(config.getParameter<int32_t>("entriesPerThread")) {
    if (!(R_ > 0.) || !std::isfinite(R_))
      throw cms::Exception("Configuration") << "FlashJet: rParam must be finite and positive, got " << R_;
    if (entriesPerThread_ <= 0)
      throw cms::Exception("Configuration") << "FlashJet: entriesPerThread must be positive, got " << entriesPerThread_;
    auto const& sd = config.getParameter<edm::ParameterSet>("softDrop");
    softDrop_.enable = sd.getParameter<bool>("enable");
    softDrop_.zcut = sd.getParameter<double>("zcut");
    softDrop_.beta = sd.getParameter<double>("beta");
    softDrop_.R0 = sd.getParameter<double>("R0");
    if (softDrop_.enable && (!(softDrop_.R0 > 0.) || !std::isfinite(softDrop_.R0)))
      throw cms::Exception("Configuration") << "FlashJet: softDrop.R0 must be finite and positive";
  }

  void FlashJetAlgo::fillPSetDescription(edm::ParameterSetDescription& desc,
                                         std::string const& algorithm,
                                         double rParam) {
    desc.add<std::string>("jetAlgorithm", algorithm)->setComment("AntiKt, Kt or CambridgeAachen");
    desc.add<double>("rParam", rParam);
    desc.add<int32_t>("entriesPerThread", 1)
        ->setComment("CPU backends: entries per work item; GPU backends: device threads per block");
    edm::ParameterSetDescription sd;
    sd.add<bool>("enable", false)->setComment("soft drop the hardest jet of every entry");
    sd.add<double>("zcut", 0.1);
    sd.add<double>("beta", 0.);
    sd.add<double>("R0", 0.8);
    desc.add<edm::ParameterSetDescription>("softDrop", sd);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
