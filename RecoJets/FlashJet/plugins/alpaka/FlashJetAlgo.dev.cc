#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoJets/FlashJet/interface/FlashJetCore.h"

#include "FlashJetAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  namespace {
    class FlashJetKernel {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    flashjet::FlashJetSoA::View view,
                                    double R,
                                    double p,
                                    double* fscratch,
                                    int32_t* iscratch) const {
        // a single work item: the event is clustered by one device thread
        for ([[maybe_unused]] int32_t item : uniform_elements(acc, 1)) {
          const int32_t n = view.metadata().size();
          const ::flashjet::Scratch s{fscratch,
                                      fscratch + n,
                                      fscratch + 2 * n,
                                      fscratch + 3 * n,
                                      fscratch + 4 * n,
                                      fscratch + 5 * n,
                                      fscratch + 6 * n,
                                      fscratch + 7 * n,
                                      fscratch + 8 * n,
                                      iscratch,
                                      iscratch + n,
                                      iscratch + 2 * n,
                                      iscratch + 3 * n,
                                      iscratch + 4 * n};
          view.nJets() = ::flashjet::clusterEvent(n,
                                                  R,
                                                  p,
                                                  view.px().data(),
                                                  view.py().data(),
                                                  view.pz().data(),
                                                  view.e().data(),
                                                  view.histP1().data(),
                                                  view.histP2().data(),
                                                  view.histChild().data(),
                                                  view.histD().data(),
                                                  view.jetIdx().data(),
                                                  view.jetPx().data(),
                                                  view.jetPy().data(),
                                                  view.jetPz().data(),
                                                  view.jetE().data(),
                                                  s);
        }
      }
    };
  }  // namespace

  void FlashJetAlgo::cluster(Queue& queue, flashjet::FlashJetDeviceCollection& collection) const {
    // nJets is expected to be initialised to 0 by the caller
    const int32_t n = collection->metadata().size();
    if (n == 0)
      return;
    auto fscratch = make_device_buffer<double[]>(queue, 9 * n);
    auto iscratch = make_device_buffer<int32_t[]>(queue, 6 * n);
    auto workDiv = make_workdiv<Acc1D>(1, 1);
    alpaka::exec<Acc1D>(queue, workDiv, FlashJetKernel{}, collection.view(), R_, p_, fscratch.data(), iscratch.data());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
