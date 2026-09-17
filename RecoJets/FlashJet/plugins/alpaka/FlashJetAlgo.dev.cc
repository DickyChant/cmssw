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
                                    flashjet::FlashJetParticleSoA::View particles,
                                    flashjet::FlashJetEntrySoA::View entries,
                                    double R,
                                    double p,
                                    FlashJetAlgo::SoftDrop softDrop,
                                    double* fscratch,
                                    int32_t* iscratch) const {
        for (int32_t b : uniform_elements(acc, entries.metadata().size())) {
          const int32_t off = entries.offset()[b];
          const int32_t n = entries.size()[b];
          // scratch is laid out like the particles, scaled per particle
          const auto s = ::flashjet::makeScratch(
              fscratch + ::flashjet::kFloatScratch * off, iscratch + ::flashjet::kIntScratch * off, n);
          const int32_t nJets = ::flashjet::clusterEvent(n,
                                                         R,
                                                         p,
                                                         particles.px().data() + off,
                                                         particles.py().data() + off,
                                                         particles.pz().data() + off,
                                                         particles.e().data() + off,
                                                         particles.histP1().data() + off,
                                                         particles.histP2().data() + off,
                                                         particles.histChild().data() + off,
                                                         particles.histD().data() + off,
                                                         particles.jetIdx().data() + off,
                                                         particles.jetPx().data() + off,
                                                         particles.jetPy().data() + off,
                                                         particles.jetPz().data() + off,
                                                         particles.jetE().data() + off,
                                                         s);
          entries.nJets()[b] = nJets;
          ::flashjet::SoftDropResult sd{0., 0., 0., 0., 0., 0., 0};
          if (softDrop.enable) {
            sd = ::flashjet::softDrop(n,
                                      nJets,
                                      softDrop.zcut,
                                      softDrop.beta,
                                      softDrop.R0,
                                      particles.px().data() + off,
                                      particles.py().data() + off,
                                      particles.pz().data() + off,
                                      particles.e().data() + off,
                                      particles.histP1().data() + off,
                                      particles.histP2().data() + off,
                                      particles.histChild().data() + off,
                                      particles.jetPx().data() + off,
                                      particles.jetPy().data() + off,
                                      s);
          }
          entries.groomedPx()[b] = sd.px;
          entries.groomedPy()[b] = sd.py;
          entries.groomedPz()[b] = sd.pz;
          entries.groomedE()[b] = sd.e;
          entries.zg()[b] = sd.zg;
          entries.rg()[b] = sd.rg;
          entries.nDropped()[b] = sd.nDropped;
        }
      }
    };
  }  // namespace

  void FlashJetAlgo::cluster(Queue& queue, flashjet::FlashJetDeviceCollection& collection) const {
    const int32_t nParticles = collection.view().particles().metadata().size();
    const int32_t nEntries = collection.view().entries().metadata().size();
    if (nParticles == 0 || nEntries == 0)
      return;
    auto fscratch = make_device_buffer<double[]>(queue, ::flashjet::kFloatScratch * nParticles);
    auto iscratch = make_device_buffer<int32_t[]>(queue, ::flashjet::kIntScratch * nParticles);
    const uint32_t items = entriesPerThread_;
    auto workDiv = make_workdiv<Acc1D>(divide_up_by(nEntries, items), items);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        FlashJetKernel{},
                        collection.view().particles(),
                        collection.view().entries(),
                        R_,
                        p_,
                        softDrop_,
                        fscratch.data(),
                        iscratch.data());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
