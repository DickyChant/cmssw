#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoJets/FlashJet/interface/FlashJetCore.h"

#include "FlashJetAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  namespace {

    // threads cooperating on one entry in the block-parallel kernel
    constexpr uint32_t kBlockThreads = 256;

    ALPAKA_FN_ACC inline ::flashjet::Scratch entryScratch(double* fscratch, int32_t* iscratch, int32_t off, int32_t n) {
      return ::flashjet::makeScratch(
          fscratch + ::flashjet::kFloatScratch * off, iscratch + ::flashjet::kIntScratch * off, n);
    }

    ALPAKA_FN_ACC inline void writeSoftDrop(flashjet::FlashJetEntrySoA::View entries,
                                            int32_t b,
                                            ::flashjet::SoftDropResult const& sd) {
      entries.groomedPx()[b] = sd.px;
      entries.groomedPy()[b] = sd.py;
      entries.groomedPz()[b] = sd.pz;
      entries.groomedE()[b] = sd.e;
      entries.zg()[b] = sd.zg;
      entries.rg()[b] = sd.rg;
      entries.nDropped()[b] = sd.nDropped;
    }

    // One entry per work item: the whole sequential recombination runs in a
    // single thread (CPU backends, where a block holds one thread anyway).
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
          const auto s = entryScratch(fscratch, iscratch, off, n);
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
          writeSoftDrop(entries, b, sd);
        }
      }
    };

    // One entry per BLOCK: the merge sequence stays serial (it is inherently
    // so), but every scan over the particles of the entry -- the argmin of the
    // candidate distances, the nearest-neighbour sweep after a merge, and the
    // rescans of the rows a merge invalidated -- is shared by the threads of
    // the block.  Cost per merge step drops from O(n) to O(n / threads).
    //
    // This is the layout of the Triton kernels this was ported from (one
    // program per event, vector operations over the particles), and it gives
    // exactly the same merge history as FlashJetKernel: every reduction breaks
    // ties towards the lowest slot index, and the distances are computed by
    // the same expressions.
    class FlashJetBlockKernel {
    public:
      ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                    flashjet::FlashJetParticleSoA::View particles,
                                    flashjet::FlashJetEntrySoA::View entries,
                                    double R,
                                    double p,
                                    FlashJetAlgo::SoftDrop softDrop,
                                    double* fscratch,
                                    int32_t* iscratch) const {
        using namespace ::flashjet::detail;

        // block-wide reduction buffers and the state the threads share
        auto& redDist = alpaka::declareSharedVar<double[kBlockThreads], __COUNTER__>(acc);
        auto& redIdx = alpaka::declareSharedVar<int32_t[kBlockThreads], __COUNTER__>(acc);
        auto& iSel = alpaka::declareSharedVar<int32_t, __COUNTER__>(acc);
        auto& jSel = alpaka::declareSharedVar<int32_t, __COUNTER__>(acc);
        auto& isPair = alpaka::declareSharedVar<int32_t, __COUNTER__>(acc);
        auto& nextId = alpaka::declareSharedVar<int32_t, __COUNTER__>(acc);
        auto& nJets = alpaka::declareSharedVar<int32_t, __COUNTER__>(acc);

        const int32_t nEntries = entries.metadata().size();
        const uint32_t tid = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
        const uint32_t threads = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
        const uint32_t blocks = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];

        for (int32_t b = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]; b < nEntries; b += blocks) {
          const int32_t off = entries.offset()[b];
          const int32_t n = entries.size()[b];
          const auto s = entryScratch(fscratch, iscratch, off, n);
          double const* inPx = particles.px().data() + off;
          double const* inPy = particles.py().data() + off;
          double const* inPz = particles.pz().data() + off;
          double const* inE = particles.e().data() + off;
          int32_t* histP1 = particles.histP1().data() + off;
          int32_t* histP2 = particles.histP2().data() + off;
          int32_t* histChild = particles.histChild().data() + off;
          double* histD = particles.histD().data() + off;
          double* jetPx = particles.jetPx().data() + off;
          double* jetPy = particles.jetPy().data() + off;
          double* jetPz = particles.jetPz().data() + off;
          double* jetE = particles.jetE().data() + off;
          const double invR2 = 1. / (R * R);

          if (n <= 0) {
            if (tid == 0) {
              entries.nJets()[b] = 0;
              writeSoftDrop(entries, b, ::flashjet::SoftDropResult{0., 0., 0., 0., 0., 0., 0});
            }
            alpaka::syncBlockThreads(acc);
            continue;
          }

          auto candidate = [&](int32_t k, int32_t nn, double d) {
            const double wk = s.w[k];
            const double wn = (nn >= 0) ? s.w[nn] : kInf;
            const double pair = ((wk < wn) ? wk : wn) * d * invR2;
            return (pair < wk) ? pair : wk;
          };

          // geometric nearest neighbour of slot k over all active slots
          auto rescan = [&](int32_t k) {
            double best = kInf;
            int32_t bj = -1;
            for (int32_t m = 0; m < n; ++m) {
              if (!s.act[m] || m == k)
                continue;
              const double d = dr2(s.rap[k], s.phi[k], s.rap[m], s.phi[m]);
              if (d <= best && (d < best || m < bj)) {
                best = d;
                bj = m;
              }
            }
            s.nnd[k] = best;
            s.nni[k] = bj;
            s.cand[k] = candidate(k, bj, best);
          };

          // block reduction of (distance, slot), ties to the lowest slot
          auto reduceMin = [&](double value, int32_t index) {
            redDist[tid] = value;
            redIdx[tid] = index;
            alpaka::syncBlockThreads(acc);
            for (uint32_t half = threads / 2; half > 0; half /= 2) {
              if (tid < half) {
                const double other = redDist[tid + half];
                const int32_t otherIdx = redIdx[tid + half];
                if (otherIdx >= 0 &&
                    (redIdx[tid] < 0 || other < redDist[tid] || (other == redDist[tid] && otherIdx < redIdx[tid]))) {
                  redDist[tid] = other;
                  redIdx[tid] = otherIdx;
                }
              }
              alpaka::syncBlockThreads(acc);
            }
          };

          for (int32_t k = tid; k < n; k += threads) {
            s.px[k] = inPx[k];
            s.py[k] = inPy[k];
            s.pz[k] = inPz[k];
            s.e[k] = inE[k];
            s.ids[k] = k;
            s.act[k] = 1;
            s.stale[k] = 0;
            double kt2;
            rapPhiKt2(s.px[k], s.py[k], s.pz[k], s.e[k], s.rap[k], s.phi[k], kt2);
            s.w[k] = weight(kt2, p);
          }
          alpaka::syncBlockThreads(acc);

          for (int32_t k = tid; k < n; k += threads)
            rescan(k);
          if (tid == 0) {
            nextId = n;
            nJets = 0;
          }
          alpaka::syncBlockThreads(acc);

          for (int32_t step = 0; step < n; ++step) {
            // the merge with the smallest candidate distance
            double bestCand = kInf;
            int32_t bestSlot = -1;
            for (int32_t k = tid; k < n; k += threads) {
              if (s.cand[k] < bestCand) {
                bestCand = s.cand[k];
                bestSlot = k;
              }
            }
            reduceMin(bestCand, bestSlot);

            // the merge itself: bookkeeping on a single thread
            if (tid == 0) {
              const int32_t i = redIdx[0];
              iSel = i;
              isPair = 0;
              jSel = i;
              if (i >= 0) {
                const double gbest = s.cand[i];
                const double wi = s.w[i];
                const int32_t j = s.nni[i];
                const double dpair = (j >= 0) ? ((wi < s.w[j]) ? wi : s.w[j]) * s.nnd[i] * invR2 : kInf;
                if (dpair < wi) {
                  isPair = 1;
                  jSel = j;
                  histP1[step] = s.ids[i];
                  histP2[step] = s.ids[j];
                  histChild[step] = nextId;
                  histD[step] = gbest;
                  s.px[i] += s.px[j];
                  s.py[i] += s.py[j];
                  s.pz[i] += s.pz[j];
                  s.e[i] += s.e[j];
                  s.act[j] = 0;
                  s.cand[j] = kInf;
                  double kt2;
                  rapPhiKt2(s.px[i], s.py[i], s.pz[i], s.e[i], s.rap[i], s.phi[i], kt2);
                  s.w[i] = weight(kt2, p);
                  s.ids[i] = nextId++;
                } else {
                  histP1[step] = s.ids[i];
                  histP2[step] = -1;
                  histChild[step] = -1;
                  histD[step] = gbest;
                  jetPx[nJets] = s.px[i];
                  jetPy[nJets] = s.py[i];
                  jetPz[nJets] = s.pz[i];
                  jetE[nJets] = s.e[i];
                  ++nJets;
                  s.act[i] = 0;
                  s.cand[i] = kInf;
                }
              }
            }
            alpaka::syncBlockThreads(acc);
            const int32_t i = iSel;
            if (i < 0)
              break;  // unreachable while any slot is live
            const int32_t j = jSel;
            const bool pair = isPair != 0;

            // one sweep over the entry: the new pseudojet's nearest neighbour,
            // the rows that improved towards it, and the rows it invalidated
            double bestD = kInf;
            int32_t bestJ = -1;
            for (int32_t k = tid; k < n; k += threads) {
              if (!s.act[k] || k == i)
                continue;
              const double d = dr2(s.rap[i], s.phi[i], s.rap[k], s.phi[k]);
              if (d < bestD || (d == bestD && k < bestJ)) {
                bestD = d;
                bestJ = k;
              }
              if (s.nni[k] == i || s.nni[k] == j) {
                s.stale[k] = 1;
                continue;
              }
              if (pair && d < s.nnd[k]) {
                s.nnd[k] = d;
                s.nni[k] = i;
                s.cand[k] = candidate(k, i, d);
              }
            }
            reduceMin(bestD, bestJ);
            if (tid == 0 && pair) {
              const int32_t bj = redIdx[0];
              const double best = (bj >= 0) ? redDist[0] : kInf;
              s.nnd[i] = best;
              s.nni[i] = bj;
              s.cand[i] = candidate(i, bj, best);
            }
            alpaka::syncBlockThreads(acc);

            for (int32_t k = tid; k < n; k += threads) {
              if (s.stale[k]) {
                rescan(k);
                s.stale[k] = 0;
              }
            }
            alpaka::syncBlockThreads(acc);
          }

          // particle -> jet and the optional soft drop: O(n) walks over the history
          if (tid == 0) {
            int32_t jet = nJets;
            for (int32_t step = n - 1; step >= 0; --step) {
              if (histP2[step] < 0) {
                s.jetOf[histP1[step]] = --jet;
              } else {
                const int32_t owner = s.jetOf[histChild[step]];
                s.jetOf[histP1[step]] = owner;
                s.jetOf[histP2[step]] = owner;
              }
            }
            entries.nJets()[b] = nJets;
            ::flashjet::SoftDropResult sd{0., 0., 0., 0., 0., 0., 0};
            if (softDrop.enable)
              sd = ::flashjet::softDrop(n,
                                        nJets,
                                        softDrop.zcut,
                                        softDrop.beta,
                                        softDrop.R0,
                                        inPx,
                                        inPy,
                                        inPz,
                                        inE,
                                        histP1,
                                        histP2,
                                        histChild,
                                        jetPx,
                                        jetPy,
                                        s);
            writeSoftDrop(entries, b, sd);
          }
          alpaka::syncBlockThreads(acc);
          for (int32_t k = tid; k < n; k += threads)
            particles.jetIdx().data()[off + k] = s.jetOf[k];
          alpaka::syncBlockThreads(acc);
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
    if constexpr (requires_single_thread_per_block_v<Acc1D>) {
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
    } else {
      // one block per entry, up to the number of entries
      auto workDiv = make_workdiv<Acc1D>(nEntries, kBlockThreads);
      alpaka::exec<Acc1D>(queue,
                          workDiv,
                          FlashJetBlockKernel{},
                          collection.view().particles(),
                          collection.view().entries(),
                          R_,
                          p_,
                          softDrop_,
                          fscratch.data(),
                          iscratch.data());
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
