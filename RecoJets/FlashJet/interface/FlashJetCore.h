#ifndef RecoJets_FlashJet_interface_FlashJetCore_h
#define RecoJets_FlashJet_interface_FlashJetCore_h

// Generalized-kt sequential recombination (anti-kt p=-1, C/A p=0, kt p=1),
// E-scheme, with FastJet's N2Plain nearest-neighbour strategy.
//
// This is a port of the FlashJet C++ kernel (flashjet/_cpu_kernel.cpp, the
// "plain" path; https://github.com/jet-universe/FlashJet) to a single
// allocation-free function that runs unchanged on the host and inside an
// alpaka kernel.  All working memory is passed in by the caller.
//
//   * each live slot caches its GEOMETRIC nearest neighbour (min dR^2); by the
//     Cacciari-Salam lemma the global d_ij minimum is realized at some slot's
//     geometric NN, so the per-slot candidate is
//         min(w_i, w_nn) * dR2_nn / R^2   vs the beam distance w_i;
//   * a merge overwrites slot i and kills slot j; only rows that pointed at i
//     or j are rescanned, so the cost is O(n^2) per event.
//
// Ties follow the FlashJet kernel: every argmin keeps the FIRST minimum
// (lowest slot index), and pair-vs-beam ties go to the beam.  An incremental
// update keeps the neighbour it has when a new distance ties it exactly, as
// FlashJet and FastJet both do; only a full rescan applies the lowest-index
// rule among equals.
//
// A pair is only a candidate below R: FastJet records no nearest neighbour at
// dist == R^2 (ClusterSequence.hh), so a pair exactly at dR = R is left to the
// beam.  Upstream FlashJet merges it instead when the softer particle has the
// lower slot index; this port follows FastJet.

#include <cmath>
#include <cstdint>
#include <limits>

#include <alpaka/alpaka.hpp>

namespace flashjet {

  namespace detail {
    inline constexpr double kTwoPi = 6.283185307179586;
    inline constexpr double kPi = 3.141592653589793;
    inline constexpr double kMaxRap = 1e5;
    inline constexpr double kTiny = 1e-300;
    inline constexpr double kInf = std::numeric_limits<double>::infinity();

    // FastJet's numerically stable rapidity; phi folded into [0, 2pi)
    ALPAKA_FN_HOST_ACC inline void rapPhiKt2(
        double px, double py, double pz, double e, double& rap, double& phi, double& kt2) {
      kt2 = px * px + py * py;
      phi = std::atan2(py, px);
      if (phi < 0.)
        phi += kTwoPi;
      double m2 = e * e - pz * pz - kt2;
      if (!(m2 > 0.))
        m2 = 0.;
      const double apz = std::fabs(pz);
      if (kt2 + m2 <= 0.) {
        rap = (pz >= 0.) ? kMaxRap + apz : -(kMaxRap + apz);
      } else {
        const double denom = (e + apz) * (e + apz);
        const double halfLog = 0.5 * std::log((kt2 + m2) / denom);
        rap = (pz >= 0.) ? -halfLog : halfLog;
      }
    }

    // d_iB = kt^(2p), integral exponents special-cased as in FlashJet
    ALPAKA_FN_HOST_ACC inline double weight(double kt2, double p) {
      const double k = (kt2 > kTiny) ? kt2 : kTiny;
      if (p == -1.)
        return 1. / k;
      if (p == 0.)
        return 1.;
      if (p == 1.)
        return k;
      return std::pow(k, p);
    }

    ALPAKA_FN_HOST_ACC inline double dr2(double rapA, double phiA, double rapB, double phiB) {
      double dphi = phiA - phiB;
      if (dphi > kPi)
        dphi -= kTwoPi;
      else if (dphi < -kPi)
        dphi += kTwoPi;
      const double drap = rapA - rapB;
      return drap * drap + dphi * dphi;
    }
  }  // namespace detail

  // Working memory for one event of n particles: kFloatScratch * n doubles
  // and kIntScratch * n int32 (see makeScratch).
  struct Scratch {
    double* px;       // [n] current pseudojet momenta (mutated)
    double* py;       // [n]
    double* pz;       // [n]
    double* e;        // [n]
    double* rap;      // [n]
    double* phi;      // [n]
    double* w;        // [n] beam distance kt^(2p)
    double* nnd;      // [n] geometric NN distance dR^2
    double* cand;     // [n] per-slot candidate distance
    double* pjPx;     // [2n] four-momentum of every pseudojet id (grooming)
    double* pjPy;     // [2n]
    double* pjPz;     // [2n]
    double* pjE;      // [2n]
    int32_t* nni;     // [n] geometric NN slot
    int32_t* ids;     // [n] pseudojet id held by the slot
    int32_t* act;     // [n] 1 if the slot is live
    int32_t* stale;   // [n] rows to rescan after a merge
    int32_t* jetOf;   // [2n] jet index of every pseudojet id (decode)
    int32_t* stepOf;  // [2n] history step that created a pseudojet id (grooming)
  };

  inline constexpr int32_t kFloatScratch = 17;
  inline constexpr int32_t kIntScratch = 8;

  ALPAKA_FN_HOST_ACC inline Scratch makeScratch(double* f, int32_t* i, int32_t n) {
    return Scratch{f,
                   f + n,
                   f + 2 * n,
                   f + 3 * n,
                   f + 4 * n,
                   f + 5 * n,
                   f + 6 * n,
                   f + 7 * n,
                   f + 8 * n,
                   f + 9 * n,
                   f + 11 * n,
                   f + 13 * n,
                   f + 15 * n,
                   i,
                   i + n,
                   i + 2 * n,
                   i + 3 * n,
                   i + 4 * n,
                   i + 6 * n};
  }

  // Clusters one event.
  //
  // inputs:  px, py, pz, e      [n]
  // outputs: histP1/P2/Child/D  [n] merge history, one row per step
  //          jetIdx             [n] jet index of each input particle
  //          jetPx..jetE        [n] jet four-momenta in beam-merge order
  // returns the number of jets
  ALPAKA_FN_HOST_ACC inline int32_t clusterEvent(int32_t n,
                                                 double R,
                                                 double p,
                                                 double const* inPx,
                                                 double const* inPy,
                                                 double const* inPz,
                                                 double const* inE,
                                                 int32_t* histP1,
                                                 int32_t* histP2,
                                                 int32_t* histChild,
                                                 double* histD,
                                                 int32_t* jetIdx,
                                                 double* jetPx,
                                                 double* jetPy,
                                                 double* jetPz,
                                                 double* jetE,
                                                 Scratch const& s) {
    using namespace detail;
    if (n <= 0)
      return 0;
    const double invR2 = 1. / (R * R);
    const double R2 = R * R;

    for (int32_t k = 0; k < n; ++k) {
      s.px[k] = inPx[k];
      s.py[k] = inPy[k];
      s.pz[k] = inPz[k];
      s.e[k] = inE[k];
      s.ids[k] = k;
      s.act[k] = 1;
      double kt2;
      rapPhiKt2(s.px[k], s.py[k], s.pz[k], s.e[k], s.rap[k], s.phi[k], kt2);
      s.w[k] = weight(kt2, p);
    }

    auto candidate = [&](int32_t k, int32_t nn, double d) -> double {
      const double wk = s.w[k];
      // a neighbour at or beyond R is no neighbour: only the beam is left
      const double wn = (nn >= 0 && d < R2) ? s.w[nn] : kInf;
      const double pair = ((wk < wn) ? wk : wn) * d * invR2;
      return (pair < wk) ? pair : wk;
    };

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

    for (int32_t k = 0; k < n; ++k)
      rescan(k);

    int32_t nextId = n;
    int32_t nJets = 0;
    for (int32_t step = 0; step < n; ++step) {
      double gbest = kInf;
      int32_t i = -1;
      for (int32_t k = 0; k < n; ++k) {
        if (s.cand[k] < gbest) {
          gbest = s.cand[k];
          i = k;
        }
      }
      if (i < 0)
        break;  // unreachable while any slot is live

      const double wi = s.w[i];
      const int32_t j = s.nni[i];
      const double dpair = (j >= 0 && s.nnd[i] < R2) ? ((wi < s.w[j]) ? wi : s.w[j]) * s.nnd[i] * invR2 : kInf;
      int32_t nStale = 0;

      if (dpair < wi) {
        // pair merge: i absorbs j, j dies
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

        // one sweep: the new pseudojet's NN, rows that went stale, and rows
        // whose NN improved to the new pseudojet
        double best = kInf;
        int32_t bj = -1;
        for (int32_t k = 0; k < n; ++k) {
          if (!s.act[k] || k == i)
            continue;
          const double d = dr2(s.rap[i], s.phi[i], s.rap[k], s.phi[k]);
          if (d < best || (d == best && k < bj)) {
            best = d;
            bj = k;
          }
          if (s.nni[k] == i || s.nni[k] == j) {
            s.stale[nStale++] = k;
            continue;
          }
          if (d < s.nnd[k]) {
            s.nnd[k] = d;
            s.nni[k] = i;
            s.cand[k] = candidate(k, i, d);
          }
        }
        s.nnd[i] = best;
        s.nni[i] = bj;
        s.cand[i] = candidate(i, bj, best);
      } else {
        // beam merge: slot i becomes jet nJets
        histP1[step] = s.ids[i];
        histP2[step] = -1;
        histChild[step] = -1;
        histD[step] = gbest;
        jetPx[nJets] = s.px[i];
        jetPy[nJets] = s.py[i];
        jetPz[nJets] = s.pz[i];
        jetE[nJets] = s.e[i];
        ++nJets;
        for (int32_t k = 0; k < n; ++k) {
          if (s.act[k] && s.nni[k] == i && k != i)
            s.stale[nStale++] = k;
        }
        s.act[i] = 0;
        s.cand[i] = kInf;
      }
      for (int32_t t = 0; t < nStale; ++t)
        rescan(s.stale[t]);
    }

    // decode particle -> jet by walking the history backwards: a beam merge
    // names the jet of its pseudojet, and a pair merge hands its child's jet
    // down to both parents
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
    for (int32_t k = 0; k < n; ++k)
      jetIdx[k] = s.jetOf[k];

    return nJets;
  }

  struct SoftDropResult {
    double px, py, pz, e;  // groomed jet
    double zg, rg;         // momentum sharing and opening angle of the accepted split (0 if none)
    int32_t nDropped;      // branches removed
  };

  // Soft drop (FastJet contrib SoftDrop, scalar_z symmetry measure, larger_pt
  // recursion) of the hardest jet of an event clustered by clusterEvent.
  // Meant for C/A histories: the declustering follows the merge tree, and a
  // split with z = min(pt1, pt2) / (pt1 + pt2) is kept when
  //   z > zcut * (dR12 / R0)^beta.
  ALPAKA_FN_HOST_ACC inline SoftDropResult softDrop(int32_t n,
                                                    int32_t nJets,
                                                    double zcut,
                                                    double beta,
                                                    double R0,
                                                    double const* inPx,
                                                    double const* inPy,
                                                    double const* inPz,
                                                    double const* inE,
                                                    int32_t const* histP1,
                                                    int32_t const* histP2,
                                                    int32_t const* histChild,
                                                    double const* jetPx,
                                                    double const* jetPy,
                                                    Scratch const& s) {
    using namespace detail;
    SoftDropResult result{0., 0., 0., 0., 0., 0., 0};
    if (n <= 0 || nJets <= 0)
      return result;

    // four-momentum of every pseudojet, summed in merge order
    for (int32_t k = 0; k < n; ++k) {
      s.pjPx[k] = inPx[k];
      s.pjPy[k] = inPy[k];
      s.pjPz[k] = inPz[k];
      s.pjE[k] = inE[k];
    }
    int32_t hardest = 0;
    for (int32_t j = 1; j < nJets; ++j) {
      if (jetPx[j] * jetPx[j] + jetPy[j] * jetPy[j] > jetPx[hardest] * jetPx[hardest] + jetPy[hardest] * jetPy[hardest])
        hardest = j;
    }
    int32_t node = -1;
    int32_t beams = 0;
    for (int32_t step = 0; step < n; ++step) {
      const int32_t a = histP1[step];
      if (histP2[step] < 0) {
        if (beams++ == hardest)
          node = a;
        continue;
      }
      const int32_t b = histP2[step];
      const int32_t c = histChild[step];
      s.pjPx[c] = s.pjPx[a] + s.pjPx[b];
      s.pjPy[c] = s.pjPy[a] + s.pjPy[b];
      s.pjPz[c] = s.pjPz[a] + s.pjPz[b];
      s.pjE[c] = s.pjE[a] + s.pjE[b];
      s.stepOf[c] = step;
    }
    if (node < 0)
      return result;

    const double invR0 = 1. / R0;
    while (node >= n) {
      const int32_t step = s.stepOf[node];
      int32_t a = histP1[step];
      int32_t b = histP2[step];
      double pt2a = s.pjPx[a] * s.pjPx[a] + s.pjPy[a] * s.pjPy[a];
      double pt2b = s.pjPx[b] * s.pjPx[b] + s.pjPy[b] * s.pjPy[b];
      if (pt2a < pt2b) {
        const int32_t t = a;
        a = b;
        b = t;
        const double u = pt2a;
        pt2a = pt2b;
        pt2b = u;
      }
      const double pta = std::sqrt(pt2a);
      const double ptb = std::sqrt(pt2b);
      const double z = ptb / (pta + ptb);
      double rapA, phiA, rapB, phiB, kt2;
      rapPhiKt2(s.pjPx[a], s.pjPy[a], s.pjPz[a], s.pjE[a], rapA, phiA, kt2);
      rapPhiKt2(s.pjPx[b], s.pjPy[b], s.pjPz[b], s.pjE[b], rapB, phiB, kt2);
      const double dR = std::sqrt(dr2(rapA, phiA, rapB, phiB));
      const double threshold = (beta == 0.) ? zcut : zcut * std::pow(dR * invR0, beta);
      if (z > threshold) {
        result.zg = z;
        result.rg = dR;
        break;
      }
      ++result.nDropped;
      node = a;
    }
    result.px = s.pjPx[node];
    result.py = s.pjPy[node];
    result.pz = s.pjPz[node];
    result.e = s.pjE[node];
    return result;
  }

}  // namespace flashjet

#endif  // RecoJets_FlashJet_interface_FlashJetCore_h
