#ifndef RecoJets_FlashJet_interface_FlashJetTiled_h
#define RecoJets_FlashJet_interface_FlashJetTiled_h

// Tiled, heap-driven variant of flashjet::clusterEvent: same generalized-kt
// sequential recombination and the same merge history, in O(n log n) instead
// of O(n^2).  It mirrors the TILED + HEAP path of the FlashJet C++ kernel
// (flashjet/_cpu_kernel.cpp), which in turn mirrors FastJet's N2MinHeapTiled:
//
//   * a pair exactly at dR = R is left to the beam, as in FastJet and in
//     FlashJetCore.h;
//   * slots live in a (rapidity, phi) grid whose cells are at least R across,
//     so every neighbour within R sits in the slot's own cell or the 8 around
//     it.  A realized merge always has dR < R -- the winning pair (a, b) has
//     w_a = min(w_a, w_b), so d_ab = w_a dR^2 / R^2 only beats d_aB = w_a when
//     dR < R -- so a slot with no neighbour within R can be given no nearest
//     neighbour at all: its candidate is then its (exact) beam distance.
//   * the winning slot comes off an indexed binary heap keyed on
//     (candidate, slot), so a step costs O(log n) instead of a scan over n.
//   * a reverse index lists the slots pointing at a given slot, so a merge
//     invalidates only the rows it must.
//
// Everything is preallocated by the caller (TiledScratch), as in
// FlashJetCore.h: no allocation, no standard containers.  Ties break exactly
// as in the plain path (lowest slot index; pair-vs-beam ties to the beam), so
// both produce the same history.

#include <cmath>
#include <cstdint>

#include "RecoJets/FlashJet/interface/FlashJetCore.h"

namespace flashjet {

  // Cells per particle (the grid is kept to O(n) cells) and the integer and
  // floating-point scratch this needs, per particle, on top of Scratch.
  inline constexpr int32_t kCellsPerParticle = 3;
  inline constexpr int32_t kTiledIntScratch = 8 + kCellsPerParticle;  // 8 per slot + the cell heads
  inline constexpr int32_t kTiledFloatScratch = 0;

  struct TiledScratch {
    int32_t* cellOf;    // [n] cell of each slot, -1 when dead
    int32_t* cellNext;  // [n] intrusive list of the slots of a cell
    int32_t* cellPrev;  // [n]
    int32_t* heapNode;  // [n] heap position -> slot
    int32_t* heapPos;   // [n] slot -> heap position, -1 when removed
    int32_t* revHead;   // [n] first slot whose nearest neighbour is this one
    int32_t* revNext;   // [n] reverse-index list
    int32_t* revPrev;   // [n]
    int32_t* cellHead;  // [kCellsPerParticle * n] first slot of each cell
  };

  ALPAKA_FN_HOST_ACC inline TiledScratch makeTiledScratch(int32_t* i, int32_t n) {
    return TiledScratch{i, i + n, i + 2 * n, i + 3 * n, i + 4 * n, i + 5 * n, i + 6 * n, i + 7 * n, i + 8 * n};
  }

  // the tiled path only pays off once there are enough slots (FlashJet uses
  // the same threshold)
  inline constexpr int32_t kTiledMin = 12;

  namespace detail {

    struct Grid {
      int32_t nRap, nPhi;
      double rapMin, invDrap, invDphi;
    };

    ALPAKA_FN_HOST_ACC inline int32_t cellIndex(Grid const& g, double rap, double phi) {
      int32_t ir = static_cast<int32_t>((rap - g.rapMin) * g.invDrap);
      ir = (ir < 0) ? 0 : ((ir >= g.nRap) ? g.nRap - 1 : ir);
      int32_t ip = 0;
      if (g.nPhi > 1) {
        ip = static_cast<int32_t>(phi * g.invDphi);
        ip = (ip < 0) ? 0 : ((ip >= g.nPhi) ? g.nPhi - 1 : ip);
      }
      return ir * g.nPhi + ip;
    }

  }  // namespace detail

  // Clusters one event; same arguments and outputs as clusterEvent, plus the
  // tiled scratch.  Falls back to clusterEvent for small events.
  ALPAKA_FN_HOST_ACC inline int32_t clusterEventTiled(int32_t n,
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
                                                      Scratch const& s,
                                                      TiledScratch const& t) {
    using namespace detail;
    if (n < kTiledMin)
      return clusterEvent(
          n, R, p, inPx, inPy, inPz, inE, histP1, histP2, histChild, histD, jetIdx, jetPx, jetPy, jetPz, jetE, s);

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
      t.revHead[k] = -1;
      t.revNext[k] = -1;
      t.revPrev[k] = -1;
      s.nni[k] = -1;
    }

    // ---- the grid: cells at least R across, O(n) of them ----
    Grid g;
    {
      double lo = s.rap[0], hi = s.rap[0];
      for (int32_t k = 1; k < n; ++k) {
        lo = (s.rap[k] < lo) ? s.rap[k] : lo;
        hi = (s.rap[k] > hi) ? s.rap[k] : hi;
      }
      const double span = hi - lo;
      double drap = R;
      // a beam-like rapidity can blow the range up; widening cells is always safe
      if (span > drap * (n + 1))
        drap = span / (n + 1);
      int32_t nRap = static_cast<int32_t>(span / drap) + 1;
      nRap = (nRap < 1) ? 1 : nRap;
      int32_t nPhi = static_cast<int32_t>(kTwoPi / R);
      if (nPhi < 3)
        nPhi = 1;  // too few columns to tile in phi
      const int64_t cap = kCellsPerParticle * static_cast<int64_t>(n);
      while (static_cast<int64_t>(nRap) * nPhi > cap) {
        if (nRap >= nPhi && nRap > 1) {
          drap += drap;
          nRap = static_cast<int32_t>(span / drap) + 1;
        } else if (nPhi > 1) {
          nPhi = (nPhi >= 6) ? nPhi / 2 : 1;
        } else {
          break;
        }
      }
      g.nRap = nRap;
      g.nPhi = nPhi;
      g.rapMin = lo;
      g.invDrap = 1. / drap;
      g.invDphi = (nPhi > 1) ? nPhi / kTwoPi : 0.;
    }
    const int32_t nCells = g.nRap * g.nPhi;
    for (int32_t c = 0; c < nCells; ++c)
      t.cellHead[c] = -1;

    auto cellInsert = [&](int32_t k, int32_t c) {
      t.cellOf[k] = c;
      t.cellPrev[k] = -1;
      t.cellNext[k] = t.cellHead[c];
      if (t.cellHead[c] >= 0)
        t.cellPrev[t.cellHead[c]] = k;
      t.cellHead[c] = k;
    };
    auto cellErase = [&](int32_t k) {
      const int32_t c = t.cellOf[k];
      if (c < 0)
        return;
      if (t.cellPrev[k] >= 0)
        t.cellNext[t.cellPrev[k]] = t.cellNext[k];
      else
        t.cellHead[c] = t.cellNext[k];
      if (t.cellNext[k] >= 0)
        t.cellPrev[t.cellNext[k]] = t.cellPrev[k];
      t.cellOf[k] = -1;
      t.cellNext[k] = -1;
      t.cellPrev[k] = -1;
    };

    for (int32_t k = 0; k < n; ++k)
      cellInsert(k, cellIndex(g, s.rap[k], s.phi[k]));

    // ---- the reverse index: which slots point at a given slot ----
    auto revDetach = [&](int32_t k) {
      const int32_t a = s.nni[k];
      if (a < 0)
        return;
      if (t.revPrev[k] >= 0)
        t.revNext[t.revPrev[k]] = t.revNext[k];
      else
        t.revHead[a] = t.revNext[k];
      if (t.revNext[k] >= 0)
        t.revPrev[t.revNext[k]] = t.revPrev[k];
      t.revPrev[k] = -1;
      t.revNext[k] = -1;
    };
    auto setNN = [&](int32_t k, int32_t nn, double d) {
      revDetach(k);
      s.nnd[k] = d;
      s.nni[k] = nn;
      if (nn >= 0) {
        t.revPrev[k] = -1;
        t.revNext[k] = t.revHead[nn];
        if (t.revHead[nn] >= 0)
          t.revPrev[t.revHead[nn]] = k;
        t.revHead[nn] = k;
      }
    };

    auto candidate = [&](int32_t k, int32_t nn, double d) {
      const double wk = s.w[k];
      const double wn = (nn >= 0 && d < R2) ? s.w[nn] : kInf;
      const double pair = ((wk < wn) ? wk : wn) * d * invR2;
      return (pair < wk) ? pair : wk;
    };

    // nearest neighbour of slot k within its 3x3 cell block, capped at R
    auto rescan = [&](int32_t k) {
      double best = kInf;
      int32_t bj = -1;
      const int32_t c = t.cellOf[k];
      const int32_t ir = c / g.nPhi, ip = c % g.nPhi;
      for (int32_t dr = -1; dr <= 1; ++dr) {
        const int32_t jr = ir + dr;
        if (jr < 0 || jr >= g.nRap)
          continue;
        for (int32_t dp = -1; dp <= 1; ++dp) {
          int32_t jp = ip + dp;
          if (g.nPhi == 1) {
            if (dp != 0)
              continue;
            jp = 0;
          } else {
            jp = (jp < 0) ? jp + g.nPhi : ((jp >= g.nPhi) ? jp - g.nPhi : jp);
          }
          for (int32_t m = t.cellHead[jr * g.nPhi + jp]; m >= 0; m = t.cellNext[m]) {
            if (m == k)
              continue;
            const double d = dr2(s.rap[k], s.phi[k], s.rap[m], s.phi[m]);
            if (!(d < R2))  // no neighbour at or beyond R, as in FastJet
              continue;
            if (d <= best && (d < best || m < bj)) {
              best = d;
              bj = m;
            }
          }
        }
      }
      setNN(k, bj, best);
      s.cand[k] = candidate(k, bj, best);
    };

    for (int32_t k = 0; k < n; ++k)
      rescan(k);

    // ---- indexed binary min-heap on (cand, slot) ----
    int32_t heapSize = n;
    auto heapLess = [&](int32_t a, int32_t b) { return s.cand[a] < s.cand[b] || (s.cand[a] == s.cand[b] && a < b); };
    auto heapPlace = [&](int32_t at, int32_t slot) {
      t.heapNode[at] = slot;
      t.heapPos[slot] = at;
    };
    auto heapUp = [&](int32_t at) {
      const int32_t slot = t.heapNode[at];
      while (at > 0) {
        const int32_t parent = (at - 1) / 2;
        if (!heapLess(slot, t.heapNode[parent]))
          break;
        heapPlace(at, t.heapNode[parent]);
        at = parent;
      }
      heapPlace(at, slot);
    };
    auto heapDown = [&](int32_t at) {
      const int32_t slot = t.heapNode[at];
      for (;;) {
        int32_t child = 2 * at + 1;
        if (child >= heapSize)
          break;
        if (child + 1 < heapSize && heapLess(t.heapNode[child + 1], t.heapNode[child]))
          ++child;
        if (!heapLess(t.heapNode[child], slot))
          break;
        heapPlace(at, t.heapNode[child]);
        at = child;
      }
      heapPlace(at, slot);
    };
    auto heapUpdate = [&](int32_t slot) {
      const int32_t at = t.heapPos[slot];
      if (at < 0)
        return;
      heapUp(at);
      if (t.heapPos[slot] == at)
        heapDown(at);
    };
    auto heapErase = [&](int32_t slot) {
      const int32_t at = t.heapPos[slot];
      if (at < 0)
        return;
      const int32_t last = --heapSize;
      const int32_t moved = t.heapNode[last];
      t.heapPos[slot] = -1;
      if (at != last) {
        heapPlace(at, moved);
        heapUpdate(moved);
      }
    };
    for (int32_t k = 0; k < n; ++k)
      heapPlace(k, k);
    for (int32_t k = n / 2; k-- > 0;)
      heapDown(k);

    int32_t nextId = n;
    int32_t nJets = 0;
    for (int32_t step = 0; step < n; ++step) {
      if (heapSize <= 0)
        break;
      const int32_t i = t.heapNode[0];
      const double gbest = s.cand[i];
      const double wi = s.w[i];
      const int32_t j = s.nni[i];
      const double dpair = (j >= 0 && s.nnd[i] < R2) ? ((wi < s.w[j]) ? wi : s.w[j]) * s.nnd[i] * invR2 : kInf;
      int32_t nStale = 0;

      // the rows a merge invalidates: those pointing at i (which moves or
      // dies) or at j (which dies)
      auto takeStale = [&](int32_t target, int32_t skip) {
        for (int32_t k = t.revHead[target]; k >= 0; k = t.revNext[k]) {
          if (k != skip)
            s.stale[nStale++] = k;
        }
      };

      if (dpair < wi) {
        histP1[step] = s.ids[i];
        histP2[step] = s.ids[j];
        histChild[step] = nextId;
        histD[step] = gbest;
        takeStale(i, j);
        takeStale(j, i);

        s.px[i] += s.px[j];
        s.py[i] += s.py[j];
        s.pz[i] += s.pz[j];
        s.e[i] += s.e[j];
        // j dies
        revDetach(j);
        s.nni[j] = -1;
        s.act[j] = 0;
        s.cand[j] = kInf;
        heapErase(j);
        cellErase(j);

        double kt2;
        rapPhiKt2(s.px[i], s.py[i], s.pz[i], s.e[i], s.rap[i], s.phi[i], kt2);
        s.w[i] = weight(kt2, p);
        s.ids[i] = nextId++;
        const int32_t cNew = cellIndex(g, s.rap[i], s.phi[i]);
        if (cNew != t.cellOf[i]) {
          cellErase(i);
          cellInsert(i, cNew);
        }

        // one sweep of the new neighbourhood: the new pseudojet's own nearest
        // neighbour, and the rows that just got closer to it
        double best = kInf;
        int32_t bj = -1;
        const int32_t ir = t.cellOf[i] / g.nPhi, ip = t.cellOf[i] % g.nPhi;
        for (int32_t dr = -1; dr <= 1; ++dr) {
          const int32_t jr = ir + dr;
          if (jr < 0 || jr >= g.nRap)
            continue;
          for (int32_t dp = -1; dp <= 1; ++dp) {
            int32_t jp = ip + dp;
            if (g.nPhi == 1) {
              if (dp != 0)
                continue;
              jp = 0;
            } else {
              jp = (jp < 0) ? jp + g.nPhi : ((jp >= g.nPhi) ? jp - g.nPhi : jp);
            }
            for (int32_t k = t.cellHead[jr * g.nPhi + jp]; k >= 0; k = t.cellNext[k]) {
              if (k == i)
                continue;
              const double d = dr2(s.rap[i], s.phi[i], s.rap[k], s.phi[k]);
              if (!(d < R2))
                continue;
              if (d <= best && (d < best || k < bj)) {
                best = d;
                bj = k;
              }
              // a row that is about to be rescanned anyway is left alone
              if (s.nni[k] == i || s.nni[k] == j)
                continue;
              if (d < s.nnd[k]) {
                setNN(k, i, d);
                s.cand[k] = candidate(k, i, d);
                heapUpdate(k);
              }
            }
          }
        }
        setNN(i, bj, best);
        s.cand[i] = candidate(i, bj, best);
        heapUpdate(i);
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
        takeStale(i, i);
        revDetach(i);
        s.nni[i] = -1;
        s.act[i] = 0;
        s.cand[i] = kInf;
        heapErase(i);
        cellErase(i);
      }

      for (int32_t k = 0; k < nStale; ++k) {
        rescan(s.stale[k]);
        heapUpdate(s.stale[k]);
      }
    }

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

}  // namespace flashjet

#endif  // RecoJets_FlashJet_interface_FlashJetTiled_h
