// Compares flashjet::clusterEvent (host build of the alpaka kernel core)
// with FastJet on random events: identical jet constituents and
// four-momenta equal up to rounding.

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include <catch2/catch_all.hpp>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include "RecoJets/FlashJet/interface/FlashJetCore.h"
#include "RecoJets/FlashJet/interface/FlashJetTiled.h"

namespace {

  struct Particles {
    std::vector<double> px, py, pz, e;
  };

  Particles makeEvent(int n, std::mt19937_64& rng) {
    std::exponential_distribution<double> ptDist(0.2);
    std::uniform_real_distribution<double> etaDist(-5., 5.);
    std::uniform_real_distribution<double> phiDist(-M_PI, M_PI);
    std::uniform_real_distribution<double> massDist(0., 0.5);
    Particles p;
    for (int k = 0; k < n; ++k) {
      const double pt = 0.1 + ptDist(rng);
      const double eta = etaDist(rng);
      const double phi = phiDist(rng);
      const double m = (k % 3 == 0) ? 0. : massDist(rng);
      const double pz = pt * std::sinh(eta);
      p.px.push_back(pt * std::cos(phi));
      p.py.push_back(pt * std::sin(phi));
      p.pz.push_back(pz);
      p.e.push_back(std::sqrt(pt * pt + pz * pz + m * m));
    }
    return p;
  }

  using Jets = std::vector<std::pair<std::vector<int>, std::array<double, 4>>>;

  void sortJets(Jets& jets) {
    for (auto& j : jets)
      std::sort(j.first.begin(), j.first.end());
    std::sort(jets.begin(), jets.end());
  }

  // the tiled strategy must give the same merge history as the plain one
  void checkTiledMatchesPlain(Particles const& in, double R, double p) {
    const int n = in.px.size();
    std::vector<int32_t> h1(n), h2(n), hc(n), jetIdx(n), iscr(flashjet::kIntScratch * n);
    std::vector<double> hd(n), jpx(n), jpy(n), jpz(n), je(n), fscr(flashjet::kFloatScratch * n);
    std::vector<int32_t> h1t(n), h2t(n), hct(n), jetIdxT(n), iscrT(flashjet::kIntScratch * n);
    std::vector<double> hdt(n), jpxt(n), jpyt(n), jpzt(n), jet(n), fscrT(flashjet::kFloatScratch * n);
    std::vector<int32_t> tiled(flashjet::kTiledIntScratch * n);
    const auto s = flashjet::makeScratch(fscr.data(), iscr.data(), n);
    const auto st = flashjet::makeScratch(fscrT.data(), iscrT.data(), n);
    const auto t = flashjet::makeTiledScratch(tiled.data(), n);
    const int nPlain = flashjet::clusterEvent(n,
                                              R,
                                              p,
                                              in.px.data(),
                                              in.py.data(),
                                              in.pz.data(),
                                              in.e.data(),
                                              h1.data(),
                                              h2.data(),
                                              hc.data(),
                                              hd.data(),
                                              jetIdx.data(),
                                              jpx.data(),
                                              jpy.data(),
                                              jpz.data(),
                                              je.data(),
                                              s);
    const int nTiled = flashjet::clusterEventTiled(n,
                                                   R,
                                                   p,
                                                   in.px.data(),
                                                   in.py.data(),
                                                   in.pz.data(),
                                                   in.e.data(),
                                                   h1t.data(),
                                                   h2t.data(),
                                                   hct.data(),
                                                   hdt.data(),
                                                   jetIdxT.data(),
                                                   jpxt.data(),
                                                   jpyt.data(),
                                                   jpzt.data(),
                                                   jet.data(),
                                                   st,
                                                   t);
    REQUIRE(nTiled == nPlain);
    REQUIRE(h1t == h1);
    REQUIRE(h2t == h2);
    REQUIRE(hct == hc);
    REQUIRE(jetIdxT == jetIdx);
    for (int k = 0; k < nPlain; ++k) {
      REQUIRE(jpxt[k] == jpx[k]);
      REQUIRE(jpyt[k] == jpy[k]);
      REQUIRE(jpzt[k] == jpz[k]);
      REQUIRE(jet[k] == je[k]);
    }
  }

  Jets runFlashJet(Particles const& in, double R, double p) {
    const int n = in.px.size();
    std::vector<int32_t> h1(n), h2(n), hc(n), jetIdx(n), iscr(flashjet::kIntScratch * n + 1);
    std::vector<double> hd(n), jpx(n), jpy(n), jpz(n), je(n), fscr(flashjet::kFloatScratch * n + 1);
    const auto s = flashjet::makeScratch(fscr.data(), iscr.data(), n);
    const int nJets = flashjet::clusterEvent(n,
                                             R,
                                             p,
                                             in.px.data(),
                                             in.py.data(),
                                             in.pz.data(),
                                             in.e.data(),
                                             h1.data(),
                                             h2.data(),
                                             hc.data(),
                                             hd.data(),
                                             jetIdx.data(),
                                             jpx.data(),
                                             jpy.data(),
                                             jpz.data(),
                                             je.data(),
                                             s);
    Jets jets(nJets);
    for (int j = 0; j < nJets; ++j)
      jets[j].second = {jpx[j], jpy[j], jpz[j], je[j]};
    for (int k = 0; k < n; ++k) {
      REQUIRE(jetIdx[k] >= 0);
      REQUIRE(jetIdx[k] < nJets);
      jets[jetIdx[k]].first.push_back(k);
    }
    sortJets(jets);
    return jets;
  }

  Jets runFastJet(Particles const& in, double R, double p) {
    std::vector<fastjet::PseudoJet> inputs;
    for (size_t k = 0; k < in.px.size(); ++k) {
      inputs.emplace_back(in.px[k], in.py[k], in.pz[k], in.e[k]);
      inputs.back().set_user_index(k);
    }
    const fastjet::JetDefinition def = (p == -1.)  ? fastjet::JetDefinition(fastjet::antikt_algorithm, R)
                                       : (p == 0.) ? fastjet::JetDefinition(fastjet::cambridge_algorithm, R)
                                       : (p == 1.) ? fastjet::JetDefinition(fastjet::kt_algorithm, R)
                                                   : fastjet::JetDefinition(fastjet::genkt_algorithm, R, p);
    fastjet::ClusterSequence cs(inputs, def);
    Jets jets;
    for (auto const& j : cs.inclusive_jets(0.)) {
      std::vector<int> idx;
      for (auto const& c : j.constituents())
        idx.push_back(c.user_index());
      jets.emplace_back(std::move(idx), std::array<double, 4>{j.px(), j.py(), j.pz(), j.e()});
    }
    sortJets(jets);
    return jets;
  }

}  // namespace

TEST_CASE("FlashJet core matches FastJet", "[FlashJet]") {
  std::mt19937_64 rng(12345);
  for (double p : {-1., 0., 1., 0.5}) {
    for (double R : {0.4, 0.8, 1.5}) {
      for (int n : {1, 2, 3, 7, 30, 150, 600, 2000}) {
        DYNAMIC_SECTION("p=" << p << " R=" << R << " n=" << n) {
          const auto event = makeEvent(n, rng);
          checkTiledMatchesPlain(event, R, p);
          const auto flash = runFlashJet(event, R, p);
          const auto fast = runFastJet(event, R, p);
          REQUIRE(flash.size() == fast.size());
          for (size_t j = 0; j < flash.size(); ++j) {
            REQUIRE(flash[j].first == fast[j].first);
            for (int c = 0; c < 4; ++c) {
              const double a = flash[j].second[c];
              const double b = fast[j].second[c];
              REQUIRE(std::abs(a - b) <= 1e-9 * std::max(1., std::abs(b)));
            }
          }
        }
      }
    }
  }
}

TEST_CASE("FlashJet soft drop matches fastjet::contrib::SoftDrop", "[FlashJet]") {
  std::mt19937_64 rng(4242);
  std::normal_distribution<double> spread(0., 0.35);
  std::exponential_distribution<double> ptDist(0.1);
  std::uniform_real_distribution<double> unit(0., 1.);
  for (double beta : {0., 1.}) {
    for (double zcut : {0.1, 0.3}) {
      for (int n : {1, 2, 5, 40, 150}) {
        for (int trial = 0; trial < 20; ++trial) {
          DYNAMIC_SECTION("beta=" << beta << " zcut=" << zcut << " n=" << n << " trial=" << trial) {
            // a jet-like spray around a random axis
            const double eta0 = 4. * unit(rng) - 2., phi0 = 2. * M_PI * unit(rng);
            Particles in;
            for (int k = 0; k < n; ++k) {
              const double pt = 0.5 + ptDist(rng), eta = eta0 + spread(rng), phi = phi0 + spread(rng);
              in.px.push_back(pt * std::cos(phi));
              in.py.push_back(pt * std::sin(phi));
              in.pz.push_back(pt * std::sinh(eta));
              in.e.push_back(pt * std::cosh(eta));
            }
            const double R = 1000., p = 0., R0 = 0.8;
            std::vector<int32_t> h1(n), h2(n), hc(n), jetIdx(n), iscr(flashjet::kIntScratch * n);
            std::vector<double> hd(n), jpx(n), jpy(n), jpz(n), je(n), fscr(flashjet::kFloatScratch * n);
            const auto s = flashjet::makeScratch(fscr.data(), iscr.data(), n);
            const int nJets = flashjet::clusterEvent(n,
                                                     R,
                                                     p,
                                                     in.px.data(),
                                                     in.py.data(),
                                                     in.pz.data(),
                                                     in.e.data(),
                                                     h1.data(),
                                                     h2.data(),
                                                     hc.data(),
                                                     hd.data(),
                                                     jetIdx.data(),
                                                     jpx.data(),
                                                     jpy.data(),
                                                     jpz.data(),
                                                     je.data(),
                                                     s);
            REQUIRE(nJets == 1);
            const auto sd = flashjet::softDrop(n,
                                               nJets,
                                               zcut,
                                               beta,
                                               R0,
                                               in.px.data(),
                                               in.py.data(),
                                               in.pz.data(),
                                               in.e.data(),
                                               h1.data(),
                                               h2.data(),
                                               hc.data(),
                                               jpx.data(),
                                               jpy.data(),
                                               s);

            std::vector<fastjet::PseudoJet> inputs;
            for (int k = 0; k < n; ++k)
              inputs.emplace_back(in.px[k], in.py[k], in.pz[k], in.e[k]);
            fastjet::ClusterSequence cs(
                inputs, fastjet::JetDefinition(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R));
            const auto jets = fastjet::sorted_by_pt(cs.inclusive_jets(0.));
            REQUIRE(jets.size() == 1);
            fastjet::contrib::SoftDrop softDrop(beta, zcut, R0);
            softDrop.set_verbose_structure(true);
            const auto groomed = softDrop(jets[0]);
            auto const& info = groomed.structure_of<fastjet::contrib::SoftDrop>();

            const double scale = std::max(1., groomed.e());
            REQUIRE(std::abs(sd.px - groomed.px()) <= 1e-9 * scale);
            REQUIRE(std::abs(sd.py - groomed.py()) <= 1e-9 * scale);
            REQUIRE(std::abs(sd.pz - groomed.pz()) <= 1e-9 * scale);
            REQUIRE(std::abs(sd.e - groomed.e()) <= 1e-9 * scale);
            REQUIRE(sd.nDropped == info.dropped_count());
            if (groomed.has_pieces()) {
              REQUIRE(std::abs(sd.zg - info.symmetry()) <= 1e-9);
              REQUIRE(std::abs(sd.rg - info.delta_R()) <= 1e-9);
            }
          }
        }
      }
    }
  }
}

// FastJet records no nearest neighbour at exactly dist == R^2, so a pair at
// dR = R goes to the beam.  Upstream FlashJet merges it when the softer
// particle has the lower slot index; this port follows FastJet.
TEST_CASE("FlashJet leaves a pair at exactly dR = R to the beam", "[FlashJet]") {
  const double R = 0.4;
  for (double p : {-1., 0., 1.}) {
    for (bool softFirst : {false, true}) {
      DYNAMIC_SECTION("p=" << p << " softFirst=" << softFirst) {
        const double pt0 = softFirst ? 10. : 100., pt1 = softFirst ? 100. : 10.;
        Particles in;
        in.px = {pt0, pt1 * std::cos(R)};
        in.py = {0., pt1 * std::sin(R)};
        in.pz = {0., 0.};
        in.e = {pt0, pt1};
        // the two really are exactly R apart in the clustering's own measure
        REQUIRE(std::atan2(in.py[1], in.px[1]) * std::atan2(in.py[1], in.px[1]) == R * R);
        const auto flash = runFlashJet(in, R, p);
        const auto fast = runFastJet(in, R, p);
        REQUIRE(flash.size() == 2);
        REQUIRE(flash.size() == fast.size());
        REQUIRE(flash == fast);
      }
    }
  }
}
