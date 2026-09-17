// Compares flashjet::clusterEvent (host build of the alpaka kernel core)
// with FastJet on random events: identical jet constituents and
// four-momenta equal up to rounding.

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include <catch2/catch_all.hpp>
#include <fastjet/ClusterSequence.hh>

#include "RecoJets/FlashJet/interface/FlashJetCore.h"

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

  Jets runFlashJet(Particles const& in, double R, double p) {
    const int n = in.px.size();
    std::vector<int32_t> h1(n), h2(n), hc(n), jetIdx(n), iscr(6 * n);
    std::vector<double> hd(n), jpx(n), jpy(n), jpz(n), je(n), fscr(9 * n);
    const flashjet::Scratch s{&fscr[0],
                              &fscr[n],
                              &fscr[2 * n],
                              &fscr[3 * n],
                              &fscr[4 * n],
                              &fscr[5 * n],
                              &fscr[6 * n],
                              &fscr[7 * n],
                              &fscr[8 * n],
                              &iscr[0],
                              &iscr[n],
                              &iscr[2 * n],
                              &iscr[3 * n],
                              &iscr[4 * n]};
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
