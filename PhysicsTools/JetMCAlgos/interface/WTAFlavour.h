#ifndef PhysicsTools_JetMCAlgos_WTAFlavour_h
#define PhysicsTools_JetMCAlgos_WTAFlavour_h

#include <cmath>
#include <stdexcept>
#include <vector>
#include "fastjet/ClusterSequence.hh"

namespace jetflavour {
  // Terminal-parton WTA label at the input shower resolution. This is soft
  // insensitive, but not a collinear-safe fixed-order flavour observable.
  // user_index identifies a physical input constituent, never a flavour ghost.
  inline int wtaWinner(const std::vector<fastjet::PseudoJet>& inputs) {
    if (inputs.empty())
      return -1;
    for (const auto& input : inputs) {
      if (!std::isfinite(input.px()) || !std::isfinite(input.py()) || !std::isfinite(input.pz()) ||
          !std::isfinite(input.E()) || input.pt() <= 0. || input.E() < std::abs(input.pz()) ||
          input.user_index() < 0)
        throw std::invalid_argument("WTA requires finite physical momenta, positive pT and input indices");
    }
    // Recluster all members of ONE existing jet. A large radius prevents a
    // second small-R acceptance cut from discarding any of its constituents.
    const fastjet::JetDefinition definition(
        fastjet::kt_algorithm, fastjet::JetDefinition::max_allowable_R, fastjet::WTA_pt_scheme);
    fastjet::ClusterSequence sequence(inputs, definition);
    auto jet = sequence.exclusive_jets(1).front();
    fastjet::PseudoJet first, second;
    while (jet.has_parents(first, second)) {
      // Read the actual WTA axis, including FastJet's own equal-pT convention.
      // It coincides with the winning parent's direction at every merge.
      jet = jet.squared_distance(first) <= jet.squared_distance(second) ? first : second;
    }
    return jet.user_index();
  }
}  // namespace jetflavour

#endif
