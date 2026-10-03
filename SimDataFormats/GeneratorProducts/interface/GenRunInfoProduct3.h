#ifndef SimDataFormats_GeneratorProducts_GenRunInfoProduct3_h
#define SimDataFormats_GeneratorProducts_GenRunInfoProduct3_h

#include <string>
#include <utility>
#include <vector>

#include <HepMC3/Data/GenRunInfoData.h>

/** \class GenRunInfoProduct3
 *
 *  The HepMC3 run information of an input file: the names of the event
 *  weights, the tools which produced the events, and the run attributes.
 *  The event weights of GenEventInfoProduct3 are in the order of weightNames().
 */

class GenRunInfoProduct3 {
public:
  GenRunInfoProduct3() = default;
  explicit GenRunInfoProduct3(HepMC3::GenRunInfoData data) : data_(std::move(data)) {}

  const HepMC3::GenRunInfoData &runInfo() const { return data_; }

  const std::vector<std::string> &weightNames() const { return data_.weight_names; }

  // two runs can be merged if they have the same run information
  bool isProductEqual(GenRunInfoProduct3 const &other) const;

private:
  HepMC3::GenRunInfoData data_;
};

#endif  // SimDataFormats_GeneratorProducts_GenRunInfoProduct3_h
