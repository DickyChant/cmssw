#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct3.h"

bool GenRunInfoProduct3::isProductEqual(GenRunInfoProduct3 const &other) const {
  const HepMC3::GenRunInfoData &a = data_;
  const HepMC3::GenRunInfoData &b = other.data_;
  return a.weight_names == b.weight_names && a.tool_name == b.tool_name && a.tool_version == b.tool_version &&
         a.tool_description == b.tool_description && a.attribute_name == b.attribute_name &&
         a.attribute_string == b.attribute_string;
}
