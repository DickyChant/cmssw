#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// Explicit output type: the release's GenParticleRefSelector actually uses
// the default GenParticleCollection output, despite its historical name.
using GenPartonRefVectorSelector = SingleObjectSelector<
    reco::GenParticleCollection,
    StringCutObjectSelector<reco::GenParticle>,
    reco::GenParticleRefVector>;
DEFINE_FWK_MODULE(GenPartonRefVectorSelector);
