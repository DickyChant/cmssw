# RecoJets/FlashJet

Prototype integration of [FlashJet](https://github.com/jet-universe/FlashJet),
batched generalized-kt jet clustering (anti-kt, kt, Cambridge/Aachen, E-scheme),
in two flavours that produce the same `reco::PFJet` / `reco::GenJet` /
`reco::BasicJet` collections as `FastjetJetProducer` (no jet areas).

## Alpaka

* `interface/FlashJetCore.h`: the clustering itself, a port of FlashJet's C++
  kernel (FastJet N2Plain nearest-neighbour strategy, double precision) written
  as one allocation-free `ALPAKA_FN_HOST_ACC` function.
* `FlashJetProducer@alpaka`: candidates (`edm::View<reco::Candidate>`) ->
  `flashjet::FlashJetDeviceCollection` (`DataFormats/FlashJet`): particle -> jet
  index, jet four-momenta and the full merge history.
* `FlashJetRecoJetProducer`: host copy of that collection -> reco jets.

One event is clustered by one device thread, so today the GPU backends only
offload the work; the speedup is meant to come from batching many clustering
problems (e.g. jet constituents for substructure) into one kernel call.

## SONIC

* `FlashJetSonicProducer`: sends the candidates of each event to the `flashjet`
  model on a Triton inference server and builds reco jets from the returned
  particle -> jet map.
* `data/models/flashjet`: Triton Python-backend model.  Requests arriving
  together are padded into one batch and clustered with the FlashJet Triton
  kernels (GPU server with torch + triton), FlashJet's C++ kernel (CPU), or the
  NumPy reference.  Point it to a FlashJet checkout with
  `test/setupFlashJetModel.sh /path/to/FlashJet`, which also builds the C++ kernel.

## Validation

* `testFlashJetCore`: `FlashJetCore.h` vs FastJet on random events
  (all algorithms, several R, up to 2000 particles): identical constituents,
  momenta equal to 1e-9.
* `test/testFlashJet_cfg.py`: FastJet vs FlashJet (alpaka and/or SONIC) on
  generator-level jets, jet by jet.
