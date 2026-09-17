# RecoJets/FlashJet

Prototype integration of [FlashJet](https://github.com/jet-universe/FlashJet),
batched generalized-kt jet clustering (anti-kt, kt, Cambridge/Aachen, E-scheme),
in two flavours that produce the same `reco::PFJet` / `reco::GenJet` /
`reco::BasicJet` collections as `FastjetJetProducer` (no jet areas).

## Alpaka

* `interface/FlashJetCore.h`: the clustering itself, a port of FlashJet's C++
  kernel (FastJet N2Plain nearest-neighbour strategy, double precision) written
  as one allocation-free `ALPAKA_FN_HOST_ACC` function.
* `flashjet::FlashJetDeviceCollection` (`DataFormats/FlashJet`) holds a batch of
  independent clustering problems ("entries"): a particle block (inputs,
  particle -> jet index, jet four-momenta, merge history) and an entry block
  (offset/size, number of jets, soft-drop results).  One kernel call processes
  all entries, one device thread per entry.
* `FlashJetProducer@alpaka`: all candidates of the event as one entry;
  `FlashJetRecoJetProducer` turns the host copy into reco jets.
* `FlashJetReclusterProducer@alpaka`: the constituents of every jet of a jet
  collection as one entry each (C/A with R = 1000 by default, like
  `fastjet::contrib::Recluster`), optionally soft-dropped on the device
  (`softDrop.enable`); `FlashJetSoftDropProducer` stores mass, pt, zg, rg and
  nDropped as `ValueMap<float>`s.  `FastjetSoftDropProducer` computes the same
  ValueMaps with FastJet contrib as a reference.

The whole-event mode is a single device thread per event; GPU gains are
expected from the batched modes (many jets per call), and across events from
SONIC, where the server batches requests from many streams.

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
* `testFlashJetCore` also checks `softDrop` against `fastjet::contrib::SoftDrop`
  (groomed four-momentum, zg, rg, dropped count).
* `test/testFlashJet_cfg.py`: FastJet vs FlashJet (alpaka and/or SONIC) on
  generator-level jets, jet by jet; `--recluster` adds the batched soft drop of
  AK8 jets vs FastJet contrib.

## Benchmarks and GPU jobs on lxplus

* `test/benchmarkFlashJet_cfg.py`: FastTimerService timing of FastJet, FlashJet
  alpaka (any backend) and SONIC on synthetic events of tunable multiplicity,
  for whole-event AK4 clustering and batched AK8 soft drop;
  `test/summarizeBenchmarks.py` tabulates the JSON output.
* `test/condor/`: HTCondor GPU jobs for the EosSubmit schedds.
  ```
  cmsenv
  test/condor/prepareCondor.sh /eos/user/X/USER/flashjet_condor /path/to/FlashJet
  test/condor/makeSonicEnv.sh /eos/user/X/USER/flashjet_condor   # once, for SONIC jobs
  # on lxplus
  module load lxbatch/eossubmit
  condor_submit /eos/user/X/USER/flashjet_condor/flashjet_gpu.sub     # validation + alpaka benchmarks
  condor_submit /eos/user/X/USER/flashjet_condor/flashjet_sonic.sub   # Triton server with FlashJet GPU kernels
  ```
  Edit `tasks_gpu.txt` / `tasks_sonic.txt` in the EOS directory to change the job list.
