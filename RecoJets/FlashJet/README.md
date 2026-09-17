# RecoJets/FlashJet

Prototype integration of [FlashJet](https://github.com/jet-universe/FlashJet) --
batched generalized-kt jet clustering (anti-kt, kt, Cambridge/Aachen, E-scheme)
written for GPUs -- into CMSSW, in two flavours that produce the same
`reco::PFJet` / `reco::GenJet` / `reco::BasicJet` collections as
`FastjetJetProducer`, plus soft-drop observables as `ValueMap<float>`.

Status: **validated, not yet competitive on GPU.** The clustering agrees with
FastJet exactly on generator-level, Phase-2 PU200 and Run-3 scouting events
(see [Validation](#validation)), but only the CPU backend is fast enough to be
interesting today (see [Performance](#performance)).

## Clustering

`interface/FlashJetCore.h` and `interface/FlashJetTiled.h` hold the algorithm
itself, as allocation-free `ALPAKA_FN_HOST_ACC` functions that take all their
working memory from the caller, so the same code runs on the host and inside an
alpaka kernel:

| | strategy | used by |
|---|---|---|
| `clusterEvent` | plain nearest neighbours, O(n^2) | GPU backends, and small entries |
| `clusterEventTiled` | linked cells on a (rapidity, phi) grid + indexed binary heap + reverse nearest-neighbour index, O(n log n) | CPU backends above 12 particles |
| `softDrop` | de-clustering of the merge history (`fastjet::contrib::SoftDrop`, scalar_z, larger_pt) | both |

Both strategies produce the same merge history, and the unit test compares them
directly.  Numerics follow FlashJet: double precision, FastJet's stable
rapidity, E-scheme sums accumulated in merge order, and every full rescan
breaking ties towards the lowest slot index.  One deviation from upstream
FlashJet is deliberate: a pair at exactly `dR = R` is left to the beam, as
FastJet does (upstream merges it when the softer particle has the lower slot
index).

## Modules

`flashjet::FlashJetDeviceCollection` (`DataFormats/FlashJet`) holds a batch of
independent clustering problems ("entries"): a particle block (inputs, particle
-> jet index, jet four-momenta, merge history) and an entry block (offset/size,
number of jets, soft-drop result).  One kernel call processes every entry.

* `FlashJetProducer@alpaka` -- all candidates of the event as one entry.
  `FlashJetRecoJetProducer` turns the host copy into reco jets.
* `FlashJetReclusterProducer@alpaka` -- the constituents of every jet of a jet
  collection as one entry each (C/A with R = 1000 by default, like
  `fastjet::contrib::Recluster`), optionally soft-dropped on the device.
  `FlashJetSoftDropProducer` writes mass, pt, zg, rg and nDropped as ValueMaps;
  `FastjetSoftDropProducer` computes the same with FastJet as a reference.
* `FlashJetSonicProducer` -- sends the candidates to the `flashjet` model on a
  Triton inference server (`data/models/flashjet`) and builds reco jets from the
  returned particle -> jet map.  The model pads the requests that arrive
  together into one batch and clusters them with FlashJet's own Triton kernels
  (GPU), its C++ kernel, or NumPy.  `test/setupFlashJetModel.sh` points it at a
  FlashJet checkout and builds the C++ kernel.

On GPU backends one *block* handles one entry, with the threads sharing every
scan over its particles; on CPU backends one work item handles one entry.

## Validation

`scram b runtests` runs both:

* `testFlashJetCore` -- the algorithm against FastJet and
  `fastjet::contrib::SoftDrop` on random events (all algorithms, several radii,
  up to 2000 particles), tiled against plain, and the `dR = R` boundary.
* `testFlashJetModules.sh` -- `cmsRun` against `FastjetJetProducer` on generated
  TTbar events, for the alpaka backend given as its argument.

`test/benchmarkFlashJet_cfg.py --compare` does the same on real events.  Exact
agreement (0 mismatches) has been established for:

| Input | Whole-event AK4 | Reclustering + soft drop |
|---|---|---|
| Generated TTbar (gen particles) | anti-kt, kt, C/A | AK8 |
| Phase-2 PU200 TTbar MiniAOD (`packedPFCandidates`, N ~ 10200) | 37933 jets | 126010 values (AK4), 27145 (AK8) |
| Run-3 PF scouting (`hltScoutingPFPacker`, N ~ 320) | 177631 jets | 43660 values (AK4), 6890 (AK8) |

on the CPU backend, the CUDA backend (H100) and SONIC (GPU server).

## Performance

`test/benchmarkFlashJet_cfg.py` (FastTimerService throughput after a warm-up)
and `test/summarizeBenchmarks.py`.  Numbers below: 4 CPU threads, GPU an H100
MIG 1g.12gb slice (about 1/7 of a card), best over 1--16 CMSSW streams.

Whole-event anti-kt R = 0.4, throughput in events/s:

| Input (N/event) | FastJet | FlashJet CPU (tiled) | FlashJet CUDA | SONIC |
|---|---|---|---|---|
| Scouting (~320) | 6476 | 5415 | 624 | 1745 |
| Synthetic (~2000) | 1630 | 1152 | 23 | 541 |
| Synthetic (~5500) | 579 | 270 | 3.2 | 86 |
| PU200 MiniAOD (~10200) | 97 | 69 | 0.7 | 23 |

Reclustering every jet of the event and soft-dropping it -- many entries per
kernel call -- in events/s:

| Workload | FastJet | FlashJet CPU (tiled) | FlashJet CUDA |
|---|---|---|---|
| AK4 jets, scouting | 5997 | 5143 | 4357 |
| AK8 jets, synthetic (~2000) | 1219 | 1105 | 1151 |
| AK8 jets, synthetic (~5500) | 331 | 266 | **455** |

The last row is the one case measured so far where the GPU wins (1.37x over
FastJet), and it scales with streams (62 / 235 / 381 / 455 at 1 / 4 / 8 / 16):
enough jets per event, each with enough constituents to fill a block.  The
SONIC column above is a separate round (its server batches across streams
rather than within the event).

Numbers on these shared MIG slices vary by up to ~35% between rounds, so
compare implementations within one job, not across rounds.

Per-event module time, the same workload:

| Input | FastJet | FlashJet CPU (tiled) | FlashJet CPU (plain, before) |
|---|---|---|---|
| Scouting | 0.33 ms | **0.275 ms** | 1.52 ms |
| Synthetic ~2000 | 1.84 ms | 2.89 ms | 32.0 ms |
| PU200 | ~15 ms | 40.8 ms | 1318 ms |

What this says:

* The **CPU backend is the usable one**: faster than FastJet on scouting-sized
  events, within 1.5--3x at higher multiplicity.
* **Whole-event clustering does not suit a GPU here.**  An event is a serial
  chain of N merges, and CMSSW hands a module one event at a time, so the
  parallelism FlashJet was built for (many events per kernel launch) is absent.
  Giving an event a whole block instead of one thread bought 3--5x; cheaper
  reductions on top of that bought nothing measurable.
* **Batching is what helps.**  Per-jet reclustering (many entries per event)
  brings the GPU level with the CPU, and ahead of it for AK8 jets at high
  multiplicity; SONIC -- which batches requests from concurrent streams
  server-side -- scales from 125 ev/s at one stream to 1745 at sixteen, though
  it still loses to FastJet on this slice.
* The GPU numbers come from 1/7 of an H100; a full card would change them, but
  not the structure of the problem.

## Benchmarks and GPU jobs on lxplus

`test/condor/` submits to the EosSubmit schedds (everything on EOS, nothing on
AFS):

```
cmsenv
test/condor/prepareCondor.sh /eos/user/X/USER/flashjet_condor /path/to/FlashJet
test/condor/makeInputs.sh    /eos/user/X/USER/flashjet_condor   # real inputs, once
test/condor/makeSonicEnv.sh  /eos/user/X/USER/flashjet_condor   # torch+triton for the server, once
# on lxplus
module load lxbatch/eossubmit
condor_submit /eos/user/X/USER/flashjet_condor/flashjet_gpu.sub    # validation + synthetic scans
condor_submit /eos/user/X/USER/flashjet_condor/flashjet_real.sub   # verification + scans on real inputs
condor_submit /eos/user/X/USER/flashjet_condor/flashjet_sonic.sub  # Triton server with the GPU kernels
python3 test/summarizeBenchmarks.py /eos/user/X/USER/flashjet_condor/results/<cluster>
```

Note the jobs require a GPU newer than Volta: CMSSW_20_1 is built with CUDA 13,
which needs compute capability >= 7.5, so the submit files exclude V100/P100.

## Known gaps

* **No jet areas or ghosts**, so this cannot replace jets that carry the
  rho-area pileup correction.
* **PUPPI weights are not applied**: the recluster producers take the daughter
  four-vectors as they are.  The FastJet reference producer does the same, so
  their agreement does not by itself prove agreement with standard weighted
  grooming.
* **Backend reproducibility**: the alpaka path is double precision, while
  FlashJet's Triton kernels (used through SONIC) are float32 and can order
  near-degenerate merges differently; SONIC also sums jet constituents in input
  order rather than merge order.
* **Only ValueMaps** are produced for grooming: no groomed jet collection and no
  subjet association.
* **HLT**: no multiplicity or latency guard, no fallback policy if a server is
  unavailable, and no tail-latency measurement.
* **Licensing**: FlashJet is GPL-3.0 and CMSSW is Apache-2.0.  The kernels here
  are a port of GPL code and would have to be relicensed, or kept in an external
  package, before this could go to `cms-sw/cmssw`.
