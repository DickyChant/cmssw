# HDF5 input milestone validation — 2026-09-11

Branch: `improve_LHE_interface`, based on official `CMSSW_20_1_X` commit
`bc23eec048cf641dc788b0617f588d6c7be83c0c`.

## Environment and limits of the result

The isolated SCRAM area uses installed `CMSSW_20_1_0_pre1`, architecture
`el9_amd64_gcc13`, and HighFive 2.10.1. Before these changes, both
`GeneratorInterface/LHEInterface` and `SimDataFormats/GeneratorProducts` were
identical between that release and the pinned 20_1 branch. The changed LHE
package was rebuilt against the installed framework; this is not a full build
of all of 20_1_X or central CI certification.

Runtime tests used the installed CERN site configuration for InputFileCatalog;
all event inputs and outputs were local temporary files. No EOS event transfer,
grid submission or remote HDF5 workflow was exercised.

## Passed checks

- Single-job SCRAM package build, including registration of both `LHESource`
  and the unchanged legacy `LH5Source`.
- `scram b -j 1 runtests_TestLHEH5Input`: passed. This includes direct-reader
  XML/HDF5 parity, malformed/unsupported input rejection, and actual EDM
  write/readback for basic input (3 events), skipping across empty/nonempty
  files with an event limit (1 event), incompatible run metadata (6 events),
  Pythia showering (3 events), and explicit lossy input (1 event).
- Persisted original event numbers and signed weights agree across backends;
  the lossy-input case verifies the event-level loss declaration.
- Empty/short inputs and 0/1/1599/1600/1601-event boundaries; reversed column
  order and dataset-named label attributes; unit-weight and weighted modes;
  fractional offsets, oversized particles/chunks, bad mothers/process IDs,
  malformed versions/labels, non-finite weights and unsupported protocols.
- Existing XML `testMerging.sh` merge/non-merge output regressions: passed.
- Actual unmodified `pylhe` writer at
  `f7889ce52b1d93dd0262656a624b53b20885d5d3`: its HDF5 output agrees with the
  corresponding XML input in the CMSSW reader test (three core events).
  Python reference dependencies were installed only in an isolated test area;
  they are not runtime dependencies of the CMSSW backend.
- The Pythia smoke sample reports no generator warnings or errors. It checks
  interface operation and persisted output, not broad physics validation.
- Source comparison confirms no changes to the legacy LH5 reader, helper or
  source implementation. No legacy-format runtime fixture was exercised.

## Memory observations

The XML/HDF5 streaming comparison was run in separate processes, retaining only
the current pair of events, with the bounded HDF5 metadata cache enabled:

| Events compared | Peak RSS (KiB) | Elapsed time |
| --- | ---: | ---: |
| 1,000 | 48,980 | 0.56 s |
| 50,000 | 66,080 | 3.50 s |
| 200,000 | 82,792 | 12.62 s |

All comparisons passed. These are measured low-RSS results, not a claim of
perfectly constant total RSS: allocator, XML and HDF5 internals still contribute.
The first full serial package build peaked at 587,512 KiB. Payload/chunk/cache
bounds are documented in [HDF5Input.md](HDF5Input.md).

The first push contains this HDF5 milestone only. Full MadGraph LO/NLO
reweighting-payload persistence, lossless LHE export and numerical systematics
closure remain the next milestone; none is claimed by these tests.
