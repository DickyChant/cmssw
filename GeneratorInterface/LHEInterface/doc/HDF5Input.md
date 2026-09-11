# XML and consolidated HDF5 input through LHESource

This additive backend keeps `LHEEventProduct` and `LHERunInfoProduct` as the
common hard-process interface. XML remains the default. The legacy `LH5Source`
and `LH5Reader` are not changed, redirected or removed.

## Configuration

```python
process.source = cms.Source(
    "LHESource",
    fileNames=cms.untracked.vstring("file:/path/events.h5"),
    inputFormat=cms.untracked.string("hdf5"),
    hdf5AllowUnsupportedMetadata=cms.untracked.bool(False),
    hdf5MaxParticlesPerEvent=cms.untracked.uint32(100000),
)
```

`inputFormat="xml"` selects the existing XML parser. File formats are selected
explicitly, not inferred from suffixes. A source's file list uses one encoding.
`skipEvents`, framework `maxEvents`, file/run transitions and EDM conversion
remain in the shared source path. `ExternalLHEProducer` retains its XML-reader
API and behavior. No new cmsdist dependency or CMake build is required.

HDF5 currently accepts local `file:` paths, including locally mounted EOS FUSE.
It does not implement direct XRootD, remote HDF5 access, automatic format
detection or an HDF5 output writer.

## Precisely supported schema

The initial allowlist is consolidated **LHEH5 2.0.0**, not every historical LH5
or future 2.x layout. `/version` must be an integer triplet. The numeric tables
are float64: `/init` (one-dimensional), `/procInfo`, `/events`, `/particles`
(two-dimensional). Columns must have unique labels in a `properties` attribute,
or the older dataset-named attribute. Column order is not fixed.

The mapping was checked against the actual `pylhe` writer at
[`f7889ce52b1d93dd0262656a624b53b20885d5d3`](https://github.com/scikit-hep/pylhe/blob/f7889ce52b1d93dd0262656a624b53b20885d5d3/src/pylhe/lheh5.py).
See the [schema description](https://scikit-hep.org/pylhe/lheh5.html).
This does not certify every producer's generator-specific extensions.

| Input | Destination / behavior |
| --- | --- |
| `init` beams, energies, PDF identifiers, weighting strategy | HEPRUP, without a hard-coded IDWTUP |
| `procInfo` cross sections/errors/maxima and process IDs | HEPRUP process records |
| `events` process ID, nominal `weight` or `NOMINAL`, scale and couplings | HEPEUP plus `originalXWGTUP` |
| `particles` IDs, statuses, mothers, colours, momenta, masses, lifetimes, spin | HEPEUP in original particle order; mothers remain one-based |
| Process/event `npLO`, `npNLO` | Existing multiplicity fields; missing/NaN is `-99` |
| Optional `events/event_num` | Existing CMSSW event number; absent is `-1`; this is an explicitly supported extension |
| `events/start`, `nparticles` | Validated particle slice; storage offsets are not physics metadata |
| `trials`, `rscale`, `fscale` | NaN means missing; finite values require explicit lossy mode |
| Extra datasets, unknown columns or metadata attributes | Strict rejection, or declared omission in explicit lossy mode |

The original signed nominal weight is retained separately from the shower-input
weight. For `abs(IDWTUP)==3`, the latter follows the existing XML reader's unit-
weight convention. Negative weights with positive IDWTUP are rejected instead
of silently changing their sign. Malformed core fields are rejected in both
strict and lossy modes.

## Metadata preservation boundary

Strict mode is the default. It must not silently discard named alternative
weights, generator headers or counterterm information. This first milestone has
no mappings for those HDF5 extensions, so their presence is rejected. A core-only
file whose optional scale/trial columns are all NaN can be read strictly.

Setting `hdf5AllowUnsupportedMetadata=True` deliberately permits core-event-only
input. It logs the omitted fields and stores a loss declaration in both run and
event comments; the latter survives merging run products across files. This is
not a lossless conversion, not full MadGraph post-hoc reweighting support, and
not permission to classify anonymous columns as PDF/scale weights. Retain the
original file when using this mode.

MadGraph payload persistence and numerical systematics validation are a separate
subsequent milestone. No NanoAOD columns or persistent product schema are added
here. Existing XML named-weight handling is unchanged.

## Memory and record validation

Read one event row and its particle hyperslab at a time. There is no requirement
that file length be divisible by 1,600, and no whole-file event array. The default
particle bound is 100,000 per event. Additional guards bound columns to 256,
process records to 100,000, uncompressed HDF5 storage chunks to 64 MiB, and the
adaptive HDF5 metadata cache to 4 MiB. These bound the main reader allocations,
not the RSS of the entire framework or arbitrary third-party HDF5 internals.

The reader checks shapes, labels, process references, finite core values,
integral/range-safe integer fields, particle offsets and mother references.
Unsupported schema versions and protocols fail explicitly. Oversized storage
chunks are rejected before hyperslab reading, since decompression can otherwise
allocate much more than the requested slice.

## Tests and associated XML correction

`testLHEH5Reader` compares XML/HDF5 records, exercises empty/short files,
0/1/1599/1600/1601-event lengths, column permutations, signed weights, skips,
multiple files and invalid inputs. The fixture generator uses test-only h5py.
`TestLHEH5Input` additionally runs `cmsRun`, reads the persisted EDM products,
checks run transitions and event limits, verifies loss declarations, and runs
a small Pythia shower test. Existing XML merging tests remain enabled.

The parity test exposed an existing XML event-number lifetime bug: a pointer
outlived a temporary transcoded string. The parser now retains the string during
checked integer conversion and rejects missing/invalid `num` attributes.
`LHESource` also copies the parsed event number into the persistent product.
Valid syntax is `<event_num num="123"/>`. This is the only intentional XML
data-handling correction in this milestone, with dedicated regression coverage.

From a configured SCRAM area with a valid site configuration:

```sh
scram b -j 1
scram b -j 1 runtests_TestLHEH5Input
```

The test prints its temporary artifact directory for inspection. No broad
CMSSW build or full-tree indexing is needed.
