# hep_to_parquet.py — format-agnostic HiForest / HIN-NanoAOD → Parquet (or HDF5)

A small RDataFrame wrapper that converts **both** HiForest ntuples and HIN
NanoAOD to a modern columnar format (Apache Parquet by default, HDF5 optional),
with **no per-format logic**.

## Why one wrapper handles both

Both formats are read through `ROOT::RDataFrame`, which exposes every jagged
branch — the forest's C-style `jtpt[nref]/F` arrays and NanoAOD's `akCs4PFJet_pt`
vector branches alike — uniformly as `ROOT::VecOps::RVec<T>`, and every per-event
quantity as a scalar. That uniformity is what makes the conversion agnostic:

* **jagged** branches  → Arrow **list** columns (`list<float>`, …) — variable length per event
* **scalar** branches  → Arrow primitive columns (`int32`, `float`, …)

The output is then readable identically by awkward-array, pandas, polars, DuckDB,
Spark, … regardless of whether it came from a forest or a nano file.

## Usage

```bash
cmsenv   # need PyROOT + pyarrow (+ h5py for --format hdf5); all in the el9_amd64_gcc12 externals

# one input -> OUTDIR/<stem>/<tree>.parquet (one parquet per TTree)
hep_to_parquet.py forest_had.root -o out
hep_to_parquet.py hin_had_doc.root -o out

# many inputs at once
hep_to_parquet.py forest_had.root forest_upc.root hin_had_doc.root hin_upc_mc_nano.root -o out

# pick trees / drop branches / cap events / choose format
hep_to_parquet.py forest_had.root -o out --trees akCs4PFJetAnalyzer/t,hiEvtAnalyzer/HiTree
hep_to_parquet.py hin_had_doc.root  -o out --exclude '^HLT_|^L1_'   # drop trigger-bit branches
hep_to_parquet.py forest_upc.root   -o out --max-events 1000
hep_to_parquet.py forest_had.root   -o out --format hdf5
```

Each output directory gets a `_manifest.json` (source, kind, per-tree row/column counts).

## Unified schema (`--unify`) — **either input, the same output**

The default mode keeps each format's native branch names. `--unify` instead emits a
**single canonical schema**, so a forest file and a nano file produce **identical object
tables and branch names** for the aligned content — run one analysis on either source:

```bash
hep_to_parquet.py forest_had.root  -o uni --unify   # -> uni/forest_had/{HIJet,Track,Event,ZDCRecHit,ZDCsum}.parquet
hep_to_parquet.py hin_had_doc.root -o uni --unify   # -> uni/hin_had_doc/{HIJet,Track,Event,ZDCRecHit,ZDCsum}.parquet  (SAME columns)
```

Canonical objects → one `<Object>.parquet` each, columns `<Object>_<field>` plus a shared
scalar `run`/`lumi`/`event` for joining:

| object | canonical columns (examples) | forest source | nano source |
|---|---|---|---|
| `HIJet`     | `HIJet_pt/eta/phi/mass`, `HIJet_chargedSum`, `HIJet_PfCHF`, `HIJet_refpt` | `akCs4PFJetAnalyzer/t` (`jtpt`…) | `akCs4PFJet_*` |
| `Track`     | `Track_pt/eta/charge/nHits/normChi2`, `Track_pfEnergy` | `PbPbTracks/trackTree` (`trkPt`…) | `Trk_*` |
| `Event`     | `Event_hiBin/hiHF/hiNtracks/hiZDC` | `hiEvtAnalyzer/HiTree` (`hiBin`…) | `GO_*` |
| `ZDCRecHit` | `ZDCRecHit_energy/time/zside/section` | `zdcanalyzer/zdcrechit` (`energy`…) | `ZDC_*` |
| `ZDCsum`    | `ZDCsum_Plus/Minus` | `zdcanalyzer/zdcrechit` (`sumPlus`…) | `ZDCsum_*` |
| `Muon`      | `Muon_pt/eta/charge`, `Muon_looseId/tightId`, `Muon_isoTrkR03` | `muonAnalyzer/MuonTree` (`recoPt`…) | `Muon_*` |
| `Electron`  | `Electron_pt/eta/charge/sieie/hoe`, `Electron_hiPFChIso03` | `ggHiNtuplizer/EventTree` (`elePt`…) | `Electron_*` |
| `Photon`    | `Photon_pt/eta/r9/hoe/sieie`, `Photon_pixelSeed` | `ggHiNtuplizer/EventTree` (`phoEt`…) | `Photon_*` |

For the **HINHAD** flavour (full HI content) every applicable object is byte-identical
forest↔nano. Differences only arise from genuine content gaps (e.g. HINUPC nano omits the
HI iso extensions `Muon_isoTrkR03` / `Electron_hiPF*`; the hadronic forest config omits
`muonAnalyzer`) — the canonical names are still correct wherever the field exists.
Forest muon `reco*`/`inner*`/`global*` are separate collections of different length; only
the per-muon `reco*` block is mapped.

(`HIJet`, not `Jet`, avoids colliding with NanoAOD's standard `Jet` collection; it also
covers the UPC `ak0PF` jets.) The mapping lives in `UNIFIED_SCHEMA` at the top of the
script — add objects/fields there. Fields absent in a given file (e.g. `HIJet_ref*` on
data) are simply omitted on **both** sides, so forest and nano stay schema-identical.

```python
import pyarrow.parquet as pq
pq.read_schema("uni/forest_had/HIJet.parquet") == pq.read_schema("uni/hin_had_doc/HIJet.parquet")  # True
```

## Reading the output

```python
import pyarrow.parquet as pq
jets = pq.read_table("out/forest_had/akCs4PFJetAnalyzer__t.parquet")
jets.column("jtpt")[0]              # -> list of the event-0 jet pTs

import awkward as ak                # jagged columns load natively as awkward arrays
ev = ak.from_parquet("out/hin_had_doc/Events.parquet")
ev.akCs4PFJet_pt                    # var * float

import pandas as pd                 # scalar (per-event) columns
go = pd.read_parquet("out/forest_had/hiEvtAnalyzer__HiTree.parquet", columns=["hiBin","hiHF"])
```

Because the forest and nano are converted on the *same events*, a cross-format
check is a join (or positional compare) on `run`/`lumi`/`event` — e.g. nano
`akCs4PFJet_pt` vs forest `akCs4PFJetAnalyzer__t.jtpt`.

## Notes / limitations

* Signed `char` branches (e.g. `trkCharge`, `trkNHits`) are recast to `int` inside
  RDataFrame — PyROOT cannot numpy-convert `RVec<char>` directly.
* Genuinely unsupported column types are skipped with a warning (and listed nowhere
  in the output): string-as-char-array provenance (`HiForestInfo`), EDM provenance
  (`MetaData`/`ParameterSets`), `pair<…>` (PDF info), nested `RVec<vector<int>>`
  (gen mother/daughter index lists).
* `AsNumpy` loads whole columns into memory; for very large files, restrict branches
  with `--exclude`/`--trees` or process in passes with `--max-events`.
