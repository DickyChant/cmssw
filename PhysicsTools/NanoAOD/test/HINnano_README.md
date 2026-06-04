# Heavy-Ion (HIN) NanoAOD flavours

Two NanoAOD flavours for Run-3 heavy-ion analyses, built on top of the standard
`PHYS` content (they only *add* tables). They port the most-used HiForest pieces
(the `HiEvtAnalyzer` event-activity ntuple and the ZDC ntuple) to NanoAOD
`FlatTable`s, and add BTVNano-style PF-candidate / track content.

| Flavour  | Use case                         | Adds |
|----------|----------------------------------|------|
| `HINUPC` | ultra-peripheral collisions      | PF cands + lost tracks (with track details), ZDC rechit ntuple, **HF / ZDC event-activity table (no centrality bin)** |
| `HINHAD` | hadronic heavy-ion (e.g. PbPb)   | same PF/track + ZDC, plus the **full centrality table incl. `hiBin`** |

For UPC the centrality *bin* is undefined, so it is omitted; the HF and ZDC sums —
the important UPC observables — are still written (`Cent_hiHF*`, `Cent_hiZDC*`,
`ZDCsum_Plus/Minus`).

## Running

```bash
# HINUPC (data, PbPb era)
cmsDriver.py hin_upc_nano -s NANO:@HINUPC --data --era Run3_pp_on_PbPb \
  --conditions 151X_dataRun3_Prompt_v1 --eventcontent NANOAOD --datatier NANOAOD \
  --filein /store/hidata/.../MINIAOD/... --fileout file:hin_upc_nano.root -n 1000

# HINHAD (hadronic PbPb)
cmsDriver.py hin_had_nano -s NANO:@HINHAD --data --era Run3_pp_on_PbPb \
  --conditions 151X_dataRun3_Prompt_v1 --eventcontent NANOAOD --datatier NANOAOD \
  --filein /store/hidata/.../MINIAOD/... --fileout file:hin_had_nano.root -n 1000
```

Optional process modifier `--procModifiers run3_nanoAOD_HIN` applies a soft
`pt > 0.5` threshold to the PF-candidate / lost-track tables to keep the output
size under control in high-multiplicity events.

## Implementation

* `python/custom_hin_cff.py` — customise functions `HINUPCCustomNanoAOD` /
  `HINHADCustomNanoAOD` and the table/producer definitions. Registered in
  `python/autoNANO.py` as `HINUPC` / `HINHAD`.
* `plugins/CentralityTableProducer.cc` — singleton `Cent` table from a
  `reco::Centrality` (`hiCentrality`) plus the optional `int` centrality bin
  (`centralityBin:HFtowers`). Same observables as `HiEvtAnalyzer`
  (`hiHF*`, `hiHFhit*`, `hiEB/hiEE/hiET`, `hiZDC*`, `hiNpix`, `hiNtracks`, ...).
  Leave `srcBin` empty to skip `hiBin` (UPC).
* `plugins/ZDCTableProducer.cc` — variable-length `ZDC` table (one row per ZDC
  rechit: `zside`, `section`, `channel`, `depth`, `energy`, `time`, ...) and a
  singleton `ZDCsum` table with the per-side EM+HAD energy sums (0n/Xn tagging).
* `Configuration/ProcessModifiers/python/run3_nanoAOD_HIN_cff.py` — the
  `run3_nanoAOD_HIN` process modifier.

### Inputs (heavy-ion MiniAOD)

* ZDC rechits are reconstructed on the fly by the central
  `RecoLocalCalo/HcalRecProducers` `ZdcHitReconstructor_Run3` from the
  `hcalDigis:ZDC` digis kept in the HI MiniAOD — no forest dependency.
* `hiCentrality` and `centralityBin:HFtowers` are taken from the input event
  (as in the forest). If they are absent, schedule the central
  `RecoHI/HiCentralityAlgos` producers upstream.

### Notes / possible extensions

`HINHAD` is a solid hadronic foundation (PF/track + ZDC + full centrality). Natural
follow-ups for hadronic analyses: event-plane / Q-vector table (`RecoHI/HiEvtPlaneAlgos`)
and HI constituent-subtracted jet tables.
