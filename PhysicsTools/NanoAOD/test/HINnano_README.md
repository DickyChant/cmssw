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

## HINHADSKIM: fully hadronic ttbar skim flavour

`HINHADSKIM` = `HINHAD` + `custom_hin_cff.addFullHadSkim`: keeps only events
with >= 5 akCs3PF jets with UParT-regressed pT (`rawPt * UParT ptcorr`)
>= 25 GeV and |eta| <= 2.4 — strictly looser than the offline full-hadronic
selection (>= 6 jets, pT > 30, |eta| < 2.1 signal region, nJet == 5 sideband).
No b-tag cut is applied so the full b-tag spectrum (0B/1B mixtures, negative
tags) survives for the data-driven (CWoLa / anomaly-detection) region
building; `requireBTags=True` re-enables >= 2 jets with normalized UParT
b-discriminant >= 0.15. An unbiased >= 5-jet precount on the bare
`patJetsAKCs3PF` keeps the b-tag inference off for most events; the skim
filters are prepended to `nanoAOD_step` (gates the scheduled nano tables) and
the NanoAOD output gets `SelectEvents = ["nanoAOD_step"]`.

The akCs3PF chain (UParT pinned to the forest's
`HeavyIonsAnalysis/Configuration/data/PbPb_AK3_2024_v6.onnx` training, so the
scores and the measured skim rate match the HiForest twin; the central
cms-data 2023 model regresses low-pT jets ~35% higher) and an
`akCs3PFJet` table — including the embedded UParT outputs as `UParT_*` columns
(new `discriminatorTags/Names` params of `HiInclusiveJetTableProducer`) — are
added automatically. The chain also runs the **negative-tag UParT**
(`UParTNeg_*` columns; negative-SV reco + flipped-sign tagger with the same
PbPb model and HI fallbacks) for the data-driven light-mistag control regions
— the BTVNano negative-tag pattern, but on the HI jets/model (the `@BTV`
flavour itself targets `slimmedJetsPuppi`, absent in HI MiniAOD). With
`doBtagging=True` the taginfo cone (`jet_radius`, IP `maxDeltaR`) now follows
the actual jet radius (it was silently 0.4 for all radii before).

Requires `HeavyIonsAnalysis/JetAnalysis` (forest branch) checked out in the
same area, like the HINHAD jet tables.

```bash
cmsDriver.py hin_fullhad_nano -s NANO:@HINHADSKIM --data --era Run3_pp_on_PbPb \
  --conditions 151X_dataRun3_Prompt_v1 --eventcontent NANOAOD --datatier NANOAOD \
  --filein /store/hidata/HIRun2025A/HIPhysicsRawPrime0/MINIAOD/... \
  --fileout file:hin_fullhad_nano.root -n -1
```

Reference: twin HiForest skim
(`HeavyIonsAnalysis/Configuration/test/forest_miniAOD_ParticleTransformer_run3_SKIM_FULLHAD_DATA.py`
on the `HIForest_TTBAR_Run3_2025_PbPb` forest branch), validated on 2025 PbPb
data (run 400059): ~8e-5 pass rate at >= 5 jets (all events 0-5% central),
~0.4 s/event.

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
