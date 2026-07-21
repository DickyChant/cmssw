# Masking pixel layers in premixing with `OverrideLayerEfficiency`

Companion recipe for the `OverrideLayerEfficiency` knob this branch adds to
`SiPixelDigitizerAlgorithm` (commit fe1d63e). Full study write-up with plots:
<https://sqian.web.cern.ch/sqian/mask_first_bpix/>

## The knob

```python
OverrideLayerEfficiency = cms.untracked.vdouble(b1, b2, b3, b4, f1, f2, f3)
```

Indices 0..3 = BPix1..BPix4, 4..6 = FPix1..FPix3. Per entry:

| value | meaning                                                          |
| ----- | ---------------------------------------------------------------- |
| `< 0` | use the `SiPixelDynamicInefficiency` payload polynomial (default) |
| `= 0` | kill the layer/disk (`pu_scale` replaced by 0)                    |
| `> 0` | replace `pu_scale` by this fixed value (`1.0` = no PU inefficiency) |

The parameter is untracked and defaults to seven `-1`s, so existing configs
are bit-identical. It is applied uniformly at all four `pu_scale` use sites
(BPix cfg path, FPix inner/outer cfg paths, DB path).

## Where to set it

| workflow                        | PSet                              |
| ------------------------------- | --------------------------------- |
| classic mixing / premix stage 1 | `process.mix.digitizers.pixel`    |
| premix stage 2 (DATAMIX)        | `process.mixData.workers.pixel`   |

## Premixing, done properly

Standard premixing applies pixel inefficiency **once, at stage 2**, on the
combined signal+PU charge map: the `premix_stage1` procModifier sets
`AddPixelInefficiency = False` (the library stores raw charge sums), and the
stage-2 pixel worker re-injects the library charge via `setSimAccumulator()`
before the normal digitization (noise, dynamic inefficiency via
`calculateInstlumiFactor(PileupSummaryInfo)`, thresholds) runs.

Two ways to mask with the override:

1. **PU-only masking (library side — what the study used).** At stage 1,
   re-enable `AddPixelInefficiency` and set the override on
   `process.mix.digitizers.pixel`. BPix1 is killed in the PU library while the
   signal, digitized at stage 2, keeps its BPix1 hits. No stage-2
   customisation needed. Caveat: `-1` entries make stage 1 apply the payload
   polynomial to the *other* layers too, which stage 2 then applies again on
   signal+PU — fine for an A/B comparison of libraries, but use
   `0.0,1,1,1,1,1,1` if only BPix1 should be touched.
2. **signal+PU masking (combined, stage 2).** Set the override on
   `process.mixData.workers.pixel` — nothing else, `AddPixelInefficiency` is
   already `True` at stage 2. Equivalent in spirit to masking in classic
   mixing.

## Scripts

| script                        | what it does                                                        |
| ----------------------------- | ------------------------------------------------------------------- |
| `run_premix_stage1.sh`        | builds baseline + BPix1-masked PREMIX libraries (parallel cmsRun)    |
| `run_premix_stage2_compare.sh`| same signal vs both libraries through `DIGI,DATAMIX,…,RECO`; optional `STAGE2_MASK` variant for combined masking |
| `verify_premix.py`            | FWLite per-layer `simSiPixelDigis` counts in the libraries           |
| `verify_stage2_compare.py`    | FWLite per-layer `siPixelClusters` + `generalTracks` counts at RECO  |

Quick start:

```bash
export RELEASE_DIR=/path/to/CMSSW_14_0_X/src     # built from this branch
export PU_FILELIST=$PWD/pu_files.txt             # MinBias GEN-SIM list
./run_premix_stage1.sh
SIGNAL_GS=/path/to/gs_Mu.root ./run_premix_stage2_compare.sh
python3 verify_premix.py premix_study/data/premix_*.root
python3 verify_stage2_compare.py premix_study/data/stage2_*.root
```

The customisations can also be attached via the module added in
`SimTracker/SiPixelDigitizer/python/customizePixelLayerEff.py`
(after `scram b python`), e.g.

```
--customise SimTracker/SiPixelDigitizer/customizePixelLayerEff.maskBPix1PU_stage1
--customise SimTracker/SiPixelDigitizer/customizePixelLayerEff.maskBPix1_stage2
```

## Motivation

The `SiPixelDynamicInefficiency_phase1_2023_v2_fix3` payload's BPix1
polynomials cross zero around PU ≈ 100–110 and are multiplied into
`columnEfficiency` without a clamp, silently killing BPix1 at digi time at
high PU. This knob makes the effect reproducible (set `0`), tunable (set a
fraction), or removable (set `1.0`) per layer, in both classic mixing and
premixing. See also the driver study repo:
<https://gitlab.cern.ch/sqian/time_event_scan> (branch `high_PU_RSS`).
