"""cmsDriver customisation hooks for the OverrideLayerEfficiency knob.

This branch (private_dev/customize_pixel_layer_eff) adds an optional,
untracked vdouble to SiPixelDigitizerAlgorithm:

    OverrideLayerEfficiency = cms.untracked.vdouble(b1, b2, b3, b4, f1, f2, f3)

with indices 0..3 = BPix1..BPix4 and 4..6 = FPix1..FPix3.  Per entry:

    < 0   use the SiPixelDynamicInefficiency payload polynomial (default)
    = 0   kill the layer/disk (pu_scale replaced by 0)
    > 0   replace pu_scale by this fixed value (1.0 = no PU inefficiency)

The PSet the knob must be set on depends on the workflow:

    classic mixing / premix stage 1 : process.mix.digitizers.pixel
    premix stage 2 (DATAMIX)        : process.mixData.workers.pixel

Stage-1 caveat: the premix_stage1 procModifier disables AddPixelInefficiency
(inefficiency is normally applied once, at stage 2, on the combined signal+PU
charge map).  To mask a layer in the PU library itself, AddPixelInefficiency
must be re-enabled at stage 1 -- the *_stage1 helpers below do both.

Study write-up: https://sqian.web.cern.ch/sqian/mask_first_bpix/
"""

import os
import FWCore.ParameterSet.Config as cms

MASK_BPIX1 = [0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
# Variant that masks BPix1 *without* applying the payload polynomial to the
# other layers at stage 1 (avoids double-application of the few-% dynamic
# inefficiency there, since stage 2 applies it again on signal+PU):
MASK_BPIX1_ONLY = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


def _valuesFromEnv(default):
    raw = os.environ.get("OVERRIDE_LAYER_EFF", "").strip()
    if not raw:
        return list(default)
    values = [float(x) for x in raw.split(",")]
    if len(values) != 7:
        raise ValueError(
            "OVERRIDE_LAYER_EFF needs 7 comma-separated entries "
            "(BPix1-4, FPix1-3), got %d: %r" % (len(values), raw))
    return values


def _stage1(process, values):
    pixel = process.mix.digitizers.pixel
    pixel.AddPixelInefficiency = cms.bool(True)  # premix_stage1 turns this off
    pixel.OverrideLayerEfficiency = cms.untracked.vdouble(*values)
    print("customizePixelLayerEff: mix.digitizers.pixel.OverrideLayerEfficiency =", values)
    return process


def _stage2(process, values):
    # AddPixelInefficiency is already True in the stage-2 pixel worker
    pixel = process.mixData.workers.pixel
    pixel.OverrideLayerEfficiency = cms.untracked.vdouble(*values)
    print("customizePixelLayerEff: mixData.workers.pixel.OverrideLayerEfficiency =", values)
    return process


# --- entry points for cmsDriver --customise ---------------------------------

def maskBPix1PU_stage1(process):
    """Premix stage 1: kill BPix1 in the PU library, payload polynomial elsewhere."""
    return _stage1(process, MASK_BPIX1)


def maskBPix1Only_stage1(process):
    """Premix stage 1: kill BPix1 only, other layers left at efficiency 1."""
    return _stage1(process, MASK_BPIX1_ONLY)


def maskBPix1_stage2(process):
    """Premix stage 2: kill BPix1 for signal+PU combined."""
    return _stage2(process, MASK_BPIX1)


def overrideFromEnv_stage1(process):
    """Premix stage 1, vector from OVERRIDE_LAYER_EFF.  The env var must also
    be set at cmsRun time -- the customise re-runs when the config is loaded."""
    return _stage1(process, _valuesFromEnv(MASK_BPIX1))


def overrideFromEnv_stage2(process):
    """Premix stage 2, vector from OVERRIDE_LAYER_EFF (see stage-1 note)."""
    return _stage2(process, _valuesFromEnv(MASK_BPIX1))
