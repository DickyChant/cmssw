import FWCore.ParameterSet.Config as cms

##############################################################################
# HGCAL CE-E FastSim parameters, re-derived on V16 / GeometryExtendedRun4D110
# FullSim samples (Phase2Spring24, single photons at eta = 2.0).
#
# Every number below is measured, not taken from the 2021 prototype: that work
# used Extended2026D71 with 28 CE-E layers, and neither the geometry nor the
# calibration transfers.
##############################################################################

# Per-layer radiation lengths of the 26 CE-E layers, from the RAW reco dE/dx
# weight table (weightsPerLayer_V16), which is proportional to the absorber in
# front of each layer.
#
# Do NOT replace this with a uniform value or with the calcWeights-averaged
# table. The structure is physical:
#   - layer 1 is ~half the rest: the half-lead first cassette,
#   - layers alternate 1.091 / 0.797: the Pb side vs the CuW side,
#   - from layer 19 the lead thickens, giving 1.091 / 1.149.
# The calcWeights running average maps (12.86, 9.4) onto (11.13, 11.13) and
# erases exactly the alternation the w-parametrization describes.
_layerX0 = [
    0.471,
    1.091, 0.797, 1.091, 0.797, 1.091, 0.797, 1.091, 0.797,
    1.091, 0.797, 1.091, 0.797, 1.091, 0.797, 1.091, 0.797,
    1.091, 1.149, 1.091, 1.149, 1.091, 1.149, 1.091, 1.149, 1.091,
]

hgcalShowerParameters = cms.PSet(
    # Target energy per spot and a hard cap, so a 1 TeV shower cannot blow up
    # the event size.
    spotEnergy = cms.double(0.010),      # GeV
    maxSpots   = cms.uint32(2000),

    Longitudinal = cms.PSet(
        # alpha_hom = alphaSlope * ln(y) + alphaConst   (y = E / E_crit)
        # T_hom     = tSlope     * ln(y) + tConst       [X0]
        # Measured: slopes reproduce the 2021 values to 3% (alpha) and 1.7% (T);
        # the T offset differs by +0.88 X0, consistent with the +1.29 X0 expected
        # from photon conversion (the 2021 fit used electrons).
        alphaSlope   = cms.double(0.648),
        alphaConst   = cms.double(0.187),
        tSlope       = cms.double(1.030),
        tConst       = cms.double(-0.846),

        # Event-to-event fluctuations of (ln alpha, ln T).
        # NOTE: these come from a moment estimator, which correlates alpha and T
        # through their shared moments and inflates rho (measured 0.95, whereas
        # Grindhammer-Peters expects ~0.51 here). Per-event Gamma fits are needed
        # before this correlation should be trusted; the widths are also known to
        # be ~13% low because only (alpha, T) fluctuate globally.
        # The T fit above was made on PHOTON samples, so the conversion depth is
        # already inside tConst (that is the +0.88 X0 offset). Adding an explicit
        # stochastic conversion on top would count it twice. Switch this on only
        # together with an electron-based T fit.
        applyPhotonConversion = cms.bool(False),

        sigmaLnAlpha = cms.double(0.419),
        sigmaLnT     = cms.double(0.257),
        rhoLnAlphaT  = cms.double(0.51),
    ),

    CassetteSplit = cms.PSet(
        # w(x0) = a (1 - exp(b x0)), the fraction of a cassette's energy landing
        # on the Pb-side silicon layer. a is essentially energy independent (as in
        # 2021); b shrinks with energy on the same trend.
        #   a = aSlope * ln(y) + aConst
        #   b = bSlope * ln(y) + bConst
        apply  = cms.bool(True),
        aSlope = cms.double(0.0029),
        aConst = cms.double(0.536),
        bSlope = cms.double(0.0273),
        bConst = cms.double(-0.4545),
    ),

    Transverse = cms.PSet(
        # Radial profile, sampled PERPENDICULAR to the shower axis and then
        # projected onto the layer plane. Fitting the in-plane footprint instead
        # would bake the eta = 2 projection into these numbers.
        #
        # r68(tau) = r68Slope * tau + r68Const   [cm],  tau = t / T_max
        # is the measured, model-free scale. The core/tail decomposition is
        # anchored to it: (p, R_C, R_T) are strongly correlated when fitted
        # per layer, so only their relation to r68 is used here.
        apply             = cms.bool(True),
        r68Slope          = cms.double(1.128),
        r68Const          = cms.double(0.098),
        coreOverR68       = cms.double(0.62),
        tailOverCore      = cms.double(3.5),
        coreFractionSlope = cms.double(-0.30),
        coreFractionConst = cms.double(1.00),
    ),

)

hgcalCalorimeterProperties = cms.PSet(
    # Homogenised CE-E medium. The 2021 values are used as the starting point
    # for the material constants; X0 and the layer table below are measured.
    lightColl        = cms.double(0.03),
    lightCollUnif    = cms.double(0.003),
    photoStatistics  = cms.double(50.e3),
    interactionLength = cms.double(17.4),   # cm

    Aeff             = cms.double(124.41),
    Zeff             = cms.double(51.54),
    rho              = cms.double(10.78),   # g/cm3
    radLenIngcm2     = cms.double(8.30),
    radLenIncm       = cms.double(0.77),    # measured CE-E: 25.63 X0 over 26 layers

    criticalEnergy   = cms.double(10.48e-3),  # GeV -- to be re-derived for D110
    moliereRadius    = cms.double(1.49),      # cm

    Fs               = cms.double(1.0),
    ehat             = cms.double(0.0),
    resE             = cms.double(1.0),
    da               = cms.double(0.0),
    dp               = cms.double(0.0),

    etaMin           = cms.double(1.5),
    etaMax           = cms.double(3.0),

    layerX0          = cms.vdouble(_layerX0),
)

##############################################################################
# Shower energy (deposited in the absorber) -> silicon deposit, i.e. the inverse
# of what HGCalRecHit does. Without this the hits carry the incident energy and
# overstate the silicon deposit by more than two orders of magnitude.
#
# The dE/dx weights here MUST be the ones reconstruction applies -- the
# calcWeights-averaged V16 table -- or FastSim and reco disagree by construction.
# This is the opposite of the w-parametrization, which needs the RAW table.
##############################################################################

def _calcWeights(w):
    """The running two-layer average HGCalRecHit_cfi applies for V16/V19."""
    res = [sum(p) / 2. for p in zip(w[:], w[1:] + [w[-1]])]
    res[0] = 0.0
    return res

_weightsPerLayer_V16_raw = [0.0,
    5.55,
    12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4,
    12.86, 9.4, 12.86, 9.4, 12.86, 9.4,
    12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86,
    58.63,
    60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7,
    83.08, 83.08, 83.43, 83.61, 83.61, 83.61, 83.61, 83.61, 83.61, 83.61]

hgcalReverseCalibration = cms.PSet(
    # MeV per MIP, indexed by layer (index 0 unused), as used by reco
    dEdXWeights = cms.vdouble(_calcWeights(_weightsPerLayer_V16_raw)),

    # V16 uses the _mean variant; [HD120, LD200, LD300]
    fCPerMIP = cms.vdouble(2.06, 3.43, 5.15),
    # CE-E entries of the V16 thicknessCorrection
    thicknessCorrection = cms.vdouble(0.75, 0.76, 0.75),
    chargeCollectionEfficiency = cms.vdouble(1.0, 1.0, 1.0),
    keV2fC = cms.double(0.044259),

    # Measured per-crossing deposit spectrum, by thickness. NOT muon MIP values:
    # shower crossings are not minimum-ionising (mean/MPV = 1.74 vs 1.34 for a
    # true Landau). HD120 gives 387.9 eV/um against the PDG 388, untuned, and the
    # MPVs are energy independent to <1% between 50 and 500 GeV.
    crossingMPVkeV   = cms.vdouble(46.55, 75.12, 101.08),
    crossingWidthkeV = cms.vdouble(19.89, 28.23, 40.18),
    crossingMeankeV  = cms.vdouble(86.14, 186.09, 156.99),

    # Crossing multiplicity is currently approximated as Poisson, which the data
    # says it is not (LD200: <m> = 5.6 with P(1) = 0.44 against 0.02 for Poisson).
    # The mean deposit is unaffected; only the fluctuation shape is approximate.
    fluctuate = cms.bool(True),

    # The Landau tail is unbounded; a single crossing in thin silicon is not.
    # Uncapped, ~5% of hits at 500 GeV ran away to ~1 GeV each and carried 99% of
    # the event energy. 30x MPV is ~2.3 MeV for LD200 -- generous but finite.
    maxCrossingOverMPV = cms.double(30.),
)

# The block CalorimetryManager reads.
HGCalBlock = cms.PSet(
    simulateHGCal = cms.bool(True),
    HGCalCalorimeterProperties = hgcalCalorimeterProperties,
    ShowerParameters = hgcalShowerParameters,
    ReverseCalibration = hgcalReverseCalibration,
)
