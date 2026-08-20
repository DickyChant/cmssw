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

# CE-H silicon, 21 layers. Only the electromagnetic tail that punches through the
# back of CE-E lands here (0.07% of a 50 GeV photon, rising with energy), but
# without these entries the Gamma is truncated at layer 26 and detector 9 stays
# empty for photons.
#
# Depth from the geometry: the CE-H layer pitch is 6.3 cm (z = 368.0, 374.3,
# 380.6, ...) against 1.757 cm X0 for the stainless absorber, so a steel-filled
# layer is ~3.6 X0. These are NOT a tuning handle for the punch-through: the
# first CE-H bin starts at the fixed CE-E total depth (25.63 X0), so it collects
# the whole Gamma tail no matter how the entries are chosen.
#
# MEASURED, 50 GeV photons vs the V16 FullSim reference (target 0.071% of the
# shower in CE-H):
#
#   sigmaLnAlpha = 0.419 (tuned)   ->  1.04%   15x too much
#   sigmaLnAlpha = 0               ->  0.169%  2.4x too much
#
# The second number matches the analytic tail Q(alpha=5.68, beta*t=15.2) = 0.17%,
# which confirms the mechanism. So the punch-through is over-predicted twice
# over: the mean-parameter Gamma tail is itself 2.4x too fat, and the log-normal
# alpha fluctuation inflates it a further 6x because Q is strongly convex in
# 1/alpha, so low-alpha events dominate the deep tail.
#
# This is a diagnostic of the longitudinal model, not of CE-H. The CE-E-only
# profile comparison could not see it -- CE-E layer fractions are dominated by
# the peak, where a fat tail costs almost nothing. Fixing it means refitting
# sigmaLnAlpha (and rhoLnAlphaT, currently +0.51) against the punch-through and
# the CE-E spread jointly, rather than damping the tail here.
_cehFirstX0 = 3.6    # first CE-H layer (the CE-E/CE-H gap is not resolved here)
_cehLayerX0 = 3.2    # per CE-H layer thereafter
_layerX0 += [_cehFirstX0] + [_cehLayerX0] * 20

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
        # A conversion-convolved refit against the FullSim absorber profile was
        # tried (analysis/refit_longitudinal_raw.py). It returns
        # 0.5704/-0.6898, 1.0655/-2.1236, reassuringly close to the 2021 ELECTRON
        # values -- but it fits the profile poorly (residual ~0.21) and is clearly
        # WORSE in the per-layer comparison: the ratio spread goes 0.337 -> 0.730
        # and the first core layer 1.68 -> 3.09. Kept the values below, which are
        # the better description empirically. The refit needs a better target than
        # a single convolved Gamma before it can replace these.
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
        # Sample the photon conversion depth explicitly, -ln(u) * 9/7 X0.
        #
        # It was first assumed that conversion was already inside tConst, since
        # the T fit was made on photons, and the explicit term was left off. The
        # per-layer comparison against FullSim says otherwise: without it the
        # shower peaks at layer 7 instead of 9-11, about 1.2 X0 early -- close to
        # the 1.29 X0 mean conversion depth. Turning it on halves the spread of
        # the per-layer ratio (0.574 -> 0.337) and fixes layer 1 (1.71 -> 0.79),
        # while leaving the total unchanged (1.017 -> 1.013).
        applyPhotonConversion = cms.bool(True),

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
    # Cross-checked on PIONS: the LD200 per-crossing MPV agrees with photons to
    # 2% (59.9 vs 61.1 keV), so these constants are a sensor property and are not
    # photon-specific. Pions do give lower crossing multiplicity (2.33 vs 3.99)
    # and much wider showers (10x more LD300 hits), but multiplicity is derived
    # per cell from the cell energy here, so it adapts on its own.
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

##############################################################################
# Hadronic showers, measured on the V16/D110 single-pion samples.
#
# These are not the electromagnetic numbers with different values: the shower
# maximum sits at layer 25 (the CE-E/CE-H boundary) against layer 9 for photons,
# CE-H carries 20% of the visible energy at 5 GeV rising to 49% at 500 GeV, and
# the silicon sampling fraction is roughly half the electromagnetic one.
##############################################################################

hgcalHadronParameters = cms.PSet(
    criticalEnergy = cms.double(10.48e-3),   # GeV, only used to form y = E/E_crit
    spotEnergy     = cms.double(0.010),
    maxSpots       = cms.uint32(2000),

    Longitudinal = cms.PSet(
        # Gamma in LAYER units over the full 47-layer stack. Layer units rather
        # than interaction lengths because CE-E and CE-H differ in material, so a
        # single lambda scale would not span both; the geometry supplies the
        # positions. Fitted on pions: residual ~0.16-0.19.
        alphaSlope   = cms.double(0.5350),
        alphaConst   = cms.double(-1.1245),
        tSlope       = cms.double(2.6730),
        tConst       = cms.double(-5.5673),
        # Hadronic showers fluctuate more than electromagnetic ones.
        sigmaLnAlpha = cms.double(0.35),
        sigmaLnT     = cms.double(0.30),
    ),

    # Visible-energy suppression relative to an EM shower of the same energy.
    # Measured 0.623 / 0.557 / 0.548 at 5 / 50 / 500 GeV -- energy dependent, so
    # it cannot be folded into a single constant.
    ehSlope = cms.double(-0.0163),
    ehConst = cms.double(0.7239),

    Transverse = cms.PSet(
        # Hadronic showers are much wider than electromagnetic ones (ten times the
        # LD300 hits at the same energy), so the radii are correspondingly larger.
        r68Slope     = cms.double(2.20),
        r68Const     = cms.double(1.20),
        coreOverR68  = cms.double(0.62),
        tailOverCore = cms.double(4.0),
        coreFraction = cms.double(0.60),
    ),
)

# The block CalorimetryManager reads.
HGCalBlock = cms.PSet(
    simulateHGCal = cms.bool(True),
    HGCalCalorimeterProperties = hgcalCalorimeterProperties,
    ShowerParameters = hgcalShowerParameters,
    ReverseCalibration = hgcalReverseCalibration,
    HadronParameters = hgcalHadronParameters,
)
