import FWCore.ParameterSet.Config as cms

from FastSimulation.SimplifiedGeometryPropagator.TrackerMaterial_cfi import TrackerMaterialBlock 

#############
### Hack to interface "old" calorimetry with "new" propagation in tracker
#############

CaloMaterialBlock = cms.PSet(
    CaloMaterial = cms.PSet(
        maxRadius = cms.untracked.double(500.),
        maxZ = cms.untracked.double(1200.),
        
        ######
        # The calorimetry
        # Positions used from old ParticlePropagator. Do not really agree with the CMS ECAL/HCAL TDR values...
        ######

        # Coverage usually provided as eta, e.g. barrel ECAL abs(eta) < 1.479
        # Use definition of pseurorapidity: theta = 2*arctan(e^-eta)
        # And theta = cos(z/sqrt(R^2+z^2))
        # Better: eta = -0.5*ln((1-cos(theta))/(1+cos(theta))) 
        BarrelLayers = cms.VPSet(
            ########### ECAL ###########
            cms.PSet(
                radius = cms.untracked.double(129.0),
                limits = cms.untracked.vdouble(0.0, 268.4),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("ECAL")
            ),
            ########### ECAL (barrel cut corner) ###########
            cms.PSet(
                radius = cms.untracked.double(152.6),
                limits = cms.untracked.vdouble(268.4, 320.9),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("ECAL")
            ),
            ########### HCAL ###########
            cms.PSet(
                radius = cms.untracked.double(177.5),
                limits = cms.untracked.vdouble(0.0, 335.0),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("HCAL")
            ),
            ########### HCAL (barrel cut corner) ###########
            cms.PSet(
                radius = cms.untracked.double(300.0),
                limits = cms.untracked.vdouble(335.0, 400.458),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("HCAL")
            ),
            ########### Acts as end of detector to speed up simulation ###########
            cms.PSet(
                radius = cms.untracked.double(400.0),
                limits = cms.untracked.vdouble(0., 1110.0),
                thickness = cms.untracked.vdouble(0.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("VFCAL")
            ),
        ),

        EndcapLayers = cms.VPSet(
            ########### PreShowerLayer1 ###########
            cms.PSet(
                z = cms.untracked.double(303.353),
                limits = cms.untracked.vdouble(45., 125.),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("PRESHOWER1")
            ),
            ########### PreShowerLayer2 ###########
            cms.PSet(
                z = cms.untracked.double(307.838),
                limits = cms.untracked.vdouble(45., 125.),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("PRESHOWER2")
            ),
            ########### ECAL ###########
            cms.PSet(
                z = cms.untracked.double(320.9),
                limits = cms.untracked.vdouble(32.0, 152.6),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("ECAL")
            ),
            ########### HCAL ###########
            cms.PSet(
                z = cms.untracked.double(400.458),
                limits = cms.untracked.vdouble(39.9, 300.),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("HCAL")
            ),
            ########### VFCAL ###########
            cms.PSet(
                z = cms.untracked.double(1110.0),
                limits = cms.untracked.vdouble(12.2, 110.9),
                thickness = cms.untracked.vdouble(1.),
                interactionModels = cms.untracked.vstring(),
                caloType = cms.untracked.string("VFCAL")
            ),
        )
    )
)
    
if hasattr(TrackerMaterialBlock.TrackerMaterial, 'magneticFieldZ'):
    CaloMaterialBlock.CaloMaterial.magneticFieldZ = TrackerMaterialBlock.TrackerMaterial.magneticFieldZ

#############
### Phase-2: HGCAL replaces the endcap ECAL, the preshower and the endcap HCAL.
### Leaving the Run-2 endcap layers in place would flag particles as onEcal()==2
### and send them into the PbWO4 + preshower path, whose geometry is null in
### Phase-2. So the whole endcap set is replaced.
###
### Geometry (V16 / GeometryExtendedRun4D110), from the DetId LUT:
###   CE-E front face  z = 322.1 cm,  CE-H back face  z = 513.3 cm
### The radial limits follow 1.5 < |eta| < 3.0 at the CE-E front face:
###   r = z / sinh(eta):  eta=3.0 -> 32.1 cm,  eta=1.5 -> 151.3 cm
#############

from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal

_hgcalEndcapLayers = cms.VPSet(
    ########### HGCAL CE-E front face ###########
    cms.PSet(
        z = cms.untracked.double(322.1),
        limits = cms.untracked.vdouble(31.5, 152.0),
        thickness = cms.untracked.vdouble(1.),
        interactionModels = cms.untracked.vstring(),
        caloType = cms.untracked.string("HGCAL")
    ),
    ########### VFCAL (end of detector; keeps the propagation terminating) ###########
    cms.PSet(
        z = cms.untracked.double(1110.0),
        limits = cms.untracked.vdouble(12.2, 110.9),
        thickness = cms.untracked.vdouble(1.),
        interactionModels = cms.untracked.vstring(),
        caloType = cms.untracked.string("VFCAL")
    ),
)

phase2_hgcal.toModify(CaloMaterialBlock.CaloMaterial, EndcapLayers = _hgcalEndcapLayers)
    
