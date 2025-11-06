import FWCore.ParameterSet.Config as cms
generator = cms.EDFilter("Pythia8PtGun",
                         PGunParameters = cms.PSet(
        MaxPt = cms.double(1000.),
        MinPt = cms.double(1.),
        ParticleID = cms.vint32(22),
        AddAntiParticle = cms.bool(False), # false in pythia6 version but photons are their own anti-particle
        MaxEta = cms.double(4.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-4.5),
        MinPhi = cms.double(-3.14159265359) ## in radians
        ),
                         Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
                         psethack = cms.string('single gamma pt 1 to 1000'),
                         firstRun = cms.untracked.uint32(1),
                         PythiaParameters = cms.PSet(parameterSets = cms.vstring())
                         )
