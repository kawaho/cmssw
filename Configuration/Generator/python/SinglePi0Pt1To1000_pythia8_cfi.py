import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8PtGun",
                         PGunParameters = cms.PSet(
        ParticleID = cms.vint32(111),
        AddAntiParticle = cms.bool(False),
        MaxEta = cms.double(4.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-4.5),
        MinPt = cms.double(1),
        MinPhi = cms.double(-3.14159265359), ## in radians
        MaxPt = cms.double(1000.0)
        ),
                         Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
                         psethack = cms.string('single pi0 E 10'),
                         firstRun = cms.untracked.uint32(1),
                         PythiaParameters = cms.PSet(parameterSets = cms.vstring())
                         )
