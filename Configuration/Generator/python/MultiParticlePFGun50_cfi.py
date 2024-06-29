import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatEtaRangeGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32([11, -11, 13, -13, 211, -211, 22, 130, -130, 111, -111]),
        MinEta = cms.double(-5),
        MaxEta = cms.double(5),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(1.0),
        MaxE = cms.double(1000.0),
        nParticles = cms.int32(50),
        exactShoot = cms.bool(False),
        randomShoot = cms.bool(True),
    ),
    AddAntiParticle = cms.bool(False),
    debug=cms.untracked.bool(True),
)
