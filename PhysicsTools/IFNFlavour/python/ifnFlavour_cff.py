import FWCore.ParameterSet.Config as cms

# Interleaved Flavour Neutralisation (IFN) parton-flavour producers and NanoAOD
# tables. These run IN PARALLEL with the stock ghost-based flavour producers so
# the two parton-flavour definitions can be compared per-jet. Adding them is
# opt-in via addIFNFlavour(process, ...).
#
# Design (NO ghost association anywhere). IFN flavour can be defined at PARTON
# level (default) or at HADRON level via the `useHadrons` switch:
#   - GEN jets, useHadrons=False (default): cluster the physics-parton collection
#     (patJetPartons:physicsPartons) directly with IFN built on anti-kt(R). Each
#     parton enters with its real momentum and pdg-decoded net flavour. Match each
#     slimmedGenJets / slimmedGenJetsAK8 jet to the nearest IFN jet (dR < R/2) and
#     assign its signed parton flavour. (JetFlavourClusteringIFN)
#   - GEN jets, useHadrons=True            : cluster the HADRON-level final state
#     (packedGenParticles, minus neutrinos) with IFN, with heavy-hadron decays
#     SWITCHED OFF -- every final-state particle descending from a b/c hadron is
#     replaced by the parent heavy hadron (net b/c flavour).
#   - RECO jets: NO reco-level clustering. The reco IFN flavour is the flavour of
#     the matched GEN jet, obtained by dR-matching reco jets to the GEN IFN
#     flavour collection. (JetFlavourIFNTableProducer with the gen collection)
#
# Branches added (when enabled):
#   Jet_partonFlavourIFN        (reco AK4, = matched gen AK4 jet's IFN flavour)
#   GenJet_partonFlavourIFN     (gen AK4)
#   GenJetAK8_partonFlavourIFN  (gen AK8)
#
# Inputs: patJetPartons (HadronAndPartonSelector) for the physics-parton
# collection and (in hadron mode) the b/c-hadron collections; packedGenParticles
# for the hadron-level final state. All set up by the stock NanoAOD jet flavour
# wiring.

_GENPARTICLES = cms.InputTag("packedGenParticles")
_BHADRONS = cms.InputTag("patJetPartons", "bHadrons")
_CHADRONS = cms.InputTag("patJetPartons", "cHadrons")
_PARTONS = cms.InputTag("patJetPartons", "algorithmicPartons")
_DOC = "parton flavour from Interleaved Flavour Neutralisation (IFN)"


def addIFNFlavour(process, addReco=True, addGen=True, addGenAK8=True, useHadrons=False):
    """Add the IFN parton-flavour producers and NanoAOD table columns in parallel
    with the existing ghost-based flavour branches.

    useHadrons=False (default): IFN flavour is defined at PARTON level
                                (physics partons fed in via `genParticles`).
    useHadrons=True            : IFN flavour is defined at HADRON level
                                (b/c hadrons as flavour carriers, decays off;
                                 hadron-level final state fed in via `genParticles`).
    """
    from PhysicsTools.NanoAOD.jets_cff import genJetTable, genJetAK8Table

    # In parton mode `genParticles` IS the physics-parton collection; in hadron
    # mode it is the hadron-level final state.
    _GP = _GENPARTICLES if useHadrons else _PARTONS

    newModules = []

    # ---- GEN AK4: cluster with IFN, match to slimmedGenJets ----
    if addGen or addReco:
        process.genJetFlavourAssociationIFN = cms.EDProducer("JetFlavourClusteringIFN",
            jets = genJetTable.src,            # slimmedGenJets
            genParticles = _GP,
            bHadrons = _BHADRONS,
            cHadrons = _CHADRONS,
            useHadrons = cms.bool(useHadrons),
            jetAlgorithm = cms.string("AntiKt"),
            rParam = cms.double(0.4),
            deltaR = cms.double(0.2),          # match IFN jet to gen jet, R/2
            alpha = cms.double(2.0),
            omega = cms.double(1.0),
        )
        newModules.append(process.genJetFlavourAssociationIFN)

    if addGen:
        process.genJetFlavourIFNTable = cms.EDProducer("JetFlavourIFNTableProducer",
            name = cms.string("GenJet"),
            branchName = cms.string("partonFlavourIFN"),
            doc = cms.string(_DOC),
            src = genJetTable.src,
            cut = genJetTable.cut,
            deltaR = cms.double(0.1),
            jetFlavourInfos = cms.InputTag("genJetFlavourAssociationIFN"),
        )
        newModules.append(process.genJetFlavourIFNTable)

    # ---- RECO AK4: NO clustering. Reco IFN flavour = matched GEN jet flavour. ----
    # The table matches reco jets to the GEN IFN flavour collection by dR. Since
    # the gen collection is keyed by slimmedGenJets, a reco<->gen dR match here
    # directly realises "reco IFN flavour = matched gen jet's IFN flavour".
    if addReco:
        process.JetFlavourIFNTable = cms.EDProducer("JetFlavourIFNTableProducer",
            name = cms.string("Jet"),
            branchName = cms.string("partonFlavourIFN"),
            doc = cms.string(_DOC + " (from matched gen jet)"),
            src = cms.InputTag("linkedObjects", "jets"),
            cut = cms.string(""),
            deltaR = cms.double(0.2),  # reco<->gen jet match, half the AK4 cone
            jetFlavourInfos = cms.InputTag("genJetFlavourAssociationIFN"),
        )
        newModules.append(process.JetFlavourIFNTable)

    # ---- GEN AK8 ----
    if addGenAK8:
        process.genJetAK8FlavourAssociationIFN = cms.EDProducer("JetFlavourClusteringIFN",
            jets = genJetAK8Table.src,         # slimmedGenJetsAK8
            genParticles = _GP,
            bHadrons = _BHADRONS,
            cHadrons = _CHADRONS,
            useHadrons = cms.bool(useHadrons),
            jetAlgorithm = cms.string("AntiKt"),
            rParam = cms.double(0.8),
            deltaR = cms.double(0.4),          # match IFN jet to gen jet, R/2
            alpha = cms.double(2.0),
            omega = cms.double(1.0),
        )
        process.genJetAK8FlavourIFNTable = cms.EDProducer("JetFlavourIFNTableProducer",
            name = cms.string("GenJetAK8"),
            branchName = cms.string("partonFlavourIFN"),
            doc = cms.string(_DOC),
            src = genJetAK8Table.src,
            cut = genJetAK8Table.cut,
            deltaR = cms.double(0.1),
            jetFlavourInfos = cms.InputTag("genJetAK8FlavourAssociationIFN"),
        )
        newModules += [process.genJetAK8FlavourAssociationIFN, process.genJetAK8FlavourIFNTable]

    # Attach to the MC jet sequence so it runs (associations precede their tables).
    if hasattr(process, "jetMC"):
        for m in newModules:
            process.jetMC += m
    else:
        process.ifnFlavourTask = cms.Task(*newModules)
        if hasattr(process, "schedule"):
            process.schedule.associate(process.ifnFlavourTask)

    return process
