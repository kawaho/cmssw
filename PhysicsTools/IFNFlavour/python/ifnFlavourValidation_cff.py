import FWCore.ParameterSet.Config as cms

# Side-by-side validation of the IFN parton flavour against the stock ghost
# parton/hadron flavour. Adds an EDAnalyzer that writes:
#   - per-IFN-category TH2D of (ghost partonFlavour vs gen jet pt)
#   - per-IFN-category TH2D of (ghost hadronFlavour vs gen jet pt)
#   - per-flavour TH2D eta-phi maps for both gen jets (coloured by ghost
#     partonFlavour) and IFN jets (coloured by IFN net flavour)
#   - a TTree with both jet collections per event for any further slicing
# All into a TFileService output file.
#
# This is purely an analyzer; it does NOT change the NanoAOD payload. Use it
# alongside addIFNFlavour() to get the IFN flavour info collection AND the IFN
# jets (BasicJetCollection + flat netFlavour vector) that the analyzer reads --
# the producer publishes both, so no reclustering happens here.


def addIFNFlavourValidation(process,
                            outputFile="ifn_validation.root",
                            genJets=None,
                            ghostFlavourInfos=None,
                            ifnFlavourInfos=None,
                            ifnJets=None,
                            ifnJetNetFlavour=None,
                            strict=False,
                            ptMin=10.0):
    """Wire up the IFNFlavourValidator EDAnalyzer + TFileService.

    Defaults target the AK4 gen path that addIFNFlavour() sets up:
      genJets          = slimmedGenJets
      ghostFlavourInfos= genJetFlavourAssociation         (ghost-based)
      ifnFlavourInfos  = genJetFlavourAssociationIFN      (IFN per-gen-jet flav)
      ifnJets          = genJetFlavourAssociationIFN:ifnJets   (IFN jets)
      ifnJetNetFlavour = genJetFlavourAssociationIFN:ifnJetNetFlavour
                         (flat vector<int>, length 7 * nIFNJets;
                          jet j's netFlavour[k] at 7*j + k)
    """
    from PhysicsTools.NanoAOD.jets_cff import genJetTable

    if genJets is None:
        genJets = genJetTable.src                                  # slimmedGenJets
    if ghostFlavourInfos is None:
        ghostFlavourInfos = cms.InputTag("genJetFlavourAssociation")
    if ifnFlavourInfos is None:
        ifnFlavourInfos = cms.InputTag("genJetFlavourAssociationIFN")
    if ifnJets is None:
        ifnJets = cms.InputTag("genJetFlavourAssociationIFN", "ifnJets")
    if ifnJetNetFlavour is None:
        ifnJetNetFlavour = cms.InputTag("genJetFlavourAssociationIFN", "ifnJetNetFlavour")

    genEvent = cms.InputTag("generator")

    process.TFileService = cms.Service("TFileService",
        fileName = cms.string(outputFile),
        closeFileFast = cms.untracked.bool(True),
    )

    process.ifnFlavourValidator = cms.EDAnalyzer("IFNFlavourValidator",
	generator = genEvent,
        genJets = genJets,
        ghostFlavourInfos = ghostFlavourInfos,
        ifnFlavourInfos = ifnFlavourInfos,
        ifnJets = ifnJets,
        ifnJetNetFlavour = ifnJetNetFlavour,
        strict = cms.bool(False),
        ptMin = cms.double(ptMin),
    )

    # The analyzer needs the IFN producer's three products already produced this
    # event; the ghost association too. Hang off a new EndPath -- jetMC has run.
    process.ifnFlavourValidation_step = cms.EndPath(process.ifnFlavourValidator)
    if hasattr(process, "schedule"):
        process.schedule.append(process.ifnFlavourValidation_step)

    return process
