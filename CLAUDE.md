# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a CMSSW (CMS Software) EDAnalyzer package for reconstructing rare quarkonia decay final states from CMS MINIAOD data:
- **JpsiJpsiPhi**: J/ψ(→μμ) + J/ψ(→μμ) + φ(→KK) → 4μ + 2K
- **JpsiJpsiUps**: J/ψ(→μμ) + J/ψ(→μμ) + Υ(→μμ) → 6μ
- **JpsiUpsPhi**: J/ψ(→μμ) + Υ(→μμ) + φ(→KK) → 4μ + 2K

The codebase was refactored in March 2026 to modularize the analysis pipeline and externalize all selection cuts to runtime configuration.

## Common Development Commands

### Building
```bash
cd $CMSSW_BASE/src
cmsenv
scram b -j 4
```

### Running Locally
```bash
# Data (with era-specific global tag)
cmsRun HeavyFlavorAnalysis/TPS-Onia2MuMu/test/ConfFile_cfg.py \
    inputFiles=file:myData.root \
    outputFile=output.root \
    era=Run2023C

# MC (via ConfFile_cfg.py VarParsing)
cmsRun HeavyFlavorAnalysis/TPS-Onia2MuMu/test/ConfFile_cfg.py \
    inputFiles=file:myMC.root \
    outputFile=mc_output.root \
    runOnMC=True \
    era=Run2022

# MC single-object efficiency ntuple (skip composite building)
cmsRun HeavyFlavorAnalysis/TPS-Onia2MuMu/test/ConfFile_cfg.py \
    inputFiles=file:myMC.root \
    outputFile=mc_singles.root \
    runOnMC=True \
    era=Run2022 \
    keepAllSingleObjectCandsInMC=True \
    skipCompositeCandBuildingWhenKeepingSingles=True

# Switch analysis mode
cmsRun HeavyFlavorAnalysis/TPS-Onia2MuMu/test/ConfFile_cfg.py \
    analysisMode=JpsiUpsPhi
```

### CRAB Submission
The `test/crabData/` directory contains CRAB job submission infrastructure (git submodule to CRAB-Tool):
```bash
cd test/crabData
git submodule update --init --recursive
python3 generate_crab_configs.py --analysis-mode JpsiJpsiPhi --use-proxy
./submit.sh
```

## Architecture

### Main Analyzer: MultiLepPAT

The core analysis is in `MultiLepPAT` (`interface/MultiLepPAT.h`, `src/MultiLepPAT.cc`). The refactored `analyze()` method follows this modular pipeline:

1. **processMCGenInfo()** - Extract gen-level particles for MC truth matching
2. **processHLTInfo()** - Record HLT trigger decisions
3. **reconstructPrimaryVertex()** - Select primary vertex with configurable strategy
4. **fillMuonBlock()** - Store muon kinematics and ID variables
5. **pairMuons()** - Build dimuon candidates (J/ψ or Υ) via kinematic vertex fit
6. **pairTracks()** - Build track-pair candidates (φ or other mesons)
7. **storeAllSingleObjectCandidatesForMC()** - Store single-object J/ψ, Υ, φ candidates and the RecoKaonTrack block (MC efficiency studies)
8. **combineCandidates()** - Combine resonances into 3-body candidates
9. **clearEventData()** - Clear event-level storage

Pairs are built at most once per event. Capability flags (`canTrySingleJpsi`, `canTryFullMuonSide`, etc.) and `hasFullCandidateInputs()` gate whether each step runs. `shouldFillCurrentEvent()` controls TTree fill based on candidate presence, MC mode, and `RequireAcceptedCandidatesForMonteCarloTree`.

### Configuration Architecture

All selection cuts are externalized via `StringCutObjectSelector` syntax, enabling runtime modification without recompilation:

- `MuonSelection` - String cut on `pat::Muon` members (e.g., `"pt > 2.5 && abs(eta) < 2.4"`)
- `TrackSelection` - String cut on `pat::PackedCandidate` for kaons/pions
- Mass windows (`JpsiMassMin/Max`, `UpsMassMin/Max`, `PhiMassMin/Max`)
- Candidate pT/eta pre-cuts (`JpsiCandPtMin`, `JpsiCandEtaMax`, `UpsCandPtMin`, `UpsCandEtaMax`, `PhiCandPtMin`, `PhiCandEtaMax`)
- Vertex probability cuts (`OniaDecayVtxProbCut`, `PriVtxProbCut`)
- Primary vertex quality (`PVNdofMin`, `PVMaxAbsZ`, `PVMaxRho`)
- Vertex fit toggles (`DoJpsiDecayVtxFit`, `DoUpsDecayVtxFit`, `DoPhiDecayVtxFit`, `DoDiOniaVtxFit`, `DoPriVtxFit`) — enable/disable individual kinematic fits

### MC Control (VarParsing CLI only, not in canned configs)

- `runOnMC` — enables `DoMonteCarloTree` (bool, default False)
- `keepAllSingleObjectCandsInMC` — store all single-object J/ψ, Υ, φ candidates and RecoKaonTrack block (default False)
- `skipCompositeCandBuildingWhenKeepingSingles` — skip `combineCandidates()` after storing singles (default False)
- `requireAcceptedCandidatesForMonteCarloTree` — only fill TTree entries with ≥1 accepted candidate (default False)

### Multi-threading (VarParsing CLI)

- `numThreads` / `numStreams` — CMSSW multi-threading (default: 1 / 0=auto)

### Analysis Modes

The `AnalysisMode` parameter (passed from config) determines the decay chain:
- `"JpsiJpsiPhi"` (default): Two J/ψ + φ
- `"JpsiJpsiUps"`: Two J/ψ + Υ
- `"JpsiUpsPhi"`: J/ψ + Υ + φ

Each mode sets which mass windows and resonance candidates to use. The code uses `AnalysisChannel` enum internally.

### Key Utility Classes

- **VertexReProducer/OniaVtxReProducer**: Re-fit primary vertices excluding signal muons
- **Onia2MuMuPAT**: EDProducer for dimuon vertexing (separate from main analyzer)
- **TriggerBooking**: Handles HLT trigger path matching (stored in large header file)

### TTree Output Structure

The main TTree (`X_One_Tree`) stores:
- Event info: run/lumi/event numbers, trigger decisions
- Primary vertex: position, errors, χ², beamspot-corrected variants
- All muons: kinematics, ID flags, impact parameters, trigger matching
- Resonance candidates (J/ψ, Υ, φ): mass, cτ, vertex prob, momentum (with uncertainties)
- Kaon daughters of φ (`Phi_K_1_*`, `Phi_K_2_*`): per-track kinematics, PV diagnostics
- `Phi_K_*_RecoKaonTrackIdx`: indices linking φ daughters to the flat RecoKaonTrack block
- Single-object MC branches (`SingleJpsi_*`, `SingleUps_*`, `SinglePhi_*`): per-candidate kinematics, GEN-matching, daughter indices
- `SinglePhi_K*_RecoKaonTrackIdx`: indices linking single-φ daughters to the RecoKaonTrack block
- RecoKaonTrack block (`RecoKaonTrack_*`): flat per-track storage of (GEN-matched quality kaons) ∪ (φ-candidate daughters), with PV diagnostics, GEN-matching, and `usedInSinglePhi` flag
- 3-body primary vertex: combined fit results
- MC gen-level (`MC_GenPart_*`): flat storage of all relevant gen particles
- Legacy MC branches (`MC_X_*`, `MC_Dau_*`): kept for backward compatibility

### Global Tag Selection

`ConfFile_cfg.py` contains dictionaries mapping era → global tag:
- Data: `Run2022C-G` → `124X_dataRun3_*`, `Run2023C-D` → `130X_dataRun3_*`, `Run2024C-I` → `150X_dataRun3_v2`
- MC: `Run2022` → `130X_mcRun3_2022_realistic_v5`, `Run2023` → `130X_mcRun3_2023_realistic_v14`

## Code Conventions

### CMSSW-Specific Patterns
- Use `EDGetToken`/`ESGetToken` for event data product access
- `StringCutObjectSelector<T>` for runtime-configurable cuts
- `RefCountedKinematicParticle` for kinematic fit objects
- Use `TransientTrackBuilder` from EventSetup for vertex fitting

### Naming Conventions
- Resonance candidates: `Jpsi_1_*`, `Jpsi_2_*`, `Ups_*`, `Phi_*` (suffixes: mass, ctau, px/py/pz, pt/eta/phi)
- Muon indices: `Jpsi_1_mu_1_Idx`, `Jpsi_1_mu_2_Idx` (links muons to resonances)
- Kaon RecoKaonTrack indices: `Phi_K_1_RecoKaonTrackIdx`, `SinglePhi_K1_RecoKaonTrackIdx` (links into `RecoKaonTrack_*` block)
- Single-object candidates: `SingleJpsi_*`, `SingleUps_*`, `SinglePhi_*` (all-store MC efficiency branches)
- Uncertainties: `*_pxErr`, `*_pyErr`, `*_pzErr`, `*_ptErr`
- RecoKaonTrack block: `RecoKaonTrack_*` (flat per-track storage, used for phi/kaon efficiency)
- Counter branches: `nRecoKaonTrack`, `nSingleJpsiCand`, `nSingleUpsCand`, `nSinglePhiCand`

### Adding New Branches

1. Declare branch pointer in `MultiLepPAT.h` (e.g., `vector<float> *newBranch`)
2. Initialize in `beginJob()`: `newBranch = new vector<float>()`
3. Fill in appropriate analysis step (e.g., `combineCandidates()`)
4. Register TTree branch in `beginJob()`
5. Clear in `clearEventData()`

### Modifying Selection Cuts

Edit `test/ConfFile_cfg.py`:
- `MuonSelection`/`TrackSelection` for object-level cuts
- Mass window parameters (GeV units)
- Vertex probability thresholds
- `AnalysisMode` to switch decay channel

No recompilation needed for these changes.

## Known Issues and Legacy Code

The following are kept for backward compatibility but may be deprecated:
- `muD0`, `muD0E`, `muDz` - Raw impact parameters (prefer vertex-corrected)
- `muChi2`, `muGlChi2`, `muNDF`, `muGlNDF` - Old track quality variables
- `muMVAMuonID` - BDT muon ID (requires `data/TMVAClassification_BDT.class.C`)
- `MC_X_*`, `MC_Dau_*` - Legacy X(3872) decay chain branches (superseded by `MC_GenPart_*`)
