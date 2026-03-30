# TPS-Onia2MuMu: MultiLepPAT EDAnalyzer

Reconstructs quarkonia+quarkonia+meson final states from CMS MINIAOD data.

## Supported Analysis Modes

| `AnalysisMode`    | Decay Chain                              | Daughters       |
|-------------------|------------------------------------------|-----------------|
| `JpsiJpsiPhi`     | J/ψ(→μμ) + J/ψ(→μμ) + φ(→KK)           | 4μ + 2K         |
| `JpsiJpsiUps`     | J/ψ(→μμ) + J/ψ(→μμ) + Υ(→μμ)           | 6μ              |
| `JpsiUpsPhi`      | J/ψ(→μμ) + Υ(→μμ)  + φ(→KK)            | 4μ + 2K         |

## Configuration (ConfFile_cfg.py)

All selection cuts are externalized and can be changed **without recompilation**:

### Runtime String Cuts (StringCutObjectSelector)

```python
MuonSelection  = cms.untracked.string("pt > 2.5 && abs(eta) < 2.4")
TrackSelection = cms.untracked.string("pt > 2.0 && abs(eta) < 2.5 && numberOfHits > 4")
```

These accept any valid expression on `pat::Muon` or `pat::PackedCandidate` members.

### Mass Windows (GeV)

| Parameter     | Default | Description                     |
|---------------|---------|---------------------------------|
| `JpsiMassMin` | 1.0    | J/ψ dimuon mass lower bound     |
| `JpsiMassMax` | 4.0    | J/ψ dimuon mass upper bound     |
| `UpsMassMin`  | 8.0    | Υ dimuon mass lower bound       |
| `UpsMassMax`  | 12.0   | Υ dimuon mass upper bound       |
| `PhiMassMin`  | 0.8    | φ(→KK) mass lower bound         |
| `PhiMassMax`  | 1.2    | φ(→KK) mass upper bound         |

### Primary Vertex Cuts

| Parameter   | Default | Description                      |
|-------------|---------|----------------------------------|
| `PVNdofMin` | 5       | Minimum ndof for good PV         |
| `PVMaxAbsZ` | 24.0    | Maximum |z| of PV (cm)           |
| `PVMaxRho`  | 2.0     | Maximum transverse distance (cm) |

### Vertex Probability Cuts

| Parameter             | Default | Description                          |
|-----------------------|---------|--------------------------------------|
| `OniaDecayVtxProbCut` | 0.001   | Minimum vertex prob for dimuon fits  |
| `PriVtxProbCut`       | 0.0     | Minimum vertex prob for primary vtx  |

## New Branches (vs. previous version)

### Momentum Uncertainties
For each resonance (`Jpsi_1`, `Jpsi_2`, `Phi`, `Pri`):
- `*_pxErr`, `*_pyErr`, `*_pzErr`: Component errors from `KinematicParametersError`
- `*_ptErr`: Propagated transverse momentum error

### MC Gen-level (flat storage)
- `MC_GenPart_pdgId`, `MC_GenPart_status`, `MC_GenPart_motherPdgId`
- `MC_GenPart_px/py/pz/mass/pt/eta/phi`

Stores all relevant gen particles (μ=13, K=321, J/ψ=443, Υ=553, φ=333) in flat vectors for flexible offline matching.

## Legacy Variables

The following branches are retained for backward compatibility but may be deprecated:

- `muD0`, `muD0E`, `muDz` — Track impact parameters (now prefer vertex-corrected)
- `muChi2`, `muGlChi2`, `muNDF`, `muGlNDF` — Old-style muon track quality
- `mufHits`, `muFirstBarrel`, `muFirstEndCap` — Tracker hit pattern details
- `muType`, `muQual`, `muTrack` — Encoded muon type/quality bitmasks
- `muIsGoodLooseMuon`, `muIsGoodLooseMuonNew` — Deprecated soft muon selectors
- `muMVAMuonID` — BDT-based muon ID (uses `TMVAClassification_BDT.class.C`)
- `MC_X_*`, `MC_Dau_*`, `MC_Granddau_*`, `MC_Grandgranddau_*` — Legacy X(3872) decay-chain branches; superseded by `MC_GenPart_*`

## Fresh CMSSW Setup

Create a new CMSSW release area, clone this package directly under `src/HeavyFlavorAnalysis/`, and build from there:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW CMSSW_15_0_15
cd CMSSW_15_0_15/src
eval "$(scramv1 runtime -sh)"

mkdir -p HeavyFlavorAnalysis
git clone --recursive git@github.com:Eric100911/TPS-Onia2MuMu.git HeavyFlavorAnalysis/TPS-Onia2MuMu

scram b -j 4 HeavyFlavorAnalysis/TPS-Onia2MuMu
```

If you cloned without `--recursive`, initialize the CRAB helper submodule afterwards:

```bash
cd $CMSSW_BASE/src/HeavyFlavorAnalysis/TPS-Onia2MuMu
git submodule update --init --recursive
```

## Repository Structure

The package is expected to live at `src/HeavyFlavorAnalysis/TPS-Onia2MuMu` inside the CMSSW release area:

```text
CMSSW_15_0_15/
└── src/
    └── HeavyFlavorAnalysis/
        └── TPS-Onia2MuMu/
            ├── interface/
            │   ├── MultiLepPAT.h              # Main analyzer header
            │   └── VertexReProducer.h         # Vertex re-fitting utility
            ├── src/
            │   ├── MultiLepPAT.cc             # Main analyzer implementation
            │   └── VertexReProducer.cc        # Vertex re-fitting implementation
            ├── data/
            │   └── TMVAClassification_BDT.class.C
            ├── python/
            │   └── onia2MuMuPAT_cfi.py        # Onia2MuMuPAT producer config
            ├── test/
            │   ├── ConfFile_cfg.py            # Default config
            │   ├── runMultiLepPAT_MCRun3_miniAOD_Run2022.py
            │   └── crabData/                  # Git submodule -> CRAB-Tool
            └── BuildFile.xml
```

## Building

From the CMSSW release area:

```bash
cd $CMSSW_BASE/src
scram b -j 4 HeavyFlavorAnalysis/TPS-Onia2MuMu
```

To rebuild the full release instead, use `scram b -j 4`.

## Running

From `$CMSSW_BASE/src`:

```bash
# Data
cmsRun HeavyFlavorAnalysis/TPS-Onia2MuMu/test/ConfFile_cfg.py \
    inputFiles=file:myData.root outputFile=output.root

# MC
cmsRun HeavyFlavorAnalysis/TPS-Onia2MuMu/test/runMultiLepPAT_MCRun3_miniAOD_Run2022.py \
    inputFiles=file:myMC.root outputFile=mc_output.root

# Switch analysis mode (e.g. J/psi+Upsilon+phi)
cmsRun HeavyFlavorAnalysis/TPS-Onia2MuMu/test/ConfFile_cfg.py \
    AnalysisMode=JpsiUpsPhi
```

## Additional Documentation

The previous README contained longer reference material that is now preserved separately:

- `doc/Physics_Overview_and_Method.md`: physics motivation, event-building strategy, phi reconstruction notes, and implementation walkthroughs preserved from the previous README
- `doc/Datasets_and_Operations.md`: preserved dataset inventories, DAS usage notes, and run-era operational references
- `doc/Historical_Project_Notes.md`: preserved dated project status, legacy branch notes, and contributor information
