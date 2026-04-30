# CMS Onia2MuMu Analysis - MultiLepPAT

## Project Overview

`MultiLepPAT` is an EDAnalyzer module for the CMS experiment designed to search for new X particles through the reconstruction of X(3872) decay processes. This analysis focuses on the four-muon plus two-pion final state topology, using kinematic vertex fitting and mass constraint techniques to precisely measure particle properties.

## Physics Goals

The physics goals of this project include:

1. **Search for associated production of X(3872) and J/ψ particles**
2. **Search for X particle decays to J/ψ and ψ(2S)**
3. **Search for X particle decays to two ψ(2S) mesons**

### Decay Topology

```
                    X [candidate resonance]
                       │
           ┌───────────┴───────────┐
           │                       │
      ψ(2S) / X(3872)             J/ψ / ψ(2S)
    [mu1+mu2+pi1+pi2]            [mu3+mu4]
           │
     ┌─────┴─────┐
     │           │
   J/ψ₁        π⁺π⁻
 [mu1+mu2]   [pi1+pi2]
```

Final state: **4 muons + 2 pions**

The X particle is an unknown resonance, potentially X(6900) or X(7100).

## Three Mass Hypotheses

To comprehensively explore potential new resonances, the program employs three mass constraint hypotheses:

| Hypothesis Number | Mass Constraint Configuration                                                          | Physical Interpretation                             |
| ----------------- | -------------------------------------------------------------------------------------- | --------------------------------------------------- |
| **Hypothesis 1**  | mumupipi constrained to ψ(2S) (3.6861 GeV)Other mumu constrained to J/ψ (3.0969 GeV)   | Search for X → ψ(2S) J/ψ process                    |
| **Hypothesis 2**  | mumupipi constrained to X(3872) (3.872 GeV)Other mumu constrained to J/ψ (3.0969 GeV)  | Search for associated production of X(3872) and J/ψ |
| **Hypothesis 3**  | mumupipi constrained to ψ(2S) (3.6861 GeV)Other mumu constrained to ψ(2S) (3.6861 GeV) | Search for X → ψ(2S) ψ(2S) process                  |

## Core Functional Modules

```
MultiLepPAT Main Module
├── Event Preprocessing Layer
│   ├── HLT Trigger Matching
│   ├── Primary Vertex Acquisition
│   └── Event-level Pre-selection
├── Particle Selection Layer
│   ├── Track Pre-selection
│   ├── Muon Pre-selection
│   └── Trigger Filter Matching
├── Particle Combination Layer
│   ├── μ⁺μ⁻ Pair Combination
│   ├── π⁺π⁻ Pair Combination
│   ├── J/ψ - ππ Combination
│   └── J/ψ - J/ψ Combination
├── Vertex Fitting Layer
│   ├── Ordinary Vertex Fitting
│   └── Mass-constrained Vertex Fitting
├── Physics Calculation Layer
│   ├── Four-momentum Calculation
│   ├── Vertex Probability Calculation
│   ├── δR Correlation Calculation
│   └── Impact Parameter Calculation
└── Data Storage Layer
    ├── TTree Branch Management
    ├── Candidate Data Filling
    └── Variable Reset
```

## Particle Selection Criteria

### Track Pre-selection

| Criterion               | Threshold  | Physical Meaning                |
| ----------------------- | ---------- | ------------------------------- |
| High Purity Flag        | highPurity | Track quality indicator         |
| Transverse Momentum     | > 0.5 GeV  | Minimum momentum requirement    |
| Pseudorapidity          | < 2.4      | Detector acceptance range       |
| Relative Momentum Error | < 0.1      | Momentum resolution requirement |
| Valid Hits              | ≥ 10       | Track detection completeness    |
| Fit Quality             | < 0.18     | Track fitting quality           |

### Muon Pre-selection

| Criterion                  | Requirement                 | Description                                          |
| -------------------------- | --------------------------- | ---------------------------------------------------- |
| Soft Muon ID               | `isSoftMuon()`              | Track-based soft muon identification                 |
| Piecewise p\_T-η Threshold | See below                   | Momentum requirements for different rapidity regions |
| Charge Balance             | 2 positive μ + 2 negative μ | Necessary condition for two J/ψ                      |

Piecewise p\_T(η) function:

```
p_T^cut(η) =
  3.5 GeV,           when |η| < 1.2
  5.47 - 1.89*|η|,  when 1.2 < |η| < 2.1
  1.5 GeV,           when 2.1 < |η| < 2.4
```

## Program Execution Flow

### Lifecycle Overview

```
Program Start
   │
   ▼
Constructor
├── Read configuration parameters
├── Initialize EDM tokens
└── Trigger name sorting optimization
   │
   ▼
beginJob()
├── Create TTree
└── Create all branches (~130 total)
   │
   ▼
Event Loop (per event)
├── Data acquisition & pre-selection
├── HLT trigger matching
├── Primary vertex acquisition
├── Track pre-selection & charge grouping
├── Muon pre-selection & charge grouping
├── J/ψ candidate pre-caching
├── Nested loop combination & fitting
├── Three mass hypothesis constrained fitting
├── Physics quantity calculation
└── TTree Fill()
   │
   ▼
endJob()
   │
   ▼
Program End
```

### analyze() Function Selection Hierarchy

```
Start
  ↓
[Event-level] Muon count >= 4?
  ↓ NO → Return
  ↓ YES
[Event-level] HLT trigger matched?
  ↓ NO → Return
  ↓ YES
[Event-level] Valid primary vertex exists?
  ↓ NO → Return
  ↓ YES
[Particle-level] Track pre-selection count >= 2?
  ↓ NO → Return
  ↓ YES
[Particle-level] Muon pre-selection count >= 4?
  ↓ NO → Return
  ↓ YES
[Combination-level] μ₁+μ₂ charge sum = 0?
  ↓ NO → Continue
  ↓ YES
[Combination-level] J/ψ₁ invariant mass in window?
  ↓ NO → Continue
  ↓ YES
[Fitting-level] J/ψ₁ vertex fit successful?
  ↓ NO → Continue
  ↓ YES
[Fitting-level] J/ψ₁ vertex probability > 0.01?
  ↓ NO → Continue
  ↓ YES
[Combination-level] π₁+π₂ combination selection
  ↓
[Fitting-level] ψ(2S) 4-track vertex fit
  ↓
[Fitting-level] ψ(2S) vertex probability > 0.005
  ↓
[Fitting-level] ψ(2S) p_T > 4 GeV
  ↓
[Constraint-level] Three mass hypothesis constrained fits
  ↓
[Calculation-level] All physics quantity calculation
  ↓
[Storage-level] TTree filling
End
```

## Physics Quantity Storage Structure

### TTree Branch Hierarchy

**TTree Name**: `X_data`

**Total Branches**: \~130

```
X_data TTree
├── Event Information Layer (4)
│   ├── runNum, evtNum, lumiNum, nGoodPrimVtx
│
├── X Candidate Layer (24)
│   ├── Hypothesis 1: X_PJ_mass, X_PJ_VtxProb, X_PJ_massErr, X_PJ_pt, X_PJ_px, ...
│   ├── Hypothesis 2: X_XJ_mass, X_XJ_VtxProb, X_XJ_massErr, X_XJ_pt, X_XJ_px, ...
│   └── Hypothesis 3: X_PP_mass, X_PP_VtxProb, X_PP_massErr, X_PP_pt, X_PP_px, ...
│
├── ψ(2S) Candidate Layer (9)
│   └── Psi2S_mass, Psi2S_VtxProb, Psi2S_massErr, Psi2S_pt, Psi2S_absEta, ...
│
├── J/ψ Candidate Layer (16)
│   ├── Jpsi1_mass, Jpsi1_VtxProb, Jpsi1_massErr, Jpsi1_pt, ...
│   └── Jpsi2_mass, Jpsi2_VtxProb, Jpsi2_massErr, Jpsi2_pt, ...
│
├── Muon Information Layer (68)
│   ├── μ₁ (16 variables): pt, eta, charge, trackIso, d0BS, d0PV, dzPV, ...
│   ├── μ₂ (16 variables)
│   ├── μ₃ (16 variables)
│   ├── μ₄ (16 variables)
│   └── Muon ID Counts: nLooseMuons, nTightMuons, nSoftMuons, nMediumMuons
│
├── Pion Information Layer (10)
│   ├── π₁: pi1_pt, pi1_px, pi1_py, pi1_pz, pi1_absEta
│   └── π₂: pi2_pt, pi2_px, pi2_py, pi2_pz, pi2_absEta
│
└── δR Correlation Layer (14)
    ├── dR_mu1_mu2, dR_mu3_mu4, dR_pi1_pi2
    ├── dR_Psi2S_Jpsi1, dR_Psi2S_Jpsi2
    └── dR_Psi2S_pi1, dR_Psi2S_pi2, ...
```

### Invalid Value Handling

All physics quantities are set to sentinel values when calculation fails or values are invalid:

| Data Type | Invalid Value | Example                      |
| --------- | ------------- | ---------------------------- |
| Float     | `-999.0`      | `mu1_pt = -999.0`            |
| Integer   | `0`           | `nGoodPrimVtx = 0`           |
| Boolean   | `false`       | `mu1_hasFilterMatch = false` |

## Code Structure

### Main Files

```
Onia2MuMu/
├── interface/
│   └── MultiLepPAT.h         # Header file, class definition & JpsiCandidate struct
├── src/
│   └── MultiLepPAT.cc        # Source file, implementation code
└── test/
    └── runMultiLepPAT_miniAOD.py   # CMSSW configuration file
```

### Key Data Structure: JpsiCandidate

Defined in `interface/MultiLepPAT.h`, used to cache J/ψ candidate fitting results and avoid redundant calculations in nested loops:

```cpp
struct JpsiCandidate {
    // Particle iterators
    edm::View<pat::Muon>::const_iterator muPlus;
    edm::View<pat::Muon>::const_iterator muMinus;
    
    // Ordinary vertex fit results
    double mass, vtxProb, massErr;
    ROOT::Math::PxPyPzMVector p4;
    RefCountedKinematicParticle kinematicParticle;
    RefCountedKinematicVertex vertex;
    
    // Mass-constrained fit results
    bool hasConstraintFit = false;
    double constraintMass, constraintVtxProb, constraintMassErr;
    RefCountedKinematicParticle constraintParticle;
    RefCountedKinematicVertex constraintVertex;
    
    // Candidate type flags
    bool isJpsiCandidate = false;
    bool isPsi2SCandidate = false;
    
    // Additional information
    reco::TransientTrack muonTT1, muonTT2;
    bool filterMatchPlus, filterMatchMinus;
};
```

## Core Optimization Strategies

1. **Trigger Name Pre-sorting**: Sorted by length ascending, shorter patterns matched first for early loop exit on average
2. **Charge Pre-grouping**: Tracks and muons grouped by charge before loops, reducing invalid combinations
3. **J/ψ Candidate Pre-caching**: All valid J/ψ candidates pre-computed outside nested loops to avoid redundant fitting
4. **Cheap Cut First**: Low-computational-cost selections (e.g., quick mass calculation) placed innermost in loops for early invalid combination rejection

## Compilation and Running

### Environment Requirements

- CMSSW release (adjust according to actual environment)
- ROOT 6+
- gcc 7+

### Compilation Steps

```bash
# Enter CMSSW working directory
cd $CMSSW_BASE/src

# Clone code (if needed)
git clone <repository-url>

# Compile
scram b
```

### Running

```bash
# Run with configuration file
cmsRun test/runMultiLepPAT_miniAOD.py
```

## Document Version

- Document Version: v1.0
- Creation Date: 2026.04.30
- Applicable Code Version: MultiLepPAT (2024 Refactored)

