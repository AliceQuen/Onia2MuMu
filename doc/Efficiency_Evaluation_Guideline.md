# Efficiency and Acceptance Evaluation Guideline

Maps the factorized efficiency scheme in `Efficiency_scheme.md` to the current
ntuple branch structure and configuration system.

---

## 1. Required MC runs

Two configurations are needed, both from `ConfFile_cfg.py`:

### Run A — Object-level (singles only)

```bash
cmsRun ConfFile_cfg.py \
    analysisMode=JpsiJpsiPhi \
    inputFiles=file:myMC.root \
    outputFile=eff_singles.root \
    runOnMC=True era=Run2022 \
    keepAllSingleObjectCandsInMC=True \
    skipCompositeCandBuildingWhenKeepingSingles=True
```

Provides: `SingleJpsi_*`, `SinglePhi_*`, `RecoKaonTrack_*`, `MC_GenPart_*`.
Used for: acceptance, muonRECO, kaonRECO, muonID, kaonID, dimuon, dikaon steps.

Use **loose object-level cuts** so that denominator populations are not cut away
before efficiency evaluation. The tight defaults (`JpsiCandPtMin=4.0`,
`PhiCandPtMin=2.0`) will bias your denominators. Either pass loose values on the
command line (not available via VarParsing for these parameters) or maintain a
loose-cut config variant with `JpsiCandPtMin=0.0`, `PhiCandPtMin=0.0`,
`JpsiMassMin=1.0`, `JpsiMassMax=4.0`, `PhiMassMin=0.8`, `PhiMassMax=1.2`.

### Run B — Full chain (singles + composite)

```bash
cmsRun ConfFile_cfg.py \
    analysisMode=JpsiJpsiPhi \
    inputFiles=file:myMC.root \
    outputFile=eff_fullchain.root \
    runOnMC=True era=Run2022 \
    keepAllSingleObjectCandsInMC=True \
    skipCompositeCandBuildingWhenKeepingSingles=False
```

Provides: all of Run A plus `Jpsi_1_*`, `Jpsi_2_*`, `Phi_*`, `DiOnia_*`, `Pri_*`,
`Phi_K_*_RecoKaonTrackIdx`, HLT/trigger branches.
Used for: HLT, four-muon vertexing, triOnia vertexing steps.

---

## 2. Branch-to-step mapping

### 2.1 Acceptance (GEN-level, no RECO needed)

| Quantity | Branch / logic |
|----------|---------------|
| GEN J/ψ in acceptance | `MC_GenPart_pdgId == 443` and `abs(MC_GenPart_pdgId[mother]) != 443` (direct J/ψ, not feed-down). Apply generator-level pT/rapidity cuts. |
| GEN φ in acceptance | `MC_GenPart_pdgId == 333` (direct φ). Same logic. |
| GEN μ from J/ψ | `MC_GenPart_pdgId == 13` or `-13`, mother is the GEN J/ψ. Apply fiducial μ cuts (pT, η). |
| GEN K from φ | `MC_GenPart_pdgId == 321`, mother is the GEN φ. Apply fiducial K cuts. |

Numerator: events where both GEN daughters pass fiducial cuts.
Denominator: events with GEN J/ψ (or φ) in the kinematic acceptance.

The `MC_GenPart_*` branches store all relevant GEN particles in flat vectors.
Use the mother-daughter linking via `MC_GenPart_motherPdgId` and the stored
particle ordering (daughter follows mother in the GEN record). The
`handleToNtupleIndex_` / `ntupleToHandleIndex_` maps (populated in
`processMCGenInfo()` but not persisted to the TTree) are not needed here —
work directly with the flat `MC_GenPart_*` vectors.

### 2.2 muonRECO — both daughter muons reconstructed and matched

| Quantity | Branch |
|----------|--------|
| RECO muons matched to GEN J/ψ daughters | For each `SingleJpsi` candidate: `SingleJpsi_mu1_genMatchIdx >= 0` AND `SingleJpsi_mu2_genMatchIdx >= 0` |
| Denominator | `SingleJpsi` candidates where both GEN daughter muons are in fiducial acceptance (cross-reference `MC_GenPart_*` via the stored daughter PDG IDs) |

The `SingleJpsi_*` branches are filled **before** the final `JpsiCandPtMin` /
`JpsiCandEtaMax` cuts in `pairMuons()`, so they capture all dimuon pairs that
survive the mass window and vertex fit. The `SingleJpsi_mu1_Idx` /
`SingleJpsi_mu2_Idx` give indices into the `mu*` branches for the daughter
muons; `SingleJpsi_mu1_genMatchIdx` / `SingleJpsi_mu2_genMatchIdx` are indices
into `MC_GenPart_*`.

The RECO-level matching for muons uses the configurable `RecoGenMuonMatchChi2Max`
(default 25.0). Muons without a valid GEN match have `genMatchIdx = -1`.

### 2.3 kaonRECO — both daughter kaons reconstructed and matched

| Quantity | Branch |
|----------|--------|
| RECO kaons matched to GEN φ daughters | For each `SinglePhi` candidate: `SinglePhi_K1_genMatchIdx >= 0` AND `SinglePhi_K2_genMatchIdx >= 0` |
| Denominator | `SinglePhi` candidates where both GEN daughter kaons are in fiducial acceptance |

The `SinglePhi_K1_RecoKaonTrackIdx` / `SinglePhi_K2_RecoKaonTrackIdx` give indices
into `RecoKaonTrack_*` for the daughter kaons, which contains the GEN-match
information (`RecoKaonTrack_genMatchIdx`).

The RECO-level matching for kaons uses `RecoGenKaonMatchChi2Max` (default 25.0).

Alternatively, work entirely through the `RecoKaonTrack` block: iterate
`RecoKaonTrack_*` entries with `genMatchIdx >= 0`, group by event, and check
whether both GEN-matched kaons from a given φ are present.

### 2.4 muonID — both matched muons pass ID

| Quantity | Branch |
|----------|--------|
| Muon ID flags | `muIsGoodTightMuon`, `muIsPatTightMuon`, `muIsPatMediumMuon`, `muIsPatSoftMuon`, `muIsGlobalMuon` |
| Per-J/ψ daughter indices | `SingleJpsi_mu1_Idx`, `SingleJpsi_mu2_Idx` — use these to look up the muon block |

The specific ID working point is analysis-dependent. A typical tight-muon
definition checks `muIsGoodTightMuon && muIsGlobalMuon`. Apply the ID selection
to both daughter muons of the J/ψ candidate.

Denominator: `SingleJpsi` candidates passing muonRECO (both muons GEN-matched).

### 2.5 kaonID — both matched kaons pass track quality

| Quantity | Branch |
|----------|--------|
| Kaon track quality | `RecoKaonTrack_passDzPV`, `RecoKaonTrack_passDxyPV`, `RecoKaonTrack_passTrackPV`, `RecoKaonTrack_fromPV` |
| Per-φ daughter RecoKaonTrack indices | `SinglePhi_K1_RecoKaonTrackIdx`, `SinglePhi_K2_RecoKaonTrackIdx` |

The `RecoKaonTrack_*` block stores the same quality flags used in `pairTracks()`
for both GEN-matched tracks and φ-candidate daughters. Look up the daughter's
`RecoKaonTrack_*` entry via the `RecoKaonTrackIdx` and check the quality flags.

Denominator: `SinglePhi` candidates passing kaonRECO (both kaons GEN-matched).

### 2.6 dimuon — valid J/ψ dimuon candidate

| Quantity | Branch |
|----------|--------|
| Vertex fit valid | `SingleJpsi_fitValid > 0` |
| Vertex fit passes prob cut | `SingleJpsi_fitPass > 0` (respects `JpsiDecayVtxProbCut`) |
| Mass window | Already applied in `pairMuons()`; `SingleJpsi_mass` is the post-fit mass |

A "valid dimuon candidate" is one where `fitValid && fitPass && massDiff` is
within acceptable range. The `SingleJpsi_massDiff` branch stores the difference
between pre-fit and post-fit mass.

Denominator: `SingleJpsi` candidates passing muonID.

### 2.7 dikaon — valid φ → K⁺K⁻ candidate

| Quantity | Branch |
|----------|--------|
| Vertex fit valid | `SinglePhi_fitValid > 0` |
| Vertex fit passes prob cut | `SinglePhi_fitPass > 0` (respects `PhiDecayVtxProbCut`) |
| Mass window | Already applied in `pairTracks()`; `SinglePhi_mass` is the post-fit mass |
| Common PV association | `SinglePhi_commonAssocPVPass` |
| Track PV compatibility | `SinglePhi_trackPVPass` |

Denominator: `SinglePhi` candidates passing kaonID.

### 2.8 HLT — trigger OR of J/ψ triggers

| Quantity | Branch |
|----------|--------|
| Event-level trigger decisions | `trigRes` (HLT results), `trigNames` (HLT path names) |
| J/ψ trigger matching | `MatchJpsiTriggerNames` (which J/ψ trigger paths matched), `muIsJpsiTrigMatch` (per-muon trigger match) |
| J/ψ filter matching | `muIsJpsiFilterMatch` |

The HLT step requires at least one of the configured J/ψ trigger paths to fire
AND have trigger-object matching to the candidate muons.

Denominator: events with at least one valid `Jpsi_1` candidate (from the full
chain ntuple, Run B).

### 2.9 four-muon vertexing — valid J/ψ J/ψ vertex

| Quantity | Branch |
|----------|--------|
| DiOnia fit valid | `DiOnia_fitValid > 0` |
| DiOnia fit passes prob cut | `DiOnia_fitPass > 0` (respects `DiOniaVtxProbCut`) |
| Common reco-vertex check | `DiOnia_commonRecVtxPass > 0` |

This step combines the two J/ψ candidates into a common vertex. It requires both
`Jpsi_1` and `Jpsi_2` branches from the full chain ntuple (Run B).

Denominator: events passing HLT with at least one valid dimuon pair for each of
the two J/ψ slots.

### 2.10 triOnia — final three-meson vertex + event-level cuts

| Quantity | Branch |
|----------|--------|
| Primary vertex fit valid | `Pri_fitValid > 0` |
| Primary vertex fit passes prob cut | `Pri_fitPass > 0` (respects `PriVtxProbCut`) |
| Common PV association | `Pri_assocPVPass > 0` |
| Track PV compatibility | `Pri_trackPVPass > 0` |
| Final mass check | Internal (`CheckFinalMass` flag) |

The code produces four parallel quality flags for the primary vertex. The
default endpoint for corrected-yield studies is `Pri_assocPVPass / four_muon_vtx`
as noted in `Efficiency_scheme.md`.

Denominator: events passing four-muon vertexing, split by φ pT bins.

---

## 3. Configuration considerations

### 3.1 Cuts that must be loose for denominator studies

The following cuts are applied **before** single-object branches are filled and
cannot be undone downstream. They must be set loosely for efficiency denominator
evaluation:

| Parameter | Tight default | Recommended loose |
|-----------|--------------|-------------------|
| `JpsiMassMin` / `JpsiMassMax` | 2.8 / 3.3 | 1.0 / 4.0 |
| `PhiMassMin` / `PhiMassMax` | 0.97 / 1.07 | 0.8 / 1.2 |
| `OniaDecayVtxProbCut` | 5e-4 | 0.0 (or keep) |
| `JpsiDecayVtxProbCut` | 5e-4 | 0.0 (or keep) |
| `PhiDecayVtxProbCut` | 5e-4 | 0.0 (or keep) |
| `MuonSelection` | `pt > 2.5 && abs(eta) < 2.4` | Match your fiducial definition |
| `TrackSelection` | `pt > 1.0 && abs(eta) < 2.5 && numberOfHits > 4` | Match your fiducial definition |

These are hardcoded in `ConfFile_cfg.py` and not exposed via `VarParsing` (the
`VarParsing` system only exposes the MC control switches, not the physics cut
values). Maintain a separate loose-cut config file for efficiency studies, or
override the values directly before `cmsRun`.

### 3.2 Cuts you can tighten in post-processing

| Parameter | Why post-processable |
|-----------|---------------------|
| `JpsiCandPtMin`, `JpsiCandEtaMax` | These are applied in `pairMuons()` as pre-cuts on the composite candidate. The `SingleJpsi_pt` branch records the value; you can apply any pT cut in analysis. |
| `PhiCandPtMin`, `PhiCandEtaMax` | Same reasoning — `SinglePhi_pt` records the value. |
| `MinTrackFromPV` | Recorded in `RecoKaonTrack_fromPV` and `SinglePhi_K*_fromPV`. |
| Track `normalizedChi2`, `highPurity` | Recorded (indirectly) via `RecoKaonTrack_passTrackPV`. |

### 3.3 MC control flags

| Flag | Effect on efficiency evaluation |
|------|--------------------------------|
| `keepAllSingleObjectCandsInMC=True` | **Required** — enables `SingleJpsi_*`, `SinglePhi_*`, `RecoKaonTrack_*` |
| `skipCompositeCandBuildingWhenKeepingSingles=True` | Optional — skips `combineCandidates()` for a smaller, faster ntuple when only per-object steps are needed |
| `requireAcceptedCandidatesForMonteCarloTree=False` | Recommended for efficiency — keeps every event even when no candidate passes |

### 3.4 Vertex fit toggles

| Toggle | Default | Notes |
|--------|---------|-------|
| `DoJpsiDecayVtxFit` | True | Needed for `SingleJpsi_fitValid` |
| `DoPhiDecayVtxFit` | True | Needed for `SinglePhi_fitValid` |
| `DoDiOniaVtxFit` | True | Needed for `DiOnia_*` branches |
| `DoPriVtxFit` | True | Needed for `Pri_*` branches |

---

## 4. Suggested analysis workflow

### Step 1: Produce efficiency ntuples

Run both Run A (singles-only) and Run B (full chain) with loose physics cuts on
the full MC sample(s). Use `maxEvents=-1` to process all events.

### Step 2: Build per-object efficiency maps

For each J/ψ and φ kinematic bin `(pT, |y|)`:

1. **Acceptance**: From `MC_GenPart_*`, count GEN J/ψ (or φ) in the bin.
   Numerator: daughters in fiducial region. Denominator: all GEN J/ψ (or φ).
   Output: `A_Jpsi(pT, |y|)`, `A_phi(pT, |y|)`.

2. **muonRECO / kaonRECO**: From `SingleJpsi_*` / `SinglePhi_*` (Run A).
   Numerator: candidates with both daughters `genMatchIdx >= 0`.
   Denominator: candidates with both GEN daughters in fiducial acceptance.
   Output: `eps_muReco(pT, |y|)`, `eps_kReco(pT, |y|)`.

3. **muonID / kaonID**: From `SingleJpsi_*` + `mu*` / `SinglePhi_*` + `RecoKaonTrack_*`.
   Numerator: RECO-passing candidates where both daughters pass ID.
   Denominator: RECO-passing candidates.
   Output: `eps_muID(pT, |y|)`, `eps_kID(pT, |y|)`.

4. **dimuon / dikaon**: From `SingleJpsi_*` / `SinglePhi_*`.
   Numerator: candidates with `fitValid && fitPass`.
   Denominator: ID-passing candidates.
   Output: `eps_mumu(pT, |y|)`, `eps_KK(pT, |y|)`.

### Step 3: Build event-level efficiency maps

From the full chain ntuple (Run B):

1. **HLT**: Count events where trigger fired + trigger-object match exists.
   Denominator: events with ≥1 valid dimuon candidate.

2. **four-muon vertexing** (`eps_4muVtx`): From `DiOnia_*`.
   Numerator: `DiOnia_fitValid && DiOnia_fitPass`.
   Denominator: HLT-passing events.

3. **triOnia** (`eps_triOnia`): From `Pri_*`, split by φ pT bins.
   Numerator: `Pri_assocPVPass` (or chosen endpoint).
   Denominator: four-muon-vertexing-passing events.

### Step 4: Apply correction

The per-event total efficiency weight is the product of all per-object and
event-level factors evaluated at the event's kinematics, times the factorized
acceptance. See `Efficiency_scheme.md` Section "Definitions" for the product
formula.

---

## 5. Cross-checks

### 5.1 Factorized vs. cumulative comparison

Compare the factorized efficiency product against a direct "all steps pass"
count on the full chain ntuple. Significant disagreement indicates correlations
between the per-object steps that the factorized approach misses.

### 5.2 RecoKaonTrack coverage

Verify that every φ-candidate daughter has a `RecoKaonTrack_*` entry:
`SinglePhi_K*_RecoKaonTrackIdx >= 0` for all `SinglePhi` candidates. A `-1`
sentinel indicates a kaon track missing from the RecoKaonTrack block (should
not happen when `KeepAllSingleObjectCandsInMC=True`).

### 5.3 Bin-edge sensitivity

Vary the fine/coarse bin thresholds (`N_min_fine`, `N_min_coarse`) and confirm
the corrected yield is stable within MC statistical uncertainties.

### 5.4 GEN-RECO daughter matching consistency

For `SinglePhi` candidates, compare the direct GEN-match info
(`SinglePhi_K1_genMatchIdx`) with the RecoKaonTrack-level match
(`RecoKaonTrack_genMatchIdx[SinglePhi_K1_RecoKaonTrackIdx]`). They must agree.
