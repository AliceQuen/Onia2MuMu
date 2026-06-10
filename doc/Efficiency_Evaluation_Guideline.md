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

### 2.1 Acceptance — $A_{J/\psi}$, $A_{\phi}$

The only step with an unconditional denominator. GEN-level, no RECO needed.

**$J/\psi$ acceptance** — $A_{J/\psi}(p_T, |y|)$:

| Role | Logic |
|------|-------|
| Denominator | GEN J/ψ in kinematic acceptance: `MC_GenPart_pdgId == 443` and `abs(MC_GenPart_pdgId[mother]) != 443` (direct J/ψ, not feed-down) |
| Numerator | Denominator J/ψ where both GEN daughter muons are in fiducial region: `MC_GenPart_pdgId == 13` or `-13` with mother = the GEN J/ψ, passing fiducial μ cuts (pT, η) |

**$\phi$ acceptance** — $A_{\phi}(p_T, |y|)$:

| Role | Logic |
|------|-------|
| Denominator | GEN φ in kinematic acceptance: `MC_GenPart_pdgId == 333` (direct φ) |
| Numerator | Denominator φ where both GEN daughter kaons are in fiducial region: `MC_GenPart_pdgId == 321` with mother = the GEN φ, passing fiducial K cuts |

The `MC_GenPart_*` branches store all relevant GEN particles in flat vectors.
Use mother-daughter linking via `MC_GenPart_motherPdgId` and the stored
particle ordering (daughter follows mother in the GEN record). The
`handleToNtupleIndex_` / `ntupleToHandleIndex_` maps (populated in
`processMCGenInfo()` but not persisted to the TTree) are not needed here —
work directly with the flat `MC_GenPart_*` vectors.

### 2.2 muonRECO — $\varepsilon_{\mu\mathrm{Reco}|J/\psi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed acceptance) | `SingleJpsi` candidates where both GEN daughter muons are in fiducial acceptance (cross-reference `MC_GenPart_*` via daughter PDG IDs) |
| Numerator (passes THIS step) | Denominator `SingleJpsi` where `SingleJpsi_mu1_genMatchIdx >= 0` AND `SingleJpsi_mu2_genMatchIdx >= 0` |

The `SingleJpsi_*` branches are filled **before** the final `JpsiCandPtMin` /
`JpsiCandEtaMax` cuts in `pairMuons()`, so they capture all dimuon pairs that
survive the mass window and vertex fit. `SingleJpsi_mu1_Idx` /
`SingleJpsi_mu2_Idx` give indices into the `mu*` branches for the daughter
muons; `SingleJpsi_mu1_genMatchIdx` / `SingleJpsi_mu2_genMatchIdx` are indices
into `MC_GenPart_*`.

The RECO-level matching for muons uses the configurable `RecoGenMuonMatchChi2Max`
(default 25.0). Muons without a valid GEN match have `genMatchIdx = -1`.

### 2.3 kaonRECO — $\varepsilon_{K\mathrm{Reco}|\phi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed acceptance) | `SinglePhi` candidates where both GEN daughter kaons are in fiducial acceptance |
| Numerator (passes THIS step) | Denominator `SinglePhi` where `SinglePhi_K1_genMatchIdx >= 0` AND `SinglePhi_K2_genMatchIdx >= 0` |

`SinglePhi_K1_RecoKaonTrackIdx` / `SinglePhi_K2_RecoKaonTrackIdx` give indices
into the `RecoKaonTrack_*` block (always populated from φ-candidate daughters).
The `RecoKaonTrack_genMatchIdx` field mirrors the inline `SinglePhi_K*_genMatchIdx`
— they must agree (see Section 5.4).

The RECO-level matching for kaons uses `RecoGenKaonMatchChi2Max` (default 25.0).
Kaons without a valid GEN match have `genMatchIdx = -1`.

Alternatively, work entirely through the `RecoKaonTrack` block: iterate
`RecoKaonTrack_*` entries with `genMatchIdx >= 0`, group by event, and check
whether both GEN-matched kaons from a given φ are present.

### 2.4 muonID — $\varepsilon_{\mu\mathrm{ID}|J/\psi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed muonRECO) | `SingleJpsi` candidates with both muons GEN-matched (muonRECO numerator) |
| Numerator (passes THIS step) | Denominator `SingleJpsi` where both daughter muons pass the chosen ID working point |

Use `SingleJpsi_mu1_Idx` / `SingleJpsi_mu2_Idx` to look up the muon block.
Available muon ID flags: `muIsGoodTightMuon`, `muIsPatTightMuon`, `muIsPatMediumMuon`,
`muIsPatSoftMuon`, `muIsGlobalMuon`.

The baseline ID working point is **soft muon**: `muIsPatSoftMuon`. The working
point is analysis-configurable and can be tightened in post-processing.

### 2.5 kaonID — $\varepsilon_{K\mathrm{ID}|\phi}$

Track-quality identification applied to the reconstructed (GEN-matched) kaon tracks.
This is the third tier of the kaon efficiency chain, after fiducial acceptance and
reconstruction.

The kaon efficiency factorizes into three sequential tiers:

| Tier | Step | Selection | Config / branches |
|------|------|-----------|-------------------|
| Fiducial acceptance | via $A_{\phi}$ (Section 2.1) | Phase-space only: $p_T$, $\eta$ cuts on GEN kaons | `TrackSelection` (`"pt > 1.0 && abs(eta) < 2.5"`) |
| kaonRECO | $\varepsilon_{K\mathrm{Reco}}^{\phi}$ (Section 2.3) | Both kaons GEN-matched | `SinglePhi_K*_genMatchIdx >= 0` |
| kaonID | $\varepsilon_{K\mathrm{ID}}^{\phi}$ (this section) | Track quality on reconstructed kaons | `TrackQuality` + `RequireRecoKaonTrackHighPurity` |

Note: pT/eta cuts are part of fiducial acceptance, NOT kaonID. The `TrackSelection`
string configures the phase-space denominator (see Section 3.1).

| Role | Branch / logic |
|------|----------------|
| Denominator (passed kaonRECO) | `SinglePhi` candidates with both kaons GEN-matched (kaonRECO numerator) |
| Numerator (passes THIS step) | Denominator `SinglePhi` where both daughter kaons pass track-quality ID |

Use `SinglePhi_K1_RecoKaonTrackIdx` / `SinglePhi_K2_RecoKaonTrackIdx` to look up
the `RecoKaonTrack_*` block (always populated from φ-candidate daughters).

The baseline track-quality criteria are defined by two config parameters (stored
in `X_Config_Tree`):

| Parameter | Default | Branch equivalents |
|-----------|---------|-------------------|
| `TrackQuality` | `"normalizedChi2 < 8 && numberOfValidHits > 4"` | `RecoKaonTrack_normalizedChi2`, `RecoKaonTrack_numberOfHits` |
| `RequireRecoKaonTrackHighPurity` | `True` | `RecoKaonTrack_isHighPurity == 1` |

Both criteria must be satisfied for each daughter kaon to count as passing kaonID.
The quality cuts are fully configurable at production time and can be varied in
offline analysis without re-ntuplizing, since the raw quality values are stored
per-track in the `RecoKaonTrack_*` block.

PV-compatibility flags (`RecoKaonTrack_passDzPV`, `RecoKaonTrack_passDxyPV`,
`RecoKaonTrack_passTrackPV`, `RecoKaonTrack_fromPV`) are reserved for the dikaon
and triOnia vertexing steps below.

### 2.6 dimuon — $\varepsilon_{\mu\mu|J/\psi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed muonID) | `SingleJpsi` candidates where both muons pass the chosen ID (muonID numerator) |
| Numerator (passes THIS step) | Denominator `SingleJpsi` where `SingleJpsi_fitValid > 0` AND `SingleJpsi_fitPass > 0` |

The baseline criterion for a valid dimuon candidate is that the kinematic vertex
fit converged and passes the vertex probability cut (`JpsiDecayVtxProbCut`). The
mass window is already applied in `pairMuons()`; `SingleJpsi_mass` is the
post-fit mass. Alternative or additional vertexing criteria may be applied
depending on the analysis working point.

### 2.7 dikaon — $\varepsilon_{KK|\phi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed kaonID) | `SinglePhi` candidates where both kaons pass track-quality ID (kaonID numerator) |
| Numerator (passes THIS step) | Denominator `SinglePhi` where `SinglePhi_fitValid > 0` AND `SinglePhi_fitPass > 0` |

The baseline criterion for a valid dikaon candidate is that the kinematic vertex
fit converged and passes the vertex probability cut (`PhiDecayVtxProbCut`). The
mass window is already applied in `pairTracks()`. Alternative or additional
vertexing criteria (e.g., `SinglePhi_trackPVPass`, `SinglePhi_commonAssocPVPass`)
may be applied depending on the analysis working point.

### 2.8 HLT — $\varepsilon_{\mathrm{HLT}}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed per-object steps) | Events with at least one valid dimuon candidate (Run B, full chain) |
| Numerator (passes THIS step) | Denominator events where at least one configured J/ψ trigger path fired AND the candidate muons have trigger-object matching |

**Event-level trigger information** (one entry per HLT path in the event):

| Branch | Content |
|--------|---------|
| `TrigRes` | Accept bit (0/1) for each HLT path |
| `TrigNames` | Name string for each HLT path |
| `MatchJpsiTriggerNames` | Deduplicated list of which configured J/ψ trigger patterns fired |
| `MatchUpsTriggerNames` | Same for Upsilon trigger patterns |
| `L1TrigRes` | L1 technical trigger bits |

**Per-muon trigger matching** (one entry per muon in the event):

| Branch | Content |
|--------|---------|
| `muIsJpsiTrigMatch` | Event-level boolean: 1 if any configured J/ψ trigger fired (same value for all muons) |
| `muIsUpsTrigMatch` | Same for Upsilon triggers |
| `muIsJpsiFilterMatch` | 1 if this muon's trigger objects match at least one configured J/ψ filter label |
| `muIsUpsFilterMatch` | Same for Upsilon filters |
| `muJpsiMatchedTriggerIndices` | Per-muon: list of trigger pattern indices (into `TriggersForJpsi`) that this muon matched |
| `muJpsiMatchedFilterIndices` | Per-muon: list of filter label indices (into `FiltersForJpsi`) that this muon matched |
| `muUpsMatchedTriggerIndices` | Per-muon: list of trigger pattern indices (into `TriggersForUpsilon`) that this muon matched |
| `muUpsMatchedFilterIndices` | Per-muon: list of filter label indices (into `FiltersForUpsilon`) that this muon matched |

**Configuration reference** (from `X_Config_Tree`):
`TriggersForJpsi`, `FiltersForJpsi`, `TriggersForUpsilon`, `FiltersForUpsilon` —
the configured trigger path substrings and filter labels against which matching
was performed.

**Trigger matching workflow:**

1. Check `MatchJpsiTriggerNames` is non-empty — at least one J/ψ trigger fired.
2. For each `Jpsi_1` / `Jpsi_2` candidate, use `Jpsi_1_mu_1_Idx` /
   `Jpsi_1_mu_2_Idx` to look up the muon indices.
3. Check `muJpsiMatchedTriggerIndices[muIdx]` — non-empty means this muon's
   trigger objects matched a configured J/ψ trigger path.
4. A candidate passes trigger matching when at least one of its daughter muons
   has a non-empty `muJpsiMatchedTriggerIndices`.

The `muIsJpsiTrigMatch` flag alone is NOT sufficient for trigger-object
matching — it only indicates the trigger fired for the event, not that a
specific muon was responsible. Use `muJpsiMatchedTriggerIndices` for the
per-muon association.

### 2.9 four-muon vertexing — $\varepsilon_{4\mu\mathrm{vtx}}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed HLT) | Events passing HLT trigger + trigger-object matching, with at least one valid dimuon pair for each of the two J/ψ slots |
| Numerator (passes THIS step) | Denominator events where `DiOnia_fitValid > 0` AND `DiOnia_fitPass > 0` |

This step combines two J/ψ candidates (`Jpsi_1` + `Jpsi_2`) into a common
four-muon vertex (Run B). The baseline criterion is that the kinematic vertex
fit converged and passes the vertex probability cut (`DiOniaVtxProbCut`).
Alternative or additional criteria (e.g., `DiOnia_commonRecVtxPass`) may be
applied depending on the analysis working point.

### 2.10 triOnia — $\varepsilon_{\mathrm{triOnia}}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed four-muon vertexing) | Events passing four-muon vertexing ($\varepsilon_{4\mu\mathrm{vtx}}$ numerator), split by φ pT bins |
| Numerator (passes THIS step) | Denominator events passing the chosen triOnia endpoint |

The code produces four parallel quality flags for the three-meson primary vertex:

| Flag | Meaning |
|------|---------|
| `Pri_fitValid > 0` | 3-body kinematic vertex fit converged |
| `Pri_fitPass > 0` | Fit passes `PriVtxProbCut` |
| `Pri_assocPVPass > 0` | Common PV association of daughter tracks |
| `Pri_trackPVPass > 0` | Track-PV compatibility for daughter tracks |

The default endpoint for corrected-yield studies is `Pri_assocPVPass` (per
`Efficiency_scheme.md`). The denominator is four-muon-vertexing-passing events,
split by φ pT bins. Alternative endpoints may use `Pri_fitPass` or
`Pri_trackPVPass` depending on the analysis.

---

## 3. Configuration considerations

### 3.1 Cuts you can tighten in post-processing

| Parameter | Why post-processable |
|-----------|---------------------|
| `JpsiCandPtMin`, `JpsiCandEtaMax` | These are applied in `pairMuons()` as pre-cuts on the composite candidate. The `SingleJpsi_pt` branch records the value; you can apply any pT cut in analysis. |
| `PhiCandPtMin`, `PhiCandEtaMax` | Same reasoning — `SinglePhi_pt` records the value. |
| Track quality (`normalizedChi2`, `numberOfValidHits`, `isHighPurity`) | Directly recorded in `RecoKaonTrack_normalizedChi2`, `RecoKaonTrack_numberOfHits`, `RecoKaonTrack_isHighPurity`. Configurable via `TrackQuality` string and `RequireRecoKaonTrackHighPurity` flag. |
| `MinTrackFromPV` | Recorded in `RecoKaonTrack_fromPV` and `SinglePhi_K*_fromPV`. |

### 3.2 MC control flags

| Flag | Effect on efficiency evaluation |
|------|--------------------------------|
| `keepAllSingleObjectCandsInMC=True` | **Required** — enables `SingleJpsi_*`, `SinglePhi_*` |
| `skipCompositeCandBuildingWhenKeepingSingles=True` | Optional — skips `combineCandidates()` for a smaller, faster ntuple when only per-object steps are needed |
| `requireAcceptedCandidatesForMonteCarloTree=False` | Recommended for efficiency — keeps every event even when no candidate passes |

### 3.3 Vertex fit toggles

| Toggle | Default | Notes |
|--------|---------|-------|
| `DoJpsiDecayVtxFit` | True | Needed for `SingleJpsi_fitValid` |
| `DoPhiDecayVtxFit` | True | Needed for `SinglePhi_fitValid` |
| `DoDiOniaVtxFit` | True | Needed for `DiOnia_*` branches |
| `DoPriVtxFit` | True | Needed for `Pri_*` branches |

---

## 4. Suggested analysis workflow

### Step 1: Produce efficiency ntuples

Run both Run A (singles-only) and Run B (full chain) on the full MC sample(s).
Use `maxEvents=-1` to process all events.

### Step 2: Build per-object efficiency maps

Each efficiency is conditional: $\varepsilon_{\mathrm{step}} = N(\text{pass step}) \;/\; N(\text{pass previous step})$.
For each J/ψ and φ kinematic bin `(pT, |y|)`:

1. **Acceptance** — $A_{J/\psi}$, $A_{\phi}$: From `MC_GenPart_*`.
   - Denominator: all GEN J/ψ (or φ) in the kinematic bin.
   - Numerator: denominator mesons where both GEN daughters are in the fiducial region.

2. **muonRECO** — $\varepsilon_{\mu\mathrm{Reco}|J/\psi}$: From `SingleJpsi_*` (Run A).
   - Denominator: `SingleJpsi` where both GEN daughters are in fiducial acceptance.
   - Numerator: denominator candidates with both `mu*_genMatchIdx >= 0`.

3. **kaonRECO** — $\varepsilon_{K\mathrm{Reco}|\phi}$: From `SinglePhi_*` (Run A).
   - Denominator: `SinglePhi` where both GEN daughters are in fiducial acceptance.
   - Numerator: denominator candidates with both `K*_genMatchIdx >= 0`.

4. **muonID** — $\varepsilon_{\mu\mathrm{ID}|J/\psi}$: From `SingleJpsi_*` + `mu*` (Run A).
   - Denominator: `SingleJpsi` passing muonRECO (both muons GEN-matched).
   - Numerator: denominator candidates where both daughter muons pass the chosen ID (baseline: `muIsPatSoftMuon`).

5. **kaonID** — $\varepsilon_{K\mathrm{ID}|\phi}$: From `SinglePhi_*` + `RecoKaonTrack_*` (Run A).
   - Denominator: `SinglePhi` passing kaonRECO (both kaons GEN-matched).
   - Numerator: denominator candidates where both daughter kaons satisfy the track-quality criteria:
     - `RecoKaonTrack_normalizedChi2[idx] < 8 && RecoKaonTrack_numberOfHits[idx] > 4`
     - `RecoKaonTrack_isHighPurity[idx] == 1`
   - The exact criteria are documented in `X_Config_Tree::TrackQuality` and `X_Config_Tree::RequireRecoKaonTrackHighPurity`.
   - Quality cuts can be varied offline without re-ntuplizing since the raw values are stored.

6. **dimuon** — $\varepsilon_{\mu\mu|J/\psi}$: From `SingleJpsi_*` (Run A).
   - Denominator: `SingleJpsi` passing muonID.
   - Numerator: denominator candidates with `SingleJpsi_fitValid && SingleJpsi_fitPass`.

7. **dikaon** — $\varepsilon_{KK|\phi}$: From `SinglePhi_*` (Run A).
   - Denominator: `SinglePhi` passing kaonID.
   - Numerator: denominator candidates with `SinglePhi_fitValid && SinglePhi_fitPass`.

### Step 3: Build event-level efficiency maps

From the full chain ntuple (Run B):

1. **HLT** — $\varepsilon_{\mathrm{HLT}}$:
   - Denominator: events with ≥1 valid dimuon candidate.
   - Numerator: denominator events where `MatchJpsiTriggerNames` is non-empty AND at least one daughter muon of the candidate has non-empty `muJpsiMatchedTriggerIndices`.

2. **four-muon vertexing** — $\varepsilon_{4\mu\mathrm{vtx}}$: From `DiOnia_*` (Run B).
   - Denominator: events passing HLT + trigger matching, with valid dimuon pairs in both J/ψ slots.
   - Numerator: denominator events with `DiOnia_fitValid && DiOnia_fitPass`.

3. **triOnia** — $\varepsilon_{\mathrm{triOnia}}$: From `Pri_*` (Run B), split by φ pT bins.
   - Denominator: events passing four-muon vertexing.
   - Numerator: denominator events passing the chosen endpoint (default: `Pri_assocPVPass`).

### Step 4: Apply correction

The per-event total efficiency weight is the product of all per-object and
event-level factors evaluated at the event's kinematics, times the factorized
acceptance. See `Efficiency_scheme.md` for the product formula.

---

## 5. Cross-checks

### 5.1 Factorized vs. cumulative comparison

Compare the factorized efficiency product against a direct "all steps pass"
count on the full chain ntuple. Significant disagreement indicates correlations
between the per-object steps that the factorized approach misses.

### 5.2 RecoKaonTrack coverage

The `RecoKaonTrack_*` block is always populated with the union of:
- GEN-matched fiducial kaons — all kaons passing phase-space `TrackSelection` (MC with `keepAllSingleObjectCandsInMC=True`), independent of track quality,
- φ-candidate daughter tracks from `KPairCand_Meson_` (always).

The `RecoKaonTrack_usedInSinglePhi` flag distinguishes the two sources
(1 = φ-candidate daughter, 0 = GEN-matched only).

Verify that every φ-candidate daughter has a `RecoKaonTrack_*` entry:
`SinglePhi_K*_RecoKaonTrackIdx >= 0` for all `SinglePhi` candidates, and
`Phi_K_*_RecoKaonTrackIdx >= 0` for all composite φ candidates. A `-1`
sentinel should never occur.

### 5.3 Bin-edge sensitivity

Vary the fine/coarse bin thresholds (`N_min_fine`, `N_min_coarse`) and confirm
the corrected yield is stable within MC statistical uncertainties.

### 5.4 GEN-RECO daughter matching consistency

For `SinglePhi` candidates, the inline GEN-match info (`SinglePhi_K1_genMatchIdx`,
`SinglePhi_K2_genMatchIdx`) should agree with the RecoKaonTrack-level match
(`RecoKaonTrack_genMatchIdx[SinglePhi_K1_RecoKaonTrackIdx]`,
`RecoKaonTrack_genMatchIdx[SinglePhi_K2_RecoKaonTrackIdx]`). They are derived
from the same `buildPhiKaonDiagnostics` call and must be identical.

For composite φ candidates (`Phi_K_1_*` / `Phi_K_2_*`), the inline
`Phi_K_1_genMatchIdx` / `Phi_K_2_genMatchIdx` can similarly be verified against
`RecoKaonTrack_genMatchIdx[Phi_K_1_RecoKaonTrackIdx]` /
`RecoKaonTrack_genMatchIdx[Phi_K_2_RecoKaonTrackIdx]`.
