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
| Denominator (passed acceptance) | GEN J/ψ in acceptance (Section 2.1 numerator): direct J/ψ (`MC_GenPart_pdgId == 443`, `abs(MC_GenPart_motherPdgId) != 443`) with both GEN daughter muons in fiducial region |
| Numerator (passes THIS step) | Denominator GEN J/ψ where both GEN daughter muons have a RECO-level match in the `mu*` branches |

This step is evaluated entirely at GEN level: start from each acceptance-passing
GEN J/ψ, identify its two GEN daughter muons, and check whether both have a
corresponding reconstructed muon.

**Finding GEN daughter muons.** For each GEN J/ψ at index `j` in `MC_GenPart_*`,
its daughter muons are GEN particles where
`abs(MC_GenPart_pdgId) == 13` AND `MC_GenPart_motherGenIdx == j`.
Each direct J/ψ should have exactly two such daughters (μ⁺ and μ⁻).

**RECO matching.** The `muGenMatchIdx` branch (populated in `fillMuonBlock()` for
every RECO muon, regardless of whether it ends up in a `SingleJpsi` candidate)
stores the `MC_GenPart_*` index of the matched GEN muon. It uses unrestricted
matching — the GEN muon's mother PDG ID is NOT required to be 443. For each GEN
daughter muon at index `g`, check whether any RECO muon has `muGenMatchIdx == g`.

**Why not use `SingleJpsi` or `SingleJpsi_mu*_genMatchIdx`?** The `SingleJpsi`
branches are populated from `muPairCand_Onia1_`, which only stores pairs that
pass the J/ψ vertex fit (`onia1FitPass == true`, line 1643 of `MultiLepPAT.cc`).
Using them as the RECO denominator would fold the dimuon vertexing efficiency
into the RECO measurement. The inline `SingleJpsi_mu*_genMatchIdx` further
restricts matches to GEN muons whose mother is PDG 443, which is unnecessarily
narrow for the RECO step. The `muGenMatchIdx` branch has neither limitation.

### 2.3 kaonRECO — $\varepsilon_{K\mathrm{Reco}|\phi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed acceptance) | GEN φ in acceptance (Section 2.1 numerator): direct φ (`MC_GenPart_pdgId == 333`) with both GEN daughter kaons in fiducial region |
| Numerator (passes THIS step) | Denominator GEN φ where both GEN daughter kaons have a RECO-level match in the `RecoKaonTrack_*` block |

Same GEN-level approach as muonRECO. For each acceptance-passing GEN φ at
index `j`, its daughter kaons are GEN particles where
`abs(MC_GenPart_pdgId) == 321` AND `MC_GenPart_motherGenIdx == j`.

**RECO matching via `RecoKaonTrack_genMatchIdx`.** The `RecoKaonTrack_*` block is
populated independently of `SinglePhi` candidates — it contains the union of all
GEN-matched fiducial kaons and all φ-candidate daughter tracks (see Section 6.3).
For each GEN daughter kaon at index `g`, check whether any entry in the
`RecoKaonTrack_*` block has `RecoKaonTrack_genMatchIdx == g`.

**Why not use `SinglePhi` or `SinglePhi_K*_genMatchIdx`?** `SinglePhi` is
populated from `KPairCand_Meson_`, which only stores pairs that pass the φ vertex
fit (`particlesToVtx(... PhiDecayVtxProbCut_)`, line 1828 of `MultiLepPAT.cc`).
The `SinglePhi_K*_genMatchIdx` values are mother-PDG-restricted (require PDG 333).
The `RecoKaonTrack_genMatchIdx` is unrestricted and available for all
reconstructed tracks, making it the correct source for the RECO denominator.

`SinglePhi_K1_RecoKaonTrackIdx` / `SinglePhi_K2_RecoKaonTrackIdx` give indices
into the `RecoKaonTrack_*` block for φ-candidate daughter tracks. The
`RecoKaonTrack_genMatchIdx` field mirrors the inline `SinglePhi_K*_genMatchIdx`
for candidates that exist — they must agree (see Section 6.5). But the
`RecoKaonTrack` block covers additional tracks beyond φ candidates, which is
essential for the RECO step.

### 2.4 muonID — $\varepsilon_{\mu\mathrm{ID}|J/\psi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed muonRECO) | GEN J/ψ in acceptance where both GEN daughter muons have a RECO match (muonRECO numerator) |
| Numerator (passes THIS step) | Denominator GEN J/ψ where both RECO-matched daughter muons pass the chosen ID working point |

This step is evaluated on the same GEN J/ψ population as muonRECO, now checking
the ID quality of the matched RECO muons. For each RECO-matched GEN daughter muon,
find the RECO muon index `r` where `muGenMatchIdx[r] == gen_mu_idx`, then check
the ID flag at that index.

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
| kaonRECO | $\varepsilon_{K\mathrm{Reco}}^{\phi}$ (Section 2.3) | Both GEN daughter kaons have RECO match | `RecoKaonTrack_genMatchIdx >= 0` |
| kaonID | $\varepsilon_{K\mathrm{ID}}^{\phi}$ (this section) | Track quality on reconstructed kaons | `TrackQuality` + `RequireRecoKaonTrackHighPurity` |

Note: pT/eta cuts are part of fiducial acceptance, NOT kaonID. The `TrackSelection`
string configures the phase-space denominator (see Section 3.1).

| Role | Branch / logic |
|------|----------------|
| Denominator (passed kaonRECO) | GEN φ in acceptance where both GEN daughter kaons have a RECO match (kaonRECO numerator) |
| Numerator (passes THIS step) | Denominator GEN φ where both RECO-matched daughter kaons pass track-quality ID |

For each RECO-matched GEN daughter kaon, find the `RecoKaonTrack_*` entry index `k`
where `RecoKaonTrack_genMatchIdx[k] == gen_k_idx`, then check the track-quality
criteria at that index.

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
| Denominator (passed muonID) | GEN J/ψ in acceptance where both RECO-matched daughter muons pass ID (muonID numerator) |
| Numerator (passes THIS step) | Denominator GEN J/ψ where the two RECO+ID muons form a valid `SingleJpsi` candidate |

This step checks whether the pair of RECO+ID muons identified in the previous steps
survive the full dimuon selection: opposite-charge pair, mass window, pT/η pre-cuts,
and kinematic vertex fit. In the code, this corresponds to whether the pair appears
in `muPairCand_Onia1_` (which requires `onia1FitPass == true`, line 1643 of
`MultiLepPAT.cc`) and is stored in the `SingleJpsi_*` branches.

**Implementation.** For each GEN J/ψ, identify the two RECO muon indices `r1`, `r2`
(from the muonID step). Then check whether any `SingleJpsi` entry `j` satisfies:
`({SingleJpsi_mu1_Idx[j], SingleJpsi_mu2_Idx[j]} == {r1, r2})` (order-independent)
AND `SingleJpsi_fitValid[j] > 0` AND `SingleJpsi_fitPass[j] > 0`.

The baseline criterion for a valid dimuon candidate is that the kinematic vertex
fit converged and passes the vertex probability cut (`JpsiDecayVtxProbCut`). The
mass window is already applied in `pairMuons()`; `SingleJpsi_mass` is the
post-fit mass. Alternative or additional vertexing criteria may be applied
depending on the analysis working point.

**Why `SingleJpsi` appears only here.** In the corrected chain, `SingleJpsi` is
used for the FIRST time in this step. The RECO and ID steps operate directly on
GEN particles and the `mu*` branches, without needing `SingleJpsi`. This avoids
folding the dimuon vertexing/mass-window efficiency into the earlier steps.

### 2.7 dikaon — $\varepsilon_{KK|\phi}$

| Role | Branch / logic |
|------|----------------|
| Denominator (passed kaonID) | GEN φ in acceptance where both RECO-matched daughter kaons pass ID (kaonID numerator) |
| Numerator (passes THIS step) | Denominator GEN φ where the two RECO+ID kaons form a valid `SinglePhi` candidate |

Analogous to dimuon. For each GEN φ, identify the two `RecoKaonTrack_*` indices
`k1`, `k2` (from the kaonID step). Then check whether any `SinglePhi` entry `j`
satisfies: the pair of `RecoKaonTrack` indices matches `{SinglePhi_K1_RecoKaonTrackIdx[j],
SinglePhi_K2_RecoKaonTrackIdx[j]}` AND `SinglePhi_fitValid[j] > 0` AND
`SinglePhi_fitPass[j] > 0`. This is the first step where `SinglePhi` appears
in the efficiency chain.

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

**Prerequisite — Extract GEN mesons.** From `MC_GenPart_*`:
- Identify direct J/ψ: `abs(MC_GenPart_pdgId) == 443` and `abs(MC_GenPart_motherPdgId) != 443`.
  Compute rapidity `y` from `(px, py, pz, mass)` since no `MC_GenPart_y` branch exists.
- Identify direct φ: `MC_GenPart_pdgId == 333` (or `-333`, though typically `333`).
- For each such meson at index `m`, find its daughter muons (for J/ψ) or kaons (for φ)
  by scanning for GEN particles `d` where `MC_GenPart_motherGenIdx[d] == m` AND
  `abs(MC_GenPart_pdgId[d]) == 13` (muon) or `321` (kaon). Each direct J/ψ→μμ
  or φ→KK should have exactly two such daughters.

1. **Acceptance** — $A_{J/\psi}$, $A_{\phi}$: From `MC_GenPart_*` only.
   - Denominator: all direct GEN J/ψ (or φ) in the kinematic bin $(p_T, |y|)$.
   - Numerator: denominator mesons where both GEN daughter muons/kaons are in the
     fiducial region (phase-space cuts on $p_T$, $\eta$ of the daughters).

2. **muonRECO** — $\varepsilon_{\mu\mathrm{Reco}|J/\psi}$: From `mu*` + `MC_GenPart_*` (Run A).
   - Denominator: GEN J/ψ passing acceptance (numerator of step 1).
   - Numerator: denominator J/ψ where both GEN daughter muons have a RECO match:
     for each GEN daughter index `g`, check `any(muGenMatchIdx[:] == g)`.
   - **Do NOT** use `SingleJpsi` or `SingleJpsi_mu*_genMatchIdx` — those require the
     vertex fit to have passed (line 1643 of `MultiLepPAT.cc`) and a J/ψ mother.

3. **kaonRECO** — $\varepsilon_{K\mathrm{Reco}|\phi}$: From `RecoKaonTrack_*` + `MC_GenPart_*` (Run A).
   - Denominator: GEN φ passing acceptance (numerator of step 1).
   - Numerator: denominator φ where both GEN daughter kaons have a RECO match:
     for each GEN daughter index `g`, check `any(RecoKaonTrack_genMatchIdx[:] == g)`.
   - **Do NOT** use `SinglePhi` or `SinglePhi_K*_genMatchIdx` — those require the
     φ vertex fit to have passed (line 1828).

4. **muonID** — $\varepsilon_{\mu\mathrm{ID}|J/\psi}$: From `mu*` (Run A).
   - Denominator: GEN J/ψ passing muonRECO.
   - Numerator: denominator J/ψ where both RECO-matched muons pass the chosen ID:
     find the RECO muon index `r` where `muGenMatchIdx[r] == g` for each GEN daughter `g`,
     then check `muIsPatSoftMuon[r] == 1` (baseline working point).

5. **kaonID** — $\varepsilon_{K\mathrm{ID}|\phi}$: From `RecoKaonTrack_*` (Run A).
   - Denominator: GEN φ passing kaonRECO.
   - Numerator: denominator φ where both RECO-matched kaons pass track-quality ID:
     find the `RecoKaonTrack_*` index `k` where `RecoKaonTrack_genMatchIdx[k] == g`,
     then check `RecoKaonTrack_normalizedChi2[k] < 8`,
     `RecoKaonTrack_numberOfHits[k] > 4`, and `RecoKaonTrack_isHighPurity[k] == 1`.
   - The exact criteria are documented in `X_Config_Tree::TrackQuality` and
     `X_Config_Tree::RequireRecoKaonTrackHighPurity`.
   - Quality cuts can be varied offline without re-ntuplizing since the raw values
     are stored.

6. **dimuon** — $\varepsilon_{\mu\mu|J/\psi}$: From `SingleJpsi_*` + `mu*` (Run A).
   - Denominator: GEN J/ψ passing muonID.
   - Numerator: denominator J/ψ where the two RECO+ID muon indices `(r1, r2)` appear
     together as a pair in `SingleJpsi` with `SingleJpsi_fitValid > 0` AND
     `SingleJpsi_fitPass > 0`. Check both orderings: `{SingleJpsi_mu1_Idx[j],
     SingleJpsi_mu2_Idx[j]} == {r1, r2}`.
   - This is the **first** step where `SingleJpsi` branches are used. The RECO and
     ID steps operate entirely from `MC_GenPart_*` and `mu*` branches.

7. **dikaon** — $\varepsilon_{KK|\phi}$: From `SinglePhi_*` + `RecoKaonTrack_*` (Run A).
   - Denominator: GEN φ passing kaonID.
   - Numerator: denominator φ where the two RECO+ID RecoKaonTrack indices `(k1, k2)`
     appear together as a pair in `SinglePhi` with `SinglePhi_fitValid > 0` AND
     `SinglePhi_fitPass > 0`. Match via `SinglePhi_K1_RecoKaonTrackIdx` /
     `SinglePhi_K2_RecoKaonTrackIdx`.

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

## 5. Vectorized analysis with Awkward Arrays

The corrected efficiency chain requires per-event GEN-to-RECO matching across
flat vector branches. This section describes how to implement the matching
efficiently using Awkward Arrays (loaded via uproot).

### 5.1 Data model

The ntuple stores flat variable-length arrays per event. Key branches:

| Branch | Type | Shape |
|--------|------|-------|
| `MC_GenPart_pdgId` | `var * int32` | (ngen,) |
| `MC_GenPart_motherPdgId` | `var * int32` | (ngen,) |
| `MC_GenPart_motherGenIdx` | `var * int32` | (ngen,) — index in `MC_GenPart_*`, -1 if none |
| `MC_GenPart_pt`, `_eta`, `_phi` | `var * float32` | (ngen,) |
| `MC_GenPart_px`, `_py`, `_pz`, `_mass` | `var * float32` | (ngen,) |
| `muGenMatchIdx` | `var * int32` | (nmu,) — index in `MC_GenPart_*`, -1 if unmatched |
| `muIsPatSoftMuon` | `var * int32` | (nmu,) |
| `RecoKaonTrack_genMatchIdx` | `var * int32` | (nk,) — index in `MC_GenPart_*`, -1 if unmatched |
| `RecoKaonTrack_normalizedChi2` | `var * float32` | (nk,) |
| `RecoKaonTrack_numberOfHits` | `var * int32` | (nk,) |
| `RecoKaonTrack_isHighPurity` | `var * int32` | (nk,) |
| `SingleJpsi_mu1_Idx`, `_mu2_Idx` | `var * int32` | (nsj,) — indices in `mu*` |
| `SingleJpsi_fitValid`, `_fitPass` | `var * int32` | (nsj,) |
| `SinglePhi_K1_RecoKaonTrackIdx`, `_K2_` | `var * int32` | (nsp,) — indices in `RecoKaonTrack_*` |
| `SinglePhi_fitValid`, `_fitPass` | `var * int32` | (nsp,) |

No `MC_GenPart_y` branch exists — compute it from `(E, pz)` where
`E = sqrt(px² + py² + pz² + mass²)`.

### 5.2 Finding GEN J/ψ and their daughter muons

```python
import awkward as ak
import numpy as np

# Load data — Run A (singles) is sufficient for per-object efficiencies
# arrays = uproot.open("eff_singles.root:mkcands/X_data").arrays(library="ak")

gen_pdg = arrays["MC_GenPart_pdgId"]
gen_mother_pdg = arrays["MC_GenPart_motherPdgId"]
gen_mother_idx = arrays["MC_GenPart_motherGenIdx"]

# Select direct J/psi: pdgId == 443, mother is not also 443
is_direct_jpsi = (abs(gen_pdg) == 443) & (abs(gen_mother_pdg) != 443)
jpsi_idx = ak.local_index(gen_pdg)[is_direct_jpsi]           # (events, njpsi)

# Find GEN daughter muons for each J/psi
# Broadcast mother match: (events, njpsi, ngen)
is_muon = (abs(gen_pdg) == 13)
is_daughter = gen_mother_idx == jpsi_idx[:, :, np.newaxis]
is_jpsi_mu = is_muon & is_daughter                           # (events, njpsi, ngen)
daughter_mu = ak.local_index(gen_pdg)[is_jpsi_mu]           # (events, njpsi, 2+)
# Each direct J/psi->mumu should have exactly 2 daughters
```

**Memory caveat.** The broadcast `gen_mother_idx == jpsi_idx[:, :, np.newaxis]`
produces an intermediate array of shape `(events, njpsi, ngen)`. For events
with many GEN particles (O(10³)) this can be large. Pre-filter the GEN table
to interesting PDG IDs (443, 553, 333, 13, 321) to reduce `ngen`, or use
Numba-accelerated loops for the daughter-finding step.

### 5.3 RECO matching (muonRECO)

```python
mu_gen_match = arrays["muGenMatchIdx"]  # (events, nmu)

# For each GEN daughter muon, check if any RECO muon matches it
# Broadcast: (events, njpsi, 2, nmu)
reco_ok = ak.any(mu_gen_match == daughter_mu, axis=-1)  # (events, njpsi, 2)
passes_reco = ak.all(reco_ok, axis=-1)                   # (events, njpsi)
```

### 5.4 ID check (muonID)

```python
mu_soft = arrays["muIsPatSoftMuon"]

# For each GEN daughter, find the RECO index r where muGenMatchIdx[r] == gen_dau
match_pos = ak.local_index(mu_gen_match)                 # (events, nmu)
is_match = mu_gen_match == daughter_mu                   # (events, njpsi, 2, nmu)
# Take the best matching RECO index per daughter
reco_for_dau = ak.max(ak.where(is_match, match_pos, -1), axis=-1)  # (events, njpsi, 2)

valid = reco_for_dau >= 0
safe_idx = ak.where(valid, reco_for_dau, 0)
id_pass = ak.where(valid, mu_soft[safe_idx] == 1, False)
passes_id = ak.all(id_pass, axis=-1)                     # (events, njpsi)
```

### 5.5 Dimuon check

```python
sj_mu1 = arrays["SingleJpsi_mu1_Idx"]       # (events, nsj)
sj_mu2 = arrays["SingleJpsi_mu2_Idx"]       # (events, nsj)
sj_fit_ok = (arrays["SingleJpsi_fitValid"] > 0) & (arrays["SingleJpsi_fitPass"] > 0)

# reco_for_dau has shape (events, njpsi, 2) — the RECO indices for each daughter
# valid has shape (events, njpsi, 2) — True if that daughter had a RECO match

# We need both daughters valid, so we can reconstruct r1, r2
r1 = reco_for_dau[:, :, 0]  # (events, njpsi)
r2 = reco_for_dau[:, :, 1]  # (events, njpsi)

# Check ordered and reversed against SingleJpsi entries
match_ordered = (sj_mu1[:, np.newaxis, :] == r1[:, :, np.newaxis]) & \
                (sj_mu2[:, np.newaxis, :] == r2[:, :, np.newaxis])
match_reversed = (sj_mu1[:, np.newaxis, :] == r2[:, :, np.newaxis]) & \
                 (sj_mu2[:, np.newaxis, :] == r1[:, :, np.newaxis])
pair_match = match_ordered | match_reversed              # (events, njpsi, nsj)

passes_dimuon = ak.any(pair_match & sj_fit_ok[:, np.newaxis, :], axis=-1)
passes_dimuon = passes_dimuon & passes_id                # only where ID also passed
```

### 5.6 Kaon-side matching

Analogous patterns using `RecoKaonTrack_*` branches:

```python
# kaonRECO
k_gen_match = arrays["RecoKaonTrack_genMatchIdx"]
reco_k_ok = ak.any(k_gen_match == daughter_k, axis=-1)
passes_k_reco = ak.all(reco_k_ok, axis=-1)

# kaonID — find RecoKaonTrack index, then check quality
k_match_pos = ak.local_index(k_gen_match)
k_is_match = k_gen_match == daughter_k
k_reco_idx = ak.max(ak.where(k_is_match, k_match_pos, -1), axis=-1)
valid_k = k_reco_idx >= 0
safe_k = ak.where(valid_k, k_reco_idx, 0)

k_chi2 = arrays["RecoKaonTrack_normalizedChi2"][safe_k]
k_nhits = arrays["RecoKaonTrack_numberOfHits"][safe_k]
k_hp = arrays["RecoKaonTrack_isHighPurity"][safe_k]
k_id_ok = ak.where(valid_k, (k_chi2 < 8) & (k_nhits > 4) & (k_hp == 1), False)
passes_k_id = ak.all(k_id_ok, axis=-1) & ak.all(valid_k, axis=-1)

# dikaon — check pair in SinglePhi
sp_k1 = arrays["SinglePhi_K1_RecoKaonTrackIdx"]
sp_k2 = arrays["SinglePhi_K2_RecoKaonTrackIdx"]
sp_fit_ok = (arrays["SinglePhi_fitValid"] > 0) & (arrays["SinglePhi_fitPass"] > 0)

kr1 = k_reco_idx[:, :, 0]
kr2 = k_reco_idx[:, :, 1]
phi_ordered = (sp_k1[:, np.newaxis, :] == kr1[:, :, np.newaxis]) & \
              (sp_k2[:, np.newaxis, :] == kr2[:, :, np.newaxis])
phi_reversed = (sp_k1[:, np.newaxis, :] == kr2[:, :, np.newaxis]) & \
               (sp_k2[:, np.newaxis, :] == kr1[:, :, np.newaxis])
phi_pair_match = phi_ordered | phi_reversed
passes_dikaon = ak.any(phi_pair_match & sp_fit_ok[:, np.newaxis, :], axis=-1)
passes_dikaon = passes_dikaon & passes_k_id
```

### 5.7 Building per-bin efficiency maps

Flatten the nested arrays to build per-object DataFrames for binning:

```python
import pandas as pd

# GEN J/psi kinematics — compute rapidity
gen_px = arrays["MC_GenPart_px"][jpsi_idx]
gen_py = arrays["MC_GenPart_py"][jpsi_idx]
gen_pz = arrays["MC_GenPart_pz"][jpsi_idx]
gen_mass = arrays["MC_GenPart_mass"][jpsi_idx]
gen_e = np.sqrt(gen_px**2 + gen_py**2 + gen_pz**2 + gen_mass**2)
gen_y = 0.5 * np.log((gen_e + gen_pz) / (gen_e - gen_pz + 1e-15))
gen_pt = arrays["MC_GenPart_pt"][jpsi_idx]

# Build flat DataFrames for histogram-based efficiency binning
df_muon_reco = pd.DataFrame({
    "pt": ak.to_numpy(ak.flatten(gen_pt)),
    "abs_y": ak.to_numpy(np.abs(ak.flatten(gen_y))),
    "pass": ak.to_numpy(ak.flatten(passes_reco)),
})
# Then: efficiency = df.groupby(pT_bins, y_bins)["pass"].mean()
```

---

## 6. Cross-checks

### 6.1 Factorized vs. cumulative comparison

Compare the factorized efficiency product against a direct "all steps pass"
count on the full chain ntuple. In the factorized approach, each step's
denominator is the previous step's numerator at GEN level; in the cumulative
approach, count GEN mesons that pass ALL steps (acceptance → RECO → ID →
vertexing) in one pass. Significant disagreement indicates correlations
between the per-object steps that the factorized approach misses.

### 6.2 muGenMatchIdx vs SingleJpsi_mu*_genMatchIdx consistency

The `muGenMatchIdx` branch (unrestricted, matches any GEN muon) and
`SingleJpsi_mu*_genMatchIdx` (requires mother PDG == 443) serve different
purposes in the corrected chain:

- `muGenMatchIdx` is used for the RECO step — it answers "does this GEN muon
  (regardless of origin) have a RECO counterpart?"
- `SingleJpsi_mu*_genMatchIdx` is used for the dimuon step — it identifies the
  specific GEN J/ψ mother.

For a GEN muon that IS from a J/ψ decay and appears in a `SingleJpsi` candidate,
the two should agree: `muGenMatchIdx[r] == SingleJpsi_mu*_genMatchIdx[j]` where
`r == SingleJpsi_mu*_Idx[j]`. Verify this equality for all `SingleJpsi` entries.
A discrepancy indicates the chi2-based matching in `fillMuonBlock()` found a
different best match than the mother-restricted matching in
`storeSingleJpsiCandidatesForMC()`.

### 6.3 RecoKaonTrack coverage

The `RecoKaonTrack_*` block is always populated with the union of:
- GEN-matched fiducial kaons — all kaons passing phase-space `TrackSelection` (MC with `keepAllSingleObjectCandsInMC=True`), independent of track quality,
- φ-candidate daughter tracks from `KPairCand_Meson_` (always).

The `RecoKaonTrack_usedInSinglePhi` flag distinguishes the two sources
(1 = φ-candidate daughter, 0 = GEN-matched only).

Verify that every φ-candidate daughter has a `RecoKaonTrack_*` entry:
`SinglePhi_K*_RecoKaonTrackIdx >= 0` for all `SinglePhi` candidates, and
`Phi_K_*_RecoKaonTrackIdx >= 0` for all composite φ candidates. A `-1`
sentinel should never occur.

### 6.4 Bin-edge sensitivity

Vary the fine/coarse bin thresholds (`N_min_fine`, `N_min_coarse`) and confirm
the corrected yield is stable within MC statistical uncertainties.

### 6.5 GEN-RECO daughter matching consistency

This cross-check is **diagnostic only** — it verifies internal consistency of
the ntuple, not the efficiency calculation itself. The efficiency chain uses
`muGenMatchIdx` and `RecoKaonTrack_genMatchIdx` directly, not the inline match
indices on the candidate branches.

For `SinglePhi` candidates, the inline GEN-match info (`SinglePhi_K1_genMatchIdx`,
`SinglePhi_K2_genMatchIdx`) should agree with the RecoKaonTrack-level match
(`RecoKaonTrack_genMatchIdx[SinglePhi_K1_RecoKaonTrackIdx]`,
`RecoKaonTrack_genMatchIdx[SinglePhi_K2_RecoKaonTrackIdx]`). They are derived
from the same `buildPhiKaonDiagnostics` call and must be identical.

For composite φ candidates (`Phi_K_1_*` / `Phi_K_2_*`), the inline
`Phi_K_1_genMatchIdx` / `Phi_K_2_genMatchIdx` can similarly be verified against
`RecoKaonTrack_genMatchIdx[Phi_K_1_RecoKaonTrackIdx]` /
`RecoKaonTrack_genMatchIdx[Phi_K_2_RecoKaonTrackIdx]`.

Similarly for muons: verify that for each `SingleJpsi` entry `j`,
`muGenMatchIdx[SingleJpsi_mu1_Idx[j]] == SingleJpsi_mu1_genMatchIdx[j]` and
`muGenMatchIdx[SingleJpsi_mu2_Idx[j]] == SingleJpsi_mu2_genMatchIdx[j]`.
Any disagreement indicates that the unrestricted chi2 matching in
`fillMuonBlock()` selected a different GEN match than the mother-restricted
matching in `storeSingleJpsiCandidatesForMC()`.
