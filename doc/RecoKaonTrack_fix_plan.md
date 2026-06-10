# Fix: Apply track quality filter to GEN-matched tracks in RecoKaonTrack

## Problem

The current Pass A in `fillRecoKaonTrackBlockForMC()` (line 1863-1864) accepts every GEN-matched track unconditionally:

```cpp
if (diagnostics.genMatchIdx >= 0) {
    selectedNonMuonTrackIdx.insert(nonMuonIdx);
}
```

This stores 82 GEN-only tracks in 50 events, many of which are soft (pt down to 0.28 GeV) or forward (|η| up to 2.87) — tracks that would never enter a phi candidate because they fail the quality cuts in `pairTracks()`. For acceptance evaluation only GEN-level is needed (no RECO block). For efficiency, only tracks within fiducial + quality acceptance matter.

## Fix (applied in v1.8)

The solution factorizes the kaon selection into three tiers:

1. **`TrackSelection`** (phase-space only): `"pt > 1.0 && abs(eta) < 2.5"` — defines the fiducial denominator for the RecoKaonTrack block.
2. **`TrackQuality`** + **`RequireRecoKaonTrackHighPurity`**: new config parameters defining the kaonID quality criteria, applied in `pairTracks()` to gate φ candidate building and stored in the config TTree for offline numerator computation.
3. **Raw quality values** (`RecoKaonTrack_normalizedChi2`, `_numberOfHits`, `_isHighPurity`): stored per-track in the RecoKaonTrack block (added in v1.7), enabling offline variation of quality cuts without re-ntuplizing.

The hardcoded quality cuts in `fillRecoKaonTrackBlock()` Step A are removed
entirely. Step A now requires only phase-space `trackSelector_` + `genMatchIdx >= 0`.
The hardcoded cuts in `pairTracks()` are replaced with the configurable
`trackQualitySelector_(*bestTrack)` and `requireRecoKaonTrackHighPurity_` flag.

## Code changes (v1.8, applied)

### `fillRecoKaonTrackBlock()` Step A (simplified)

**File:** `src/MultiLepPAT.cc`, the Pass A loop.

The entry criteria are simplified to phase-space + GEN-match only (no quality cuts):

```cpp
    for (unsigned int nonMuonIdx = 0; nonMuonIdx < nonMuonTrack_.size(); ++nonMuonIdx) {
        const auto& cand = *nonMuonTrack_[nonMuonIdx];
        const auto diagnostics = buildPhiKaonDiagnostics(cand, thePrimaryV_, genParticles);
        diagnosticsCache[nonMuonIdx] = diagnostics;

        const bool passesTrackSel = trackSelector_(cand);
        if (diagnostics.genMatchIdx >= 0 && passesTrackSel) {
            selectedNonMuonTrackIdx.insert(nonMuonIdx);
        }
    }
```

### `pairTracks()` quality cuts (configurable)

Hardcoded `normalizedChi2 <= 8.0` and `quality(reco::Track::highPurity)` replaced with:

```cpp
if (!trackQualitySelector_(*cand.bestTrack()) ||
    (requireRecoKaonTrackHighPurity_ &&
     !cand.bestTrack()->quality(reco::Track::highPurity))) continue;
```

### New config parameters

Added to constructor, config TTree, and all Python config files:
- `TrackQuality` (string, default `"normalizedChi2 < 8 && numberOfValidHits > 4"`)
- `RequireRecoKaonTrackHighPurity` (bool, default `True`)

## Three-tier selection summary

| Tier | Where applied | Cuts | Config |
|------|--------------|------|--------|
| Fiducial acceptance | `fillRecoKaonTrackBlock()` Step A | `trackSelector_` (p_T, η phase-space only) | `TrackSelection` |
| Reconstruction | GEN-matching in `buildPhiKaonDiagnostics()` | `genMatchIdx >= 0` (GEN kaon from φ mother) | `RecoGenKaonMatchChi2Max` |
| KaonID | `pairTracks()` quality gate + offline | `trackQualitySelector_` + `requireRecoKaonTrackHighPurity_` | `TrackQuality` + `RequireRecoKaonTrackHighPurity` |

## Expected impact

Compared to the pre-v1.8 architecture (where Step A applied hardcoded quality cuts):

- **RecoKaonTrack block expands**: Now includes all fiducial GEN-matched kaons,
  not just those passing quality cuts. This enables proper denominator definition
  for kaonID efficiency.
- **φ candidates unchanged**: With default quality settings, `pairTracks()`
  applies the same cuts as before, so φ candidate yields are identical.
- **Offline flexibility**: Raw quality values stored per-track allow kaonID
  working point to be varied without re-ntuplizing.

## Verification

1. **Build**: `scram b -j 4`
2. **Run on MC**: `cmsRun test/ConfFile_cfg.py AnalysisMode=JpsiJpsiPhi inputFiles=file:test_muonPackedMatch_mcRun2022_numEvent2000.root outputFile=test_recoKaonTrack.root`
3. **Check output**: GEN-only count should drop significantly; phi-only count unchanged; "neither" category stays zero.
