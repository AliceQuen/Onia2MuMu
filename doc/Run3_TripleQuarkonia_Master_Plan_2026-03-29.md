# Run 3 Triple-Quarkonia Analysis Master Plan

Date: 2026-03-29
Updated: 2026-04-13

Analysts:
- Chi Wang
- Xing Cheng (`Endymion2288`)

Physics channels in scope:
- `JpsiJpsiUpsilon -> 6mu`
- `JpsiJpsiPhi -> 4mu + 2K`
- `JpsiUpsilonPhi -> 4mu + 2K`

Working milestones requested by the analysis team:
- `2026-04-01 12:00 Europe/Zurich`: report slide decks ready
- `2026-04-02`: individual analysis reports delivered
- `2026-04-29`: internal flagship-analysis milestone for first official efficiencies and first component-separation result
- `2026-05-20`: bachelor thesis nearly complete in writing, even if the full final analysis result is not yet frozen
- `2026-08-31`: flagship-analysis freeze target for review-grade summer maturity
- `November 2026`: projected pre-approval month, exact week to be fixed with the review chain

## 1. Executive Summary

This program is already beyond the "idea stage". There is now clear evidence of active reconstruction code, MC production, ntuple production, validation studies, fit-oriented workspaces, and channel-specific analysis paths across local EOS areas and GitHub repositories. The most mature common analysis backbone is the refactored `TPS-Onia2MuMu` package inside `CMSSW_15_0_15_JpsiJpsiPhi_refactor`, with additional dedicated work in `MuonPackedPFCandMatch`, `ParticleCand`, `ntuples-toHTCondor`, the `JpsiUpsPhi` production-and-hadd area, and the gen-level `JpsiJpsiPhi_MiniAOD_GenAnalyzer`.

As of `2026-04-13`, the shared `TPS-Onia2MuMu` EDAnalyzer version for the flagship `JpsiJpsiPhi` path is already pinned down and submitted to CRAB for full `2022-2025` processing. The immediate work is therefore no longer to freeze the flagship analyzer itself, but to monitor the running campaign, keep bookkeeping synchronized, and turn the pinned baseline into validation, efficiency, and component-separation products.

The key strategic conclusion is that the updated schedule is realistic only if the work is tiered against the real deadlines now in hand:
- `2026-04-01 12:00 Europe/Zurich` and `2026-04-02` are status-and-direction deadlines, not result deadlines.
- `2026-05-20` is a thesis-writing deadline, not a full-analysis freeze.
- `November 2026` is the review deadline for a pre-approval-grade flagship result.

The recommended channel strategy is therefore:
- `JpsiJpsiPhi` should remain the flagship near-term channel for the `2026-04-29` and summer freeze milestones.
- `JpsiJpsiUpsilon` should be kept as the MC-rich and six-muon-methods secondary channel.
- `JpsiUpsilonPhi` should be treated as an actively maintained secondary channel with stronger production-assembly maturity than was apparent in the first pass, especially on the 2022-2024 ntuple-merging and trigger-inspection side.

The most important CMS-style adjustment is to define the near-term deliverable as:
- baseline corrected yields or fiducial cross sections for the lead channel
- plus component separation by likelihood or template fit
- plus efficiency and closure validation

rather than promising a fully mature "unfolding" program for all three channels at once. In CMS language, the immediate task is better described as `component separation plus acceptance-times-efficiency correction`, not a full detector unfolding, unless differential observables are explicitly introduced.

### 1.1 Updated Hard Deadlines

- `2026-04-01 12:00 Europe/Zurich`: slide decks completed for the next-day reporting cycle
- `2026-04-02`: analyst reports delivered, with the thesis mid-term reporting step carried out smoothly
- `2026-04-29`: internal technical milestone for first official efficiency objects and first flagship component-separation result
- `2026-05-20`: bachelor thesis text nearly complete, with methods, datasets, workflow, and preliminary physics material stable
- `2026-08-31`: flagship-analysis freeze target for the first review-grade package
- `November 2026`: projected pre-approval month

### 1.2 Live Progress Tracking Boards

#### Chi Wang Tracker

Last updated:
- `2026-04-13`

Weekly note:
- flagship `TPS-Onia2MuMu` baseline is pinned and the full `2022-2025` `JpsiJpsiPhi` CRAB campaign is already in flight; the active work has moved to bookkeeping, validation, and baseline-policy decisions.

Artifacts or links:
- `test/crabData/`
- `test/crabScript_New/`
- `doc/multileppat_vertex_batch_integration_plan.md`
- `doc/study_multileppat_vertexing_full_mc_match_plan.md`

This week:
- [ ] monitor and summarize CRAB task status from `test/crabData/` for the pinned `2022-2025` `JpsiJpsiPhi` campaign
- [ ] assemble one code, data, and MC bookkeeping table shared by all three channels
- [ ] summarize object-validation status and decide whether the packed-candidate study changes the nominal baseline or remains a validation appendix

In progress:
- [ ] integrate the downstream `JpsiJpsiPhi` selector-comparison and mass-study workflow around `multileppat_vertex_batch`

Blocked:
- [ ]

Done:
- [x] pinned the flagship `TPS-Onia2MuMu` EDAnalyzer version and submitted it to CRAB for full `2022-2025` processing
- [x] completed the `2026-04-01` to `2026-04-02` report sprint, including the thesis mid-term report

#### Xing Cheng Tracker

Last updated:
- `2026-04-13`

Weekly note:
- `JpsiJpsiPhi` downstream fit and selector studies are now converging onto a package-level workflow instead of notebook-local logic.

Artifacts or links:
- `doc/multileppat_vertex_batch_integration_plan.md`
- `doc/study_mass_spectra_roofit_run2023d.ipynb`
- `doc/study_mass_spectra_iminuit_run2023d.ipynb`

This week:
- [ ] compare the flagship selector choices in the integrated RooFit and `iminuit` workflows
- [ ] identify the nominal observable set and fit-stability criteria for the April flagship result
- [ ] keep the secondary-channel fit and production notes synchronized with the shared bookkeeping table

In progress:
- [ ] compact flagship selector-comparison workflow for `JpsiJpsiPhi`

Blocked:
- [ ]

Done:
- [x] completed the `2026-04-02` reporting checkpoint

#### Shared Tracker

- [x] `2026-04-01 12:00 Europe/Zurich`: slide decks frozen
- [x] `2026-04-02`: reports delivered
- [x] flagship-channel declaration recorded in writing
- [ ] one shared bookkeeping table for all channels exists
- [ ] one MC campaign-definition cross-check exists between GitHub production repos and the shared CMSSW analysis
- [ ] thesis chapter outline frozen

## 2. Evidence-Based Current Status

### 2.1 Shared reconstruction and ntupling backbone already exists

The strongest common code base is:
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/TPS-Onia2MuMu/`

What is already present there:
- multi-channel `AnalysisMode` support for `JpsiJpsiPhi`, `JpsiJpsiUps`, and `JpsiUpsPhi`
- refactored runtime-configurable object and mass cuts
- dataset and operations notes
- validation test planning
- active CRAB support for recent run eras, with `test/crabData/` as the main CRAB task-state directory and `test/crabScript_New/` as the year-split generation and helper layer
- a pinned `TPS-Onia2MuMu` EDAnalyzer version already submitted to CRAB for full `2022-2025` flagship processing
- integrated downstream `JpsiJpsiPhi` mass-study and selector-comparison interfaces in `test/multileppat_vertex_batch/`
- new `doc/` content covering validation tests, mass-study integration, and vertex/truth study design

This means the shared analyzer baseline is now fixed for large-scale flagship running. The active development has shifted to CRAB monitoring, bookkeeping, validation studies, and downstream fit preparation rather than core EDAnalyzer churn.

### 2.2 Muon-to-packed-candidate validation is already a dedicated subproject

There is a separate repository and CMSSW package:
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/MuonPackedPFCandMatch/`
- GitHub remote: `Eric100911/MuonPackedPFCandMatch`

Already demonstrated:
- dedicated analyzer, notebook, and technical note
- smoke-tested on both MC and data
- one-event and O(10^3)-event debug outputs already present locally
- matching diagnostics for several matching definitions

This is important because it reduces one of the major detector-object risks for low-pt multi-muon analyses. It is not yet the final physics result, but it is exactly the kind of object-level validation that CMS internal review expects.

### 2.3 Strong local evidence of data and MC production

`JpsiJpsiPhi` local status:
- `test/crabData/` already contains `JpsiJpsiPhi` CRAB task directories across `Run2022C-G`, `Run2023C-D`, `Run2024C-I`, and `Run2025C-G` for `ParkingDoubleMuonLowMass[0-7]`
- `test/crabScript_New/` now provides the single-template, year-split CRAB config generation layer for those run eras
- `rootNtuple/ParkingDoubleMuonLowMass[0-7]/crab3_noTriVtx_*_Run2025C...G_MINIAOD` directories already exist
- GEN-level, HepMC, miniAOD, and ntuple directories exist under `MC_samples/`
- local test outputs exist for:
  - `test_JpsiJpsiPhi.root`
  - `test_JpsiJpsiUps.root`
  - `mymultilep_JpsiUpsPhi_dataRun2023D.root`
  - `mymultilep_JpsiUpsPhi_mcRun2022_TPS.root`
- `JpsiJpsiPhi` downstream mass studies are now organized around `multileppat_vertex_batch`, with package-level selector comparison and fit-preparation interfaces already smoke-tested on `Run2023D`

`JpsiJpsiUpsilon` local status:
- very substantial MC preparation under `JpsiJpsiUps/MC_samples/`
- DPS, SPS, and TPS directories exist across LHE, HepMC, GENSIM, miniAOD, and root ntuple stages
- dedicated TPS miniAOD and ntuple outputs already exist
- `ParticleCand` repo is active on branch `6mu-roofit`
- `ntuples-toHTCondor` repo exists and is versioned
- many six-muon preselection and candidate-level files already exist

`JpsiUpsilonPhi` local status:
- `rootNtuple/ParkingDoubleMuonLowMass[0-7]/crab3_*_Run2023C/D_MINIAOD` directories already exist
- CMSSW no-triple-vertex area exists
- `rootNtuple/` contains `hadd_2022.sh`, `hadd_2023.sh`, `hadd_2023_half.sh`, and `hadd_2024.sh`
- `rootNtuple/` also contains `gen_hadd_2022.sh`, `gen_hadd_2023.sh`, and `gen_hadd_2024.sh` for campaign-level assembly command generation
- `viewTrigger.C` and `viewTrigger.h` exist for trigger-content inspection on merged ntuples
- local output files such as `mymultilep.root`, `mymultilep_Run2023_looseMassRange.root`, and `sample_JUP_2023.root` exist

Interpretation:
- `JpsiJpsiPhi` looks strongest on common framework maturity, already-running `2022-2025` CRAB processing, and the new integrated downstream mass-study path.
- `JpsiJpsiUpsilon` looks strongest on MC preparation and exploratory six-muon downstream analysis.
- `JpsiUpsilonPhi` is stronger on run-era production assembly than the first pass suggested, especially on the 2022-2024 merge workflow and trigger inspection side, though it is still less integrated than `JpsiJpsiPhi` in the shared refactored analysis framework.

### 2.4 GitHub-side work split is already visible

Shared reconstruction repo:
- `Eric100911/TPS-Onia2MuMu`
- local branch: `refactor`
- recent commits include:
  - fixing 2025 muon JSON handling
  - user proxy handling
  - major refactor completion
  - CRAB helper integration

Authorship signal from the local repo history:
- strong contributions from both `Eric100911` and `Endymion`
- merge of `Endymion2288/Dev-J-J-P` into the shared repo is already visible in the local history
- `Endymion` commits in the shared code base include adding `fromPV` and `pvQuality`, a new test program, and follow-up repository repair
- this is the clearest sign that the main reconstruction pipeline is already collaborative

Muon packed-candidate validation repo:
- `Eric100911/MuonPackedPFCandMatch`
- currently driven mainly by Chi Wang

Six-muon downstream candidate repo:
- `Eric100911/ParticleCand`
- active branch `6mu-roofit`
- recent commits indicate:
  - RooFit restructuring
  - candidate storage in ROOT
  - branch cleanup and six-muon analysis evolution

Xing Cheng public workspace repos:
- `Endymion2288/JpsiJpsiPhi-workspace`
- `Endymion2288/JpsiJpsiUps-workspace`
- `Endymion2288/JpsiUpsPhi-workspace`

Recent commit messages show work on:
- drawing pipeline
- first-edition fit pipeline
- adding `Lxy` pieces and a "Stefanos cut"
- bootstrapping `JpsiJpsiUps` workflow from the `JpsiUpsPhi` work

Additional workflow evidence consistent with Xing's production role:
- `Endymion2288/JpsiUpsPhi-workspace` is much larger than the other two public workspaces
- the matching local `JpsiUpsPhi/rootNtuple/` area contains year-specific `gen_hadd` and `hadd` workflows for 2022, 2023, and 2024
- trigger-inspection helpers exist in the same area and read branches such as `MatchUpsTriggerNames`
- the `gen_hadd_2024.sh` script auto-generates merge commands from CRAB-style directory patterns, which is stronger evidence of campaign-level production handling than a simple hand-written one-off merge
- the concrete `hadd_2023_half.sh` and `hadd_2024.sh` scripts enumerate run-era and stream-specific merges across all `ParkingDoubleMuonLowMass[0-7]` streams

Interpretation:
- Xing's visible work is not only fit and plot development; there is also strong evidence of production-assembly and ntuple-merging workflow maturity, especially on the `JpsiUpsPhi` side
- Chi's visible GitHub work is strongest on CMSSW analyzer, object validation, production tooling, and MC-side infrastructure
- this is a good division of labor and should be made explicit instead of leaving it implicit

### 2.5 Xing's newer MC-production and gen-level repositories materially strengthen the status picture

Additional relevant repositories inspected:
- `Endymion2288/Full_MC_Production`
- `Endymion2288/JUPMCAnalyzer`
- `Endymion2288/JJP_MINIAOD_Analysis`
- `Endymion2288/TPS-Onia2MuMu`

`Endymion2288/Full_MC_Production` now provides the strongest explicit MC workflow evidence seen so far:
- a DAGMan generator for full-chain production from `LHE -> shower -> mix -> GEN-SIM -> RAW -> RECO -> MiniAOD -> ntuple`
- campaign definitions for both `JJP` and `JUP`, with SPS, DPS, and TPS topologies
- a universal `run_chain.sh` wrapper that handles LHE-pool resolution, Pythia showering, event mixing, CMSSW production, and EOS transfer
- explicit `ntuple_jjp_cfg.py` and `ntuple_jup_cfg.py` configs showing integration with JJP and JUP ntuple analyzers
- clear deployment assumptions across `CMSSW_12_4_14_patch3`, `CMSSW_14_0_18`, HTCondor, EOS, and containerized el9 execution for the ntuple stage

Recent commit evidence for `Full_MC_Production`:
- `2025-12-23`: `feat: first edition of LHE2MINIAOD`
- `2025-12-27`: `fix: successfully test`

`Endymion2288/JUPMCAnalyzer` shows a dedicated JUP MiniAOD truth-analysis layer:
- FWLite-based gen-correlation analysis for `SPS`, `DPS_1`, `DPS_2`, `DPS_3`, and `TPS`
- xrootd-based batch reading of campaign outputs
- plotting scripts and HTCondor submission flow
- recent commit `2026-01-14`: `feat: use HTCondor`

`Endymion2288/JJP_MINIAOD_Analysis` shows the corresponding JJP truth-analysis layer:
- MiniAOD truth-analysis for `SPS`, `DPS_1`, `DPS_2`, and `TPS`
- explicit iteration on topology definitions and parton-relation logic
- HTCondor support and multiprocessing updates
- recent commits include:
  - `2025-12-27`: `feat: change eta to y and use multiprocessing`
  - `2026-01-09`: `feat: strictly limit the parton process`
  - `2026-01-09`: `fix: do not limit parton relation in TPS`
  - `2026-01-14`: `feat: use HTCondor`

`Endymion2288/TPS-Onia2MuMu` remains useful as context:
- it documents the shared origin of the onia-analysis code and explicitly lists Chi and Xing as contributors
- its README is older and more `JpsiJpsiUpsilon`-oriented, but it still confirms that Xing's work is tied to the common analyzer lineage rather than being only external post-processing

Interpretation:
- Xing's role is now clearly broader than fit support or ad hoc merging; it spans full MC production design, truth-level topology studies, HTCondor execution, and downstream MiniAOD-based validation
- this materially strengthens the credibility of the JJP and JUP MC programs, especially for near-term efficiency, closure, and SPS/DPS/TPS template studies
- the immediate collaboration need is no longer to discover whether MC infrastructure exists, but to align the definitions and outputs of these repos with the shared `TPS-Onia2MuMu` analysis baseline

### 2.6 Gen-level study infrastructure already exists for JpsiJpsiPhi

Repo:
- `Eric100911/JpsiJpsiPhi_MiniAOD_GenAnalyzer`

The README clearly documents:
- gen-level correlation studies from MiniAOD
- SPS, DPS, and TPS mode handling
- plotting scripts
- HTCondor setup
- xrootd-based batch running

This is a major asset for:
- template construction
- truth-category closure tests
- acceptance studies
- DPS/TPS interpretation

### 2.7 The main gaps today

The project is advanced, but several items are still not yet in a review-ready state:

- There is still no single frozen baseline selection and ntuple schema for all three channels, even though the flagship `TPS-Onia2MuMu` EDAnalyzer version is now pinned and already in CRAB production.
- The local `TPS-Onia2MuMu` repo remains active, but the moving pieces are now concentrated in docs, notebook workflows, CRAB helper layers, and validation studies rather than the core flagship EDAnalyzer itself.
- `JpsiJpsiPhi` mass-spectrum studies are now integrated around `multileppat_vertex_batch`, but the MC truth batch path is not yet fully migrated to the same candidate skeleton.
- `MuonPackedPFCandMatch` still needs an explicit decision note stating whether it changes the nominal baseline or remains a validation appendix.
- `ParticleCand` contains many untracked files, macros, and ROOT outputs, which means the six-muon downstream chain is not yet fully reproducible from a clean checkout.
- `JpsiJpsiPhi` shows 2025 ntuple production locally, but `JpsiUpsilonPhi` does not yet show comparable 2024-2025 top-on in the same way.
- The new GitHub MC repositories still need a formal cross-check against the shared refactored analyzer so that campaign labels, topology definitions, cuts, and ntuple schemas mean the same thing everywhere.
- No unified response or efficiency production campaign is yet visible across all three channels.
- No explicit systematic-uncertainty registry is visible yet.
- No single "master note" yet ties together the three channels, the MC definitions, the fit models, the efficiency scheme, and the review path.

## 3. Recommended Physics Strategy

### 3.1 Measurement definition

For each channel, the minimal physics product should be:
- a fiducial signal yield or upper limit
- an acceptance-times-efficiency corrected yield
- a component-separated interpretation in terms of SPS, DPS, and possibly TPS

Near-term, the analysis should not over-promise a fully differential unfolding. Instead:
- use `component separation fit` for SPS/DPS/TPS-like contributions
- use acceptance and efficiency corrections for fiducial measurements
- reserve full detector unfolding for later only if a differential observable becomes central

### 3.2 Channel prioritization

Recommended prioritization as of `2026-03-29`:

Tier A, flagship:
- `JpsiJpsiPhi`

Tier B, secondary but active:
- `JpsiJpsiUpsilon`: strongest on MC depth and six-muon downstream development
- `JpsiUpsilonPhi`: strongest on run-era production assembly, hadding workflow, and trigger-inspection maturity

Reasoning:
- `JpsiJpsiPhi` still has the best evidence of an integrated path from gen study to 2025 data ntuple production and the shared refactor work.
- `JpsiJpsiUpsilon` remains the most mature secondary channel from the point of view of MC breadth and six-muon-specific downstream analysis.
- `JpsiUpsilonPhi` should no longer be treated as merely exploratory; it has meaningful workflow maturity, but should not outrun `JpsiJpsiPhi` on shared-framework stabilization unless the staffing split makes that efficient.

### 3.3 Thesis-level versus note-level deliverables

For the bachelor-thesis-level milestone on `2026-05-20`, the cleanest deliverable is:
- a nearly complete thesis text covering motivation, datasets, code architecture, object definitions, MC strategy, bookkeeping, and workflow status
- preliminary flagship-channel plots and, if available in time, first corrected or partially corrected flagship results
- supporting status plots and technical progress for the secondary channels

For the later note and pre-approval path:
- extend to full systematics
- incorporate additional certified data
- freeze a flagship result by late summer
- upgrade secondary channels from exploratory or support status to publication-quality only if the statistics and fit stability justify it

## 4. CMS-Style Analysis Workstreams

### 4.1 Workstream A: Scope, bookkeeping, and reproducibility

Tasks:
- freeze channel definitions, trigger menus, run eras, JSONs, global tags, and branching fractions
- create one bookkeeping table per channel
- tag the code states used for each production and result milestone
- enforce one script or command block per official plot and table

Deliverables:
- reproducible run cards
- dataset ledger
- git tags or frozen commit hashes
- one master note section for bookkeeping

Why CMS cares:
- review discussions fail quickly if the exact code, inputs, and conditions are ambiguous

### 4.2 Workstream B: Object and reconstruction validation

Tasks:
- finalize muon ID working points per channel
- finalize kaon track selection for phi reconstruction
- validate trigger matching behavior for all channels and run eras
- validate vertex quality, duplicate candidate handling, and no-triple-vertex behavior
- finalize the role of the packed-candidate study in the baseline object definition

Deliverables:
- channel-by-channel control plots for muon pt, eta, ID, and vertex probability
- dedicated note subsection explaining muon and track choices
- explicit decision on whether the packed-candidate study changes any baseline selection or only documents robustness

Why CMS cares:
- rare multi-object analyses live or die on object-definition trust

### 4.3 Workstream C: Data processing and ntuple production

Tasks:
- unify CRAB or condor submission patterns across channels
- keep `test/crabData/` as the main CRAB job-file and task-state area, with `test/crabScript_New/` used only as the config-generation and helper layer
- monitor the pinned `2022-2025` flagship CRAB campaign instead of reopening the shared EDAnalyzer baseline
- ensure each channel has a clean ntuple naming convention
- close run-era coverage gaps
- explicitly separate test outputs from production outputs

Deliverables:
- one standard directory layout
- one production checklist per run era and stream
- one `test/crabData/` campaign inventory with status summary for the flagship channel
- ntuple inventory table with event counts and integrated luminosities

### 4.4 Workstream D: MC truth definition and efficiency production

Tasks:
- define truth labels for SPS, DPS, and TPS hypotheses per channel
- standardize generator-level, miniAOD-level, and reconstructed-level category definitions
- build acceptance maps
- build reconstruction and selection efficiencies
- turn the ongoing `JpsiJpsiPhi` vertexing-versus-full-MC-matchedness study into a formal selector-validation input
- perform closure tests from truth to reconstructed categories

Deliverables:
- per-channel efficiency tables
- selector-versus-truth comparison tables for the flagship `JpsiJpsiPhi` baseline choice
- response matrices or template-transfer objects
- closure plots

Why CMS cares:
- correction factors without closure studies are not review-safe

### 4.5 Workstream E: Signal extraction and component separation

Tasks:
- define observables for separating:
  - combinatorial background
  - non-prompt or displaced backgrounds where relevant
  - SPS-like
  - DPS-like
  - TPS-like contributions
- build the first simultaneous fit model
- use `multileppat_vertex_batch` as the standard `JpsiJpsiPhi` mass-study frontend rather than notebook-local selection code
- compare `all6_same_recVtx`, `Pri_fitValid`, and related selector choices in both RooFit and `iminuit` before freezing the nominal flagship selection
- compare cut-based and fit-based stability
- document low-statistics protections and fallback models

Candidate observables:
- invariant masses
- vertex probability
- decay length or `Lxy`
- pair or system rapidity separations
- opening angles or `DeltaPhi`
- per-resonance or system pt

Deliverables:
- first fit workspace
- fit validation on toy MC or bootstrap resamples
- sideband and alternative-model studies

### 4.6 Workstream F: Cross section extraction

Tasks:
- define the fiducial phase space clearly
- compute corrected yields
- convert to fiducial cross sections
- define the conventions for `sigma_eff,DPS` and `sigma_eff,TPS`
- verify that the chosen formulas and symmetry factors are documented clearly

Deliverables:
- cross-section formula sheet
- one note subsection with all definitions
- cross-check with toy examples and units

### 4.7 Workstream G: Systematics

Required sources to consider:
- luminosity
- trigger efficiency
- muon reconstruction and identification
- track and kaon reconstruction
- vertexing and fit-quality cuts
- MC model dependence
- generator composition for SPS, DPS, TPS
- fit model choice
- mass-shape parameterization
- background parameterization
- limited MC statistics
- data period dependence
- pileup if applicable
- duplicate removal and event cleaning

Deliverables:
- systematic registry with owner, method, and status
- impact table
- nuisance implementation list for the final fit if needed

### 4.8 Workstream H: Documentation and review preparation

Tasks:
- maintain one live analysis note from now onward
- write methods in parallel with coding, not after the plots are done
- keep a standard validation appendix
- prepare weekly one-page progress snapshots

Deliverables:
- analysis note draft
- plot index
- reproducibility appendix
- pre-approval backup deck

Why CMS cares:
- pre-approval timing is determined as much by documentation maturity as by physics maturity

## 5. Role Split Recommended Right Now

### 5.1 Chi Wang should lead

- protect the pinned CMSSW analyzer baseline for `TPS-Onia2MuMu` and monitor the CRAB production campaign
- data and MC ntuple production bookkeeping
- trigger, JSON, global-tag, and run-era consistency
- packed-candidate and object-validation closure
- MC production chain consistency

### 5.2 Xing Cheng should lead

- full-MC DAG and campaign-definition documentation for `JJP` and `JUP`
- truth-level MiniAOD analyzers and topology-validation studies for `JJP` and `JUP`
- fit framework
- drawing and post-processing pipeline
- `JpsiUpsPhi` campaign assembly, ntuple-merging, and trigger-workflow documentation
- `Lxy` and selection optimization studies
- channel-level signal and background shape comparisons
- component-separation notebooks or scripts

### 5.3 Joint ownership

- final selection freeze
- efficiency parameterization
- systematic choices
- thesis milestone sign-off
- analysis note writing

### 5.4 Reserved Progress-Tracking Scope

This document should be updated directly by both analysts at least once per week.

Required editable blocks:
- the live trackers in Section `1.2`
- milestone status updates in Section `7`
- week-completion notes in Section `9`

Minimum update rule:
- each of Chi and Xing should mark at least one `Done`, one `In progress`, and one `Blocked` item whenever the tracker is updated during active analysis periods

## 6. What Should Be Done Immediately

These are the highest-priority next actions as of `2026-04-13`, after the `2026-04-01 12:00 Europe/Zurich` and `2026-04-02` reporting checkpoints were completed.

### 6.1 Record the completed report sprint and keep its decisions live

Action:
- record in this plan that the slides and reports were delivered on schedule
- carry the thesis mid-term report conclusions forward as the current scope reference

Required content:
- one dated checkpoint note for what was reported, completed, and deferred
- one explicit link between the report conclusions and the current April work plan

Reason:
- this prevents the early-April scope decisions from being re-litigated every week

### 6.2 Keep the lead and secondary scopes frozen in writing

Action:
- keep `JpsiJpsiPhi` explicitly recorded as the flagship channel
- keep `JpsiJpsiUpsilon` and `JpsiUpsilonPhi` explicitly recorded as active secondary channels with different strengths

Output:
- one written scope statement and one responsibilities table

Reason:
- without a written scope split, the thesis and summer plans will keep fighting each other

### 6.3 Protect the pinned reconstruction baseline during CRAB running

Action:
- keep the pinned `TPS-Onia2MuMu` EDAnalyzer version unchanged for the full `2022-2025` flagship CRAB campaign
- treat `test/crabData/` as the main task-state area and `test/crabScript_New/` as the generation and helper layer
- record the frozen commit hash, campaign naming, and run-era coverage in the bookkeeping table

Output:
- frozen flagship baseline reference and campaign summary for April production

Reason:
- efficiency studies and downstream validation must not be launched on a moving target

### 6.4 Build one bookkeeping table for every channel

Must include:
- datasets
- run eras
- JSONs
- global tags
- trigger paths
- MC campaign names and topology definitions
- source repository or script that defines each MC campaign
- ntuple output paths
- CRAB task directories and status summary for the flagship channel
- integrated luminosity target
- current status

Reason:
- this will immediately expose missing 2024-2025-2026 coverage and stop accidental double-counting
- it will also reveal where GitHub MC definitions and the shared analyzer baseline are already aligned, and where they are not

### 6.5 Run a three-channel validation matrix

Minimum matrix:
- `JpsiJpsiPhi` on lead MC and one representative data file
- `JpsiJpsiUpsilon` on lead MC and one representative data file
- `JpsiUpsilonPhi` on lead MC or best available proxy and one representative data file

Checks:
- branch filling
- candidate multiplicities
- trigger-match flags
- fit success rate
- runtime and output size

Reason:
- you need an apples-to-apples readiness comparison before promising one-month delivery across channels

### 6.6 Turn the MC into efficiency products, not just event samples

Action:
- define and produce the first official truth-to-reco efficiency objects:
  - acceptance
  - reconstruction efficiency
  - selection efficiency
  - trigger efficiency if not factorized separately
- use the ongoing vertex/truth study to justify the nominal flagship selector choice before treating it as final

Reason:
- this is the core of the April milestone

### 6.7 Start the first fit workspace now

Action:
- Xing leads a first compact fit for the flagship channel using the currently strongest observable set
- use `multileppat_vertex_batch` as the standard frontend for the flagship `JpsiJpsiPhi` selector comparison and fit-preparation path
- do not re-implement candidate gathering or selector logic notebook-by-notebook

Expected first target:
- background versus signal separation
- then SPS-like versus DPS-like separation
- TPS-like only if statistics and template stability justify it

Reason:
- fits always take longer than expected; waiting for "perfect inputs" will miss the milestone

### 6.8 Start thesis-writing structure immediately

Action:
- prepare the thesis chapter skeleton now, before the physics result is final

Must include:
- motivation and context
- dataset and trigger bookkeeping
- software and CMSSW workflow
- MC strategy
- efficiency plan
- current preliminary results section with clearly marked placeholders

Reason:
- the `2026-05-20` thesis deadline is driven by writing maturity, not only by analysis maturity

### 6.9 Decide what the packed-candidate study changes in the baseline

Action:
- write a one-page decision note:
  - does `MuonPackedPFCandMatch` only validate robustness?
  - or does it change the baseline muon or candidate treatment?

Reason:
- this study is valuable, but it must either become analysis policy or remain a validation appendix

### 6.10 Start the master note immediately

Action:
- keep one living note file with:
  - scope
  - dataset table
  - object definitions
  - fit strategy
  - efficiency plan
  - systematics registry

Reason:
- this is the single biggest lever for making the November pre-approval target plausible

## 7. Evaluation of the Updated Milestones

### 7.1 Milestone 1: `2026-04-01 12:00 Europe/Zurich`

Observed outcome:
- slide decks were ready for the Thursday report on schedule

Assessment:
- `Completed`

Checkpoint note:
- this milestone was handled as a status-and-direction checkpoint rather than a result freeze
- that interpretation was correct for the state of the analysis

Follow-up:
- keep the reported scope split, blockers, and April work plan visible in the live tracker sections rather than treating the slide deck as a dead artifact

### 7.2 Milestone 2: `2026-04-02`

Observed outcome:
- individual analysis reports were delivered

Assessment:
- `Completed`

Checkpoint note:
- the thesis mid-term reporting step was carried out smoothly
- the report cycle is now closed, and the remaining work is execution rather than report preparation

Follow-up:
- make sure the April bookkeeping, CRAB monitoring, validation, and fit tasks stay aligned with what was presented in the report cycle

### 7.3 Milestone 3: `2026-04-29`

Requested outcome:
- first official flagship efficiencies and first component-separation result

Assessment:
- `Green/Yellow` for one flagship channel
- `Red/Yellow` for all three channels at equal maturity

CMS-style interpretation:
- realistic if this means a first corrected flagship measurement with closure tests and a working template fit
- unrealistic if this means full three-channel precision extraction with finalized systematics and harmonized all-era inputs

Recommendation:
- define April 29 as:
  - pinned flagship EDAnalyzer baseline kept fixed while the downstream nominal selector is finalized
  - first efficiency maps
  - first component-separation fit
  - first corrected yield or fiducial cross section for `JpsiJpsiPhi`
  - progress-only status for the secondary channels

### 7.4 Milestone 4: `2026-05-20`

Requested outcome:
- bachelor thesis nearly complete in writing

Assessment:
- `Green/Yellow`

Why it is plausible:
- the thesis does not need the full final analysis result to be frozen by that date
- much of the thesis can and should be completed from now through early May:
  - motivation
  - datasets and triggers
  - software structure
  - analysis methodology
  - MC workflow
  - efficiency strategy
  - current status and preliminary plots

Why it is risky:
- it becomes risky if writing is delayed until the analysis is "done"
- it also becomes risky if all three channels are forced into equal narrative weight

Recommendation:
- by May 20, the thesis should be missing only:
  - final polished result wording
  - final systematic tables
  - possibly final conclusion emphasis

### 7.5 Milestone 5: `2026-08-31`

Requested outcome:
- flagship-analysis freeze for the first review-grade package

Assessment:
- `Yellow`

Why it is plausible:
- this gives the project the whole summer to mature the flagship channel
- it is a much healthier freeze point than trying to jump directly from a thesis deadline to a November pre-approval

Why it is risky:
- scope creep from the secondary channels can still destabilize the flagship package
- summer systematics and fit-robustness work may expand

Recommendation:
- treat late August as the real physics freeze, with autumn reserved for review hardening

### 7.6 Milestone 6: `November 2026`

Requested outcome:
- pre-approval

Assessment:
- `Yellow`

It is plausible only if all of the following happen:
- the flagship result is frozen by late August
- the note is mostly written by September
- secondary channels do not create late-stage scope creep
- review plots, closure tests, and systematics tables are already in place well before November

CMS-style warning:
- pre-approval is rarely blocked by only one missing plot; it is usually delayed by a combination of:
  - moving selection
  - unclear correction scheme
  - incomplete note
  - missing robustness checks

Recommendation:
- work backward from the pre-approval target and treat late August as the real freeze month

## 8. Key Risks and Mitigations

| Risk | Why it matters | Mitigation |
| --- | --- | --- |
| Scope creep across three channels | Three rare channels can easily become three half-finished analyses | Freeze `JpsiJpsiPhi` as flagship now |
| Moving code baseline | Efficiency and systematic studies become invalid if the code keeps shifting | Keep the pinned flagship baseline fixed and record campaign hashes and CRAB task mapping |
| Fragmented downstream workflow | Six-muon and fit studies may become irreproducible | Consolidate exploratory macros into tracked scripts |
| 2026 data dependence | Availability and certification are not fully under analyst control | Keep 2026 top-on as stretch until certified |
| Fit instability at low yields | DPS/TPS separation may be statistically fragile | Start with simplest stable model and add complexity only after closure |
| Late systematic explosion | Common in rare processes once the fit is frozen | Build the systematics registry in parallel starting now |
| Weak documentation pace | Pre-approval can slip even with good plots | Write note sections every week |

## 9. Week-by-Week Timeline

### Phase 0: Report Sprint and Scope Freeze

| Week | Dates | Main focus | Exit criterion |
| --- | --- | --- | --- |
| 01 | 2026-03-30 to 2026-04-05 | Completed report sprint: prepare slides by `2026-04-01 12:00 Europe/Zurich`, deliver reports on `2026-04-02`, freeze flagship and secondary scopes, and capture the thesis mid-term checkpoint | Slides delivered, reports delivered, written scope split exists, thesis mid-term report completed smoothly |

### Phase I: April Technical Consolidation

| Week | Dates | Main focus | Exit criterion |
| --- | --- | --- | --- |
| 02 | 2026-04-06 to 2026-04-12 | Protect the pinned `TPS-Onia2MuMu` baseline, launch and track the flagship CRAB campaign, and update trackers and bookkeeping | Pinned flagship baseline and submitted CRAB campaign are documented |
| 03 | 2026-04-13 to 2026-04-19 | Run representative MC and data validation for all three channels, with integrated `JpsiJpsiPhi` selector-comparison and vertex/truth study inputs | Validation matrix and first readiness comparison exist |
| 04 | 2026-04-20 to 2026-04-26 | Produce first official flagship MC efficiency objects, truth categories, and selector-versus-truth summaries; start the flagship fit workspace | First acceptance-times-efficiency tables, selector-validation summaries, and first fit skeleton exist |
| 05 | 2026-04-27 to 2026-05-03 | Deliver the April internal flagship milestone centered on `JpsiJpsiPhi` | First corrected flagship status package exists |

### Phase II: Thesis Write-Up Sprint

| Week | Dates | Main focus | Exit criterion |
| --- | --- | --- | --- |
| 06 | 2026-05-04 to 2026-05-10 | Write thesis chapters on motivation, datasets, triggers, software, and workflow; freeze chapter outline | Non-result thesis chapters are mostly drafted |
| 07 | 2026-05-11 to 2026-05-17 | Add MC strategy, efficiency plan, bookkeeping tables, and preliminary plots to the thesis | Thesis contains stable technical content and placeholders only where needed |
| 08 | 2026-05-18 to 2026-05-24 | Reach the `2026-05-20` thesis-near-final milestone and record the remaining result-dependent gaps | Thesis is nearly complete in writing, with only final-result polish left open |

### Phase III: Summer Flagship Maturation

| Week | Dates | Main focus | Exit criterion |
| --- | --- | --- | --- |
| 09 | 2026-05-25 to 2026-05-31 | Run flagship closure tests and compare alternative nominal selections | Closure plots and nominal-cut shortlist exist |
| 10 | 2026-06-01 to 2026-06-07 | Freeze the first-generation flagship fit observable set | Nominal fit observables are chosen |
| 11 | 2026-06-08 to 2026-06-14 | Evaluate trigger and muon-ID systematics | Object-level trigger and muon-ID prescriptions exist |
| 12 | 2026-06-15 to 2026-06-21 | Evaluate track, phi, and vertexing systematics for `4mu+2K` channels | Channel-specific object-systematics note exists |
| 13 | 2026-06-22 to 2026-06-28 | Evaluate MC-model and SPS/DPS/TPS composition systematics | Generator and composition uncertainty set exists |
| 14 | 2026-06-29 to 2026-07-05 | Process and document missing 2024-2025 ntuple coverage where needed | Updated run-era coverage table exists |
| 15 | 2026-07-06 to 2026-07-12 | Extend the secondary-channel fit framework to `JpsiJpsiUpsilon` | Secondary six-muon fit runs end to end |
| 16 | 2026-07-13 to 2026-07-19 | Consolidate the six-muon downstream chain into tracked and reproducible scripts | Clean `JpsiJpsiUpsilon` downstream workflow exists |
| 17 | 2026-07-20 to 2026-07-26 | Consolidate and document the `JpsiUpsilonPhi` campaign-assembly and trigger workflow | `JpsiUpsilonPhi` production-and-merge workflow is documented and reproducible |
| 18 | 2026-07-27 to 2026-08-02 | Decide whether each secondary channel is full-scope or support-scope for the note | Formal scope decision is recorded |
| 19 | 2026-08-03 to 2026-08-09 | Rerun the flagship nominal fit on the best available expanded inputs | Updated nominal flagship result exists |
| 20 | 2026-08-10 to 2026-08-16 | Build the combined flagship systematic summary table | Full flagship systematic impact table exists |
| 21 | 2026-08-17 to 2026-08-23 | Bring methods, efficiencies, and systematics sections of the note to review-draft quality | Note draft is coherent for the flagship channel |
| 22 | 2026-08-24 to 2026-08-30 | Final pre-freeze sanity checks on scope, fit stability, and bookkeeping | Summer freeze package is ready for the next week |
| 23 | 2026-08-31 to 2026-09-06 | Reach the late-summer flagship-analysis freeze milestone | Review-grade flagship package is frozen |

### Phase IV: Autumn Review Hardening and Pre-Approval Preparation

| Week | Dates | Main focus | Exit criterion |
| --- | --- | --- | --- |
| 24 | 2026-09-07 to 2026-09-13 | Produce the full review plotbook and closure appendix | Review plot index and appendix exist |
| 25 | 2026-09-14 to 2026-09-20 | Stress-test robustness with alternative sidebands, windows, and category definitions | Robustness summary exists |
| 26 | 2026-09-21 to 2026-09-27 | Hold the first internal rehearsal and collect feedback | First rehearsal completed |
| 27 | 2026-09-28 to 2026-10-04 | Address rehearsal feedback and tighten note wording, formulas, and definitions | Note and slides version 2 exist |
| 28 | 2026-10-05 to 2026-10-11 | Finalize backup studies and the likely committee-question material | Backup deck exists |
| 29 | 2026-10-12 to 2026-10-18 | Bring the note to near-final form and align all figure captions, tables, and definitions | Near-final note exists |
| 30 | 2026-10-19 to 2026-10-25 | Internal sign-off round and final technical corrections | Internal comments are mostly closed |
| 31 | 2026-10-26 to 2026-11-01 | Lock final figures, tables, and supporting scripts | Final figure list is frozen |
| 32 | 2026-11-02 to 2026-11-08 | Prepare committee-facing rehearsal and timing | Committee-style rehearsal held |
| 33 | 2026-11-09 to 2026-11-15 | First candidate pre-approval window | First November submission window available |
| 34 | 2026-11-16 to 2026-11-22 | Second candidate pre-approval window | Main November target window available |
| 35 | 2026-11-23 to 2026-11-29 | Buffer week for final actions, if needed | November target month is still protected |

## 10. Success Criteria by Milestone

### 10.1 By `2026-04-02`

Status:
- partially completed checkpoint

Completed:
- slides were ready by `2026-04-01 12:00 Europe/Zurich`
- reports were delivered by `2026-04-02`
- the thesis mid-term report was carried out smoothly
- the flagship versus secondary-channel scope split is recorded in this plan

Carry-over into Phase I:
- shared bookkeeping table draft

### 10.2 By `2026-04-29`

Required:
- pinned flagship EDAnalyzer baseline kept fixed and linked to the running `2022-2025` CRAB campaign
- nominal flagship downstream selector fixed
- efficiency objects available
- first closure test complete
- first fit-based component separation complete
- one corrected flagship number or limit produced

Nice to have:
- first secondary-channel fit status
- first note skeleton

### 10.3 By `2026-05-20`

Required:
- thesis text nearly complete in writing
- methods, datasets, and software workflow sections stable
- MC and efficiency strategy documented
- preliminary flagship plots included

Nice to have:
- first corrected flagship number or limit included in thesis draft
- clearer secondary-channel progress appendix

### 10.4 By `2026-08-31`

Required:
- full flagship systematic table
- updated data bookkeeping including available 2025 and any usable 2026 inputs
- stable flagship note draft
- secondary channels clearly classified as either full-scope or supporting-scope

### 10.5 By `November 2026`

Required:
- frozen flagship result
- review plotbook and closure appendix
- completed note draft
- pre-approval talk and backup slides

## 11. Final Recommendation

The collaboration should treat the early-April `report-and-scope sprint` as completed, the rest of April as a `pinned-baseline execution and validation sprint`, the next seven weeks as a `thesis-and-flagship-method sprint`, the summer as a `flagship-analysis freeze sprint`, and autumn as `review hardening`. That is the most realistic way to satisfy the fixed writing deadlines while still protecting a credible November 2026 pre-approval target.

The current evidence already supports optimism:
- reconstruction code is real and active
- the flagship `TPS-Onia2MuMu` baseline is already pinned and submitted to CRAB for full `2022-2025` processing
- data and MC production already exist
- fit-oriented work has started
- `JpsiJpsiPhi` downstream mass studies now have a package-level workflow through `multileppat_vertex_batch`
- Xing's visible role now includes explicit full-chain JJP/JUP MC production and MiniAOD truth-analysis repositories, not only downstream fitting or merge scripts
- the Chi/Xing skill split is complementary

The main danger is not lack of work. It is lack of scope control.

If scope is controlled now, this can become a credible CMS-style analysis program. If scope is not controlled now, the project risks reaching May with a stressed thesis schedule and reaching autumn with many interesting intermediate products but no review-safe flagship result.

## 12. Evidence Used For This Snapshot

Local paths inspected:
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/TPS-Onia2MuMu/`
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/MuonPackedPFCandMatch/`
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/TPS-Onia2MuMu/test/crabData/`
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/TPS-Onia2MuMu/test/crabScript_New/`
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/TPS-Onia2MuMu/doc/multileppat_vertex_batch_integration_plan.md`
- `CMSSW_15_0_15_JpsiJpsiPhi_refactor/src/HeavyFlavorAnalysis/TPS-Onia2MuMu/doc/study_multileppat_vertexing_full_mc_match_plan.md`
- `JpsiJpsiPhi/rootNtuple/`
- `JpsiJpsiPhi/MC_samples/`
- `JpsiJpsiUps/MC_samples/`
- `JpsiJpsiUps/ParticleCand/`
- `JpsiJpsiUps/condor_job/ntuples-toHTCondor/`
- `JpsiUpsPhi/rootNtuple/`
- `JpsiUpsPhi/rootNtuple/hadd_2022.sh`
- `JpsiUpsPhi/rootNtuple/hadd_2023.sh`
- `JpsiUpsPhi/rootNtuple/hadd_2024.sh`
- `JpsiUpsPhi/rootNtuple/gen_hadd_2022.sh`
- `JpsiUpsPhi/rootNtuple/gen_hadd_2023.sh`
- `JpsiUpsPhi/rootNtuple/gen_hadd_2024.sh`
- `JpsiUpsPhi/rootNtuple/viewTrigger.C`
- `JpsiUpsPhi/rootNtuple/viewTrigger.h`

GitHub repositories inspected:
- `Eric100911/TPS-Onia2MuMu`
- `Eric100911/MuonPackedPFCandMatch`
- `Eric100911/ParticleCand`
- `Eric100911/ntuples-toHTCondor`
- `Eric100911/JpsiJpsiPhi_MiniAOD_GenAnalyzer`
- `Endymion2288/TPS-Onia2MuMu`
- `Endymion2288/JpsiJpsiPhi-workspace`
- `Endymion2288/JpsiJpsiUps-workspace`
- `Endymion2288/JpsiUpsPhi-workspace`
- `Endymion2288/Full_MC_Production`
- `Endymion2288/JUPMCAnalyzer`
- `Endymion2288/JJP_MINIAOD_Analysis`

Selected git-history evidence cross-checked locally:
- merge commit of `Endymion2288/Dev-J-J-P` into `TPS-Onia2MuMu`
- `Endymion` commits titled `feat: Add the fromPV and the pvQuality`
- `Endymion` commits titled `feat: new test program`
- `Endymion` commits titled `fix: repair git object`

Selected GitHub file evidence inspected:
- `Endymion2288/Full_MC_Production: README.md`
- `Endymion2288/Full_MC_Production: dag_generator.py`
- `Endymion2288/Full_MC_Production: processing/run_chain.sh`
- `Endymion2288/Full_MC_Production: common/cmssw_configs/ntuple_jjp_cfg.py`
- `Endymion2288/Full_MC_Production: common/cmssw_configs/ntuple_jup_cfg.py`
- `Endymion2288/JUPMCAnalyzer: README.md`
- `Endymion2288/JUPMCAnalyzer: run_all_modes.sh`
- `Endymion2288/JJP_MINIAOD_Analysis: README.md`
- `Endymion2288/JJP_MINIAOD_Analysis: run_all_modes.sh`
- `Endymion2288/TPS-Onia2MuMu: README.md`

Selected GitHub commit evidence inspected:
- `Endymion2288/Full_MC_Production`: `feat: first edition of LHE2MINIAOD`
- `Endymion2288/Full_MC_Production`: `fix: successfully test`
- `Endymion2288/JUPMCAnalyzer`: `feat: first edition for analyze the JUP MC`
- `Endymion2288/JUPMCAnalyzer`: `feat: use HTCondor`
- `Endymion2288/JJP_MINIAOD_Analysis`: `feat: change eta to y and use multiprocessing`
- `Endymion2288/JJP_MINIAOD_Analysis`: `feat: strictly limit the parton process`
- `Endymion2288/JJP_MINIAOD_Analysis`: `fix: do not limit parton relation in TPS`
- `Endymion2288/JJP_MINIAOD_Analysis`: `feat: use HTCondor`
