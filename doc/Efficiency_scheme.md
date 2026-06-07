### Efficiency and acceptance scheme

| Step                                        | Symbol                                                  | Level                    | Map axes                                                 | Short definition                                             |
| ------------------------------------------- | ------------------------------------------------------- | ------------------------ | -------------------------------------------------------- | ------------------------------------------------------------ |
| Acceptance                                  | $A_{J/\psi}(\vec{x}_{J_i})$, $A_{\phi}(\vec{x}_{\phi})$ | per $J/\psi$ or $\phi$   | $(p_T^{\mathrm{meson}}, |y^{\mathrm{meson}}|)$           | decay daughters in fiducial region given meson in acceptance |
| muonRECOinJpsi                              | $\epsilon_{\mu\mathrm{Reco}}^{J/\psi}$                  | per $J/\psi$             | $(p_T^{J/\psi}, y^{J/\psi})$                             | both daughter muons reconstructed and matched                |
| kaonRECOinPhi                               | $\epsilon_{K\mathrm{Reco}}^{\phi}$                      | per $\phi$               | $(p_T^{\phi}, y^{\phi})$                                 | both kaons reconstructed and matched                         |
| muonID                                      | $\epsilon_{\mu\mathrm{ID}}^{J/\psi}$                    | per $J/\psi$             | $(p_T^{J/\psi}, y^{J/\psi})$                             | both matched muons pass ID                                   |
| kaonID                                      | $\epsilon_{K\mathrm{ID}}^{\phi}$                        | per $\phi$               | $(p_T^{\phi}, y^{\phi})$                                 | both matched kaons pass track/kaon selection                 |
| dimuon                                      | $\epsilon_{\mu\mu}^{J/\psi}$                            | per $J/\psi$             | $(p_T^{J/\psi}, y^{J/\psi})$                             | valid dimuon candidate                                       |
| dikaon                                      | $\epsilon_{KK}^{\phi}$                                  | per $\phi$               | $(p_T^{\phi}, y^{\phi})$                                 | valid $K^+K^-$ candidate                                     |
| HLT                                         | $\epsilon_{\mathrm{HLT}}$                               | event                    | nominally event-level                                    | trigger OR and trigger-object matching                       |
| fourMuonVertexing                           | $\epsilon_{4\mu\mathrm{vtx}}$                           | per $J/\psi J/\psi$ pair | pair-level axes                                          | valid four-muon vertex                                       |
| triOniaVertexingAndOtherEventLevelSelection | $\epsilon_{\mathrm{triOnia}}$                           | event                    | $(p_T^{J/\psi_1},p_T^{J/\psi_2})$, split by $p_T^{\phi}$ | three-meson vertex and remaining event-level cuts            |

### Current implementation notes

The parquet efficiency workflow stores object-level reconstruction and ID steps
as conditional chains, then switches to event-level flags. The nominal event
chain is:

```text
s_cand -> hlt_event -> four_muon_vtx
```

The primary-vertex diagnostics are parallel branches with denominator
`four_muon_vtx`, not a single cumulative final step:

```text
Pri_fitValid
Pri_fitPass
Pri_assocPVPass
Pri_trackPVPass
```

For corrected-yield studies the default correction is factorized rather than a
single cumulative `correlated_3d` lookup. It multiplies per-object acceptance and
conditional efficiency maps in `(pT, |y|)` by event-level maps for HLT,
four-muon vertexing, and the final tri-onia/PV endpoint. The default endpoint is
`Pri_assocPVPass / four_muon_vtx`.

Factorized map lookup uses a tiered fallback:

```text
fine bin -> coarse bin -> inclusive bin
```

The nominal central value uses observed bin efficiencies after fallback
selection. Jeffreys symmetric binomial uncertainties are stored for each bin and
used for MC-stat diagnostics. The default thresholds are `N_min_fine = 30` and
`N_min_coarse = 50`.

The stacked J/psi derived plots are diagnostics for the two J/psi objects
combined into one per-object map. Their parquet products use direct
`(pT, |y|)` binning (`y_axis == "abs_y"`):

```text
stacked_jpsi_acceptance_maps.parquet: fiducial acceptance
stacked_jpsi_efficiency_maps.parquet: muonRECO, muonID, dimuon
```

When only these plots are needed, use `build_derived_efficiency.py
--plot-scope stacked-jpsi` to skip cumulative, conditional, pair-level, and
systematics plots.

The correction workflow is:

```bash
python3 rebuild_efficiency_maps.py \
  --input-dir <merged_efficiency_dir> \
  --output-dir <fresh_rebuilt_dir> \
  --samples JJP_DPS1 JJP_DPS2_CS JJP_DPS2_G JJP_SPS_CS JJP_SPS_G

python3 build_derived_efficiency.py --input-dir <fresh_rebuilt_dir>

python3 -m efficiency_workflow.build_factorized_maps \
  --input-dir <fresh_rebuilt_dir> \
  --samples JJP_DPS1 JJP_DPS2_CS JJP_DPS2_G JJP_SPS_CS JJP_SPS_G

python3 compute_efficiency_corrected_yield.py \
  --data-input <jjp_data_selected.root> \
  --efficiency-dir <fresh_rebuilt_dir> \
  --plot-dir <yield_plot_dir> \
  --corrected-root <jjp_data_effcorr_selected.root> \
  -o <yield_summary.json>

python3 fit_splot.py \
  --channel JJP --dataset data \
  -i <jjp_data_effcorr_selected.root> \
  -o <jjp_data_effcorr_splot.root> \
  --effcorr-weight-branch effcorr_weight

python3 plot_weighted_distributions.py \
  --channel JJP --dataset data \
  -i <jjp_data_effcorr_splot.root> \
  -o <effcorr_splot_dynamics> \
  -w signal_effcorr_sw
```

### Definitions

Let

$$
J_i \equiv J/\psi_i,\qquad i=1,2,
$$

and

$$
\vec{x}*{J_i}=(p_T^{J_i}, |y^{J_i}|),\qquad
\vec{x}*{\phi}=(p_T^\phi, |y^\phi|).
$$

The general acceptance is defined as the fraction of generated events whose decay products are within the fiducial region, which is defined as:

$$
  A_{\mathrm{tot}} = \frac{
    N_{\mathrm{gen}}^{\mathrm{fid}}
  }{
    N_{\mathrm{gen}}(J_1, J_2, \phi\in\Omega_{\mathrm{meson}})
  }
$$

Equivalently, one write:

$$
A_{\mathrm{tot}} = \frac{
    N_{\mathrm{gen}}
      \left( 
        \mu^\pm\in\Omega_\mu,
        K^\pm\in\Omega_K
        \middle | 
        J_1, J_2, \phi\in\Omega_{\mathrm{meson}}
      \right)
  }{
    N_{\mathrm{gen}}
      \left(
        J_1, J_2, \phi\in\Omega_{\mathrm{meson}}
      \right)
  }
$$

A cross check is to factorize the acceptance by the $J/\psi$ and $\phi$ meson level:

$$
A_{J/\psi}(\vec{x}_{J_i})
=
\frac{
N_{\mathrm{gen}}(\mu_1\in\Omega_\mu, \mu_2\in\Omega_\mu | J_i\in\Omega_{\mathrm{meson}})
}{
N_{\mathrm{gen}}(J_i\in\Omega_{\mathrm{meson}})
}
$$

$$
A_{\phi}(\vec{x}_{\phi})
=
\frac{
N_{\mathrm{gen}}(K_1\in\Omega_K, K_2\in\Omega_K | \phi\in\Omega_{\mathrm{meson}})
}{
N_{\mathrm{gen}}(\phi\in\Omega_{\mathrm{meson}})
}
$$

And we approximate with:

$$
A_{\mathrm{tot}} (\vec{x}_{J_1}, \vec{x}_{J_2}, \vec{x}_{\phi}) \approx A_{J/\psi}(\vec{x}_{J_1}) A_{J/\psi}(\vec{x}_{J_2}) A_{\phi}(\vec{x}_{\phi}).
$$

This approximation is exact only if the acceptance of each meson is independent of the kinematics of the other mesons. Since unpolarized nature of the $J/\psi$'s are assumed, one may expect that the acceptance of each $J/\psi$ is not strongly dependent on the kinematics of the other $J/\psi$ or the $\phi$.

We further define the efficiency of each step as the fraction of events passing the step selection among the events passing the previous step selection, with the same kinematic binning as the acceptance.

For each $J/\psi$,

$$
\epsilon_{\mu\mathrm{Reco}|J_i}(\vec{x}_{J_i})
=

\frac{
N(S_{\mathrm{fid}}\cap\mu_1^{\mathrm{reco}}\cap\mu_2^{\mathrm{reco}})
}{
N(S_{\mathrm{fid}})
},
$$

$$
\epsilon_{\mu\mathrm{ID}|J_i}(\vec{x}_{J_i})
=

\frac{
N(S_{\mu\mathrm{Reco}}\cap\mu_1^{\mathrm{ID}}\cap\mu_2^{\mathrm{ID}})
}{
N(S_{\mu\mathrm{Reco}})
},
$$

$$
\epsilon_{\mu\mu|J_i}(\vec{x}_{J_i})
=

\frac{
N(S_{\mu\mathrm{ID}}\cap J_i^{\mathrm{reco}})
}{
N(S_{\mu\mathrm{ID}})
}.
$$

For the (\phi),

$$
\epsilon_{K\mathrm{Reco}|\phi}(\vec{x}_{\phi})
=

\frac{
N(S_{\mathrm{fid}}\cap K_1^{\mathrm{reco}}\cap K_2^{\mathrm{reco}})
}{
N(S_{\mathrm{fid}})
},
$$

$$
\epsilon_{K\mathrm{ID}|\phi}(\vec{x}_{\phi})
=

\frac{
N(S_{K\mathrm{Reco}}\cap K_1^{\mathrm{ID}}\cap K_2^{\mathrm{ID}})
}{
N(S_{K\mathrm{Reco}})
},
$$

$$
\epsilon_{KK|\phi}(\vec{x}_{\phi})
=

\frac{
N(S_{K\mathrm{ID}}\cap \phi^{\mathrm{reco}})
}{
N(S_{K\mathrm{ID}})
}.
$$

The trigger efficiency is

$$
\epsilon_{\mathrm{HLT}}
=

\frac{
N(S_{\mathrm{cand}}\cap \mathrm{HLT})
}{
N(S_{\mathrm{cand}})
}.
$$

The four-muon vertexing efficiency is defined at the (J/\psi J/\psi)-pair level:

$$
\epsilon_{4\mu\mathrm{vtx}}(\vec{x}_{JJ})
=

\frac{
N(S_{\mathrm{HLT}}\cap V_{4\mu})
}{
N(S_{\mathrm{HLT}})
}.
$$

The final event-level efficiency is

$$
\epsilon_{\mathrm{triOnia}}^{(k)}(p_T^{J_1},p_T^{J_2})
=

\frac{
N(S_{4\mu\mathrm{vtx}}\cap V_{J/\psi J/\psi\phi}\cap S_{\mathrm{other}};,p_T^\phi\in I_k)
}{
N(S_{4\mu\mathrm{vtx}};,p_T^\phi\in I_k)
}.
$$

Here (I_k) is the (k)-th (\phi)-(p_T) interval used for plotting.

The total fiducial efficiency is

$$
\begin{aligned}
\epsilon_{\mathrm{fid}}
={}&
\prod_{i=1}^{2}
\epsilon_{\mu\mathrm{Reco}|J_i}
\epsilon_{\mu\mathrm{ID}|J_i}
\epsilon_{\mu\mu|J_i}
\
&\times
\epsilon_{K\mathrm{Reco}|\phi}
\epsilon_{K\mathrm{ID}|\phi}
\epsilon_{KK|\phi}
\
&\times
\epsilon_{\mathrm{HLT}}
\epsilon_{4\mu\mathrm{vtx}}
\epsilon_{\mathrm{triOnia}} .
\end{aligned}
$$

And the full correction factor is

$$
A\epsilon
=

A_{\mathrm{tot}}\epsilon_{\mathrm{fid}}.
$$

### 
