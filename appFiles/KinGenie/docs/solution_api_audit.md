<!-- AI-generated audit and reproduction instructions. -->
# Solution simulation and fitting API audit

Checked on 2026-11-10 against pyKinGenie **1.0.19**, corresponding to upstream
[`40622da`](https://github.com/osvalB/pykingenie/tree/40622da1afd2343c88c17156cae06750df5e899f).
The application dependency is pinned to this release.

## Equations and parameter mappings

The panel's one-site, induced-fit and conformational-selection signals agree with
independently integrated mass-action equations. Slow and fast equilibrium
relaxation rates agree with the eigenvalues of independently constructed
Jacobians for both two-step models during the audit. These mathematical checks
are not included in the app test suite; equation correctness belongs to the
Python package's tests.

| Input | Expected meaning in pyKinGenie |
| --- | --- |
| Protein and ligand concentrations | Total concentrations after mixing, in μM; free concentrations evolve by mass conservation. |
| One-site `Kd`, `koff` | μM and s⁻¹; the solver uses `kon = koff / Kd`. |
| `kon`, `koff` | Binding and dissociation rates, μM⁻¹ s⁻¹ and s⁻¹. |
| Induced-fit `kc`, `krev` | Forward intermediate → induced complex and reverse induced → intermediate rates, s⁻¹. |
| Conformational-selection `kc`, `krev` | Inactive → active and active → inactive protein rates, s⁻¹. |
| Species signal inputs | Signal coefficients in A.U./μM, multiplied by species concentrations; these are not complete signal amplitudes. |
| Initial state | No complex; conformational selection starts with the free protein conformations at equilibrium. |
| Conformational-selection free-protein signal | The app assigns the same coefficient to both free protein conformations. |
| Fit signal checkboxes | Checked coefficients are fitted globally within a sample; unchecked species have coefficient zero. Fitted coefficients have nonnegative bounds. |
| Unequal induced-fit complex coefficients | `ESint_equals_ES = FALSE`; both complex coefficients are fitted when the bound-complex checkbox is selected. |
| One-site fitting `t0` | The model evaluates `t + t0`; positive `t0` advances the reaction. The induced-fit app dialog does not offer this option. |

The one-site helper's generic arguments `a_total` and `signal_a` consistently
receive protein concentration and protein signal coefficient in both simulation
and fitting. Its symmetric mass-action equation does not require swapping them.

The relaxation helper API uses `ke = kc`, `kr = krev` for **both** models.
This differs from the notation of the induced-fit mechanism in
[Paul and Weikl (2016)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005067),
where paper `kr` is forward and paper `ke` is reverse. The example defaults now
use `kc = 10`, `krev = 1` to match the rates in Figure 3.

The app currently offers conformational selection for simulation, but does not
expose pyKinGenie's conformational-selection fitting API. Adding that fit model
would be a separate feature.

## Application corrections

- Save both calculated relaxation-rate arrays in the simulation results so the
  relaxation plot receives its data.
- Use named solution-solver and relaxation-helper arguments to make rate and
  concentration mappings explicit.
- Correct the induced-fit example's forward/reverse defaults.
- Label total concentrations and species coefficient units, and explain initial
  states and fitting checkbox behavior.
- Reject invalid solution dilution counts, zero time steps/dilution factors,
  nonfinite values, and nonpositive parameters used for reversible equilibrium
  models before building grids. Zero-concentration simulation controls remain
  supported.
- Clear the previous exponential model's relaxation-plot state after completing
  either binding-model fit. Reject fitting with no signal-producing species or
  nonpositive/nonfinite total concentrations in the actual fitting objects.

## Remaining upstream optimizer finding

The fitting equations and app keyword mappings agree, but successful return from
the optimizer does not guarantee a good fit. With six noiseless traces, uniformly
sampled from 0 to 8 seconds (61 points), fitting all species coefficients can
converge to a poor solution even though an exact solution exists by construction.
This reproduces by calling the Python fitting classes directly, without Shiny.

The direct Python investigation used pyKinGenie 1.0.19 and compared predictions
with the noiseless input data after fitting. No Python-only test or diagnostic
script is included in this repository. The investigated data and settings were:

- Protein concentrations: 0.4 and 0.8 μM, each paired with ligands 0.3, 0.9, 2.7 μM.
- Free protein, free ligand and complex coefficients: 0.2, 0.4, 1.7 A.U./μM.
- One-site: `Kd = 0.4 μM`, `koff = 0.3 s⁻¹`; fixed zero time offset.
- Induced fit: `kon = 1.3 μM⁻¹ s⁻¹`, `koff = 0.3 s⁻¹`, `kc = 0.7 s⁻¹`,
  `krev = 1.1 s⁻¹`; intermediate coefficient 1.7, induced coefficient 2.3;
  fixed zero time offset and unequal complex coefficients.

One direct Python run returned approximately **0.50 A.U.** maximum error for
one-site and **0.092 A.U.** for induced fit, with an `ODEintWarning` during the
induced-fit search. Numerical optimizer paths and exact errors can vary by
runtime. This is a fitting robustness finding, not evidence that the differential
equations are wrong. No pyKinGenie source change is included in this app branch.

The app integration tests cover actual Shiny callbacks, input validation,
coefficient selections, result dimensions, result state and simulation plot data.
Simulation outputs are compared with the package API to check the app's parameter
forwarding. They do not reimplement ODEs, test Jacobians or assert optimizer
accuracy; those belong in the Python package.
