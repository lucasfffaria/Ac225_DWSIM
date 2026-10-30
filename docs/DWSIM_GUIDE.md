# Using the two DWSIM workspaces

Both workspaces were executed with DWSIM 10.2.4.0. The automation scripts run in a separate Python process, so an open DWSIM window does not need to be closed.

## Connected process flowsheet

Open `dwsim/Ac225_Connected_Purification.dwxml` in DWSIM and run the flowsheet calculation. The custom-operation scripts and the density correlation are embedded in the file.

| Operation | Function and main inputs |
|---|---|
| P-100 PRIMARY BOUNDARY | Primary fraction: 80 mL of 10 M HNO₃; `cooling_d`, `r226`, `r227`, `eta_primary` |
| M-101 DILUTION | Mixing with the DI-WATER stream; `duration_min` |
| C-201 DGA2 | Feed and reagent balance, product cut and waste; `flow_mL_min`, `eta_Ac`, `carry_HCl_M`, column geometry |
| R-201 ACID ADJUSTMENT | Reaction of the carried-over acidity with the NH3-ADDITION stream; `duration_min` |
| C-202 LN | Balance and collection of 10 mL of 0.7 M HNO₃; `flow_mL_min`, `eta_Ac`, column geometry |
| T-301 PROTOCOL AND QUEUE | Remaining protocol time and additional queue; `protocol_allowance_min`, `extra_queue_h` |
| EVIDENCE: CITRATE AG50 pH | Independent study: `load_pH`, `pH_halfwidth`; the whole interval must lie within pH 1.5–2.5 |

Each stream stores the isotope activities as extra properties (`A225_Bq`, `A226_Bq`, `A227_Bq`, `Ra225_Bq`) together with `age_d`, `batch_volume_mL` and `batch_basis_s`. The operation results report durations, balance residuals, activity concentrations and pressure-drop estimates. The script interpreter may also display intermediate numerical variables; use the descriptive names given in the Supporting Information and in `results/connected_unit_outputs.json`.

The mass flow shown for a stream uses an accounting basis of one batch per 3600 s: multiply kg/s by 3600 to obtain kg per batch. The actual column flow is the `flow_mL_min` parameter; changing only the basis of a stream is not equivalent to changing the pump. The inputs represent sequential additions aggregated per batch.

To explore the documented scenarios, change `flow_mL_min` in both columns and `extra_queue_h` in T-301, then recalculate. Recoveries are held constant in these cases; they do not predict recovery as a function of flow. The pressure-drop comparison also depends on `particle_um` and `voidage`. The 50 kPa criterion used in the article is illustrative, not an equipment specification.

The acidity of the DGA cut (`carry_HCl_M`) and the ammonia feed must be changed consistently. R-201 computes the balance with the amount of ammonia actually fed; it is not an automatic pH controller. The Python function `process_model.recipe(carry_M=...)` determines the 1 M NH₃ addition for the 0.05 M free-acid target and was used for the cases in `results/unreported_acid_carryover_sensitivity.csv`. The article assumes that HCl is the only acid carried over in this cut; other compositions require revising this submodel, the reagents and the reaction balance.

The property package of the conventional streams is Raoult's Law with an imposed liquid phase, and the density is replaced by the embedded aqueous correlation. Only density, composition, mass and explicitly modelled volume are used. The remaining standard thermodynamic fields, including energy and phase equilibrium, do not support any result. The isotopes are separate metadata without assigned fluid properties.

## Baseline inventory workspace

Open `dwsim/Ac225_Isotope_Timing.dwxml` to explore the two irradiation cases, ²²⁵Ra retention, the robust scenario and the documented reagent accounting. It contains six independent custom operations. Outputs prefixed with `parent_free_` always refer to the case without parent feeding, even when the operation also reports a delivered activity with Ra present.

## Complete reproduction

From the repository root in PowerShell:

```powershell
.\reproduce.ps1 -Python '.\.venv\Scripts\python.exe'
```

The command reruns both analyses, both DWSIM workspaces, the 25 flow/queue cases and the 101 evaluations of the pH model. Logs are written to `results/reproduction_logs/`. See the main `README.md` for the environment setup.

The files in `data/` record data provenance and assumptions. The analysis scripts contain the article-specific scenarios, so editing an isolated JSON file does not change all of them automatically.
