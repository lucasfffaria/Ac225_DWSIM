# Actinium-225 purification scheduling with DWSIM: supporting material

Companion repository for the article:

> Lucas Faria da Silva. *Scheduling ²²⁵Ac purification under isotope and delivery uncertainty using DWSIM process modeling.* (The full reference will be added upon publication.)

The repository contains the Supporting Information and the files needed to inspect, rerun and authenticate the calculations reported in the article: both native DWSIM workspaces, the scripts embedded in their custom operations, the CPython model and analyses, input-provenance files, numerical results and verification reports.

## Repository contents

| Path | Contents |
|---|---|
| `Supporting_Information_Ac225_DWSIM.docx` | Supporting Information (Sections S1–S10) |
| `dwsim/` | The two native DWSIM 10.2.4.0 workspaces |
| `code/` | CPython model and analyses, DWSIM automation runners and the custom-operation scripts |
| `data/` | Literature inputs with their page or table locations, and explicitly labelled assumptions |
| `results/` | Numerical outputs, verification reports and reproduction logs |
| `docs/DWSIM_GUIDE.md` | How to open and explore the two DWSIM workspaces |
| `requirements.txt` | Pinned Python packages |
| `reproduce.ps1` | One-command reproduction on Windows |
| `CHECKSUMS.sha256` | SHA-256 checksum of every file in the repository |

### DWSIM workspaces (`dwsim/`)

| File | Description | SI section |
|---|---|---|
| `Ac225_Isotope_Timing.dwxml` | Baseline workspace with six independent custom operations: proton and nitrogen isotope inventories (parent-free and ²²⁵Ra scenario), the robust upper input corner and the documented liquid ledger. It has no material streams and no property package. | S6 |
| `Ac225_Connected_Purification.dwxml` | Connected DGA2–LN purification subflowsheet with seven custom operations, sixteen material streams and twenty-one connections, an embedded aqueous density correlation, operation-level isotope clocks and an independent citrate/AG50 pH evidence node. | S7–S10 |

### Code (`code/`)

| File | Role | SI section |
|---|---|---|
| `model.py` | Analytical inventory model: purity, feasibility windows, robust start window, ²²⁵Ra retention and matrix form | S1–S2 |
| `run_analysis.py` | Isotope timing, Monte Carlo propagation, robust scheduling, recovery check, liquid ledger, published benchmark and numerical verification | S1–S5 |
| `liquid_properties.py` | Laliberté (2009) density and viscosity, solution preparation, acid conditioning and Ergun pressure drop | S8–S9 |
| `process_model.py` | Connected-process recipe, stage inventories and hydraulic reference cases | S7, S9 |
| `run_process_analysis.py` | Independent density check, component balance, hydraulic screen, link between processing time and the robust schedule, pH-recovery sensitivity | S7–S10 |
| `run_dwsim.py` | Builds, solves, saves, reloads and recalculates the baseline workspace and compares it with CPython | S6 |
| `run_connected_dwsim.py` | Builds, solves, saves and reloads the connected subflowsheet, runs the native flow/queue and pH sweeps and checks every stage against a matrix-exponential calculation | S7, S9, S10 |
| `dwsim_inventory_unit.py`, `dwsim_reagent_unit.py`, `dwsim_process_unit.py` | IronPython-compatible custom-operation scripts, also embedded in the `.dwxml` files | S6–S7 |

### Input data (`data/`)

| File | Contents | Source details |
|---|---|---|
| `source_inputs.json` | Half-lives, isotope activity ratios with uncertainties (Tan et al. 2025, Table 4), recovery observations (Tan et al. 2025, Table 2) and Monte Carlo settings | Article Table 1; SI S4 and Table S3 |
| `reagent_recipe.csv` | Documented fresh liquid additions (Tan et al. 2025, Sections 2.3–2.4) | SI S5 and Table S4 |
| `process_input_provenance.csv` | Every process input labelled as a reported value, correlation input, endpoint reconciliation or scenario assumption | SI S7–S10 |
| `liquid_property_parameters.json` | Laliberté (2009) coefficients and validity limits, transcribed from thermo 0.6.1 | SI S8 |
| `nitric_density_check.csv` | 37 density differences at 298.15 K (Hnedkovsky et al. 2025, SI Table S1) used for the independent density check | SI S8 |
| `pH_recovery_literature.csv` | Citrate/AG50 actinium recovery versus load pH (Fitzsimmons et al. 2019a) | SI S10 |

### Results (`results/`)

| File | Contents | SI section |
|---|---|---|
| `verification.json` | Numerical verification: closed form versus matrix exponential and ODE integration, ledger closure, semigroup identity, recovery cancellation and robust corner check | S3.1, Table S1 |
| `dwsim_verification.json` | Baseline native run: DWSIM version, object count, solver errors, save/reload check, comparison with CPython and SHA-256 of the workspace | S6 |
| `connected_dwsim_verification.json` | Connected native run: stream and operation counts, solver outcome, reload check, matrix-oracle comparison, global mass residual and SHA-256 of the workspace | S7 |
| `summary.json` | Isotope-timing, scheduling, recovery, ledger and benchmark results | S1–S5 |
| `process_summary.json` | Engineering values reported in the article | S7–S10 |
| `connected_stream_table.csv`, `connected_unit_outputs.json`, `connected_edges.json` | Stream inventories (complete Table S6), custom-operation outputs and native connections | S7 |
| `component_balance.csv`, `process_stage_inventory.csv` | Reaction-adjusted component balance and stage isotope inventories | S7 |
| `dwsim_process_sweep.csv`, `dwsim_pH_sweep.csv` | Native flow/queue sweep (25 cases) and pH sweep (101 cases) | S9–S10 |
| `hydraulic_flow_screen.csv`, `process_robust_schedule.csv` | CPython hydraulic scan and link between processing time and the robust schedule | S9 |
| `density_independent_check.csv` | Independent nitric-acid density comparison | S8 |
| `pH_control_sensitivity.csv`, `recovery_factorization_check.csv`, `unreported_acid_carryover_sensitivity.csv` | Recovery-response, recovery-allocation and acid-carryover sensitivities | S7, S10 |
| `time_curves.csv`, `parent_scenarios.csv`, `robust_schedule_curves.csv`, `robust_schedule_sensitivity.csv`, `correlation_sensitivity.csv`, `time_savings.csv`, `source_recovery_check.csv`, `stage_time_ledger.csv`, `mc_isotope_ratio_draws.npz` | Isotope-timing, scheduling, recovery and ledger outputs | S1–S5 |
| `reproduction_logs/` | Console output of the four runs that produced the files above | – |

## Software

Verified environment: 64-bit Windows 11, CPython 3.11.9 and DWSIM 10.2.4.0 (Automation3 interface), with the package versions pinned in `requirements.txt`: NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.5, thermo 0.6.1, chemicals 1.5.2, fluids 1.3.1 and, on Windows only, pythonnet 3.1.0. DWSIM and pythonnet are needed only for the two native runners.

## Reproducing the results

**Windows (complete reproduction, including DWSIM).** From the repository root in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\reproduce.ps1 -Python '.\.venv\Scripts\python.exe'
```

If script execution is blocked, use `powershell -ExecutionPolicy Bypass -File .\reproduce.ps1 -Python '.\.venv\Scripts\python.exe'`.

The script runs `run_analysis.py`, `run_process_analysis.py`, `run_dwsim.py` and `run_connected_dwsim.py` and stops at the first error. Console output is written to `results/reproduction_logs/`. DWSIM is located through the `DWSIM_BIN` environment variable or the `-DwsimBin` parameter; the default is `%LOCALAPPDATA%\DWSIM`. Use `-SkipDwsim` to rerun only the platform-independent analyses while keeping the native results included in this repository. The automation runs in a separate process and does not interfere with an open DWSIM window.

**Any operating system (platform-independent analyses only):**

```bash
python -m venv .venv
. .venv/bin/activate        # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
python code/run_analysis.py
python code/run_process_analysis.py
```

Rerunning the DWSIM runners rewrites both `.dwxml` files. DWSIM stores the save time inside each file, so regenerated workspaces have new SHA-256 values; the hashes recorded in the verification reports authenticate the files in this repository.

## Verifying the files

`CHECKSUMS.sha256` lists the SHA-256 checksum of every file in the repository. To check them on Linux, macOS or Git Bash:

```bash
sha256sum -c CHECKSUMS.sha256
```

In PowerShell (no output means that every file matches):

```powershell
foreach ($line in Get-Content .\CHECKSUMS.sha256) {
    $hash, $file = $line -split ' [ *]', 2
    if ((Get-FileHash -Algorithm SHA256 -LiteralPath $file).Hash.ToLower() -ne $hash) { "MISMATCH  $file" }
}
```

The two DWSIM workspaces are also authenticated by the hashes written by the runs that created them: `dwxml_sha256` in `results/dwsim_verification.json` and `sha256` in `results/connected_dwsim_verification.json`.

Reference values for a rerun (`results/summary.json`, proton case): 99% actinium-only delivery interval of 15.07–28.31 d, robust start window of 12.65–21.19 d and guaranteed normalized delivered activity of 0.1760607724. The numerical-verification values of Table S1 are in `results/verification.json`.

## Scope and limitations

- The baseline workspace evaluates isotope inventories and liquid accounting only; it has no material streams and no thermodynamic property package.
- The connected workspace starts at the primary actinium fraction and covers the downstream DGA2–LN sequence. The citrate/AG50 pH node is an independent evidence example and is not applied to the DGA2–LN recoveries.
- The calculations do not predict chromatographic selectivity, complete-facility waste, equipment pressure ratings or clinical release acceptability.
- The 99% threshold refers to actinium-only activity purity, not to a radionuclidic or radiochemical purity specification.
- No new radiochemical measurements were made. Literature values are transcribed with their source locations; the original publications are not redistributed here.

Detailed assumptions and provenance are given in the article, in the Supporting Information and in `data/process_input_provenance.csv`.
