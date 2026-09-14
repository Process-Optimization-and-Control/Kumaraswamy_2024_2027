# FOCAPO - CPC Conference 2027

Code, results and figures for the FOCAPO - CPC 2027 submission.

This folder runs a closed-loop simulation of a Nonlinear Model Predictive Controller (NMPC)
regulating a $CO_2$ pipeline-injection network. The manipulated variables are the four wellhead
choke openings and the mainline pump speed; the controlled variables are the four well flow rates.
The NMPC uses a first-principles dynamic model of the mass and momentum balances, while the
"plant" is a dynamic multiphase flow simulator, LedaFlow, driven from Julia through JavaScript at each sample. Simplified state estimation and an offset-free bias correction are used. The controller is tested against a sequence of step changes in the well flow-rate setpoints.

Two cases are run: one with the offset-free correction and one without.

The code for extracting LedaFlow results and for running LedaFlow via JavaScript is based on Marie
Sunde's master thesis. This work extends it and runs a closed-loop simulation with several step
changes in the wellhead flow-rate setpoints.

## Requirements

- **LedaFlow** (tested with v2.11.271.018)
- **Ipopt with the HSL `ma27` linear solver.** Both models set `linear_solver = "ma27"`.
- **Julia packages:** Refer to `Manifest.toml` and `Project.toml`  
- A **LedaFlow case** whose caseid matches the `lf_case_id` set in the case file. The two cases in
  this repository point at two different LedaFlow cases.

## The two cases

| Case | Offset-free bias | Purpose |
| --- | --- | --- |
| `case_final` | on — well flow and constrained-node pressure bias | Main result |
| `case_final_no_offset` | off | Baseline showing the steady-state offset from plant-model mismatch |

The two case files are otherwise identical, differing only in `CASE_NAME`, `lf_case_id`,
`well_flow_bias_bool` and `pressure_bias_bool`. 

## Running a case

Open `case_final.ipynb` (or `case_final_no_offset.ipynb`) and run the cells in order. The notebook
solves the steady-state model, initialises LedaFlow at that operating point, and then alternates
NMPC and LedaFlow one sample at a time for the whole simulation.

Afterwards, to plot the results:

1. Run `plot_lf_results.ipynb` once for the case, to dump every LedaFlow trend logger to
   `<case>/full_lf_output.csv`. The include for `extract_full_output.jl` is commented out at the top of
   `plot_lf_results.ipynb`; uncomment it for the first run.
2. Run `plot_lf_results.ipynb` to produce the closed-loop figures in `<case>/lf_plots/`.

## Simulation settings

Set in the case files. Both cases use the same values.

| Setting | Value |
| --- | --- |
| Sample time, `DT` | 300 s |
| Prediction horizon, `NP` | 3600 s (12 samples) |
| Control horizon, `NM` | 4 samples |
| Closed-loop simulation time, `NT` | 20 h (240 samples) |
| LedaFlow initialisation before the controller is switched on | 2000 s |

The LedaFlow clock therefore leads the closed-loop clock by 2000 s; `plot_lf_results.jl` calls this
offset `controller_on_time` and shifts the setpoint trajectory by it.

## Network

An arc-node network. 

| Arcs | Element |
| --- | --- |
| 1 | Mainline pump |
| 2–40 | Mainline, 39 km |
| 41–49 | Line 1, 9 km |
| 54–62 | Line 2, 9 km |
| 50, 52, 63, 65 | Wellhead chokes 1–4 |
| 51, 53, 64, 66 | Wellbores 1–4, 1.2 km each |

Node 41 splits the mainline into Lines 1 and 2; node 50 splits Line 1 to chokes 1 and 2; node 63
splits Line 2 to chokes 3 and 4. Nodes 52, 54, 65 and 67 are the bottom-hole nodes.

Soft constraints: pressure at nodes 2 (pump outlet), 52, 54, 65 and 67 (bottom hole), and the
cavitation index across each choke.

## Repository layout

| File / folder | Contents |
| --- | --- |
| `case_<name>.ipynb` | Driver notebook for case `<name>`: runs the whole closed-loop simulation and saves the results |
| `case_<name>.jl` | Configuration for case `<name>`: network, horizons, bounds, setpoints and objective function |
| `case_<name>/` | Results for case `<name>` — see below |

Contents of a results folder:

| File / folder | Contents |
| --- | --- |
| `mpc_pred/pred_t<t>.jld2` | The full predicted horizon from the NMPC at closed-loop time `<t>` seconds, one file per sample |
| `mpc_status.txt` | Ipopt termination status at each sample |
| `full_lf_output.csv` | Every LedaFlow trend logger over the whole run, written by `extract_full_output.jl` |
| `steady_state_trends.csv` | LedaFlow trends from the 2000 s initialisation |
| `trends.csv` | The subset of trends read back at each sample time; overwritten every mpc run |
| `lf_plots/` | Closed-loop figures produced by `plot_lf_results.jl` |
| `leda*.log` | LedaFlow softshell output |

Julia support files shared by every case:

- `generic_param.jl` — shared parameters: CO2 properties, the CoolProp property fits, the
  collocation matrix and the model variable lists.
- `model_ss.jl` — steady-state version of the NMPC internal model. Solved once per run to give a
  consistent operating point for LedaFlow and a starting guess for the first NMPC solve.
- `model_mpc.jl` — the first-principles dynamic model used inside the NMPC. Built once and reused
  for every sample.
- `run_mpc.jl` — assigns the initial condition, warm start, bias terms and parameter estimates to
  the NMPC, solves it, and returns the input to be applied.
- `run_sim.jl` — writes the MVs to LedaFlow, advances the simulation by one sample, reads the
  measurements back, and updates the state estimate, bias terms and parameter estimates.
- `lf_softshell.jl` — locates the LedaFlow `softsh` executable 
- `lf_ss.jl` — applies the steady-state operating point to LedaFlow and runs it to steady state
  (2000 s) before the closed loop starts.
- `extract_full_output.jl` — dumps every LedaFlow trend logger to `full_lf_output.csv`. Run once
  per case, after the simulation.

Generated JavaScript. LedaFlow is driven through its softshell, which takes a script, 
so these files are **rewritten by the Julia code on every call** and hold whatever values were
last written. They are outputs, not inputs — editing them has no effect.

- `ledaflow_ss.js` — written by `lf_ss.jl`: sets the steady-state pump speed and choke openings,
  then runs LedaFlow to steady state.
- `ledaflow_write.js` — written by `run_sim.jl`: writes the choke ramp and pump speed for the
  current sample and advances the simulation by `DT`.
- `ledaflow_extract.js` — written by `extract_full_output.jl`: dumps all trend loggers to CSV.

Plotting:

- `plot_lf_results.jl`, `plot_lf_results.ipynb` — closed-loop plant results over the whole run.
  The notebook includes the case file and then the `.jl` file.
- `plot_pred.jl`, `plot_pred.ipynb` — the NMPC prediction at a single chosen sample, read back from
  `mpc_pred/`. The x-axis spans one prediction horizon, not the whole run.
- `pressure_profile_plot.ipynb`, `pressure_profile_plot_single.ipynb` — spatial pressure profiles
  along the pipeline, from LedaFlow profile exports. Stacked panels and single-axis versions of the
  same figure. Specific to `case_final`.

Reference:

- `LF_full_output_headers.txt` — the full list of column names in `full_lf_output.csv`, for looking
  up trend logger names.
