# Optimal Inventory Control using MPC

This directory contains the code used to develop the cases in the article: 

Kumaraswamy, A., Turan, E. M., & Jäschke, J. (2025). Optimal inventory control for bottleneck isolation in general processes. In Journal of Process Control (Vol. 155). https://doi.org/10.1016/j.jprocont.2025.103557

It builds on the following work: 

Turan, E. M., Skogestad, S., & Jäschke, J. (2024). Model predictive control for bottleneck isolation with unmeasured faults. In IFAC-PapersOnLine (Vol. 58, Issue 14, pp. 767–774). https://doi.org/10.1016/j.ifacol.2024.08.430

This directory contains several networks of liquid tanks where an optimal inventory control MPC formulation has been applied. 

## Requirements

Julia 1.12.3, with:

| Purpose | Packages |
|---|---|
| Optimisation | `JuMP`, `HiGHS`, `Ipopt` |
| Weight algorithm | `Graphs`, `SimpleWeightedGraphs`, `GraphPlot`, `Compose` |
| Plotting | `Plots`, `LaTeXStrings`, `DelimitedFiles` |
| Misc | `Distributions`, `LinearAlgebra` |

## Running a case

Open the notebook for the case and run the cells in order:

1. **Include packages**
2. **Include relevant files** — loads the case definition, the MPC model, the weight algorithm and the plotting code
3. **Assign weights & build model** — runs the inventory weight algorithm and constructs the JuMP model once
4. **Run MPC** — closed-loop simulation; re-fixes initial conditions, disturbances and the previous input at each step, then re-solves
5. **Plot and save results** — writes figures, CSV profiles and a parameter record to `<case_name>_Results/`

Cell 3 must be re-run before cell 4 whenever you restart a simulation, so the initial conditions are reset.

## Repository layout

### Shared code

| File | Purpose |
|---|---|
| `MPC Model.jl` | Builds the JuMP model — variables, dynamics, constraints and the objective functions |
| `MPC assign weights.jl` | Algorithm assigning weights to each inventory based on its distance to the consumer |
| `Plot.jl` | Plots and saves flow/level profiles, CSV data and the parameter record |
| `Plot_RO.jl` | Same as `Plot.jl`, restricted to the storage tanks and using hourly units (Specific to Reverse Osmosis Case) |

### Case definitions

Each case is one `.jl` file holding every parameter for that scenario — network structure, MPC tuning, constraints and disturbances — paired with a notebook of the same name that runs it.

### Results

Each run writes `<case_name>_Results/` containing:

- `<case>_Flow.png`, `<case>_Inventory.png` — all flows / all levels
- `<case>_F0.png` … , `<case>_h1.png` … — individual profiles
- `<case>_Flow_profile.csv`, `<case>_Inventory_profile.csv` — raw trajectories
- `<case>_parameters.txt` — full parameter record for the run

## Case index

| Notebook | `case_name` | Tanks | Flows | Obj. | Distinguishing feature |
|---|---|:--:|:--:|:--:|---|
| `3 Tank` | `3Tk` | 3 | 4 | 1 | Three tanks in series; leak in tank 3 from t > 30 |
| `3 Tank Recycle` | `3Tk_Recycle` | 3 | 5 | 1 | Recycle stream; minimum-flow constraint active |
| `4 Tank - Merger - Case 1` | `4Tk_Merger_Case1` | 4 | 6 | 4 | Two producers merging; input regularisation, penalty 1e-3 |
| `4 Tank - Merger - Case 2` | `4Tk_Merger_Case2` | 4 | 6 | 4 | As Case 1 with penalty 1e-2 |
| `4 Tank - Split - Case 0` | `4Tk_Split_Case0` | 4 | 6 | 1 | Split to two consumers; baseline |
| `4 Tank - Split - Case 1` | `4Tk_Split_Case1` | 4 | 6 | 1 | Minimum-flow constraints on all flows |
| `4 Tank - Split - Case 2` | `4Tk_Split_Case2` | 4 | 6 | 1 | Leak in tank 3 from t > 30 |
| `5 Tank - Split Join - Case 0` | `5Tk_Split_Join_Case0` | 5 | 8 | 1 | Split and recombination; baseline |
| `5 Tank - Split Join - Case 1` | `5Tk_Split_Join_Case1` | 5 | 8 | 4 | Input regularisation |
| `5 Tank - Split Join - Case 2` | `5Tk_Split_Join_Case2` | 5 | 8 | 5 | Equality constraints on the parallel flows 5, 6, 7 |
| `Mini Network` | `Mini Network Scenario` | 8 | 13 | 4 | Larger network, 120 min horizon; minimum-flow and level-rate constraints |
| `Reverse Osmosis` | `Reverse Osmosis` | 7 | 10 | 7 | Process case; fixed flow ratios, 5 storage tanks and 2 process nodes |

## Objective functions

Selected per case with `obj_mpc` in the case file.

| `obj_mpc` | Objective |
|:--:|---|
| 0 | Unreachable flow setpoint |
| 1 | Weighted economic objective — consumer flows plus weighted inventory |
| 4 | As 1, with input regularisation penalising changes in flow rate |
| 5 | As 1, with equality constraints between selected flows added to the objective |
| 7 | As 1, with inventory terms and level constraints restricted to storage tanks |

Objectives 4, 5 and 7 are tuned via `input_reg_penalty`, `equality_con_param`
and `storage_tk_idx` respectively.

## Model formulation

The plant is a network of inventories linked by flows, described by an
incidence matrix `M`, with producer nodes `P` and consumer nodes `C`.
Levels follow `h(t+1) = h(t) + (M / a) f(t) + Bd`, where `a` is the tank
cross-sectional area and `Bd` a level disturbance such as a leak.

Constraints:

- Maximum flow `fmax` — hard constraint, time-varying
- Maximum and minimum level `hmax` / `hmin` — soft constraint, via slack variables
- Minimum flow `fmin` — soft, enabled with `fmin_binary`
- Maximum level change `deltah` — soft, enabled with `deltah_binary`
- Flows held constant beyond the control horizon
- Optional case-specific constraints via `special_constraints`

Constraint violations are penalised by `penalty`, scaled per constraint type by
`penalty_mult = [level, minimum flow, level rate]`.
