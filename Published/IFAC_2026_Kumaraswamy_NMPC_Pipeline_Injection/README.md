# NMPC for a $CO_2$ Pipeline-Injection Network

Code and results for the IFAC World Congress 2026 submission on the nonlinear model predictive control (NMPC) of a $CO_2$ pipeline-injection network.

The study asks two questions:

1. Can a single NMPC layer coordinate a mainline pump and four wellhead chokes so that individual
   well flow rates follow prioritised setpoints, while respecting pipeline pressure, pressure
   rate-of-change and choke cavitation limits?
2. What does adding an explicit pump-energy term to the objective cost in tracking or throughput
   performance, and what does it save in pump work?

Everything is written in Julia. The models are built in JuMP and solved with Ipopt (using the HSL
`ma57` linear solver); $CO_2$ properties come from CoolProp.

---

## The physical system

```
                                                    ┌─ choke 1 ─ well 1 ─┐
                              ┌── Line 1 (9 km) ────┤                    ├── reservoir
 P_IN ─ pump ─ mainline ──────┤   node 50           └─ choke 2 ─ well 2 ─┘
        arc 1  (39 km)        │
                node 41       │                     ┌─ choke 3 ─ well 3 ─┐
                              └── Line 2 (9 km) ────┤                    ├── reservoir
                                  node 63           └─ choke 4 ─ well 4 ─┘
```

The network is modelled using an arc-node structure: 

| Element | Arcs | Geometry |
|---|---|---|
| Mainline pump | 1 | Pump Curve, 3000 rpm reference, 85 % efficiency |
| Mainline | 2–40 | 39 × 1 km, ID 0.4572 m, −0.44° slope |
| Line 1 | 41–49 | 9 × 1 km, ID 0.2984 m |
| Line 2 | 54–62 | 9 × 1 km, ID 0.2984 m |
| Choke valves 1–4 | 50, 52, 63, 65 |  |
| Wellbores 1–4 | 51, 53, 64, 66 | 1.2 km, ID 0.1571 m, vertical (−90°) |

Boundary conditions: supply pressure `P_IN = 100 bar` at the pump inlet, and a linear reservoir
inflow `w = injectivity · (BHP − P_RES)` at each bottom hole with `P_RES = 120 bar`.

## The models

All three models share the same first-principles equations and are built from the same parameter
files, so they stay consistent with each other.

- **Mass balance** per pipe segment, `V·dρ/dt = w_in − w_out`, integrated with the trapezoidal rule.
- **Momentum balance** at segment inlet and outlet: Darcy friction (Haaland explicit correlation)
- **Fluid properties**: the compressibility factor `Z` and viscosity are quadratic least-squares
  fits of CoolProp data over 50–300 bar at a fixed 283.15 K, and density follows from the real-gas
  law. This keeps the model smooth and differentiable for Ipopt.
- **Pump**: Pump curve used to determine pressure drop and flow at various speeds
- **Chokes**: orifice equation `w = Cv · opening · √(ρ·Δp)`, with a non-negative pressure drop so
  flow cannot reverse. Cavitation index `CI = (p_in − p_vap)/(p_in − p_out)`.

| File | Role |
|---|---|
| `model_ss.jl` | Steady-state feasibility solve. No time dimension. Chokes and total throughput are fixed, leaving pump speed as the only unknown. Provides both the initial condition and the initial guess for everything downstream. |
| `model_mpc.jl` | The NMPC internal model over the prediction horizon, with MV bounds, MV rate-of-change limits and soft pressure / pressure-ROC / cavitation constraints. Built once and reused for every sample. |
| `model_sim.jl` | The "plant". Same equations, solved as a square feasibility problem over one sample with the inputs held at the move the NMPC just returned. |

## Controller setup

| Quantity | Value |
|---|---|
| Sample time `DT` | 300 s |
| Prediction horizon `NP` | 3600 s (12 steps) |
| Closed-loop run `NT` | 18 000 s (5 h, 60 samples) |
| MVs | pump speed (3000–3700 rpm, ≤ 0.03 rpm/s) and 4 choke openings (0.2–1.0, ≤ 1e-4 s⁻¹) |
| Differential state | segment average pressure `p_average` (density, and hence inventory, follows from it) |
| Soft constraints | node pressure 180–250 bar at the pump outlet and the four bottom holes; pressure ROC ≤ 0.005 bar/s; cavitation index ≥ 1.7 |

Constraints are enforced through slack variables penalised in the objective, so the problem stays feasible under transients.
Input movement is additionally regularised with an L2 penalty.

## The four cases

Wells are deliberately given unequal priority, weights **5 / 20 / 80 / 200** for wells 1–4, so the
controller has to choose which well to sacrifice when the network cannot serve everything.

| Case | Directory | Objective |
|---|---|---|
| **A** | `01_wellFC_prioritisation` | Track well flow setpoints (125 → 130 kg/s step at t = 2 h). No pump-energy term. |
| **B** | `02_wellFC_prioritisation_opt` | As A, plus a pump-work penalty. |
| **C** | `03_max_throughput` | Maximise well flow rates until a constraint becomes active. No pump-energy term. |
| **D** | `04_max_throughput_opt` | As C, plus a pump-work penalty. |

Total pump work over the 5 h run reported in `pump work calc.ipynb`

## Repository layout

```
generic_param.jl              $CO_2$ properties, CoolProp property fits, model variable name lists
0X_<case>.jl                  per-case: topology, geometry, tuning, setpoints, objective function
0X_<case>.ipynb               per-case driver notebook (steady state → build → closed-loop → plot)

model_ss.jl                   steady-state model
model_mpc.jl                  NMPC internal model
model_sim.jl                  simulation ("plant") model
run_mpc.jl                    one NMPC solve: warm start, solve, shift guess, extract optimal move
run_sim.jl                    advance the plant by one sample under the applied move
plot.jl                       per-case figures and JSON export of every model variable

IFAC_plots_wellFC.ipynb       paper figures, cases A/B
IFAC_plots_max_throughput.ipynb  paper figures, cases C/D
IFAC_plots_pump_speed_subplot.ipynb  pump-speed comparison across all four cases
pump work calc.ipynb          total pump energy per case
plot_property_fit.ipynb       quality of the CoolProp property fits

0X_<case>/                    results for that case:
    results.jld2                closed-loop trajectory of every model variable
    <case>_<var>.json           the same, one file per variable
    <case>_pltN.png             standard per-case plots
    mpc_status.txt              Ipopt termination status, one line per sample
    sim_status.txt              ditto for the simulation model
    mpc_pred/pred_t<t>.jld2     full NMPC prediction dump at each sample
```

## Running a case

Open `0X_<case>.ipynb` and run it top to bottom. Each notebook:

1. includes the shared `.jl` files and its own case parameter file,
2. creates the results directories,
3. solves the steady-state model for the initial condition and initial guess,
4. builds the NMPC and simulation models once,
5. loops for `NT/DT` samples, alternating `run_optimiser` and `run_simulator`,
6. saves `results.jld2` and calls `plot_save`.

Run all four cases before the `IFAC_plots_*` notebooks, which load `results.jld2` from each case
directory to build the comparison figures.


## Conventions

- Pressures are in **bar**, flows in **kg/s**, lengths in **m**, pump speed in **rpm**, pump work in
  **MWh** per sample interval.
- Arc indices are used everywhere as the primary key: `p_in[k, t]`, `w_in[k, t]` and so on are
  defined for all 66 arcs, and the pipe/pump/valve subsets pick out which equations apply.
- Time index `t = 1` is always the current (fixed) condition; `t = 2` is the first move the
  controller is free to choose, which is why every objective and MV constraint runs from `t = 2`.
