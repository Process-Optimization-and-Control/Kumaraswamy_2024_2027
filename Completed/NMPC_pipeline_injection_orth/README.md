# NMPC for a CO₂ Pipeline-Injection Network — Orthogonal Collocation

This repository reproduces the case study of the **IFAC 2026 conference paper** on nonlinear MPC
for a CO₂ pipeline-injection network. The plant model, the network, all physical parameters, and all four case definitions are the same as in the paper.

**The one modelling difference is the time discretisation: the dynamics here are discretised by
orthogonal collocation on finite elements instead of trapezoidal integration.** Everything
downstream of that — the objective functions, the constraint set, the closed-loop driver, the
plotting — is unchanged.

**Note: For throughput maximisation cases, w_out from the wellbores are chosen as the variable for maximization. 
This is the flow directly being injected into the reservoir and essentially, the flow we seek to maximize. Maximizing wellhead flow rates would allow the controller to unpack inventory in the wellbores. **

## Orthogonal Collocation Implementation: 

Each sample of length `DT` is one finite element carrying `D` points along a collocation
dimension `d`:

| `d`     | meaning                                                                       |
|---------|-------------------------------------------------------------------------------|
| `1`     | start of the element — holds the initial condition, no model equation written  |
| `2:D`   | the collocation roots and the element end point; all model equations live here |
| `D`     | end of the element — the value reported as the result for that sample          |

`D = 4` (three Radau IIA collocation points at 0.1551, 0.6449 and 1.0 on the unit element).
`Adot` in [generic_param.jl](generic_param.jl) is the corresponding collocation derivative matrix,
and the mass balance is written as

```julia
Adot' * rho[k, 1:D, t] .== (DT / V[k]) * (w_in[k, 2:D, t] - w_out[k, 2:D, t])
```

Consequences that show up throughout the code:

* Every model variable gains a middle index `d`, so variables are `[k, d, t]`.
* The inputs (`choke_vlv_op`, `speed_pump`) carry **no** collocation dimension — one value per
  finite element, held constant across it (zero-order hold).
* Objectives, soft constraints and plots all read states at `d = D`.

### Parameter Changes 

`NP = 2400 s` (8 elements of 300 s) rather than the `3600 s` quoted in the paper. The collocation
formulation carries roughly three times as many variables per element, and the shorter horizon
keeps the solve time workable. This is intentional — see the case files.

## Repository layout

### Model and driver files

| File | Role |
|------|------|
| [generic_param.jl](generic_param.jl) | CO₂ properties (CoolProp fits), collocation parameters `D` and `Adot`, model variable list |
| [model_ss.jl](model_ss.jl) | Steady-state model. Square feasibility problem; supplies the initial condition and the first warm start |
| [model_mpc.jl](model_mpc.jl) | NMPC internal model over `T_mpc` finite elements. Built once, reused every sample |
| [model_sim.jl](model_sim.jl) | Plant model. Same equations, one finite element, inputs fixed to the NMPC move |
| [run_mpc.jl](run_mpc.jl) | Solves the NMPC for one sample: warm start, fix IC, attach objective, return first move |
| [run_sim.jl](run_sim.jl) | Advances the plant by one element and returns the state at `d = D` |
| [plot.jl](plot.jl) | Writes the per-case PNG figures and JSON series into the case folder |

### Case definitions

Each case file holds the network, parameters, horizon, bounds, setpoints and the NMPC objective.

| Case | File | Objective |
|------|------|-----------|
| A | [case1a.jl](case1a.jl) | Track well flow setpoints |
| B | [case1b.jl](case1b.jl) | Track well flow setpoints **+ pump work penalty** |
| C | [case1c.jl](case1c.jl) | Maximise well flows (weighted linear reward, wells compete) |
| D | [case1d.jl](case1d.jl) | Maximise well flows **+ pump work penalty** |

### Notebooks

| Notebook | Purpose |
|----------|---------|
| [case1a.ipynb](case1a.ipynb) … [case1d.ipynb](case1d.ipynb) | Run one closed-loop case end to end and write `<case>/results.jld2` |
| [IFAC_plots_wellFC.ipynb](IFAC_plots_wellFC.ipynb) | Paper figures for cases A and B |
| [IFAC_plots_max_throughput.ipynb](IFAC_plots_max_throughput.ipynb) | Paper figures for cases C and D |
| [IFAC_plots_pump_speed_subplot.ipynb](IFAC_plots_pump_speed_subplot.ipynb) | Pump speed comparison across all four cases |
| [pump work calc.ipynb](pump%20work%20calc.ipynb) | Totals the pump work series for each case |
| [xx To reorganise/](xx%20To%20reorganise/) | Open-loop step test on the simulation model alone; not part of the paper results |
