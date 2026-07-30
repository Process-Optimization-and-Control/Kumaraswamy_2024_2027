# NPC Workshop 2026

This folder contains code, results, and figures submitted for the Nordic Process Control Workshop 2026. Here we build PID control frameworks in Julia and implement the calculated outputs in LedaFlow. LedaFlow models a $CO_2$ pipeline-injection network with 4 choke valves and a pump. Controllers manipulate the choke valves and pump speed.

LedaFlow is triggered to run using javascript at user-defined intervals. A number of control schemes are built and tested here.

## Repository layout

| File / folder | Contents |
| --- | --- |
| `ctrl_scheme_xx.jl` | The control scheme for case `xx`: setpoint schedule, tuning parameters and the PI calculations. Each file defines one `run_pid_controller` function. |
| `run_pid_scheme_xx.ipynb` | Driver for case `xx`. Runs LedaFlow to steady state, then advances the closed loop one time step at a time, calling the controller at each step. |
| `casexx_trends.csv` | Tabulated format of the results from LedaFlow for case `xx`. |
| `plot_lf_results_casexx.jl` / `.ipynb` | Code for generating the results. The notebook simply includes the `.jl` file and calls `plot_lf_case_results("casexx")`. |
| `casexx/` | Figures with the closed-loop simulation results. |

Support files shared by every case:

- `run_ledaflow_sim_script.jl` — writes the MVs to LedaFlow, advances the simulation by one time step, and reads the measurements back.
- `run_ledaflow_ss.jl` — resets the case and runs it to steady state (100 s) before the closed loop starts.
- `extract_full_output.jl` — writes all trend loggers to a CSV file.
- `lf_softshell.jl` — locates the LedaFlow `softsh` executable (override with the `LEDAFLOW_SOFTSH` environment variable).

## Control structure

All four schemes share the same regulatory layer:

- **Wells 2–4:** wellhead flow control via the choke valve opening.
- **Well 1:** choke on manual at a fixed opening of 1.0.
- **Mainline:** pressure control at 500 m via the pump speed.

Every controller is a PI in velocity form with the proportional term on the measurement. Outputs are clamped to the valve limits (0.2–1.0) and pump speed limits (2800–3900 rpm), and the velocity form integrates from the previous *clamped* output, which provides anti-windup.

Schemes B, C and D add two override layers on top of this:

- a **high-BHP override** on each of wells 2–4 (min-select on the choke opening), which pinches the choke back when the downhole pressure reaches its limit; and
- a **minimum mainline flow override** (max-select on the pump speed), which pushes the pump speed up when the mainline flow falls too low.

## Tuning

All gains come from the SIMC rules applied to the FOPDT models identified in `simc_tuning_calc.ipynb`:

$$K_c = \frac{\tau_1 / k}{\tau_c + \theta}, \qquad \tau_I = \min\left(\tau_1,\; 4(\tau_c + \theta)\right)$$

| Loop | MV | CV | $\tau_c$ | $K_c$ | $\tau_I$ (s) |
| --- | --- | --- | --- | --- | --- |
| Well flow control (wells 2–4) | Choke opening | Wellhead flow | 120 s | 4.27079e-5 | 0.1 |
| Mainline pressure control | Pump speed | Pressure at 500 m | 500 s | 29.5426 | 1039 |
| High-BHP override (wells 2–4) | Choke opening | Well BHP | 120 s | 0.00307629 | 7.2 |
| Minimum mainline flow override | Pump speed | Mainline flow at 19500 m | 500 s | 0.931889 | 62 |

Cases A, C and D all use this tuning. Case B keeps everything except the well flow controllers, which are deliberately spread apart in speed:

| Well flow controller | $\tau_c$ | $K_c$ | $\tau_I$ (s) |
| --- | --- | --- | --- |
| Well 2 | 200 s | 2.57855e-5 | 0.1 |
| Well 3 | 120 s | 4.27079e-5 | 0.1 |
| Well 4 | 90 s | 5.66496e-5 | 0.1 |

## The cases

Each case is run with a time step of $dt$ = 30 s. Note that the controller clock starts at 0 when the closed loop begins, whereas the trend files carry the preceding 100 s steady-state run, so setpoint overlays in the plotting scripts are offset by `t_ss = 100`.

- **Case A — baseline, no overrides** (`nt` = 7200 s). The plain regulatory layer only: well flow control and mainline pressure control, with neither override layer active. Well BHP and mainline flow are measured but unused. A step in the pressure controller's setpoint is evaluated (197.8 → 210.0 bar). This case shows what the regulatory layer does on its own, with nothing protecting the downhole pressure or the minimum mainline flow.

- **Case B — staggered well tuning** (`nt` = 20000 s). Identical to Case D in every respect — same setpoint schedule, same override layers, same override setpoints, same pump tuning — except that the three well flow controllers are given different $\tau_c$ values (200 / 120 / 90 s). Because the wells share the mainline header, giving all three loops the same speed in Case D makes them respond to a common disturbance at the same rate and fight each other through the header. Separating them in speed is intended to break that interaction.

- **Case C — isolated high-BHP override activation** (`nt` = 7200 s). Tuning is identical to Case D, but the scenario is deliberately minimal: the pump pressure and all well flow setpoints are held at nominal, so there is no slow pressure ramp and the run stays short. The event is instead a step in well 4's BHP *override* setpoint, from 248 to 243 bar at $t$ = 500 s, dropping it below the natural BHP so the min-select takes the choke on its own. This isolates the override handover from every other transient.

- **Case D — reference override scheme** (`nt` = 20000 s). The full scheme with both override layers and a uniform well tuning of $\tau_c$ = 120 s. The setpoint schedule moves one setpoint at a time and holds each phase long enough to settle:

  | Phase | Time (s) | Event |
  | --- | --- | --- |
  | 0 | 0 – 500 | Baseline: pump 197.8 bar, all wells 125 kg/s |
  | 1 | 500 – 4500 | Pump pressure 197.8 → 205 bar |
  | 2 | 4500 – 8500 | Well 4 flow 125 → 135 kg/s (BHP override engages near 248 bar) |
  | 3 | 8500 – 12500 | Well 4 flow 135 → 125 kg/s (override releases, loop recovers) |
  | 4 | 12500 – 20000 | Pump pressure 205 → 190 bar (mainline min-flow override engages) |

## Model identification (step tests)

This folder also includes the step tests used for model identification. A step test is carried out in LedaFlow for the pump speed and the choke valve. The results for these step-tests can be found in `case_speed_step_trends.csv` and `case_valve_step_trends.csv`. These results are used to fit a first-order-plus-dead-time (FOPDT) model and calculate tuning parameters using the SIMC tuning rules. These results can be found in `simc_tuning_calc.ipynb`, with the fitted responses saved as the `fit_*.png` and `*_step_test_*.png` figures.

| Step test | $k$ | $\tau_1$ (s) | $\theta$ (s) |
| --- | --- | --- | --- |
| Choke opening → wellhead flow | 19.2 | 0.1 | 1.9 |
| Choke opening → well BHP | 19.2 | 7.2 | 1.9 |
| Pump speed → mainline pressure at 500 m | 0.0702 | 1040 | 1.0 |
| Pump speed → mainline flow at 19500 m | 0.124 | 62.0 | 37.0 |

