# Interactive demo — run it cell by cell in VS Code, no terminal needed.
#
#   • Put your cursor inside a cell (the blocks between `##` lines)
#   • Press Alt+Enter (Option+Enter on Mac) to run that cell
#   • Plots appear in the VS Code plot pane
#
# Plots.jl is loaded from your global Julia environment, so it is not a
# dependency of the package itself.

##  Load the package (first run takes a few seconds)
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using GithubRepoDemo
using Plots

##  Parameters (SI units)
A    = 1.0    # tank cross-section [m²]
k    = 0.01   # outflow coefficient [m^2.5/s]
q_in = 0.02   # inflow [m³/s]
h0   = 0.5    # initial level [m]

h_ss = steady_state_level(q_in; k)
println("Steady-state level: $(h_ss) m")

##  Simulate and plot
sim = simulate_tank(h0, q_in; A, k, dt = 1.0, t_end = 3000.0)

plot(sim.t ./ 60, sim.h; label = "level h(t)", lw = 2,
     xlabel = "time [min]", ylabel = "level [m]", title = "Gravity-drained tank")
hline!([h_ss]; label = "steady state", ls = :dash)

##  Try it: change q_in above, re-run the cells, and watch the steady state move.
