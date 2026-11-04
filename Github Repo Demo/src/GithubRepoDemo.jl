"""
A deliberately tiny process model used to practise the GitHub workflow
(issues → branches → pull requests → automated tests).

Model: a gravity-drained tank with inflow `q_in` and outflow `k * sqrt(h)`
(Torricelli's law):

    A * dh/dt = q_in - k * sqrt(h)

Units are SI throughout: h [m], A [m²], q_in [m³/s], k [m^2.5/s].
"""
module GithubRepoDemo

export steady_state_level, simulate_tank

"""
    steady_state_level(q_in; k)

Level [m] at which outflow balances inflow: `h_ss = (q_in / k)^2`.
"""
function steady_state_level(q_in::Real; k::Real)
    q_in >= 0 || throw(ArgumentError("inflow q_in must be non-negative"))
    k > 0 || throw(ArgumentError("outflow coefficient k must be positive"))
    return (q_in / k)^2
end

"""
    simulate_tank(h0, q_in; A, k, dt = 1.0, t_end = 100.0)

Integrate the tank level with explicit Euler. Returns a named tuple
`(t, h)` of vectors. The level is clipped at zero (an empty tank cannot
drain further).
"""
function simulate_tank(h0::Real, q_in::Real; A::Real, k::Real,
                       dt::Real = 1.0, t_end::Real = 100.0)
    h0 >= 0 || throw(ArgumentError("initial level h0 must be non-negative"))
    A > 0 || throw(ArgumentError("tank area A must be positive"))
    dt > 0 || throw(ArgumentError("time step dt must be positive"))
    steady_state_level(q_in; k)  # reuses the q_in and k checks

    t = collect(0.0:dt:t_end)
    h = similar(t)
    h[1] = h0
    for n in 1:length(t)-1
        dhdt = (q_in - k * sqrt(h[n])) / A
        h[n+1] = max(h[n] + dt * dhdt, 0.0)
    end
    return (t = t, h = h)
end

end
