# co2_span_wagner.jl
#
# Pure CO2 property table (Span–Wagner EOS via CoolProp) and a pump
# temperature-rise check for the LedaFlow case.
#
# Setup (once):
#   julia> using Pkg; Pkg.add(["CoolProp", "DelimitedFiles"])
#
# Run:
#   julia co2_span_wagner.jl
#
# Outputs:
#   co2_sw_table.csv   – long-format table (one row per P,T point)
#   console            – pump check: isentropic and real outlet temperature

using CoolProp
using DelimitedFiles
using Printf

const FLUID = "CO2"          # CoolProp HEOS backend = Span–Wagner for CO2

# ---------------------------------------------------------------------------
# 1. Property table
# ---------------------------------------------------------------------------

# Grid – adjust to suit. Fine spacing matters near the critical point
# (31.0 °C, 73.8 bar).
P_bar = collect(1.0:1.0:250.0)          # bara
T_C   = collect(-20.0:0.5:60.0)         # °C

const HEADER = ["P_bara", "T_C", "phase",
                "rho_kg_m3", "drhodP_kg_m3_bar", "drhodT_kg_m3_C",
                "h_J_kg", "s_J_kgK", "cp_J_kgK", "cv_J_kgK",
                "mu_Pa_s", "k_W_mK", "sos_m_s", "beta_1_K"]

# Single-phase properties at (P [Pa], T [K]).
function props(p, t)
    q(k) = PropsSI(k, "P", p, "T", t, FLUID)
    rho    = q("D")
    drhodP = PropsSI("d(D)/d(P)|T", "P", p, "T", t, FLUID)   # kg/m3/Pa
    drhodT = PropsSI("d(D)/d(T)|P", "P", p, "T", t, FLUID)   # kg/m3/K
    return (rho    = rho,
            drhodP = drhodP * 1e5,          # -> kg/m3/bar
            drhodT = drhodT,                # K and °C steps are identical
            h      = q("H"),
            s      = q("S"),
            cp     = q("C"),
            cv     = q("O"),
            mu     = q("V"),
            k      = q("L"),
            sos    = q("A"),
            beta   = -drhodT / rho)
end

phase_name(p, t) = PhaseSI("P", p, "T", t, FLUID)

function build_table(P_bar, T_C)
    rows = Vector{Vector{Any}}()
    nfail = 0
    for pb in P_bar, tc in T_C
        p, t = pb * 1e5, tc + 273.15
        try
            r = props(p, t)
            push!(rows, Any[pb, tc, phase_name(p, t),
                            r.rho, r.drhodP, r.drhodT,
                            r.h, r.s, r.cp, r.cv,
                            r.mu, r.k, r.sos, r.beta])
        catch
            # Points exactly on the saturation line can fail; flag them.
            nfail += 1
            push!(rows, Any[pb, tc, "FAILED", fill(NaN, length(HEADER) - 3)...])
        end
    end
    return rows, nfail
end

function write_table(path, rows)
    open(path, "w") do io
        println(io, join(HEADER, ","))
        for r in rows
            println(io, join((x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
                              for x in r), ","))
        end
    end
end

# ---------------------------------------------------------------------------
# 2. Pump check
# ---------------------------------------------------------------------------

"""
    pump_check(P1_bar, T1_C, P2_bar; eta, mdot)

Outlet temperature for a pump on pure CO2 (Span–Wagner):
  * isentropic outlet
  * real outlet, textbook loss:  h2 = h1 + (h2s - h1)/eta
  * LedaFlow-style loss on top of isentropic compression:
      P = Q·dp/eta,  E = P·(1/eta - 1)
"""
function pump_check(P1_bar, T1_C, P2_bar; eta = 0.85, mdot = 500.0)
    p1, t1, p2 = P1_bar * 1e5, T1_C + 273.15, P2_bar * 1e5

    h1  = PropsSI("H", "P", p1, "T", t1, FLUID)
    s1  = PropsSI("S", "P", p1, "T", t1, FLUID)
    rho = PropsSI("D", "P", p1, "T", t1, FLUID)

    h2s = PropsSI("H", "P", p2, "S", s1, FLUID)
    T2s = PropsSI("T", "P", p2, "S", s1, FLUID)

    h2  = h1 + (h2s - h1) / eta
    T2  = PropsSI("T", "P", p2, "H", h2, FLUID)

    Q   = mdot / rho                       # m3/s at inlet
    Pw  = Q * (p2 - p1) / eta              # W, shaft power (LedaFlow def.)
    E   = Pw * (1 / eta - 1)               # W, LedaFlow "enthalpy added"
    T2L = PropsSI("T", "P", p2, "H", h2s + E / mdot, FLUID)

    @printf("\nPump check (Span–Wagner, pure CO2)\n")
    @printf("  Inlet          : %.1f bara, %.2f °C, rho = %.1f kg/m3\n", P1_bar, T1_C, rho)
    @printf("  Volume flow    : %.0f m3/d   head = %.0f m\n",
            Q * 86400, (p2 - p1) / (rho * 9.81))
    @printf("  Isentropic     : T2s = %.2f °C   (dT = %.2f K)\n", T2s - 273.15, T2s - t1)
    @printf("  Real, eta=%.2f : T2  = %.2f °C   (dT = %.2f K)\n", eta, T2 - 273.15, T2 - t1)
    @printf("  LedaFlow-style : T2  = %.2f °C   (dT = %.2f K)\n", T2L - 273.15, T2L - t1)
    @printf("  Shaft power    : %.2f MW,  dissipation E = %.2f MW\n", Pw / 1e6, E / 1e6)
    return (T2s = T2s - 273.15, T2 = T2 - 273.15, T2_ledaflow = T2L - 273.15)
end

# ---------------------------------------------------------------------------
# 3. Spot values to compare against your LedaFlow PVT table
# ---------------------------------------------------------------------------

function spot_values(points)
    @printf("\nSpot values (compare with your PVT table)\n")
    @printf("  %7s %6s %9s %11s %9s %9s\n",
            "P[bar]", "T[°C]", "rho", "drho/dT", "cp", "beta")
    for (pb, tc) in points
        r = props(pb * 1e5, tc + 273.15)
        @printf("  %7.1f %6.1f %9.1f %11.3f %9.0f %9.5f\n",
                pb, tc, r.rho, r.drhodT, r.cp, r.beta)
    end
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


spot_values([(100.0, 8.0), (100.0, 15.0), (150.0, 12.0),
                 (192.0, 17.0), (197.0, 17.0)])
pump_check(100.0, 8.0, 197.0; eta = 0.85, mdot = 500.0)

println("\nBuilding table ($(length(P_bar)) x $(length(T_C)) points)...")
rows, nfail = build_table(P_bar, T_C)
write_table("co2_sw_table.csv", rows)
println("Wrote co2_sw_table.csv  ($(length(rows)) rows, $nfail failed points)")
