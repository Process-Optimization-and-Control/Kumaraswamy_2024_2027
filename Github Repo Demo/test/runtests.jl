# Run by `Pkg.test()` and by GitHub Actions on every pull request.
using Test
using GithubRepoDemo

@testset "GithubRepoDemo" begin

    @testset "steady state matches hand calculation" begin
        # q_in = 0.02, k = 0.01  →  h_ss = (0.02 / 0.01)^2 = 4 m
        @test steady_state_level(0.02; k = 0.01) ≈ 4.0
        @test steady_state_level(0.0; k = 0.01) == 0.0
    end

    @testset "simulation settles at steady state" begin
        sim = simulate_tank(0.5, 0.02; A = 1.0, k = 0.01, dt = 1.0, t_end = 5000.0)
        @test sim.h[end] ≈ steady_state_level(0.02; k = 0.01) rtol = 1e-3
    end

    @testset "starting at steady state stays there" begin
        h_ss = steady_state_level(0.02; k = 0.01)
        sim = simulate_tank(h_ss, 0.02; A = 1.0, k = 0.01)
        @test all(h -> isapprox(h, h_ss), sim.h)
    end

    @testset "draining tank never goes negative" begin
        sim = simulate_tank(1.0, 0.0; A = 0.1, k = 0.05, dt = 1.0, t_end = 100.0)
        @test all(>=(0), sim.h)
        @test sim.h[end] == 0.0
    end

    @testset "invalid inputs are rejected" begin
        @test_throws ArgumentError steady_state_level(-1.0; k = 0.01)
        @test_throws ArgumentError steady_state_level(0.02; k = 0.0)
        @test_throws ArgumentError simulate_tank(-1.0, 0.02; A = 1.0, k = 0.01)
        @test_throws ArgumentError simulate_tank(1.0, 0.02; A = 0.0, k = 0.01)
    end

    @testset "simulation settles at steady state" begin
        sim = simulate_tank(0.5, 0.02; A = 1.0, k = 0.01, dt = 1.0, t_end = 5000.0)
        @test sim.h[end] ≈ steady_state_level(0.02; k = 0.01) rtol = 1e-3
    end

end
