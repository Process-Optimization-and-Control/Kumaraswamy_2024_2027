using Test
using GithubRepoDemo

@testset "resource allocation" begin
    result = solve_allocation(
        80;
        capacities = [50, 60, 40],
        costs = [3, 1, 2],
    )

    @test sum(result.supplied) ≈ 80
    @test result.supplied ≈ [0, 60, 20]
    @test result.total_cost ≈ 100
end

@testset "input validation" begin
    @test_throws ArgumentError solve_allocation(
        10;
        capacities = [5],
        costs = [1, 2],
    )
    @test_throws ArgumentError solve_allocation(
        10;
        capacities = [5],
        costs = [1],
    )
    @test_throws ArgumentError solve_allocation(
        -1;
        capacities = [5],
        costs = [1],
    )
end
