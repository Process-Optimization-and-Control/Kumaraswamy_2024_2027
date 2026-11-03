module GithubRepoDemo

using HiGHS
using JuMP
import MathOptInterface as MOI

export solve_allocation

"""
    solve_allocation(demand; capacities, costs)

Minimize the operating cost of supplying a process demand from several
candidate sources. Each source has a maximum capacity and a unit cost.
"""
function solve_allocation(
    demand::Real;
    capacities::AbstractVector{<:Real},
    costs::AbstractVector{<:Real},
)
    demand > 0 || throw(ArgumentError("demand must be positive"))
    length(capacities) == length(costs) ||
        throw(ArgumentError("capacities and costs must have the same length"))
    !isempty(capacities) || throw(ArgumentError("at least one source is required"))
    all(>=(0), capacities) || throw(ArgumentError("capacities must be non-negative"))
    all(>=(0), costs) || throw(ArgumentError("costs must be non-negative"))
    sum(capacities) >= demand ||
        throw(ArgumentError("total capacity must meet the demand"))

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= supplied[i = 1:length(capacities)] <= capacities[i])
    @constraint(model, sum(supplied) >= demand)
    @objective(model, Min, sum(costs[i] * supplied[i] for i in eachindex(costs)))
    optimize!(model)

    termination_status(model) == MOI.OPTIMAL ||
        error("optimization did not reach an optimal solution")

    return (
        supplied = value.(supplied),
        total_cost = objective_value(model),
        demand = Float64(demand),
    )
end

end
