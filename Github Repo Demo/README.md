# GitHub Repo Demo: process resource allocation

This small Julia project demonstrates a reproducible optimization experiment
and a practical GitHub workflow for PhD students. It chooses how much material
to purchase from several sources while meeting a process demand at minimum
cost.

The model is deliberately small:

```text
minimize  sum(cost[i] * supplied[i])
subject to sum(supplied) >= demand
           0 <= supplied[i] <= capacity[i]
```

It uses [JuMP](https://jump.dev/) to express the model and
[HiGHS](https://highs.dev/) as an open-source linear optimization solver.

## Run it

From this directory, install the project dependencies and run the tests:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

You can also try the model interactively:

```julia
using GithubRepoDemo

solve_allocation(80; capacities=[50, 60, 40], costs=[3, 1, 2])
```

The expected allocation is `[0, 60, 20]`: use the cheapest source first, then
the next cheapest source, while meeting the demand of `80`.

## Suggested GitHub exercise

1. Open an issue proposing a new constraint, such as a minimum purchase from a
   local supplier.
2. Create a branch for the change.
3. Ask Copilot to explain the current model and suggest tests before editing it.
4. Implement the change and add a test that would fail without it.
5. Open a pull request using the repository checklist.
6. Review the result for units, assumptions, solver status, and reproducibility.

This workflow mirrors a small research change: define the question, make the
assumption explicit, test it, and preserve a reviewable record.

## What this demonstrates

- `Project.toml` makes Julia dependencies explicit.
- `test/runtests.jl` protects model behavior and input validation.
- GitHub Actions runs the same tests automatically on every change.
- Issues and pull requests make research decisions reviewable.
- Copilot can accelerate implementation, but the researcher remains
  responsible for checking equations, units, solver status, and results.
