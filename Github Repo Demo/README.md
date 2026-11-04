# GitHub Repo Demo

A tiny Julia package for practising GitHub's research-code workflow:
**issues → branches → pull requests → automated tests (GitHub Actions)**.

The code models a gravity-drained tank, `A·dh/dt = q_in − k·√h`. It is small
enough to check by hand, so the focus stays on the workflow. Everything below
can be done from **VS Code** and the **GitHub website**, with no terminal.

```
Github Repo Demo/
├── src/GithubRepoDemo.jl   # the model: steady_state_level, simulate_tank
├── test/runtests.jl        # tests, run automatically by GitHub Actions
├── demo.jl                 # interactive demo, run cell by cell
└── Project.toml            # package metadata; no runtime dependencies
```

Repository-level GitHub files live in `.github/` at the repo root:
`workflows/julia-demo.yml` (CI), `ISSUE_TEMPLATE/` and `pull_request_template.md`.

## 1. Try the demo

Open `demo.jl`, put your cursor in the first cell (cells are separated by
`##` lines) and press **Alt+Enter**. Repeat for each cell. The plot opens in
the VS Code plot pane. (The demo uses Plots.jl from your global Julia
environment.)

The tests in `test/runtests.jl` run automatically on GitHub whenever you
open or update a pull request, so you don't need to run them yourself.

## 2. Exercises

### A: Issue → branch → pull request

1. On GitHub, go to **Issues → New issue → Feature request** and propose
   *"Add a maximum tank height `h_max` (overflow)"*. Note the issue number,
   e.g. `#12`.
2. In the issue's right sidebar, click **Development → Create a branch**.
   Choose *Checkout locally*. GitHub names the branch after the issue.
3. In VS Code, switch to that branch (click the branch name in the bottom-left
   status bar).
4. Add an `h_max` keyword to `simulate_tank`, clipping the level with
   `min(..., h_max)`, and add a `@testset` to `test/runtests.jl` that checks
   the level never exceeds `h_max`.
5. Commit and push from the **Source Control** panel (Ctrl+Shift+G).
6. On GitHub, click **Compare & pull request**. The template asks for
   `Closes #12`, so merging the PR will close the issue automatically.
7. Wait for the **Julia 1.10** and **Julia 1** checks to turn green, then merge.

### B: Watch CI catch a mistake

1. Create a branch called `break-something` (status bar → *Create new branch*).
2. In `steady_state_level`, change `(q_in / k)^2` to `q_in / k`, which is
   a plausible-looking units mistake.
3. Commit, push, and open a **draft** pull request.
4. The checks turn red ❌. Click **Details** to read which test failed and why.
5. Revert the change, push again, and watch the checks turn green ✅.
   Close the PR without merging and delete the branch.

### C: Run a workflow by hand

Go to the **Actions** tab → *Julia demo tests* → **Run workflow**. This
re-runs the tests on GitHub's servers, for example after a package update.
(The button appears once the workflow file is on `main`.)

### D (optional): Protect `main`

**Settings → Branches → Add rule** for `main`, then tick *Require a pull request
before merging* and *Require status checks to pass*, choosing both Julia
checks. After that, nobody can merge red code into `main`. On an organisation
repo you need admin rights for this.

## What each piece teaches

| Piece | Why it matters for research code |
|---|---|
| Issues | A record of *why* a change was made, linked to the code that made it |
| Branches | Try ideas without breaking the version your results came from |
| Pull requests | A reviewable diff, a discussion and a checklist before code lands |
| GitHub Actions | Tests run on a clean machine, on two Julia versions, every time |
| Tests | Hand-checkable cases (e.g. `h_ss = (0.02/0.01)² = 4 m`) catch silent errors |
| No committed `Manifest.toml` | Avoids "works on my Julia version only" failures in CI |
