# Copilot instructions for this repository

- Preserve the existing research code and avoid unrelated formatting changes.
- For the Julia demo, keep it dependency-free, Julia 1.10-compatible, and
  add tests as `@testset` blocks in `Github Repo Demo/test/runtests.jl`.
- Explain optimization assumptions, units, bounds, and solver status in
  proposed changes.
- Add a focused test for every behavior change and run the relevant Julia
  tests before claiming success.
- Never treat a plausible numerical result as proof of correctness; compare
  it with hand-checkable cases or independent calculations.
- Do not add secrets, generated solver artifacts, or large data files.
