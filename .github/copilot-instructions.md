# Copilot instructions for this repository

- Preserve the existing research code and avoid unrelated formatting changes.
- For the Julia demo, use Julia 1.11-compatible syntax and the declared
  JuMP/HiGHS dependencies.
- Explain optimization assumptions, units, bounds, and solver status in
  proposed changes.
- Add a focused test for every behavior change and run the relevant Julia
  tests before claiming success.
- Never treat a plausible numerical result as proof of correctness; compare
  it with hand-checkable cases or independent calculations.
- Do not add secrets, generated solver artifacts, or large data files.
