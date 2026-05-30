# Patch order for banded_optimizationV1

`banded_optimizationV1` must be layered after `lapack_peel_off`.

Recommended order for this branch:

1. baseline BLAS / sparse-library setup
2. large `lapack_surgery`
3. `lapack_peel_off`
4. `banded_optimizationV1`

Why this matters:

- `banded_optimizationV1` assumes the LAPACK-facing solve flow has already been
  reorganized by `lapack_peel_off`.
- The frozen `LINK3` policy depends on the peeled-off factorization/solve
  structure to decide when to:
  - stay on `DPBTRF`
  - skip dense fallback
  - rescue to `SuperLU`
  - preserve explicit bailout semantics
- The frozen `SOLVE_GMN` fallback also assumes the post-peel-off dense solve
  flow exists cleanly enough to intercept `DGETRF` failure and route `RMM` to
  `SuperLU`.

Practical guidance:

- If you are replaying patches one by one into an older tree, do not try to
  land `banded_optimizationV1` on a pre-`lapack_peel_off` source layout.
- If merge conflicts appear around `LINK3`, `SOLVE_GMN`, or matrix-export
  utilities, first confirm that `lapack_peel_off` has already been applied.
