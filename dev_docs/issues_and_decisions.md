# Issues and decisions for banded_optimizationV1

## Main decision

The frozen runtime policy is:

- use true banded factorization when the matrix is a good banded candidate
- skip banded for compact-band cases that are too expensive or almost dense
- rescue to `SuperLU` where validation or robustness requires it

## Why this was needed

The raw banded path was not enough to keep the validation suite stable after
the newer LAPACK work. The key trouble spots were:

- huge effective bandwidth cases such as the Raasch MITC4/MITC8 hooks
- constraint-heavy `RMM` solves for `RBE3`
- explicit bailout decks that encode historical validation semantics
- output paths that suppressed zero tables expected by the validator

## Important case stories

### MITC4 / MITC8 Raasch

These showed that `KLL` can be sparse while still becoming a terrible compact
band candidate. In those cases the half-band approaches the matrix size, so the
banded storage behaves almost like dense storage. The frozen policy therefore
lets `LINK3` jump to `SuperLU` instead of forcing a bad banded solve.

### `RBE3 10`

This was not really a `KLL` problem. The failure was in `RMM` inside
`SOLVE_GMN`, where dense `DGETRF` was less robust than the sparse path. The
frozen V1 state keeps the banded workflow, but rescues `RMM` to `SuperLU` when
dense factorization fails.

### `BANDED BAILOUT` and `SPARSE BAILOUT -1`

The validator is not checking one universal numerical rule here. It is checking
named historical behaviors. The frozen V1 policy preserves those semantics so
the validation deck family still means what the suite expects it to mean.

## Explicit non-goal in V1

Skyline was tested as an intermediate fallback, but it caused regressions in
the production validation run. It was removed from the active runtime path and
left for future research only.
