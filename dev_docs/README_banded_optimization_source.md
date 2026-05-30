# Banded Optimization Patch

This folder mirrors the `optimization_rcm_v2` banded-storage diagnostic patch from `MYSTRANSolver-18.0.0`.

Source commit:

```text
ca2ee287c48a9bbbd79c1b50a5d2f1d18f87ec2d
Add banded storage optimization diagnostics
```

## Contents

- `mystran/Source/...`: full changed source files, with `! --- BANDED_optimizisation -begin-- !` / `! --- BANDED_optimizisation -end-- !` markers around the added banded optimization blocks.
- `mystran/dev_docs/banded_optimization_proposal_v1.md`: staged proposal for banded storage optimization.
- `mystran/dev_docs/banded_path_audit_v1.md`: audit of the current banded path.
- `snippets_banded_optimization.md`: commit patch/snippets for manual review.

## Scope

Included:

- Stage 1 banded storage estimator for compact band, CSR, and skyline/profile memory estimates.
- Visible solver dispatch diagnostics for `BANDED_ONLY` vs sparse behavior.
- `WINAMEM` modernization so `WINAMEM <= 0` disables the historical Windows XP memory cap.
- Banded allocation request diagnostics for `ABAND`, `BBAND`, and `RFAC`.
- RCM diagnostic text confirming `GRID_SEQ/INV_GRID_SEQ` ordering is shared before DOF numbering, so `K` and `M` stay aligned.

Not included:

- No `SUBSP` dispatch experiment.
- No automatic ARPACK/SuperLU/BANDED reroute for SUBSPACE.
- No silent SuperLU fallback for `SOLLIB=BANDED`.
- No hijack BLAS changes.

## Verification Used

- `cmake --build build --target mystran --config Release -- -j8`
- Static banded smoke deck:
  `codex_mod/rcm_add/mystran/examples/midas_static_24_cquad4_moremesh_rcm_sorted_banded.dat`

