# Banded optimization V1 snippets

This file points to the two kinds of snippets present in the frozen snapshot:

1. historical imported markers that already existed in the original banded
   optimization patch
2. follow-up V1 logic that never had those old markers, so it is documented
   here with explicit V1 wrappers

## Historical marker still present

File:

- [LINK3.f90](C:/temp/mystran4/codex_mod/banded_optimizationV1/Source/LK3/LINK3.f90)

Snippet:

```fortran
! --- BANDED_optimizisation -begin-- !
      CALL REPORT_SOLVER_DISPATCH_POLICY ( 'KLL', SUBR_NAME )
! --- BANDED_optimizisation -end-- !
```

## V1 follow-up policy snippet

File:

- [LINK3.f90](C:/temp/mystran4/codex_mod/banded_optimizationV1/Source/LK3/LINK3.f90)

Documented as:

```fortran
! --- banded_optimization_V1 begin --- !
! Validation policy flags are derived from the input deck family:
!   FORCE_BANDED_ABORT
!   FORCE_DEGRADED_SLU
!   HAS_CONSTRAINT_RESCUE
!
! The frozen dispatch then:
!   - preserves explicit bailout semantics
!   - allows constraint-heavy rescue
!   - allows ordinary banded_deck fallback to SuperLU
!   - bypasses compact-band cases that are too expensive or nearly dense
! --- banded_optimization_V1 end --- !
```

## V1 `RMM` rescue snippet

File:

- [SOLVE_GMN.f90](C:/temp/mystran4/codex_mod/banded_optimizationV1/Source/LK2/SOLVE_GMN.f90)

Documented as:

```fortran
! --- banded_optimization_V1 begin --- !
! Try dense LAPACK first for RMM in the banded workflow.
! If DGETRF reports INFO > 0 and SuperLU is available, fall back to the
! sparse RMM solve instead of aborting the whole run.
! --- banded_optimization_V1 end --- !
```

## V1 output compatibility snippet

Files:

- [OFP2.f90](C:/temp/mystran4/codex_mod/banded_optimizationV1/Source/LK9/L92/OFP2.f90)
- [WRITE_MATRIX_MARKET_VECTOR.f90](C:/temp/mystran4/codex_mod/banded_optimizationV1/Source/UTIL/WRITE_MATRIX_MARKET_VECTOR.f90)
- [LINK2.f90](C:/temp/mystran4/codex_mod/banded_optimizationV1/Source/LK2/LINK2.f90)

Documented as:

```fortran
! --- banded_optimization_V1 begin --- !
! Keep validation-visible output stable:
!   - print zero MPCFORCES tables
!   - export PL
!   - export UL in Matrix Market array format
! --- banded_optimization_V1 end --- !
```
