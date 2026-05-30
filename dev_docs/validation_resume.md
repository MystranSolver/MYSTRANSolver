# Validation resume for banded_optimizationV1

Validation target:

- `C:\temp\mystran4\MYSTRAN_Validation-main\test_banded.py`

Frozen result:

- `0/2605 failed -> PASS`

Key repaired themes:

1. Constraint and rigid-element rescue

`LINK3` now allows sparse `KLL` rescue for constraint-heavy decks where a pure
banded path is not robust enough.

Examples:

- `vic/9/S30 node-surface coupling RBE3 10.bdf`
- `vic/9/S30 node-surface coupling RBE3 15.bdf`
- `vic/10/NAS S30 all zero MPC forces.bdf`
- `vic/10/NAS S30 AUTOSPC reaction forces.bdf`
- `vic/12/S30 Shell pinned support multiple nodes.bdf`

2. Explicit bailout semantics

The validation suite carries historical expectations for named bailout families.
The frozen V1 policy preserves those expectations instead of forcing one global
numerical behavior onto every bailout deck.

Families handled explicitly:

- `BANDED BAILOUT`
- `SPARSE BAILOUT -1`

3. General non-SPD banded decks

Ordinary `banded_deck` models that are not explicit bailout families may rescue
to `SuperLU` when banded or dense fallback cannot produce a valid factorization.

Examples:

- `vic/1/NAS S30 shell laminate coupling.bdf`
- `vic/1/NAS S30 shell laminate coupling_2.bdf`
- `vic/5/S30 mitc8 extreme high x low E.bdf`
- `vic/6/NAS S30 quad4 mitc4+ flat.bdf`
- `vic/7/S30 mitc4 macneal plate bending patch test.bdf`

4. `SOLVE_GMN` robustness

Problematic `RMM` solves now fall back from dense LAPACK to `SuperLU` in
banded runs, which fixed the hard `RBE3` validation failures.

5. Output and Matrix Market support

- zero-valued `MPCFORCES` tables are now printed
- `PL` export is available
- `UL` is exported in Matrix Market `array` form

Observed banded vs sparse split from `passed_banded` unique decks:

- total counted: `272`
- true banded KLL path: `260` = `95.59%`
- `SuperLU` fallback KLL path: `12` = `4.41%`
