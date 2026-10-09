# MYSTRAN performance changes and fixes, October 2026

This work builds on the reorganized source of PR #307 (`Bruno02468:cleanup` at 1b13f0b5, itself
on `dev` 8d50d390). The commits are listed in Appendix A. Each was checked against stock on
the same machine, with the method in section 4. In what follows:
- **Stock:** `dev` 8d50d390, built here with the same compiler and flags as this PR. The shipped
  19.0.0 binary is faster than that build (15.2 s against 18.1 s on the 128 × 128 plate), so
  against the release the plate speedup is about 1.8×.
- **This PR:** all the commits of this pull request.

## 1. Results at a glance

### Public decks with a reference solver

Times are wall time on one thread, one solve at a time.

**Triangular plate series** (`MYSTRAN_Resources`: static, CTRIA3, PLOAD2):

| Plate | Grids | Stock | This PR | "NASTRAN" (2022 workbook) |
|---|---:|---:|---:|---:|
| 128 × 128 | 19,235 | 18.1 s | 8.9 s | 8.06 s |
| 256 × 256 | 76,364 | 98.8 s | **42.2 s** | 64.6 s |

- **Results:** both sizes print the same displacements and stresses as stock.
- **The "NASTRAN" column** is copied from the workbook that ships with these decks. The workbook
  gives no version, hardware or thread count, so treat it as a rough reference.

**pCRM9 wing** (TU Delft, SOL 103, 30 modes, 2,043 grids). Its published files include the
reference run's F06 and log: MSC Nastran 2020, `smp=2` on an Intel Xeon E5-2640 v4 at 2.40 GHz.

| | Wall time | Modes 1 to 21 against the reference |
|---|---:|---|
| MSC Nastran 2020 (published log) | 11.65 s | reference |
| MYSTRAN, this PR, 1 thread | **7.1 s** | within 0.34% |
| MYSTRAN stock | 55.9 s | wrong: 30 negative eigenvalues (see below) |

| Mode | MYSTRAN (Hz) | Reference (Hz) | Difference |
|---:|---:|---:|---:|
| 1 | 0.6881 | 0.6878 | +0.04% |
| 2 | 1.8255 | 1.8248 | +0.04% |
| 3 | 2.1689 | 2.1678 | +0.05% |
| 5 | 5.5321 | 5.5298 | +0.04% |
| 10 | 15.5078 | 15.5595 | −0.33% |
| 16 | 25.6025 | 25.6893 | −0.34% |
| 21 | 35.5521 | 35.4964 | +0.16% |

- **The MYSTRAN copy of the deck** replaces each CBEAM with a CBAR and each PBEAML BAR with a
  PBARL BAR, since MYSTRAN's CBEAM needs a PBEAM. Nothing else is changed.
- **Stock returns 30 negative eigenvalues** of nearly the same size (0.54 to 0.67 Hz as printed)
  instead of the modes. One CONM2 with an offset in two directions gets a mass matrix with the
  wrong sign on a product of inertia, and the mass matrix is then indefinite. The CONM2 fix
  (section 2.8) corrects this.
- **Modes 22 to 30** are 0.5% to 12% high. They are probably local rib and stringer modes, where
  the CBAR substitution matters. Not investigated further.
- **The two machines differ:** a 2024 laptop core against two cores of a 2016 server, and the
  reference time includes license checkout. Read the times as comparable, not as a ranking.

### Rigid elements

The plate below, with 900 spiders added: each spider is a grid above a 3 × 3 patch of plate grids,
joined to it by an RBE2 or an RBE3 (353,886 G-set DOF, 5,400 M-set DOF). Displacements, SPC forces
and MPC forces are printed. `Build_Test_Cases/statics/make_spider_plate_decks.py` writes both
decks.

| Deck | Stock | This PR | Results against stock |
|---|---:|---:|---|
| Plate + 900 RBE2 | 507.6 s | **22.2 s (23×)** | displacements of all 58,981 grids within 9e-7 |
| Plate + 900 RBE3 | 127.2 s | **20.3 s (6.3×)** | within 8e-7 |

The differences are in the last printed digit.

The GMN solve now grows linearly with the number of rigid elements. On the same plate with RBE3
spiders, it went from 4.0 s to 0.2 s with 900 RBE3s, from 56.4 s to 0.6 s with 3,481, and from
130.2 s to 1.1 s with 6,241 (section 2.6).

### A larger private model

An aircraft shell model of 216,846 DOF, and a flat plate of 240 × 240 CQUAD4. Stock here is the
reorganized branch of PR #307 as it is; the I/O changes of section 2.3 were included, the rigid-
element changes were not yet (this model has 14 RBE3s). "16 threads" is the SuperLU_MT build with
`OMP_NUM_THREADS=16`.

| Model | Run | Threads | Stock | With these changes |
|---|---|---|---|---|
| Aircraft shell model, 216,846 DOF | static | 1 | 106.3 s | **41.0 s (2.6×)** |
| | | 16 | 42.4 s | **25.0 s (1.7×)** |
| Same model | 10 modes | 1 | 121.2 s | **51.2 s (2.4×)** |
| | | 16 | 61.5 s | **44.5 s (1.4×)** |
| Flat plate, 240 × 240 CQUAD4, 348,486 DOF | static | 1 | 45.8 s | **18.9 s (2.4×)** |
| | | 16 | 29.9 s | **15.8 s (1.9×)** |

- **Results:** the shell model's displacements agree with stock to 1e-7 (F06 precision), and all
  10 eigenvalues to 7 digits.
- **Where the time went** (shell model, static, seconds saved on 1 / 16 threads):

  | Change | 1 thread | 16 threads |
  |---|---|---|
  | Symmetric ordering in the SuperLU factorization | 46.5 | 0 |
  | STF3 reallocation in memory | 6.8 | 6.7 |
  | Matrix files between links | 6.1 | 6.3 |
  | Bulk Data reading | 2.3 | 2.3 |
  | DOF tables | 2.1 | 2.1 |
  | **Measured total** | **65.3** | **17.4** |

## 2. The changes

### 2.1 SuperLU: symmetric mode for symmetric matrices

**Before.** `SYM_MAT_DECOMP_SUPRLU` factors KLL, KOO and the shifted matrix of the eigensolvers,
and also RMM, which is not symmetric. The C driver ordered all of them with METIS on A'·A and full
partial pivoting, the setting for general unsymmetric matrices.

**Now.**
- Symmetric matrices are ordered with METIS on A'+A (`MMD_AT_PLUS_A` without METIS), in
  `SymmetricMode` with `DiagPivotThresh = 0.001`.
- It is still LU with threshold pivoting, so it is correct for any nonsingular matrix.
- RMM keeps the general path.
- SuperLU_MT keeps its general path.

On the shell model the factor went from 220.6 M to 107.0 M nonzeros.

**Checked.**
- SOL 105, Euler column: 8.224677 against 8.224670 by hand.
- Inverse power iteration: 0.1614938 against 0.1615 by hand.
- SOL 31 with SUPORT: fixed-interface modes equal the clamped SOL 103 modes.
- Static with RBE2 and RBE3: identical to stock.

A free-free Lanczos run with its shift at zero, which stopped on stock with ERROR 981 (zero pivot),
now runs and gives the right modes.

### 2.2 ESP: STF3 reallocated in memory

After the element loop, every stiffness term went to a scratch file and back one record per term,
and STF3 came back at its estimated size anyway, because LTERM was reset only after the
reallocation. It is now an in-memory copy at the right size: 11.4 s to 4.7 s on the shell model,
and STF3 from 355 MB to 179 MB.

### 2.3 Matrix files and DOF tables between links

**Before.**
- LINK1L (KGG) and the 18 LINK2 matrix files held one unformatted record per nonzero: 11.2 M
  records for KGG on the shell model.
- The DOF tables were written one integer per record.

**Now.**
- After the header, a matrix file holds a marker record (−1), then the rows, columns and values
  as one record each. READ_MATRIX_1 still reads the old form.
- Each DOF table is one record. READ_DOF_TABLES does not read the old form, which matters only
  for a restart across versions; restart decks stop with ERROR 903 on stock as well (section 5).

### 2.4 Bulk Data reading

TO_UPPER called INDEX over the alphabet for every character of every line, in both Bulk Data
passes. It now tests the character code, and the `$` is found with one INDEX. Reading went from
3.9 s to 1.6 s on the shell model.

### 2.5 Rigid elements and MPCs: three sparse kernels

Three kernels took time proportional to the matrix size squared:
- **MATMULT_SSS and MATMULT_SSS_NTERM** (KNM·GMN, GMNᵀ·KMM·GMN) visited every column of B for every
  row of A. They now multiply row by row (Gustavson's method). Every term of C is the same sum,
  added in the same order, so C is bit for bit the same.
- **SPARSE_CRS_SPARSE_CCS** looked at every term of A for each column. It is now a counting
  sort, with the same output.
- **SOLVE_GMN_SOLVER** searched RMN for each of its columns. It now puts RMN into CCS once.

On the 900-RBE2 plate, "REDUCE KGG TO KNN" went from 439.7 s to 0.7 s, and "REDUCE MGG TO MNN"
from 77.0 s to under 0.5 s.

### 2.6 GMN solve, block by block

RMM falls apart into independent blocks, since each rigid element or MPC couples only its own
dependent DOFs. The new module BLOCK_DIAGONAL_SOLVE finds the blocks and factors each block of up
to 200 DOFs with LAPACK (larger ones together with SuperLU). Each column of RMN is solved only in
the blocks it touches.
- The GMN solve now grows linearly with the number of rigid elements (section 1).
- Results change at round-off level: 68,249 numeric F06 lines of the 900-RBE3 plate agree within
  1.2e-14.
- A new F06 line gives the number of blocks and the largest.

### 2.7 Nearly singular models now stop on the SuperLU path too

The ratio of matrix diagonal to factor diagonal (PARAM MAXRATIO, BAILOUT) was checked only on the
LAPACK path. On SuperLU only an exactly zero pivot was caught, so a nearly singular model ran to
the end with very large displacements and no error.
- SYM_MAT_DECOMP_SUPRLU now gets the diagonal of U and calls BAILOUT_CHECK for KLL and KOO, as
  the LAPACK path does.
- Rows after an off-diagonal pivot are left out of the check, and the F06 says how many.
- **This changes behaviour.** Three of the 27 regression decks now stop with ERROR 982 or 983:
  - two that ran to displacements of 5E+40 and 4E+28;
  - one with a free out-of-plane mechanism under an in-plane load.

  PARAM BAILOUT −1 lets them run as before. The LAPACK path stops on all three too.
- **Test decks** in `Build_Test_Cases/statics`:
  - `mechanism_square_crod.bdf` (exact zero pivot, ERROR 981 as before);
  - `mechanism_square_crod_turned.bdf` (stock: displacements of 1.7E+12 and a normal end; now
    ERROR 982);
  - `mechanism_square_crod_weak.bdf` (stock: no message; now ERROR 983).

### 2.8 Bug fixes

- **CONM2, products of inertia and offsets.**
  - RCONM2 holds I21, I31, I32 as on the entry: products of inertia, whose negatives are the
    tensor terms. CONM2_PROC_2 moved them to the grid with the tensor sign, and the CID rotations
    rotated the products instead of the tensor. An offset with two nonzero components gave an
    indefinite mass matrix.
  - The grid point weight generator had the same sign error.
  - `conm2_mass_check.py` compares MGG with the exact matrix: stock fails all five checks, this PR
    passes them.
- **CONM2 with CID = −1** (c.g. and inertias in basic) stopped with ERROR 1822; it is now read.
- **RBE3 rank check.** The rank of the retained REFC system was judged at machine precision, so
  identical RBE3s passed or failed ERROR 1957 depending on their orientation. It uses a relative
  tolerance of 1.0E-10 now. On the 24 rotated copies of one singular RBE3 in
  `rbe3_rank_singular.bdf`, stock rejects 10 and this PR all 24.
- **CTETRA with 10 grids listed clockwise** was rejected (all 195,924 TETRA10s of NASA's HIRENASD
  model, for example). Such elements are reordered once in LINK1. Test decks:
  `tet10_orientation_rh.bdf` and `tet10_orientation_lh.bdf`.
- **Adaptive Lanczos.**
  - The NEV cap left no room for ARPACK's Krylov space, so a search that reached it stopped with
    ERROR 9776. The cap is now (DOFs with mass) − 4.
  - When the search stops at its limit before covering the requested range, the F06 now says so
    (WARNING 9105). Test decks: `Build_Test_Cases/modes/chain8_*.bdf`, an 8-mass spring chain.
- **Bulk Data entry names** were compared by prefix, so FORCE1 was read as FORCE, PBEAML as PBEAM,
  PCOMPG as PCOMP and SUPORT1 as SUPORT, and MYSTRAN's own PCOMP1 was never reached.
  - Names are now compared exactly. An entry whose name only starts with one MYSTRAN reads stops
    with ERROR 1705.
  - FORCE1, FORCE2, MOMENT1 and MOMENT2 are now read.
  - After one bad FORCE or MOMENT, every later one was dropped (a SAVEd error count); fixed.
  - Test: `force12_check.py`.
- **Every stop for fatal errors** now ends the F06 with "*** PROCESSING TERMINATED: n FATAL
  ERROR(S), LISTED ABOVE".
- **Small fixes.**
  - CMake: SuperLU's docs and tests were not really switched off, a METIS include path broke
    builds in paths with spaces, two `.f03` files with `#ifdef` were not preprocessed, and running
    CMake again switched EMBEDDED BLAS to SYSTEM.
  - The OUTPUT4 trailer record (pyNastran could not read MYSTRAN OUTPUT4).
  - The SOL 31 METHOD check never fired.
  - A debug print to the screen.
  - A PARAM POST message.
  - BAILOUT_CHECK looked up a grid name for every row, even when it printed nothing.

## 3. Test setup

- **Machine:** Windows 11, Intel i9-14900HX (8 performance and 16 efficiency cores), 64 GB.
  Identical runs vary by about ±10%.
- **Toolchain:** MSYS2 MinGW64, gfortran 15.1, CMake 4.0. Release builds with embedded
  Reference-LAPACK and BLAS.

## 4. How results were compared

- **Regression set:** 27 decks.
  - `Build_Test_Cases/statics` and `/buckling`;
  - 16 more: shells, bars, RBE3, PCOMP, MITC4+, CQUAD8; SOL 101, 103 and 105; 100 to 6,000
    grids.
- **Tools** (`Build_Test_Cases/scripts`):
  - `cmpdisp.py` compares every grid displacement of two F06 files.
  - `cmpf06.py` compares every numeric line, each number scaled by the largest magnitude within
    40 lines.
- **Outcome against stock:**
  - Every deck that prints displacements gives the same displacements, or within 1e-13.
  - New information lines: the MAXRATIO check and the RMM blocks.
  - Three decks now stop (section 2.7).
  - "Occurs at grid" tie-breaks may differ for equal maxima.

## 5. Found and not fixed

1. **SuperLU_MT symmetric mode crashes** on a 15-RBE deck; SuperLU_MT keeps its general path.
2. **One record per nonzero** is still used by LINK1E, LINK1R, LINK1J and the scratch files read by
   READ_MATRIX_2.
3. **One right-hand side per SuperLU solve** inside column loops (SOLVE_GOA, SOLVE_DLR,
   SOLVE_PHIZL1). Blocking the columns would let `dgstrs` use BLAS-3.
4. **HIRENASD** (NASA, 1,064,802 DOF) reads and assembles with these changes, once its three
   CTETRAs with only some midside grids are made 4-node (MYSTRAN takes 4 or 10 grids). It then
   fails in SuperLU (ERROR 981, memory allocation, with memory still free), which looks like
   SuperLU's 32-bit index limit on factor terms.
5. **Restart decks stop with ERROR 903**, because the file names are built from RESTART_FILNAM
   before it is set.
6. **Executive Control is case sensitive:** a lowercase `sol 101` / `cend` stops with ERROR 1011.

## 6. Next steps

- **Threading:** only the factorization is threaded. On 16 threads most of the remaining time is
  single-threaded work; OpenMP over the element loop is the largest candidate.
- **Lanczos:** several right-hand sides per solve, or a threaded triangular solve.
- **SuperLU with 64-bit indices**, for models above about 1M DOF.

## Appendix A: commits (oldest first)

| Commit | Subject |
|---|---|
| 9dac8448 | SuperLU: factor symmetric matrices in symmetric mode (METIS on A'+A) |
| ee021aa6 | CMake: turn SuperLU examples/tests/docs off for real; build in paths with spaces |
| 95b332e6 | LOADC: require METHOD for SOL 31 (GEN CB MODEL) as intended |
| b703b4f0 | OP2_GRID_OUTPUT: remove leftover debug print to the screen |
| 445289a6 | OUTPUT4: write one value in the matrix trailer record |
| 143c0c48 | BD_PARAM: POST warning gives the same allowed values in the F06 as in the ERR |
| d2a25691 | RBE3: judge the rank of the retained REFC system with a relative tolerance |
| d1fab977 | Bulk Data reading: TO_UPPER by character code; one INDEX for the $ comment |
| a3337da5 | ESP: reallocate STF3 at NTERM through an in-memory copy, not a scratch file |
| 68d90b5e | CMake: preprocess STARTUP_SPECIALS.f03 and PARAM_CARDS.f03 |
| 615103dc | CMake: keep the BLAS+LAPACK mode when CMake is run again |
| 8cf56c5d | SPARSE_CRS_SPARSE_CCS: counting sort instead of a search of A for each column |
| 14c8d0e5 | MATMULT_SSS, MATMULT_SSS_NTERM: multiply row by row instead of row by column |
| 82dd4ef6 | SOLVE_GMN_SOLVER: take the columns of RMN from a CCS copy |
| 6bacd1aa | BAILOUT_CHECK: look up the grid and component only for a message |
| 110f4bff | SuperLU: check the ratio of matrix diagonal to factor diagonal (MAXRATIO) |
| a920d935 | Build_Test_Cases/scripts: cmpf06.py and cmpdisp.py, F06 comparison used for the regression checks |
| 4500ae43 | DOF tables (LINK1C): one record per table instead of one per entry |
| da18821c | Matrix files between links: three records instead of one per term |
| abdf8daa | SOLVE_GMN_SOLVER: solve RMM block by block (new module BLOCK_DIAGONAL_SOLVE) |
| 72e27417 | CTETRA with 10 grids: accept either corner order |
| d91ab7ed | CONM2: products of inertia and offsets give the right mass matrix |
| 11b3e1e7 | CONM2 with CID = -1: the c.g. coordinates and the inertias are in basic |
| 8799a0d2 | Adaptive Lanczos: cap NEV where ARPACK can run |
| 447bedb4 | Adaptive Lanczos: WARNING 9105 in the F06 when the search stops at its limit short of the range |
| b2f158c6 | Bulk Data names compared exactly; FORCE1, FORCE2, MOMENT1, MOMENT2 |
| 56810791 | OUTA_HERE: end the F06 with a PROCESSING TERMINATED line on every stop for fatal errors |
| d0ba93a7 | Build_Test_Cases/statics: make_spider_plate_decks.py, the rigid-element timing decks |

## Appendix B: reproducing

- **Build** in a path without spaces, from an MSYS2 MinGW64 shell:

  ```
  cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DMYSTRAN_BLAS_LAPACK=EMBEDDED -S . -B build
  cmake --build build -j8
  ```

  Add `-DUSE_SUPERLU_MT=ON` for SuperLU_MT, and set `OMP_NUM_THREADS`.
- **Test decks and checks:**

  ```
  python Build_Test_Cases/statics/make_spider_plate_decks.py
  python Build_Test_Cases/statics/make_rbe3_rank_decks.py
  python Build_Test_Cases/statics/conm2_mass_check.py --mystran <exe>
  python Build_Test_Cases/statics/force12_check.py --mystran <exe>
  ```

- **Compare two runs:**

  ```
  python Build_Test_Cases/scripts/cmpf06.py a.F06 b.F06
  python Build_Test_Cases/scripts/cmpdisp.py a.F06 b.F06
  ```
