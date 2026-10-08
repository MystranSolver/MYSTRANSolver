
/*
 * -- SuperLU routine (version 6.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * March 26, 2023  Add 64-bit indexing and METIS ordering
 */

#include "slu_ddefs.h"

#define HANDLE_SIZE 8

/* kind of integer to hold a pointer.  Use 64-bit. */
typedef long long int fptr;

typedef struct {
  SuperMatrix *L;
  SuperMatrix *U;
  int *perm_c;
  int *perm_r;
} factors_t;

/*!
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs triangular solve
 *      = 3, free all the storage in the end
 *
 * f_factors (input/output) fptr*
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 */
/*
 * Matrix kind for the next factorization (iopt == 1), set by the Fortran
 * caller through c_fortran_dgssv_symmetric_:
 *   0 = general matrix: column ordering on A'*A, partial pivoting
 *       (diag_pivot_thresh = 1.0). This is the default.
 *   1 = symmetric matrix (stiffness, shifted stiffness): ordering on A'+A,
 *       SuperLU's symmetric mode and threshold pivoting that prefers the
 *       diagonal (DiagPivotThresh = 0.001). Off-diagonal pivots are still
 *       taken when a diagonal is too small, so the factorization stays
 *       correct for any nonsingular matrix; the fill is about half.
 */
static int slu_symmetric = 0;

void c_fortran_dgssv_symmetric_(int *flag) { slu_symmetric = (*flag != 0); }

/*
 * Diagonal of U of the factored matrix in f_factors, by column of A:
 * udiag[i] is the pivot of the column that is column i+1 of A.
 * exact_d[i] = 1 if that pivot is D(i) of A = L*D*L' for a symmetric A: the
 * pivot row is A's diagonal (row i+1) and no off-diagonal pivot reaches it.
 * An off-diagonal pivot at one position changes every later position it
 * updates (U(k,j) != 0), so those are marked 0 through the structure of U.
 * *available = 1 (0 from drivers that do not provide it).
 */
void c_fortran_dgssv_udiag_(fptr *f_factors, int *n, double *udiag,
                            int *exact_d, int *available) {
  factors_t *LUfactors = (factors_t *)*f_factors;
  SCformat *Lstore = (SCformat *)LUfactors->L->Store;
  NCformat *Ustore = (NCformat *)LUfactors->U->Store;
  double *Lval = (double *)Lstore->nzval;
  int *perm_c = LUfactors->perm_c;
  int *perm_r = LUfactors->perm_r;
  double *dpos;
  int *iperm_c, *clean;
  int_t p;
  int i, j, k, ok, fsupc;

  if (!(dpos = doubleMalloc(*n)))
    ABORT("Malloc fails for dpos[].");
  if (!(iperm_c = int32Malloc(*n)) || !(clean = int32Malloc(*n)))
    ABORT("Malloc fails for iperm_c[] or clean[].");
  for (i = 0; i < *n; ++i)
    iperm_c[perm_c[i]] = i;
  /* The diagonal of U is in the diagonal block of each supernode of L; the
     rest of U is in Ustore by column, with rows as pivot positions. */
  for (k = 0; k <= Lstore->nsuper; ++k) {
    fsupc = L_FST_SUPC(k);
    for (j = fsupc; j < L_FST_SUPC(k + 1); ++j) {
      dpos[j] = Lval[L_NZ_START(j) + (j - fsupc)];
      ok = (perm_r[iperm_c[j]] == j);
      for (i = fsupc; i < j && ok; ++i)
        ok = clean[i];
      for (p = U_NZ_START(j); p < U_NZ_START(j + 1) && ok; ++p)
        ok = clean[U_SUB(p)];
      clean[j] = ok;
    }
  }
  for (i = 0; i < *n; ++i) {
    udiag[i] = dpos[perm_c[i]];
    exact_d[i] = clean[perm_c[i]];
  }
  SUPERLU_FREE(dpos);
  SUPERLU_FREE(iperm_c);
  SUPERLU_FREE(clean);
  *available = 1;
}

void c_fortran_dgssv_(int *iopt, int *n, int_t *nnz, int *nrhs, double *values,
                      int_t *rowind, int_t *colptr, double *b, int *ldb,
                      fptr *f_factors, /* a handle containing the address
                                          pointing to the factored matrices */
                      int_t *info) {
  SuperMatrix A, AC, B;
  SuperMatrix *L, *U;
  int *perm_r; /* row permutations from partial pivoting */
  int *perm_c; /* column permutation vector */
  int *etree;  /* column elimination tree */
  SCformat *Lstore;
  NCformat *Ustore;
  int i, panel_size, permc_spec, relax;
  trans_t trans;
  mem_usage_t mem_usage;
  superlu_options_t options;
  SuperLUStat_t stat;
  factors_t *LUfactors;
  GlobalLU_t Glu; /* Not needed on return. */
  int_t *rowind0; /* counter 1-based indexing from Fortran arrays. */
  int_t *colptr0;

  trans = NOTRANS;

  if (*iopt == 1) { /* LU decomposition */

    /* Set the default input options. */
    set_default_options(&options);

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Adjust to 0-based indexing */
    if (!(rowind0 = intMalloc(*nnz)))
      ABORT("Malloc fails for rowind0[].");
    if (!(colptr0 = intMalloc(*n + 1)))
      ABORT("Malloc fails for colptr0[].");
    for (i = 0; i < *nnz; ++i)
      rowind0[i] = rowind[i] - 1;
    for (i = 0; i <= *n; ++i)
      colptr0[i] = colptr[i] - 1;

    dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind0, colptr0, SLU_NC,
                           SLU_D, SLU_GE);
    L = (SuperMatrix *)SUPERLU_MALLOC(sizeof(SuperMatrix));
    U = (SuperMatrix *)SUPERLU_MALLOC(sizeof(SuperMatrix));
    if (!(perm_r = int32Malloc(*n)))
      ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = int32Malloc(*n)))
      ABORT("Malloc fails for perm_c[].");
    if (!(etree = int32Malloc(*n)))
      ABORT("Malloc fails for etree[].");
/*
 * Get column permutation vector perm_c[], according to permc_spec:
 *   permc_spec = 0: natural ordering
 *   permc_spec = 1: minimum degree on structure of A'*A
 *   permc_spec = 2: minimum degree on structure of A'+A
 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
 *   permc_spec = 6: METIS ordering on structure of A'*A
 *   METIS_AT_PLUS_A: METIS ordering on structure of A'+A
 */
    if (slu_symmetric) {
      options.SymmetricMode = YES;
      options.DiagPivotThresh = 0.001;
    }
#if (HAVE_METIS)
    printf(slu_symmetric ? "USING METIS ORDERING (SYMMETRIC)\r\n"
                         : "USING METIS ORDERING\r\n");
    permc_spec = slu_symmetric ? METIS_AT_PLUS_A : METIS_ATA;
#else
    permc_spec = slu_symmetric ? MMD_AT_PLUS_A : options.ColPerm;
#endif
    //	permc_spec = 0;
    // printf("before get_perm_c: permc_spec %d, *n %d\n", permc_spec, *n);
    get_perm_c(permc_spec, &A, perm_c);
    // printf("after get_perm_c: permc_spec %d\n", permc_spec);

    sp_preorder(&options, &A, perm_c, etree, &AC);

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    dgstrf(&options, &AC, relax, panel_size, etree, NULL, 0, perm_c, perm_r, L,
           U, &Glu, &stat, info);

    if (*info == 0) {
      Lstore = (SCformat *)L->Store;
      Ustore = (NCformat *)U->Store;
      printf("No of nonzeros in factor L = %lld\n", (long long)Lstore->nnz);
      printf("No of nonzeros in factor U = %lld\n", (long long)Ustore->nnz);
      printf("No of nonzeros in L+U = %lld\n",
             (long long)Lstore->nnz + Ustore->nnz);
      dQuerySpace(L, U, &mem_usage);
      printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu / 1e6,
             mem_usage.total_needed / 1e6);
    } else {
      printf("dgstrf() error returns INFO= %lld\n", (long long)*info);
      if (*info <= *n) { /* factorization completes */
        dQuerySpace(L, U, &mem_usage);
        printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu / 1e6,
               mem_usage.total_needed / 1e6);
      }
    }

    /* Save the LU factors in the factors handle */
    LUfactors = (factors_t *)SUPERLU_MALLOC(sizeof(factors_t));
    LUfactors->L = L;
    LUfactors->U = U;
    LUfactors->perm_c = perm_c;
    LUfactors->perm_r = perm_r;
    *f_factors = (fptr)LUfactors;

    /* Free un-wanted storage */
    SUPERLU_FREE(etree);
    Destroy_SuperMatrix_Store(&A);
    Destroy_CompCol_Permuted(&AC);
    SUPERLU_FREE(rowind0);
    SUPERLU_FREE(colptr0);
    StatFree(&stat);

  } else if (*iopt == 2) { /* Triangular solve */
    int iinfo;

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Extract the LU factors in the factors handle */
    LUfactors = (factors_t *)*f_factors;
    L = LUfactors->L;
    U = LUfactors->U;
    perm_c = LUfactors->perm_c;
    perm_r = LUfactors->perm_r;

    dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);

    /* Solve the system A*X=B, overwriting B with X. */
    dgstrs(trans, L, U, perm_c, perm_r, &B, &stat, &iinfo);
    *info = iinfo;

    Destroy_SuperMatrix_Store(&B);
    StatFree(&stat);

  } else if (*iopt == 3) { /* Free storage */
    /* Free the LU factors in the factors handle */
    LUfactors = (factors_t *)*f_factors;
    SUPERLU_FREE(LUfactors->perm_r);
    SUPERLU_FREE(LUfactors->perm_c);
    Destroy_SuperNode_Matrix(LUfactors->L);
    Destroy_CompCol_Matrix(LUfactors->U);
    SUPERLU_FREE(LUfactors->L);
    SUPERLU_FREE(LUfactors->U);
    SUPERLU_FREE(LUfactors);
  } else {
    fprintf(stderr, "Invalid iopt=%d passed to c_fortran_dgssv()\n", *iopt);
    exit(-1);
  }
}
