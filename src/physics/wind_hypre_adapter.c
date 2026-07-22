/* Exact 15-point HICAR wind operator -> HYPRE ParCSR adapter.
 *
 * This deliberately assembles and solves the calibrated operator on the host.
 * It is a correctness baseline for HICAR's NCCL layout, whose CPU-only I/O
 * rank requires GPU-aware MPICH to be disabled.  HICAR's OpenACC vectors
 * remain the source of truth and the native true-residual check remains
 * mandatory after every HYPRE solve.
 */
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_utilities.h"

typedef struct {
  int xs, ys, zs, xm, ym, zm;
  HYPRE_BigInt first;
} hicar_block_t;

typedef struct {
  MPI_Comm comm;
  int nprocs, rank;
  int xs, ys, zs, xm, ym, zm, mx, my, mz;
  int ci_s, ci_e, ck_s, ck_e, cj_s, cj_e;
  HYPRE_BigInt first, last;
  hicar_block_t *blocks;
  unsigned char *zero_rows;
  HYPRE_Int *mgr_point_marker;
  HYPRE_IJMatrix ij_a;
  HYPRE_ParCSRMatrix a;
  HYPRE_Solver mgr;
  HYPRE_Solver amg;
  HYPRE_Solver fgmres;
  int solver_ready;
  int valid;
} hicar_hypre_t;

static hicar_hypre_t state = {0};
static int hypre_initialized = 0;

int hicar_hypre_initialize(void)
{
  int ierr = 0;
  if (hypre_initialized) return 0;
  fprintf(stderr, "HICAR HYPRE: initializing runtime\n"); fflush(stderr);
  ierr |= HYPRE_Initialize();
  fprintf(stderr, "HICAR HYPRE: HYPRE_Initialize rc=%d\n", ierr); fflush(stderr);
  /* HICAR's NCCL topology disables GPU-aware MPICH because its per-node I/O
   * rank has no CUDA context.  HYPRE device vectors would nevertheless use
   * MPI device transfers during AMG/FGMRES setup and crash on this stack.
   * The bridge already stages HICAR fields on the host, so keep the complete
   * HYPRE solve on host memory until a compatible all-GPU MPI topology exists. */
  if (!ierr) ierr |= HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
  fprintf(stderr, "HICAR HYPRE: set host memory rc=%d\n", ierr); fflush(stderr);
  if (!ierr) ierr |= HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
  fprintf(stderr, "HICAR HYPRE: set host execution rc=%d\n", ierr); fflush(stderr);
  if (!ierr) hypre_initialized = 1;
  return ierr;
}

static void clear_state(void)
{
  if (state.fgmres) HYPRE_ParCSRFlexGMRESDestroy(state.fgmres);
  if (state.mgr) HYPRE_MGRDestroy(state.mgr);
  if (state.amg) HYPRE_BoomerAMGDestroy(state.amg);
  if (state.ij_a) HYPRE_IJMatrixDestroy(state.ij_a);
  free(state.blocks);
  free(state.zero_rows);
  free(state.mgr_point_marker);
  memset(&state, 0, sizeof(state));
}

void hicar_hypre_destroy(void) { clear_state(); }

static HYPRE_BigInt local_index(const hicar_block_t *b, int i, int k, int j)
{
  return b->first + ((HYPRE_BigInt)(j - b->ys) * b->zm + (k - b->zs)) * b->xm + (i - b->xs);
}

static HYPRE_BigInt global_index(int i, int k, int j)
{
  int r;
  for (r = 0; r < state.nprocs; ++r) {
    const hicar_block_t *b = &state.blocks[r];
    if (i >= b->xs && i < b->xs + b->xm &&
        j >= b->ys && j < b->ys + b->ym &&
        k >= b->zs && k < b->zs + b->zm)
      return local_index(b, i, k, j);
  }
  return -1;
}

static int hicar_row_is_identity(int i, int k, int j)
{
  return i == 0 || i == state.mx - 1 || j == 0 || j == state.my - 1 ||
         k == 0 || k == state.mz - 1;
}

/* Coefficient storage is Fortran contiguous: i is the fastest index. */
static size_t coeff_index(int i, int k, int j)
{
  return ((size_t)(j - state.cj_s) * (size_t)(state.ck_e - state.ck_s + 1) +
          (size_t)(k - state.ck_s)) * (size_t)(state.ci_e - state.ci_s + 1) +
         (size_t)(i - state.ci_s);
}

int hicar_hypre_build(MPI_Fint comm_f, int xs, int ys, int zs, int xm, int ym, int zm,
                      int mx, int my, int mz, int ci_s, int ci_e, int ck_s, int ck_e, int cj_s, int cj_e,
                      const float *a, const float *b, const float *c, const float *d,
                      const float *e, const float *f, const float *g, const float *h,
                      const float *ii, const float *jj, const float *kcoef, const float *l,
                      const float *m, const float *n, const float *o)
{
  int meta[6], *allmeta = NULL, r, ierr = 0, fgmres_print_level = 0;
  HYPRE_Int mgr_num_cpoints[1] = {1}, mgr_cpoint = 0;
  HYPRE_Int *mgr_cpoints[1] = {&mgr_cpoint};
  int local_bad_diag = 0, global_bad_diag = 0;
  int local_zero_rows = 0, global_zero_rows = 0;
  int local_zero_diag_offdiag = 0, global_zero_diag_offdiag = 0;
  double local_min_abs_diag = HUGE_VAL, global_min_abs_diag = 0.0;
  HYPRE_BigInt *alloffsets = NULL;
  HYPRE_BigInt nlocal = (HYPRE_BigInt)xm * ym * zm;
  HYPRE_BigInt offset = 0;
  HYPRE_Int *ncols = NULL;
  HYPRE_BigInt *rows = NULL, *cols = NULL;
  HYPRE_Real *vals = NULL;
  HYPRE_BigInt q = 0, entry = 0;
  int li, lk, lj;
  const char *print_level_env;
  /* Match wind_iterative::spmv exactly: A, B, C, D, E, F, G, H, I,
   * J, K, L, M, N, O.  The final four terms are mixed j/k neighbors,
   * not i/j/k corner neighbors. */
  const int di[15] = {0,0,0, 1,-1,0,0, 1,-1, 1,-1, 0,0,0,0};
  const int dk[15] = {0,1,-1, 0,0,0,0, 1, 1,-1,-1, 1,1,-1,-1};
  const int dj[15] = {0,0,0, 0,0,1,-1, 0, 0, 0, 0, 1,-1,1,-1};

  clear_state();
  fprintf(stderr, "HICAR HYPRE: entering matrix build\n"); fflush(stderr);
  if (hicar_hypre_initialize()) return HYPRE_GetError();
  state.comm = MPI_Comm_f2c(comm_f);
  if (state.comm == MPI_COMM_NULL) return -1;
  MPI_Comm_rank(state.comm, &state.rank);
  MPI_Comm_size(state.comm, &state.nprocs);
  fprintf(stderr, "HICAR HYPRE rank %d/%d: communicator ready\n", state.rank, state.nprocs); fflush(stderr);
  state.xs=xs; state.ys=ys; state.zs=zs; state.xm=xm; state.ym=ym; state.zm=zm;
  state.mx=mx; state.my=my; state.mz=mz;
  state.ci_s=ci_s; state.ci_e=ci_e; state.ck_s=ck_s; state.ck_e=ck_e; state.cj_s=cj_s; state.cj_e=cj_e;
  MPI_Exscan(&nlocal, &offset, 1, HYPRE_MPI_BIG_INT, MPI_SUM, state.comm);
  if (state.rank == 0) offset = 0;
  state.first = offset; state.last = offset + nlocal - 1;

  state.blocks = calloc((size_t)state.nprocs, sizeof(*state.blocks));
  state.zero_rows = calloc((size_t)nlocal, sizeof(*state.zero_rows));
  allmeta = calloc((size_t)state.nprocs * 6, sizeof(*allmeta));
  alloffsets = calloc((size_t)state.nprocs, sizeof(*alloffsets));
  if (!state.blocks || !state.zero_rows || !allmeta || !alloffsets) { clear_state(); free(allmeta); free(alloffsets); return -2; }
  meta[0]=xs; meta[1]=ys; meta[2]=zs; meta[3]=xm; meta[4]=ym; meta[5]=zm;
  MPI_Allgather(meta, 6, MPI_INT, allmeta, 6, MPI_INT, state.comm);
  for (r=0; r<state.nprocs; ++r) {
    hicar_block_t *blk=&state.blocks[r];
    blk->xs=allmeta[6*r]; blk->ys=allmeta[6*r+1]; blk->zs=allmeta[6*r+2];
    blk->xm=allmeta[6*r+3]; blk->ym=allmeta[6*r+4]; blk->zm=allmeta[6*r+5];
    blk->first=0;
  }
  MPI_Allgather(&offset, 1, HYPRE_MPI_BIG_INT, alloffsets, 1, HYPRE_MPI_BIG_INT, state.comm);
  for (r=0; r<state.nprocs; ++r) state.blocks[r].first=alloffsets[r];
  free(allmeta);
  free(alloffsets);

  ncols = calloc((size_t)nlocal, sizeof(*ncols));
  rows  = calloc((size_t)nlocal, sizeof(*rows));
  /* 15 entries for every row is a safe upper bound; identity rows use one. */
  cols  = calloc((size_t)nlocal * 15, sizeof(*cols));
  vals  = calloc((size_t)nlocal * 15, sizeof(*vals));
  if (!ncols || !rows || !cols || !vals) { free(ncols); free(rows); free(cols); free(vals); clear_state(); return -3; }

  for (lj=0; lj<ym; ++lj) for (lk=0; lk<zm; ++lk) for (li=0; li<xm; ++li) {
    int i=xs+li, k=zs+lk, j=ys+lj, s;
    size_t p = 0;
    HYPRE_Real row_abs_sum = 0.0;
    float cv[15] = {0};
    HYPRE_BigInt row = local_index(&state.blocks[state.rank], i,k,j);
    if (!hicar_row_is_identity(i,k,j)) {
      if (i < state.ci_s || i > state.ci_e || k < state.ck_s || k > state.ck_e ||
          j < state.cj_s || j > state.cj_e) { ierr = -5; break; }
      p=coeff_index(i,k,j);
      cv[0]=a[p]; cv[1]=b[p]; cv[2]=c[p]; cv[3]=d[p]; cv[4]=e[p]; cv[5]=f[p]; cv[6]=g[p];
      cv[7]=h[p]; cv[8]=ii[p]; cv[9]=jj[p]; cv[10]=kcoef[p]; cv[11]=l[p]; cv[12]=m[p]; cv[13]=n[p]; cv[14]=o[p];
      for (s=0; s<15; ++s) row_abs_sum += fabsf(cv[s]);
      if (row_abs_sum == 0.0) { state.zero_rows[q] = 1; ++local_zero_rows; }
      if (!isfinite(cv[0]) || fabsf(cv[0]) <= 1.0e-30f) ++local_bad_diag;
      if (fabsf(cv[0]) <= 1.0e-30f && row_abs_sum > 0.0) ++local_zero_diag_offdiag;
      if (isfinite(cv[0]) && fabsf(cv[0]) < local_min_abs_diag) local_min_abs_diag = fabsf(cv[0]);
    }
    rows[q]=row;
    if (hicar_row_is_identity(i,k,j) || state.zero_rows[q]) {
      ncols[q]=1; cols[entry]=row; vals[entry++]=1.0;
    } else {
      ncols[q]=15;
      for (s=0; s<15; ++s) {
        HYPRE_BigInt col=global_index(i+di[s], k+dk[s], j+dj[s]);
        if (col < 0) { ierr=-4; break; }
        cols[entry]=col; vals[entry]=(HYPRE_Real)cv[s]; ++entry;
      }
    }
    if (ierr) break;
    ++q;
  }
  MPI_Allreduce(&local_bad_diag, &global_bad_diag, 1, MPI_INT, MPI_SUM, state.comm);
  MPI_Allreduce(&local_zero_rows, &global_zero_rows, 1, MPI_INT, MPI_SUM, state.comm);
  MPI_Allreduce(&local_zero_diag_offdiag, &global_zero_diag_offdiag, 1, MPI_INT, MPI_SUM, state.comm);
  MPI_Allreduce(&local_min_abs_diag, &global_min_abs_diag, 1, MPI_DOUBLE, MPI_MIN, state.comm);
  if (state.rank == 0) {
    fprintf(stderr, "HICAR HYPRE matrix audit: bad_diagonal_rows=%d zero_rows=%d zero_diag_with_offdiag=%d min_abs_diagonal=%g\n",
            global_bad_diag, global_zero_rows, global_zero_diag_offdiag, global_min_abs_diag); fflush(stderr);
  }
  if (!ierr) ierr |= HYPRE_IJMatrixCreate(state.comm, state.first, state.last, state.first, state.last, &state.ij_a);
  if (!ierr) ierr |= HYPRE_IJMatrixSetObjectType(state.ij_a, HYPRE_PARCSR);
  if (!ierr) ierr |= HYPRE_IJMatrixInitialize_v2(state.ij_a, HYPRE_MEMORY_HOST);
  if (!ierr) ierr |= HYPRE_IJMatrixSetValues(state.ij_a, (HYPRE_Int)nlocal, ncols, rows, cols, vals);
  if (!ierr) ierr |= HYPRE_IJMatrixAssemble(state.ij_a);
  if (!ierr) { fprintf(stderr, "HICAR HYPRE rank %d: matrix assembled\n", state.rank); fflush(stderr); }
  if (!ierr) ierr |= HYPRE_IJMatrixGetObject(state.ij_a, (void **)&state.a);
  if (!ierr) { fprintf(stderr, "HICAR HYPRE rank %d: matrix retained on host\n", state.rank); fflush(stderr); }
  /* Use an explicit 2x2x2 spatial C/F split before algebraic coarsening.
   * This supplies a distributed geometric coarse space to the nonsymmetric
   * projection operator, instead of asking a local smoother to represent its
   * long horizontal modes.  Markers use HYPRE's scalar block id 0 for C and
   * -1 for F, and are retained until MGR is destroyed. */
  state.mgr_point_marker = calloc((size_t)nlocal, sizeof(*state.mgr_point_marker));
  if (!state.mgr_point_marker) ierr = -6;
  if (!ierr) {
    for (lj=0; lj<ym; ++lj) for (lk=0; lk<zm; ++lk) for (li=0; li<xm; ++li) {
      int i=xs+li, k=zs+lk, j=ys+lj;
      HYPRE_BigInt qmark=((HYPRE_BigInt)lj*zm + lk)*xm + li;
      state.mgr_point_marker[qmark] = ((i & 1) == 0 && (j & 1) == 0 && (k & 1) == 0) ? 0 : -1;
    }
  }
  if (!ierr) ierr |= HYPRE_MGRCreate(&state.mgr);
  if (!ierr) ierr |= HYPRE_MGRSetCpointsByPointMarkerArray(state.mgr, 1, 1,
      mgr_num_cpoints, mgr_cpoints, state.mgr_point_marker);
  if (!ierr) ierr |= HYPRE_MGRSetNonCpointsToFpoints(state.mgr, 1);
  if (!ierr) ierr |= HYPRE_MGRSetMaxCoarseLevels(state.mgr, 1);
  if (!ierr) ierr |= HYPRE_MGRSetRelaxType(state.mgr, 18);
  if (!ierr) ierr |= HYPRE_MGRSetFRelaxMethod(state.mgr, 0);
  if (!ierr) ierr |= HYPRE_MGRSetNumRelaxSweeps(state.mgr, 1);
  if (!ierr) ierr |= HYPRE_MGRSetRestrictType(state.mgr, 3);
  if (!ierr) ierr |= HYPRE_MGRSetInterpType(state.mgr, 3);
  if (!ierr) ierr |= HYPRE_MGRSetMaxIter(state.mgr, 1);
  if (!ierr) ierr |= HYPRE_MGRSetTol(state.mgr, 0.0);
  if (!ierr) ierr |= HYPRE_MGRSetPrintLevel(state.mgr, 0);
  if (!ierr) ierr |= HYPRE_BoomerAMGCreate(&state.amg);
  /* The terrain-following projection operator is generally nonsymmetric.
   * Use AMG for one flexible preconditioning cycle, and let FGMRES own the
   * outer convergence test rather than applying AMG as a stationary solver. */
  if (!ierr) ierr |= HYPRE_BoomerAMGSetTol(state.amg, 0.0);
  if (!ierr) ierr |= HYPRE_BoomerAMGSetMaxIter(state.amg, 1);
  if (!ierr) ierr |= HYPRE_BoomerAMGSetPrintLevel(state.amg, 0);
  /* The default Falgout hierarchy is prohibitively expensive on the 3.75 M
   * cell Swiss operator.  HMIS plus extended+i interpolation is HYPRE's
   * scalable nonsymmetric-3D policy; bound interpolation density to keep
   * coarse operators and setup memory controlled. */
  if (!ierr) ierr |= HYPRE_BoomerAMGSetCoarsenType(state.amg, 10);
  if (!ierr) ierr |= HYPRE_BoomerAMGSetInterpType(state.amg, 6);
  if (!ierr) ierr |= HYPRE_BoomerAMGSetPMaxElmts(state.amg, 4);
  if (!ierr) ierr |= HYPRE_BoomerAMGSetStrongThreshold(state.amg, 0.25);
  if (!ierr) ierr |= HYPRE_MGRSetCoarseSolver(state.mgr,
      HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, state.amg);
  if (!ierr) ierr |= HYPRE_ParCSRFlexGMRESCreate(state.comm, &state.fgmres);
  if (!ierr) ierr |= HYPRE_ParCSRFlexGMRESSetKDim(state.fgmres, 50);
  if (!ierr) ierr |= HYPRE_ParCSRFlexGMRESSetTol(state.fgmres, 1.0e-5);
  if (!ierr) ierr |= HYPRE_ParCSRFlexGMRESSetMaxIter(state.fgmres, 1000);
  print_level_env = getenv("HICAR_HYPRE_PRINT_LEVEL");
  if (print_level_env) fgmres_print_level = atoi(print_level_env);
  if (fgmres_print_level < 0) fgmres_print_level = 0;
  if (fgmres_print_level > 2) fgmres_print_level = 2;
  if (!ierr) ierr |= HYPRE_ParCSRFlexGMRESSetLogging(state.fgmres, fgmres_print_level ? 1 : 0);
  if (!ierr) ierr |= HYPRE_ParCSRFlexGMRESSetPrintLevel(state.fgmres, fgmres_print_level);
  /* The calibrated matrix has now been assembled with its true local
   * coefficient bounds.  Apply MGR's explicit spatial coarse correction and
   * its HMIS AMG coarse solve through flexible FGMRES. */
  if (!ierr) ierr |= HYPRE_ParCSRFlexGMRESSetPrecond(state.fgmres,
      HYPRE_MGRSolve, HYPRE_MGRSetup, state.mgr);
  free(ncols); free(rows); free(cols); free(vals);
  if (ierr) { clear_state(); return ierr; }
  state.valid=1;
  fprintf(stderr, "HICAR HYPRE rank %d: matrix build complete\n", state.rank); fflush(stderr);
  return 0;
}

int hicar_hypre_solve(const double *rhs, double *x, int max_iter, double tol, int *iterations, double *residual)
{
  HYPRE_IJVector ij_b=NULL, ij_x=NULL;
  HYPRE_ParVector b=NULL, sol=NULL;
  HYPRE_BigInt nlocal=state.last-state.first+1, q;
  HYPRE_BigInt *ids=NULL;
  HYPRE_Real *bvals=NULL, *xvals=NULL;
  int ierr=0, stage=0, local_gauge_rhs=0, global_gauge_rhs=0;
  double setup_start;
  if (!state.valid) return -10;
  HYPRE_ClearAllErrors();
  ids=malloc((size_t)nlocal*sizeof(*ids)); bvals=malloc((size_t)nlocal*sizeof(*bvals)); xvals=malloc((size_t)nlocal*sizeof(*xvals));
  if (!ids || !bvals || !xvals) { free(ids); free(bvals); free(xvals); return -11; }
  for(q=0;q<nlocal;++q){
    ids[q]=state.first+q;
    if (state.zero_rows[q]) {
      if (fabs(rhs[q]) > 1.0e-12) ++local_gauge_rhs;
      bvals[q]=0.0; xvals[q]=0.0;
    } else {
      bvals[q]=rhs[q]; xvals[q]=x[q];
    }
  }
  MPI_Allreduce(&local_gauge_rhs, &global_gauge_rhs, 1, MPI_INT, MPI_SUM, state.comm);
  if (global_gauge_rhs) {
    if (state.rank == 0) fprintf(stderr, "HICAR HYPRE: zero-row gauge has %d non-zero RHS entries\n", global_gauge_rhs);
    free(ids); free(bvals); free(xvals); return -12;
  }
  #define HICAR_HYPRE_STEP(ID, CALL) do { stage=(ID); ierr=(CALL); if(ierr) goto done; } while(0)
  HICAR_HYPRE_STEP(1, HYPRE_IJVectorCreate(state.comm,state.first,state.last,&ij_b));
  HICAR_HYPRE_STEP(2, HYPRE_IJVectorSetObjectType(ij_b,HYPRE_PARCSR));
  HICAR_HYPRE_STEP(3, HYPRE_IJVectorInitialize_v2(ij_b,HYPRE_MEMORY_HOST));
  HICAR_HYPRE_STEP(4, HYPRE_IJVectorSetValues(ij_b,(HYPRE_Int)nlocal,ids,bvals));
  HICAR_HYPRE_STEP(5, HYPRE_IJVectorAssemble(ij_b));
  HICAR_HYPRE_STEP(6, HYPRE_IJVectorGetObject(ij_b,(void**)&b));
  HICAR_HYPRE_STEP(8, HYPRE_IJVectorCreate(state.comm,state.first,state.last,&ij_x));
  HICAR_HYPRE_STEP(9, HYPRE_IJVectorSetObjectType(ij_x,HYPRE_PARCSR));
  HICAR_HYPRE_STEP(10, HYPRE_IJVectorInitialize_v2(ij_x,HYPRE_MEMORY_HOST));
  HICAR_HYPRE_STEP(11, HYPRE_IJVectorSetValues(ij_x,(HYPRE_Int)nlocal,ids,xvals));
  HICAR_HYPRE_STEP(12, HYPRE_IJVectorAssemble(ij_x));
  HICAR_HYPRE_STEP(13, HYPRE_IJVectorGetObject(ij_x,(void**)&sol));
  HICAR_HYPRE_STEP(15, HYPRE_ParCSRFlexGMRESSetMaxIter(state.fgmres,max_iter));
  HICAR_HYPRE_STEP(16, HYPRE_ParCSRFlexGMRESSetTol(state.fgmres,tol));
  /* HYPRE's error state is process-global.  Surface a stale error before the
   * collective setup rather than attributing it to the solver call. */
  if (HYPRE_GetError()) {
    ierr = HYPRE_GetError();
    stage = 16;
    goto done;
  }
  if (!state.solver_ready) {
    fprintf(stderr, "HICAR HYPRE rank %d: entering FGMRES setup\n", state.rank); fflush(stderr);
    setup_start = MPI_Wtime();
    HICAR_HYPRE_STEP(17, HYPRE_ParCSRFlexGMRESSetup(state.fgmres,state.a,b,sol));
    state.solver_ready = 1;
    if (state.rank == 0) {
      fprintf(stderr, "HICAR HYPRE: FGMRES/AMG setup complete in %.3f s\n", MPI_Wtime() - setup_start);
      fflush(stderr);
    }
  }
  stage=18;
  ierr=HYPRE_ParCSRFlexGMRESSolve(state.fgmres,state.a,b,sol);
  /* A non-converged Krylov solve still has actionable iteration and residual
   * diagnostics.  Fetch them before returning its non-zero status. */
  HYPRE_ParCSRFlexGMRESGetNumIterations(state.fgmres,iterations);
  HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm(state.fgmres,residual);
  if (ierr) goto done;
  HICAR_HYPRE_STEP(21, HYPRE_IJVectorGetValues(ij_x,(HYPRE_Int)nlocal,ids,xvals));
  if(!ierr) for(q=0;q<nlocal;++q) x[q]=xvals[q];
done:
  if (ierr) {
    char error_description[256];
    HYPRE_DescribeError(ierr, error_description);
    fprintf(stderr, "HICAR HYPRE solve stage %d failed with code %d (argument %d): %s\n",
            stage, ierr, HYPRE_GetErrorArg(), error_description);
  }
  HYPRE_IJVectorDestroy(ij_b); HYPRE_IJVectorDestroy(ij_x);
  free(ids); free(bvals); free(xvals);
  #undef HICAR_HYPRE_STEP
  return ierr;
}

/* Apply the assembled ParCSR matrix to an owned local vector.  This is used
 * by HICAR's matrix-equivalence gate before accepting a HYPRE solve. */
int hicar_hypre_apply(const double *x, double *y)
{
  HYPRE_IJVector ij_x=NULL, ij_y=NULL;
  HYPRE_ParVector par_x=NULL, par_y=NULL;
  HYPRE_BigInt nlocal=state.last-state.first+1, q;
  HYPRE_BigInt *ids=NULL;
  HYPRE_Real *xvals=NULL, *yvals=NULL;
  int ierr=0;
  if (!state.valid) return -10;
  HYPRE_ClearAllErrors();
  ids=malloc((size_t)nlocal*sizeof(*ids)); xvals=malloc((size_t)nlocal*sizeof(*xvals)); yvals=calloc((size_t)nlocal,sizeof(*yvals));
  if (!ids || !xvals || !yvals) { free(ids); free(xvals); free(yvals); return -11; }
  for(q=0;q<nlocal;++q) { ids[q]=state.first+q; xvals[q]=x[q]; }
  #define HICAR_HYPRE_APPLY(CALL) do { ierr=(CALL); if(ierr) goto done; } while(0)
  HICAR_HYPRE_APPLY(HYPRE_IJVectorCreate(state.comm,state.first,state.last,&ij_x));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorSetObjectType(ij_x,HYPRE_PARCSR));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorInitialize_v2(ij_x,HYPRE_MEMORY_HOST));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorSetValues(ij_x,(HYPRE_Int)nlocal,ids,xvals));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorAssemble(ij_x));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorGetObject(ij_x,(void**)&par_x));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorCreate(state.comm,state.first,state.last,&ij_y));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorSetObjectType(ij_y,HYPRE_PARCSR));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorInitialize_v2(ij_y,HYPRE_MEMORY_HOST));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorSetValues(ij_y,(HYPRE_Int)nlocal,ids,yvals));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorAssemble(ij_y));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorGetObject(ij_y,(void**)&par_y));
  HICAR_HYPRE_APPLY(HYPRE_ParCSRMatrixMatvec(1.0,state.a,par_x,0.0,par_y));
  HICAR_HYPRE_APPLY(HYPRE_IJVectorGetValues(ij_y,(HYPRE_Int)nlocal,ids,yvals));
  for(q=0;q<nlocal;++q) y[q]=yvals[q];
done:
  if (ierr) fprintf(stderr, "HICAR HYPRE matrix apply failed with code %d\n", ierr);
  HYPRE_IJVectorDestroy(ij_x); HYPRE_IJVectorDestroy(ij_y);
  free(ids); free(xvals); free(yvals);
  #undef HICAR_HYPRE_APPLY
  return ierr;
}
