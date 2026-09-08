/*
   test_solve.c -- functional test for the SPOOLES solvers.

   Factors and solves a 3D grid Laplacian and checks the answer against the
   known exact solution: once through the serial factor and solve out of
   spooles.a, then repeatedly through FrontMtx_MT_factorInpMtx() and
   FrontMtx_MT_solve() out of spoolesMT.a.  That is the pipeline CalculiX
   and every other consumer drives.

   It exists because SPOOLES 2.2's MT factor and solve contain data races
   that return silently wrong results -- no crash, no diagnostic, just a
   solution that is wrong by O(1).  A "does the .a file exist" test cannot
   see that.

   Note that SPOOLES' own drivers cannot serve as a test.  They print a
   "maxabs error" and return 0 regardless of its value: MT/drivers/testGridMT
   on a racy library reports 3.3e+03 where 1e-12 is correct and still exits
   successfully.  None of the 156 shipped drivers compares its result against
   a threshold, and their 176 do_* wrappers are csh scripts carrying absolute
   paths from the author's 1999 machine.  They are diagnostics, not tests --
   which is a large part of why this bug survived since 1999.

   On the iteration count.  Measured per-iteration failure rates against the
   unpatched library, 20^3 grid:

       aarch64   4 threads  16.5%        8 threads  20.8%

       x86_64    2 threads   0/400       3 threads   1/400
                 4 threads   1/400       8 threads  11/1500  (0.73%)
                16 threads   1/500

   Three lessons.  x86's store ordering makes these races ~40x rarer than on
   aarch64, not harmless -- which is why they survived from 1999.  They need
   more runnable threads than cores in order to interleave, so this must not
   cap the thread count to the CPU count: on a 2-core CI runner that would
   detect nothing at all.  And detection peaks around 8 threads and falls off
   again by 16, where the oversubscription is heavy enough to serialize the
   very interleavings it is trying to provoke.  Hence 8.

   Oversubscription is not free: SPOOLES' inter-phase barriers are busy-wait
   spins with no yield, so 500 iterations at 8 threads on a 3-core conda-forge
   runner costs about 26 minutes of wall clock and 81 minutes of CPU, against
   27 s locally.  That is a deliberate trade.  The Azure job timeout is 360
   minutes and half an hour is an ordinary build-and-test, so spending ~26 of
   those minutes to raise the chance of catching an unpatched library from
   roughly 0.73% per iteration to 1 - 0.9927^500 ~= 97% is worth it for a defect
   whose entire character is that it is rare and silent.  The last argument is
   a wall-clock safety cap, not a tuning knob: it exists so that a
   pathologically slow runner bails out rather than burning the whole timeout.

   The recipe *also* checks deterministically that Utilities/SPOOLESatomic.h
   is installed.  That header exists only if the MT patch series applied, so
   the cheap regression -- a rerender or version bump silently dropping the
   patches -- is caught with certainty and at no cost, on every platform,
   including the cross-compiled and Windows ones where this loop cannot run.
   The two are complementary: the header check catches the patches going
   missing, this loop catches a genuinely new race.

   Where this runs.  conda-forge.yml cross-compiles both ARM targets
   (linux_aarch64 from linux_64, osx_arm64 from osx_64).  conda-build does
   *not* skip the test phase when cross-compiling -- it runs the commands
   anyway -- but the compiler package pulled in by test/requires is built for
   the host, so on those targets it is an executable the build machine cannot
   run ("Bad CPU type in executable").  The compile-and-run commands are
   therefore guarded by [build_platform == target_platform], leaving them on
   linux-64 and osx-64: the two platforms where the races are rarest.
   Building the ARM targets natively would make this decisive.

   Usage: test_solve [niter] [nthread] [gridsize]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "misc.h"
#include "FrontMtx.h"
#include "SolveMap.h"
#include "SymbFac.h"
#include "MT.h"      /* FrontMtx_MT_factorInpMtx, FrontMtx_MT_solve */

/*
   Tolerance.  The correct answer comes out at 1e-10..1e-8 (500 iterations at
   8 threads, worst 1.8e-08); a race drops a Schur-complement update, so the
   wrong answers are O(1) relative -- the smallest observed was 17.8, and they
   run to 1e+06.  1e-4 sits four orders above the noise and five below the
   signal, so this cannot fail spuriously on a well-conditioned run.
*/
#define TOL 1.0e-4

static FILE *quiet ;    /* SPOOLES chatters unconditionally; send it here */

/*
   nthread == 0 selects the serial factor and solve (spooles.a); anything
   larger selects the multithreaded pair (spoolesMT.a).
*/
static double
one_solve(int n, int nthread, int seed)
{
double        cpus[20], droptol = 0.0, tau = 100.0, err ;
int           error = -1, lookahead = 0, maxsize = 64, maxzeros = 0 ;
int           neqns, nrhs = 1, stats[20] ;
int           type = SPOOLES_REAL, symmetryflag = SPOOLES_SYMMETRIC ;
int           sparsityflag = 0, pivotingflag = 1 ;
Chv           *rootchv ;
ChvManager    *chvmanager ;
DV            *cumopsDV ;
DenseMtx      *mtxB, *mtxX, *mtxZ ;
ETree         *frontETree ;
FrontMtx      *frontmtx ;
IV            *frontOwnersIV ;
IVL           *symbfacIVL ;
InpMtx        *mtxA ;
SolveMap      *solvemap ;
SubMtxManager *mtxmanager ;
FILE          *msgFile = quiet ;

/* build the system: A, the exact solution X, and the rhs B = A*X */
mkNDlinsys(n, n, n, maxzeros, maxsize, type, symmetryflag, nrhs, seed,
           0, msgFile, &frontETree, &symbfacIVL, &mtxA, &mtxX, &mtxB) ;

/* map fronts to threads (unused in the serial case) */
frontOwnersIV = NULL ;
if ( nthread > 0 ) {
   cumopsDV = DV_new() ;
   DV_init(cumopsDV, nthread, NULL) ;
   DV_fill(cumopsDV, 0.0) ;
   frontOwnersIV = ETree_subtreeSubsetMap(frontETree, type,
                                          symmetryflag, cumopsDV) ;
   DV_free(cumopsDV) ;
}

frontmtx   = FrontMtx_new() ;
mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(mtxmanager, LOCK_IN_PROCESS, 0) ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL, type, symmetryflag,
              sparsityflag, pivotingflag, LOCK_IN_PROCESS, 0, NULL,
              mtxmanager, 0, msgFile) ;

/* factorization */
IVzero(20, stats) ;
DVzero(20, cpus) ;
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, LOCK_IN_PROCESS, 1) ;
if ( nthread > 0 ) {
   rootchv = FrontMtx_MT_factorInpMtx(frontmtx, mtxA, tau, droptol,
                                      chvmanager, frontOwnersIV, lookahead,
                                      &error, cpus, stats, 0, msgFile) ;
} else {
   rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, droptol,
                                   chvmanager, &error, cpus, stats,
                                   0, msgFile) ;
}
if ( rootchv != NULL || error >= 0 ) {
   fprintf(stderr, "\n factorization failed, error = %d\n", error) ;
   exit(1) ;
}
ChvManager_free(chvmanager) ;
FrontMtx_postProcess(frontmtx, 0, msgFile) ;

/* solve */
neqns = mtxB->nrow ;
mtxZ  = DenseMtx_new() ;
DenseMtx_init(mtxZ, type, 0, 0, neqns, nrhs, 1, neqns) ;
DenseMtx_zero(mtxZ) ;
solvemap = NULL ;
DVzero(20, cpus) ;
if ( nthread > 0 ) {
   solvemap = SolveMap_new() ;
   SolveMap_ddMap(solvemap, SPOOLES_SYMMETRIC,
                  FrontMtx_upperBlockIVL(frontmtx), NULL,
                  nthread, frontOwnersIV, frontmtx->tree, seed, 0, msgFile) ;
   FrontMtx_MT_solve(frontmtx, mtxZ, mtxB, mtxmanager, solvemap,
                     cpus, 0, msgFile) ;
} else {
   FrontMtx_solve(frontmtx, mtxZ, mtxB, mtxmanager, cpus, 0, msgFile) ;
}

/* compare against the exact solution */
DenseMtx_sub(mtxZ, mtxX) ;
err = DenseMtx_maxabs(mtxZ) ;

FrontMtx_free(frontmtx) ;
DenseMtx_free(mtxX) ;
DenseMtx_free(mtxB) ;
DenseMtx_free(mtxZ) ;
if ( frontOwnersIV != NULL ) IV_free(frontOwnersIV) ;
IVL_free(symbfacIVL) ;
InpMtx_free(mtxA) ;
ETree_free(frontETree) ;
SubMtxManager_free(mtxmanager) ;
if ( solvemap != NULL ) SolveMap_free(solvemap) ;

return err ;
}

int
main(int argc, char *argv[])
{
int      i, nbad = 0, ncpu, niter, nthread, n, budget, ndone ;
time_t   t0 ;
double   err, worst = 0.0 ;

if ( (quiet = fopen("/dev/null", "w")) == NULL ) {
   fprintf(stderr, "cannot open /dev/null\n") ;
   return 1 ;
}

niter   = (argc > 1) ? atoi(argv[1]) : 100 ;
nthread = (argc > 2) ? atoi(argv[2]) : 4 ;
n       = (argc > 3) ? atoi(argv[3]) : 20 ;
budget  = (argc > 4) ? atoi(argv[4]) : 0 ;   /* seconds; 0 = no limit */

/*
   Deliberately NOT capped to the online CPU count.  Measured detection rate
   against the unpatched library on x86_64 is 0/400 at two threads, 1/400 at
   three and four, and 6/1000 at eight: the races need more runnable threads
   than cores to interleave, so capping to a 2-core CI runner would make this
   loop useless.  But SPOOLES' inter-phase barriers are busy-wait spins with
   no yield, so oversubscription is expensive -- 8 threads on a 3-core runner
   took 26 minutes of wall clock for what takes 27 s locally.  Hence the wall
   clock budget below rather than a fixed iteration count: it bounds what CI
   pays without pinning the thread count to something that cannot detect
   anything.
*/
ncpu = (int) sysconf(_SC_NPROCESSORS_ONLN) ;
if ( nthread < 2 ) nthread = 2 ;

printf("SPOOLES solve test: %dx%dx%d grid Laplacian\n", n, n, n) ;
fflush(stdout) ;

/* First the serial factor and solve, out of spooles.a.  Deterministic, so
   once is enough -- but until now nothing tested it at all, and the I2Ohash
   patch touches that library. */
err = one_solve(n, 0, 10101) ;
printf("  serial factor + solve: maxabs error %.4e -- %s\n",
       err, (err < TOL) ? "ok" : "WRONG") ;
fflush(stdout) ;
if ( !(err < TOL) ) {
   printf("FAIL: the serial solver returned a wrong answer.\n") ;
   return 1 ;
}

/* Then the multithreaded pair, out of spoolesMT.a, many times. */
printf("  MT factor + solve: up to %d iterations at %d threads"
       " (%d online CPUs", niter, nthread, ncpu) ;
if ( budget > 0 ) printf(", %d s budget", budget) ;
printf(")\n") ;
fflush(stdout) ;
t0 = time(NULL) ;

ndone = 0 ;
for ( i = 0 ; i < niter ; i++ ) {
   if ( budget > 0 && i > 0 && (long)(time(NULL) - t0) >= (long) budget ) {
      printf("  stopping at %d iterations: %d s budget reached\n", i, budget) ;
      break ;
   }
   ndone++ ;
   err = one_solve(n, nthread, 10101 + i) ;
   if ( !(err < TOL) ) {          /* catches NaN as well */
      nbad++ ;
      if ( nbad <= 10 ) {         /* keep CI logs readable */
         printf("  iteration %d: maxabs error %.4e  <-- WRONG\n", i, err) ;
         if ( nbad == 10 ) printf("  ... further failures not listed\n") ;
         fflush(stdout) ;
      }
   }
   if ( err > worst ) worst = err ;
}

printf("%d/%d iterations wrong, worst maxabs error %.4e (tolerance %.1e)\n",
       nbad, ndone, worst, TOL) ;

if ( nbad > 0 ) {
   printf("FAIL: the multithreaded solver returned wrong answers.\n") ;
   return 1 ;
}
printf("PASS\n") ;
return 0 ;
}
