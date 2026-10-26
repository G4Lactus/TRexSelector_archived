# TRexSelector changelog

## TRexSelector C++ 1.0.0 / R 2.0.0 / Python 1.0.0

Release candidate on `dev_beyond` (2026-09-24); the three packages carry
their own version numbers from here on: the C++ SDK 1.0.0 (CMake project
version, `find_package(TRexSelector 1.0)`), the R package TRexSelectorNeo
2.0.0, the Python package TRexSelectorNeo 1.0.0.

### Breaking changes since 2.0.0 (June 2026)

- L-loop strategies: `ONDEMAND` and `PERMUTATION_ONDEMAND` are `SEEDED` and
  `PERMUTATION_SEEDED` (the seed-derived, nothing-stored strategies); the old strings
  are rejected by `trex_control()` and the Python enum has no alias.
- `SolverTypeForTRex` was renumbered when the CD family and the Tikhonov solvers
  were added; serialized integer codes from before July are not portable.
- CD solver family renamed: `TCCD`, `TCENET`, `TCIENET` are `TCDLASSO`, `TCDENET`,
  `TCDIENET`; the general-Tikhonov mode of `TCIENET` was extracted into
  `TCDTIKNET`, which is the Tikhonov selector's default solver. Checkpoints written
  by the old classes are unreadable by the new ones (formats incompatible).
- The CD family localizes every dummy crossing with certificates by default
  (`cd_certified_localization = TRUE`); the pre-August probe bisection is the
  LEGACY fallback and opt-in.
- `TRexGVSSelector`, `gvs_type = "IEN"`: the default solver is the native pathwise
  `TIENET`; `TIENET_AUG` and `TCDIENET` are the alternatives.
- `TRexTikhonovSelector` defaults: cluster-aware dummies (`cluster_dummies = TRUE`)
  and FOLDED dummy coupling (`fold_dummy_coupling = TRUE`); the `PILOT_MAX` /
  `PILOT_KNEE` lambda_2 rules exist next to the CV tuners (default `CV_1SE_IEN_CCD`).
- The SDK headers install under `<prefix>/include/TRexSelector/` (manual include
  paths change; exported targets do not), see the 2026-09-24 entry.
- The C++ tree has one folder per dummy mechanism: `ad/` (augmented dummies, the
  classical family), `vd/` (virtual dummies), `sd/` (sparse dummies), in both
  `tsolvers/` and `trex_selector_methods/`, and the namespaces follow:
  `trex::tsolvers::ad::linear_model::lars_based::TLARS_Solver`,
  `trex::trex_selector_methods::ad::trex_core::TRexSelector`. Include paths and
  namespaces of every C++ user change; the R and Python surfaces do not. See the
  2026-09-24 entry.
- The selectors reject non-finite entries and constant columns of X (constant
  columns were silently unselectable before), see `ScalingMode::NONE` below.
- Every T-solver is constructed as `(X, D, y, [groups | K,] options)` with one
  options structure; the positional constructors and the per-solver R / Python
  keyword arguments are gone, see below.
- Floating-point classification is bit-based and the Release build no longer uses
  `-ffast-math`: results are IEEE 754 conformant; runs from fast-math builds are
  not bit-reproducible against these.
- The T-Rex selectors mark their centered design `precentered`: on a saturating
  design the solvers stop at n - 1 active variables (the rank of the centered
  design) instead of taking an n-th, numerically collinear step. T-solver
  checkpoints written before 2026-09-27 are unreadable (the base state and T-OLS
  gained fields). See the 2026-09-27 entry.
- T-NCGMP in its line-search variant (the default) stops by the T-MP rules (the
  row bound on the active atoms, 8 x min(max_actives, row bound) steps) instead of running without
  an active-set cap for 20 x (p + L) steps; past saturation its paths are shorter
  and equal T-MP's. The augmented EN wrappers count the intercept in their
  degrees of freedom for lambda2 > 0 (one more than before). See the 2026-09-27
  entry.
- `SolverOptions::tie_seed` defaults to 0 (C++, R `tsolver_options()`, Python
  `SolverOptions`): a solver without an explicit seed breaks exact ties
  reproducibly; a negative seed keeps the former std::random_device behavior. The
  T-Rex selectors are unchanged (seeded selectors derive per-experiment seeds,
  unseeded ones pass -1).

### 2026-10-26

#### T-solvers: the augmented wrappers forward seed and tolerance to the inner solver; small items

- `TSolver_Base::setTieSeed()` and `setTolerance()` are virtual. `TENETAug_Solver` and
  `TIENETAug_Solver` hid both with forwarding members, so a call through a `TSolver_Base`
  pointer (R's `$set_tie_seed()` and `$set_tolerance()`, which hold the solver as a base
  pointer) seeded the wrapper's unused engine and set the wrapper's unused tolerance. The
  overrides forward to the inner solver. Python holds the concrete type and was not affected.
- The wrappers also built the inner solver's options without `tie_seed`, so the inner
  tie-break engine was always seeded 0 and a seeded selector's per-experiment seed (or -1 for
  an unseeded run) stopped at the wrapper. `tie_seed` is forwarded at construction. A
  selection can only move where the inner path meets an exact dummy tie (identical columns),
  which continuous designs do not produce.
- Tests: `TieSeedDispatch.AugmentedWrappersForwardThroughBasePointer` (two bitwise identical
  dummies and a budget of one; the survivor set through a base pointer or through the options
  equals the one set on the concrete type, and the seed decides it) and the same check in R
  (`test-tsolvers.R`, through `$set_tie_seed()` and `tsolver_options(tie_seed)`). The full
  C++ suite (660 tests) and the R tsolvers testthat files (547 expectations, against a
  scratch install) pass.
- `test_ml_methods_memory.cpp` included `<sys/resource.h>` and called `getrusage` unguarded.
  The Windows CI job builds every test, so the next push would have broken the MSVC job. On
  Windows the file holds one skipped test.
- The dead `.gitignore` line (`TODO before main merge.txt`) is gone. The briefs of the
  `tsolvers/ad.hpp` and `tsolvers/tsolvers.hpp` umbrella headers are one sentence each.
- R vendored copy re-synced.

#### R documentation of the TRexSelector result fields

- The field texts of `T_stop`, `L`, `voting_grid`, `fdp_hat_mat`, `phi_mat`, `r_mat` and `phi_prime` in
  `trex_selector.R` and `TRexSelector.Rd` described the wrong quantities ("probability mass function
  evaluations", "empirical CDF computations", "grid of voters over trials and variables"). They now name the
  relative occurrences after each T step, the number of selected variables and the FDP estimate per T step and
  voting level, the voting levels, the calibrated occurrences (at most 1, possibly negative), the stopping time and
  the dummy multiplier. Documentation only.

### 2026-10-24

#### FDP estimate: `Phi_prime` keeps negative values (lower clamp removed)

- `TRexSelector::computePhiPrime` clamped `Phi_prime` to [0, 1]. The upper clamp
  (Phi' > 1 -> 1) prevents a negative false-discovery mass 1 - Phi' and can only
  raise FDP_hat. The lower clamp (Phi' < 0 -> 0) capped a selected variable's
  false-discovery mass at one, where the T-Rex paper and the R reference let it
  exceed one, and made FDP_hat anti-conservative in long T-loops. In demo 09 of
  the Examples (n 300, p 1,000, 10 unit coefficients, Normal predictors and noise,
  K 20), TLARS at SNR 5 reached an FDR of 0.530 at tFDR 0.5 and 0.415 at tFDR 0.4
  (no stagnation stop, T-loop to T ~ 38), and FDP_hat at fixed T and v -> 1 fell
  below the realized FDP from T = 10 on (0.408 against 0.444 at T = 20). Replaying
  the final selection on the same 200 data sets with the unclipped Phi_prime gives
  0.462 and 0.362. Only the upper clamp remains, together with the [0, 1] clamp of
  `computeFDPHat`, which guards against the negative FDP estimates seen in
  trex_da (they come from Phi' > 1).
- The same change in `TRexSD::Phi_prime_fun_` and `TRexVD::Phi_prime_fun_`.
  trex_da, trex_gvs, trex_screening, trex_tikhonov and trex_spca reuse
  `computePhiPrime` and follow. FDP_hat can only rise, so the L-loop may calibrate
  a larger L, the T-loop may stop earlier and selections can only shrink in size.
  Results move where a selected variable had Phi' < 0, which the replay found at
  tFDR >= 0.3 and not at tFDR 0.1 in the tested cells (TLARS and TOMP, SNR 1
  and 5); the L-loop and T-loop effects need full runs.
- Tests: `TRexCoreTest.PhiPrimeKeepsNegativeValues` builds a Phi history with a
  negative Phi_prime and checks it and FDP_hat against the paper's formula. The
  Phi_prime bound checks of the core, Tikhonov, GVS and VD tests and of the R
  testthat files keep only Phi' <= 1 (`TRexVDTest.PhiPrimeIsCappedAtOne` also
  checks FDP_hat in [0, 1]). The full C++ suite (659 tests), the Python
  selector tests (172) and the R testthat files of the selector, DA, GVS,
  screening and scaling paths pass.
- R vendored copy re-synced (`bundle_for_cran.sh --check` passes).

### 2026-10-22

#### TRex-SPCA accepts the elastic net only; flaky low-memory CV budget

- `TRexSPCASelector` is defined on the elastic net, but the R control
  `trex_spca_control(gvs_type = "IEN")` and the Python docstring offered the
  informed elastic net, and the C++ core silently ran EN instead. An IEN request
  is now an error: the C++ constructor throws `std::invalid_argument` when
  `gvs_ctrl.gvs_type` is not EN (Python raises `ValueError` in `select()`), and
  `trex_spca_control()` stops for any `gvs_type` other than `"EN"`. The IEN
  example of `trex_spca_control()` is gone, the Python docstrings (binding,
  stub, wrapper) say EN only, and the Python test that ran "IEN" as EN now
  checks the error. `cv_dfmax_ratio` stays unexposed in the SPCA controls: it
  only acts on the IEN-geometry tuners, which the EN track rejects.
- `MLMemory.ElasticNetCVLowMemory`: budget 0.3 x X instead of 0.2. One
  low-memory fit needs 0.095 x X in every run, but the six sequential fits of a
  5-fold CV measured 0.12-0.23 x X for identical input (whether the allocator
  reuses the freed paths varies), so about one run in four failed. A copy of X
  would still add 1.0 x X.
- Tests: `test_trex_spca` gains `Exceptions_IENRejected` (16 tests pass),
  `test_ml_methods_memory` passes 12 of 12 repeated runs of the CV case. The R
  test that ran "IEN" as EN now checks that a hand-built IEN control list fails
  in the C++ constructor, next to a check of `trex_spca_control()`. The Python
  SPCA tests (23), the R SPCA testthat file (25 blocks) and the full C++ suite
  (658 tests) pass.
- R vendored copy re-synced (`bundle_for_cran.sh --check` passes).

### 2026-10-21

#### Settable support cap of the IEN-geometry lambda_2 tuners (`cv_dfmax_ratio`)

- The IEN-geometry tuner `ienet_cv_ccd` (`CV_1SE_IEN_CCD` / `CV_MIN_IEN_CCD`)
  stops each lambda_1 path of its reference scan once more than
  floor(`dfmax_ratio` x n_train_min) coefficients are nonzero. `TRexGVSSelector`
  and `TRexTikhonovSelector` called it without the ratio, so the cap was fixed at
  0.8 of the smallest fold's training size. The new control field
  `cv_dfmax_ratio` (default 0.8) sets it: `<= 0` disables the cap, larger values
  allow larger supports at a higher tuner cost, and a non-finite value throws
  `std::invalid_argument` at construction. The selectors now pass the tuner's
  own defaults for `max_iter` (100000) and `tol` (1e-4) explicitly, so the
  default 0.8 gives the same results as before.
- C++: `TRexGVSControlParameter::cv_dfmax_ratio` and
  `TRexTikhonovControlParameter::cv_dfmax_ratio`. Python: the same field on both
  control classes (stubs updated by hand). R: argument `cv_dfmax_ratio` of
  `trex_gvs_control()` and `trex_tikhonov_control()`, which reject a non-finite
  value as well; Rd files updated by hand.
- `lambda2_cv_diagnostics` gains `dfmax_ratio` (the ratio the IEN tuner was
  called with, as passed; 0 when no IEN-geometry tuner ran; `support_cap == -1`
  still marks a disabled cap) and `profile_at_support_cap` (the profiled
  lambda_1 optimum of the chosen lambda_2 sits on the truncation edge, and the
  reference path of that lambda_2 stopped on the support cap rather than on a
  deviance rule or the end of the path). Both are exposed in Python
  (`Lambda2CVDiagnostics`) and in the R `lambda2_cv_diagnostics` list.
- Boundary-pick warning: when the edge is the support cap it says so ("here the
  support cap (more than N nonzero coefficients, cv_dfmax_ratio = r)") and adds
  the remedy "raise the support cap (cv_dfmax_ratio)"; the other wording is
  unchanged.
- Tests: `test_trex_lambda2_guard` (GVS default and explicit 0.8 bit-identical,
  0.95 raises the cap from 43 to 51 at n = 60 with 10 folds, the ratio rides in
  the diagnostics and the GVS result, non-finite values throw; the same
  forwarding for the Tikhonov selector), `test_ml_methods_model_selection` (a
  binding cap is flagged and named in the warning, a disabled cap is not; the
  message contract), Python and R control round trips and diagnostics fields.
  All pass: the full C++ suite (657 tests), the Python GVS/Tikhonov tests (40)
  and the R GVS/Tikhonov testthat files (29 blocks).
- R vendored copy re-synced (`bundle_for_cran.sh --check` passes).

### 2026-10-20

#### Free-space check before memory-mapped files are created

- Memory-mapped files are created sparse, so a full disk only showed at a later
  mapped write and ended the process with SIGBUS (no error, no cleanup). The
  new `utils/disk_space.hpp` compares the bytes a mapping may write with the
  free space of its volume and throws `std::runtime_error` with the required
  and free GiB, the directory and a remedy instead; a volume that cannot be
  queried passes.
- `MemmapManager::initialize()` (the dummy matrices of every T-Rex selector with
  `use_memory_mapping`) checks the worst case before it creates anything: all
  files filled up to the dummy cap, K files (one if shared) of
  n x `max_dummy_multiplier` x p doubles (`requiredBytes()`). The message names
  `max_dummy_multiplier`, K and the temporary directory (TMPDIR; R sets
  `tempdir()` at load; `TempDirManager::setCustomTempDir` in C++).
- `MemoryMappedMatrix` in ReadWrite mode checks the growth of a new or enlarged
  file at its own location before creating it, so user-created designs (R
  `mem_matrix`, Python `MemoryMappedMatrix`) are covered as well.
- `TempDirManager::create` takes the description and the remedy for this
  message, `TempDirManager::baseDir()` returns the directory it creates in; the
  warm-start checkpoints keep their 5 GiB minimum with the new message.
- Test `test_utils_disk_space` (7 cases: helpers, matrix, temp dir and manager
  refuse an unfulfillable size before creating a file or directory, small sizes
  still work). R vendored copy re-synced (`bundle_for_cran.sh --check` passes).
  No result changes.

### 2026-10-08

#### Scientific references in the class documentation; master bibliography; citation files

- `references.bib` at the repository root is the master list of every paper the
  library credits (49 entries: the T-Rex method papers, the origin paper of each
  solver algorithm, the clustering, hashing, cross-validation and ridge
  references, Eigen). `tools/references/render_references.py` renders an entry
  as one plain-text line (initials and surnames, title, venue, volume, pages,
  year, DOI) wrapped for a Doxygen block, a roxygen block or a pybind11
  docstring; `references_map.json` fixes which keys each class cites, the origin
  papers first and the T-Rex paper (the greedy-pursuits paper for the greedy
  solvers) for the terminating form. No citation machinery in the sources: no
  Doxygen `\cite`, no Rdpack, no Sphinx roles.
- C++: every selector (`TRexSelector`, DA, GVS, screening, biobank screening,
  SPCA, Tikhonov), `MemoryMappedMatrix`, all 21 T-solvers, the clustering
  classes (SLINK, NN-chain, generic linkage, the dispatcher, the Lance-Williams
  coefficients, the LSH distance and update policies), the elastic-net, IEN,
  ridge and Tikhonov cross-validation tuners, the EBIC criterion, PCA, the SVD
  solver, the ridge solver and the dummy `Distribution` carry a `References:`
  paragraph in their class block. The hand-written lists of the GVS, TENET,
  TIENET, TIKNET, TMP, TGP, TACGP, SLINK, NN-chain and generic-linkage headers
  and the `@note Reference` lines are replaced by it; TCDLASSO no longer points
  to an internal note that is not in the repository, and two demo pointers are
  gone. Vendored copy in the R package re-synced (`bundle_for_cran.sh --check`
  passes).
- R: `@references` on the seven selectors, `mem_matrix` / `mem_matrix_from`, the
  21 T-solver classes, the six ML classes, the three clustering functions and
  `trex_dummy_distribution`; 40 Rd files regenerated, NAMESPACE unchanged.
  `inst/CITATION` names TRexSelectorNeo, keeps "TRexSelector Team" as the package
  author and lists the eleven T-Rex method papers with text versions rendered
  from the master file (`tools/references/make_inst_citation.py`);
  `citation("TRexSelectorNeo")` prints them.
- Python: the selectors, `DummyDistribution`, `MemoryMappedMatrix`, the T-solver
  classes, the ML classes, the clustering enums and `agglomerative_cluster` end
  their docstring with a numpy-style "References" section. The T-solver and ML
  binding helpers now set a class docstring (the classes had none, only
  `__init__` had one); stubs regenerated with pybind11-stubgen 2.5.5 against a
  pybind11 3.1.0 build (an older pybind11 found by CMake drops the typed enum
  `__eq__` / `__ne__` overloads of the committed stubs; point `pybind11_DIR` at
  3.1.0 when regenerating).
- Repository: `CITATION.cff` at the root for GitHub's "Cite this repository"
  (preferred citation the T-Rex paper, the other method papers and Eigen as
  references); the READMEs acknowledge Eigen and point to the citation files,
  the Python and R READMEs link the journal DOI of the T-Rex paper.

### 2026-09-28

#### OpenMP optional in every build; one OpenMP runtime in the R package on macOS

- C++ and Python: OpenMP is optional, as in R. `find_package(OpenMP)` without
  REQUIRED; when found, `trex_utils` and `trex_tsolvers` link `OpenMP::OpenMP_CXX`
  as before, otherwise CMake reports "OpenMP not found: building single-threaded"
  and the stubs in `utils/openmp/utils_openmp.hpp` take over (pragmas ignored,
  `omp_get_max_threads()` 1). `-DTREX_WITH_OPENMP=OFF` (default ON) forces the
  single-threaded build; Python takes it as `-C cmake.define.TREX_WITH_OPENMP=OFF`.
  A single-threaded build adds `-Wno-unknown-pragmas` (MSVC `/wd4068`), so GCC's
  one warning per inactive pragma is gone. The Homebrew LLVM libomp pin applies
  only to Homebrew's clang and when that libomp exists; MSVC keeps `-openmp:llvm`.
- Installed package config: `find_dependency(OpenMP)` only when the library was
  built with OpenMP; `TRexSelector_WITH_OPENMP` tells consumers which. A
  single-threaded install exports no OpenMP link. The same for BLAS/LAPACK:
  `find_dependency(BLAS/LAPACK)` only for a `TREX_USE_BLAS_LAPACK=ON` build
  (`TRexSelector_WITH_BLAS_LAPACK`); before, every install required them.
- R on macOS: CRAN's R has no OpenMP flags, and configure linked Homebrew's
  libomp by absolute path while `-lomp` resolved to R's bundled runtime, so two
  OpenMP runtimes were mapped. configure now probes (1) R's flags, (2) on macOS
  R's `lib/libomp.dylib` with an omp.h from Homebrew or from the mac.r-project.org
  OpenMP package in /usr/local/include (`-Xpreprocessor -fopenmp`,
  `-L${R_HOME}/lib -lomp`, the runtime data.table and RcppArmadillo load),
  (3) Homebrew's runtime alone, (4) sequential. `TREX_OPENMP=no` forces the
  sequential build and adds `-Wno-unknown-pragmas`; the automatic fallback does
  not, since R CMD check lists that flag as non-portable. configure regenerated
  with autoconf 2.73 (AC_INIT now names version 2.0.0, not 0.2.0); Makevars.win
  unchanged.
- Tests: two OpenMP stub tests failed without OpenMP (status string, the wtime
  stub returns 0 on its first call) and are fixed; `OpenMPTest.SingleThreadedStubs`
  checks the one-thread stubs. CI: `test-no-openmp` (ubuntu-24.04) builds and
  tests with `TREX_WITH_OPENMP=OFF` and fails if any target compiles with -fopenmp,
  then runs the install smoke test with OpenMP disabled for the consumer
  (`TREX_SMOKE_CONSUMER_ARGS`, new in run_install_smoke.sh), so the config of a
  single-threaded install is proven not to ask for OpenMP.
- Verified: C++ 646 of 646 with OpenMP (22 skipped) and without (24 skipped, the
  three OpenMP thread-hygiene tests; no build warnings); R 1,407 with R's runtime
  (12 threads, one libomp mapped) and with TREX_OPENMP=no (1 thread); Python 635
  with and without OpenMP; vendor check OK; install smoke tests OK without
  OpenMP, and without OpenMP, BLAS and LAPACK (consumer with all three disabled).

#### T-solvers: the held heap bytes as `memoryBytes()`, R `$get_memory_bytes()`, Python `getMemoryBytes()`

- Found in the Examples memory demo (demo_ts_perf_03: p = L = 1,000, n up to
  2,000, T = 10): the peak resident set of a child process is no measure there.
  libomp resets the thread count in every forked child on macOS (the
  "single-threaded" children ran on 12 threads), and a spawned R or Python
  interpreter serves the sub-MiB solver allocations from pages it already holds
  (R read about 0 MiB per solver). The library now reports the memory itself, so
  the demo can run in one process per language.
- `TSolver_Base::memoryBytes()` (virtual, const, noexcept) sums the heap bytes of
  the state a solver holds: the y copy, residual, correlations, Cholesky factor,
  sparse path and running coefficients, diagnostics, index sets, the writer's
  means and norms and every reservation (vectors by capacity, Eigen objects by
  size). The classes with own containers add theirs (TLARS, TIENET, TIKNET, TOMP,
  TGP, TACGP, TMP, TOOLS, TIKOMP, TCDLASSO, TCDIENET, TCDTIKNET and the augmented
  wrappers); T-LASSO, T-Stepwise, T-Stagewise, T-ENET, T-NCGMP, T-AFS and T-CDENET
  inherit. Not counted: the borrowed X and D (readers and writers alike), the
  temporaries inside a step and the object itself (sizeof, about 3.3 kB with the
  tie-break engine). The coefficient lookup map is an estimate (bucket array plus
  four words per entry). A call is O(members + steps); nothing runs on the hot
  path, no member is added and checkpoints are unchanged.
- The report counts the up-front reservations as they are: `dummies_at_step_` and
  the lambda sequences (LARS, CD) hold 8 min(max_actives, row bound) + 1 entries
  (T-AFS divides by rho), 70 % of T-LARS at n = 2,000 (256 of 367 kB) and 94 % of
  T-AFS with rho 0.1 (1,279 of 1,361 kB). T-CDLASSO keeps its Gram rows and its
  continuation cache after a step; they are counted too. The augmented EN
  wrappers count their owned augmented copies (8 (n + p + L)(p + L + 1) bytes for
  TENETAug), the diagnostic mirror, the inner solver object and
  `inner_->memoryBytes()`; the inner X and D view the owned copies, so nothing is
  counted twice.
- `solver_utils/memory_bytes.hpp`: `bytesOf()` for vectors (capacity, nested
  heap), dense and sparse Eigen objects, strings and the hash-map estimate; the
  products go through `checked_bytes`.
- Bindings: R `TSolver_Base$get_memory_bytes()` returns a double (a byte count
  overflows R integers at 2 GiB), Python `getMemoryBytes()` exists on every
  solver. RcppExports, the Rd pages and the T-solver stubs are regenerated
  (pybind11 3.1.0, pybind11-stubgen 2.5.5); the vendored copy is re-synced.
- Against the allocator (macOS malloc zones, the 15 Examples demo solvers as
  precentered readers, p = L = 1,000, T = 10): `memoryBytes()` is 88-92 % of the
  bytes in use at n = 50 and 93-98 % at n = 2,000, the rest being the solver
  object and the allocator's block rounding. The peak inside a step is 1.04-1.67x
  the held state at n = 50 and 1.02-1.29x at n = 2,000.
- Tests: C++ 645 of 645. `test_tsolvers_memory` (52): exact hand counts at
  construction for T-LARS, T-OMP, T-OLS, T-AFS (rho 0.5), T-CDLASSO and TENETAug,
  growth bounds after a step, doubling n moves only the two n-vectors, a writer
  holds 16 (p + L) bytes more than a reader, every solver reports a positive and
  growing value, the checkpoint is byte-identical. `test_tsolvers_memory_peak`
  (21): a counting allocator (malloc zones on macOS, `__libc_*` on glibc, skipped
  elsewhere and under sanitizers) bounds peak / held per family; measured up to
  1.45 (LARS), 1.20 (OMP, AFS) and 4.3 (CD, whose certified localization builds
  the restricted Gram for an eigen-solve). Verified on Linux (GCC 13, glibc 2.39)
  as well. R 1407 expectations, Python 635.

### 2026-09-27

#### T-solvers: selection read from the coefficient path, the fixed default tie seed, IEEE Examples

- `cdiagnostics::selected_variables()`: the selection of a fitted solver as the
  T-Rex selector reads it, the variables with |beta_j| > machine epsilon at the
  last recorded step, split into predictors and dummies (active ones in
  active-set order, then the frozen ones). `print_selection()` and
  `print_selection_quality()` use it; they listed the active set, which
  under-reported T-Stagewise (one demo run: 32 active, 44 with a coefficient,
  four true variables near 2.7 among the frozen ones). Frozen variables are
  printed on a line of their own; for the other solvers the output is unchanged.
  The TSTAGEWISE_Solver documentation states the frozen coefficients.
- `SolverOptions::tie_seed` defaults to 0: direct solver runs are reproducible
  also where exact ties occur (T-LASSO at T = 50 changed in 1 of 500 demo runs
  between repetitions with the random_device default). A negative seed keeps the
  nondeterministic seeding; the T-Rex selectors pass their seeds explicitly.
- Cycling at saturation (T-Stagewise through removals and re-entries up to the
  step limit, T-LASSO more mildly) is part of these solvers and is left as is.
- The Examples compile their demos without `-ffast-math`, as the library: the
  header-inline library code they instantiate stays IEEE 754 conformant, and demo
  numbers are bit-comparable with the library builds.
- Solver identities over full paths (Examples cpp/tsolvers/README.md, measured on
  demo_ts_perf_01, 100 runs to T = 100): T-OMP = T-AFS (rho 1) = T-NCGMP (fully
  corrective) bit for bit; T-Stepwise the same steps as T-OMP (coefficients within
  1.1e-12); T-MP = T-NCGMP (line search) in every step (coefficients within
  2.1e-14); T-LASSO = T-LARS bit for bit up to the first T-LASSO drop, T-Stagewise
  within 2.7e-13 up to its first NNLS removal.
- Tests: C++ 572 of 572 (new: `SelectedVariablesReadTheCoefficientPath`,
  `TStagewiseSelection.FrozenCoefficientsAreSelected`,
  `DefaultTieSeedIsReproducible`), R and Python check the default seed.

#### T-solvers: T-NCGMP stopping rules, augmented EN degrees of freedom, install timestamps

- T-NCGMP line search is T-MP step for step (Locatello et al. 2017, Alg. 4 with
  L = 1) but ran without an active-set cap for 20 x (p + L) steps, so at T = 100
  in demo_ts_perf_01 it ended with about 106 selected variables against 54 for
  T-MP. It now stops by T-MP's rules (effective_n_ distinct atoms,
  8 x min(max_actives_, effective_n_) steps); the fully corrective variant keeps
  T-OMP's. On the demo setup (100 runs, T up to 100) the line-search actions equal
  T-MP's in every run, the coefficients to 2e-14 (T-NCGMP divides the step by
  `||x_j||^2`, 1 up to rounding on normalized columns).
  `NCGMPVariantsMatchMPAndOMP` pins the full paths of both variants; the
  precentered contract test no longer exempts T-NCGMP from the n - 1 cap.
- TENETAug and TIENETAug reported the degrees of freedom of their inner solver,
  which runs intercept-free on the augmented system for lambda2 > 0, so the
  intercept was missing even with intercept = true. The wrappers add it
  (intercept or precentered), as the pathwise TENET and TIENET count it; the
  TENET / TENETAug and TIENET / TIENETAug parity tests compare the DoF.
- Install timestamps: `file(INSTALL)` gives every installed header the
  modification time of its source. A consumer object compiled after a header
  edit but before the install looked newer than the changed header, was not
  rebuilt, and the changed class layout corrupted its heap (seen in a scratch
  benchmark after the T-OMP / T-OLS members were added). `cmake --install` now
  records the SHA-256 of the header tree before pruning it and stamps new or
  changed files with the install time; unchanged headers keep theirs (verified:
  fresh install 106 files, repeat install 0, one altered header 1; the install
  smoke test passes).
- On a saturating design TCDLASSO can hit its sweep budget in `solveAtLambda`,
  after which its working active set runs ahead of the last recorded step
  (support 16, working set 18, n = 20); the recorded path and its DoF are
  consistent. Noted, not changed.
- R glue: the comment sweep had missed `R/TRexSelectorNeo/src/rcpp_*.cpp` and
  `TRexSelectorNeo_types.h`; 23 plan references, work-package tags and date
  stamps reformulated as the technical statement.

#### T-solvers: T-OLS scoring in one pass per step, the T-OMP refit, and `precentered`

- Found in the Examples performance demo (demo_ts_perf_01: n 100, L = p = 300,
  500 runs): T-OLS took 1.25 ms to T = 1 and 82 ms to T = 100, against 0.17 and
  2.50 ms for T-OMP. The scoring re-orthogonalized every candidate against the
  whole active set at every step (`h_j = X_A^T x_j` and a triangular solve per
  candidate), `O((p + L) n k^2)` per path; 97 % of the T-OLS time.
- T-OLS keeps the squared projected norms `||P_perp x_j||^2` of all candidates and
  downdates them by `(q^T x_j)^2` for each new direction `q = X_A R^{-1} e_k` of the
  active span, fused into the correlation pass: one `O(n (p + L))` pass per step,
  the modified Gram-Schmidt form of OLS (Blumensath & Davies, OMP vs. OLS) without
  an orthogonalized copy of the dictionary. 82 -> 2.9 ms to T = 100, 18.5 ->
  1.2 ms to T = 20, 1.25 -> 0.28 ms to T = 1 (1.3-1.6x T-OMP). The path equals the
  from-scratch score bit for bit up to n - 1 active variables (500 of 500 demo
  runs; pinned by `TOOLSScoring.DowndatedNormsMatchFromScratchProjection`); the
  norms are serialized, so a resumed solver continues bit for bit.
- T-OMP: `updateActiveSet()` no longer copies `R_` as a backup before
  `updateR()`, which never mutates it (T-AFS and TIKOMP already skipped the
  copy), and the refit caches `X_A^T y` in active-set order. 2.50 -> 2.13 ms to
  T = 100 at demo conditions (per core 2.30 -> 1.83 ms); paths bitwise unchanged
  (T-OMP, TOOLS and TIKOMP share the refit).
- `SolverOptions::precentered` (R `tsolver_options(precentered = )`, Python
  `SolverOptions(precentered=)`): the caller centered X, D and y. The solver
  counts the intercept as with `intercept = true` — the active set is capped at
  n - 1, the degrees of freedom and Cp include the intercept — without touching
  the data. Before, a reader on centered data allowed n active variables; at
  saturation the n-th entered as a noise pick through the machine-epsilon
  collinearity check and blew up the refit (100 demo runs at T = 100: T-LARS 98
  runs, coefficients up to 4.4x; T-OMP / T-AFS 76 runs, 7.6x; T-OLS 55 runs,
  coefficients 5e5). The T-Rex dispatch (`buildRunnerConfig`) and the GVS K-loop
  set it in every ScalingMode except NONE; the augmented EN wrappers pass it to
  the inner solver only on the lambda2 = 0 collapse (the augmented rows are not
  centered). The elastic nets with lambda2 > 0 and the line-search T-NCGMP keep
  running without the cap. Contract test for all 18 AD solvers:
  `PrecenteredCountsInterceptLikeIntercept` (same path, DoF and Cp as
  `intercept = true` on exactly centered data).
- With `precentered` left false every other solver path is bitwise unchanged
  (12 solver configurations, 100 demo runs, actions, DoF and coefficients).
- Tests: C++ 569 of 569, R 282 expectations, Python 611.

### 2026-09-26

#### CRAN readiness: the as-cran result, and the Linux installed size

- `R CMD check --as-cran` on macOS (R 4.6): Status 1 NOTE, the "new
  submission" one every first submission gets; tarball 1.7 MB, installed
  11.8 MB, examples (also `--run-donttest`), tests, both manuals and the
  vignette OK; install 323 s, vignette rebuild 158 s.
- The 292 MB installed size seen in the Ubuntu run was debug information,
  not code: R compiles with `-g` on every platform, macOS keeps the DWARF
  in the object files (the shared object is 6.8 MB), ELF platforms link
  the DWARF of the 48 Eigen-templated translation units into the shared
  object. The generated `src/Makevars.in` now strips it after linking on
  Linux (`all: $(SHLIB) strippedLib`, the CRAN-accepted post-link idiom of
  rstan and prophet; no compiler flag is overridden; a no-op elsewhere).
  Measured in the provisioned Ubuntu image: shared object 5.6 MB,
  installed package 7.1 MB. `Makevars.win` is unchanged; the templates in
  `bundle_for_cran.sh` are plain ASCII.
- Still to do before a submission: a precomputed vignette (2.6 min of
  simulations at check time), win-builder for Windows and R-devel, and
  the release from the public main repository.

### 2026-09-25

#### Linux verification of f390bfc (OrbStack): one budget, one locale

- The three Ubuntu jobs on the committed f390bfc: C++ 548 of 549, the R
  suite, `R CMD build`, `R CMD check` (the three known environment NOTEs),
  Python 610, the install smoke, and the cp312 manylinux_2_28 x86_64 wheel
  with its in-container tests — all green, except:
- `MLMemory.RidgeCV` measured 4.50 X on Ubuntu (Eigen 3.4) against the
  4 X budget set on macOS (Eigen 5.0): the Jacobi SVD keeps two internal
  copies of the fold under 3.4 (the scaled copy and the QR preconditioner)
  and one under 5.0. The budget is version-dependent now (4 X / 5 X).
- `R CMD check` reported a WARNING in the syntax check because a bare
  `ubuntu:24.04` image has no `en_US.UTF-8` locale; the stage script
  installs `locales`, generates it and runs the check under it.

#### ML methods: zero-copy harness and the read-only binding boundary (Tier A)

- Ground rule (Fabian): no binding copies X, D or a memory-mapped matrix;
  full copies of the design are avoided wherever they are not needed, and
  functionality is preserved at 100 %. Plan and tiers:
  notes/ZERO_COPY_PLAN.md; inventory: notes/R_PACKAGE_ASSESSMENT.md.
- Harness. `test_ml_methods_reference.cpp` runs every ML core on three
  fixed shapes and compares every output byte for byte against a reference
  recorded outside the repository (`TREX_ML_REF_DIR`; skipped when unset,
  recorded on first use); `tools/zero_copy_reference.{R,py}` do the same
  for the R and Python surfaces (`record` / `compare`). `test_ml_methods_memory.cpp`
  measures the peak-memory growth of each core in units of X on a wide
  500 x 8000 design and pins the current budgets: ridge CV 4.55 X, elastic
  net 1.09 X, elastic net CV 2.17 X, SVD 0.35 X, ridge 0.33 X, PCA 0.35 X,
  scalers 0.03 X (one process per test; ctest does that).
- R binding (Tier A, bit-identical): the eight ML entry points that took
  `const Eigen::MatrixXd&` — and therefore copied X in full on every call
  (RcppEigen converts only to a mutable map) — take the R matrix and build a
  read-only `Eigen::Map<const>` that binds to the cores' `Ref<const>`:
  `ridge_cv_fit`, `enet_fit`, `enet_fit_grid`, `enet_predict`, `enet_cv_fit`,
  `svd_compute`, `pca_transform`, `ridge_solve`. New memory-mapped entry
  points through the const map (read-only files accepted): `enet_fit_mmap`,
  `enet_fit_grid_mmap`, `enet_predict_mmap`, `enet_cv_fit_mmap`,
  `svd_compute_mmap`, `ridge_solve_mmap`, `pca_transform_mmap`, and
  `pca_create_mmap` / `pca_restore_mmap` (readwrite: the file is centered in
  place). The R6 classes route `mmap_matrix` inputs to them instead of
  reading the whole file through `as.matrix()`; `pca_create` pins X in the
  XPtr protection slot. R and Python references identical; new R tests
  compare every mmap path with the in-memory path bit for bit.
- Python (Tier A): `agglomerative_cluster` takes `Ref<const MatrixXd>`
  instead of `const MatrixXd&` (no copy of a Fortran-ordered array).

#### ML cores: copy fusion (Tier B, bit-identical)

- Ridge CV: the gathered training fold is standardized in place instead
  of into a second matrix (the two passes round exactly like the former
  fused expression): peak 4.55 X -> 3.64 X on the wide design; the Jacobi
  SVD's own copies remain.
- Elastic net: the gathered CV fold becomes the fit's working copy
  (`fit_owned` / `preprocess_owned`, centred and standardized in place)
  instead of being copied a second time: elastic-net CV 2.17 X -> 1.37 X.
  A fitted `enet_gaussian` releases its centred working copy after the
  path is solved (`predict` and the accessors read the p x K path only),
  so a fitted object holds p x K, not n x p.
- Ridge: the Cholesky factorization runs in place in the Gram / kernel
  matrix (`LLT<Ref<MatrixXd>>`): one p x p (or n x n) instead of two; the
  rare LDLT fallback rebuilds the matrix from X (a second pass, no copy).
- PCA: `PCA::readOnly(X, M)` fits a const matrix or a read-only
  memory-mapped file without writing (a named factory, so `PCA(X, M)` keeps
  centering and normalizing in place); `restore()` may be called before
  `fit()` (the preprocessing happens at construction; it used to throw),
  is idempotent, and checks that it is given the constructor's matrix; the
  destructor writes X back when `restore()` was never called, and the duty
  travels with a move (the T-solver rule; the Python docstring that
  claimed this is true now). X stays pinned in the R XPtr protection slot
  and by Python's keep_alive.
- Scalers: `fit()` takes a read-only reference, so the statistics can be
  fitted on a read-only memory-mapped file in R and Python (the in-place
  transforms still need a writable mapping).
- Clustering (R): the transpose copy on the in-memory path stays, on
  purpose — the engine clusters columns and the R surface clusters rows,
  and a lazy transpose would make every object a strided row; the comment
  that blamed undefined behaviour on the lazy transpose described an
  NN-chain bug the core fixed in June. The memory-mapped path reads the
  file through the const map (read-only mappings accepted).
- Verified: the 21 C++ reference cases, the R and Python references and
  all three suites are identical / green after every item.

#### ML cores: low-memory engines (Tier C, opt-in; the defaults stay bit-identical)

- Elastic net and its CV (`enet_gaussian(bool low_memory)`,
  `enet_cv_ccd(bool low_memory)`; R `MLElasticNet$new(low_memory = TRUE)`,
  `MLElasticNetCV$new(low_memory = TRUE)`; Python
  `ElasticNet(low_memory=True)`, `ElasticNetCV(low_memory=True)`): the
  coordinate descent runs on the raw X through the caller's view — a
  read-only memory-mapped file included — with the centring carried as a
  scalar residual offset (folded into the residual once per sweep, which
  keeps the training residual mean-zero and the dot products free of the
  cancellation a large column mean would cause) and the standardisation as
  a per-column scale; the CV folds are 0/1 masks and the held-out rows are
  predicted row by row. No allocation of the size of X: 1.09 X -> 0.09 X
  (path) and 1.37 X -> 0.15 X (CV) on the wide design. Same estimator,
  same grid, same deviance rules; the two engines agree to rounding (about
  1e-12 in the coefficients; the tests pin 1e-7), so a CV pick can differ
  only on a flat curve. The default engine is unchanged and bit-identical.
- Ridge CV (`ridge_cv_svd(bool low_memory)`; R `MLRidgeCV$new(low_memory
  = TRUE)`; Python `RidgeCV(low_memory=True)`): the dual (kernel) fold
  engine. The n_tr x n_tr kernel of the standardised training fold and the
  n_test x n_tr cross kernel are accumulated from column blocks of X
  standardised on the fly (read in place), eigen-decomposed, and the ridge
  path is evaluated on the held-out rows as K_te U diag(1/(s^2+lambda)) U^T y
  — no gathered fold, no p-dimensional coefficients, the SVD route's rank
  rule. 3.64 X -> 0.19 X on the wide design; memory O(n^2) instead of
  O(n p), the route for p >> n. Agreement with the SVD route to 1e-7 on
  the CV curve and identical picks in the tests; the kernel carries the
  squared condition number, so on ill-conditioned tall designs the default
  SVD route is the more accurate one.
- SVD dispatch thresholds exposed (`SVDSolver(min_dim_thresh, block_size)`
  in Python, `MLSVD$new(min_dim_thresh, block_size)` in R; defaults
  unchanged): lowering `min_dim_thresh` selects the randomized, block-wise
  path — two passes over X, no copy — for smaller designs; it is the route
  for a large memory-mapped X and is accurate for a decaying spectrum (the
  tests pin 1e-4 on the retained singular values), not for a flat one.
- Ridge `solve` on a memory-mapped X needs no new engine: the Gram or
  kernel matrix is accumulated from the map (Tier A's `ridge_solve_mmap`).
- Peak-memory budgets (`test_ml_methods_memory.cpp`), wide 500 x 8000:
  ridge CV 3.64 X (kernel engine 0.19 X), elastic net 1.09 X (raw engine
  0.09 X), elastic-net CV 1.37 X (raw engine 0.15 X), SVD 0.35 X, ridge
  0.31 X, PCA 0.35 X, scalers 0.03 X.

#### R surface catch-up and small utilities (5.4, 5.5)

- T-solvers (R6 base): `get_residuals()`, `get_beta_path()` (step 0
  included, Python's layout), `get_cp_with_sigma()`, `solver_type_to_string()`,
  `get_dropped_indices()`, `get_active_predictor_indices()`,
  `get_active_dummy_indices()`, `get_num_actives()`, `get_dummy_start_index()`,
  `is_connected()`, `set_tolerance()`, `set_tie_seed()`,
  `set_exchangeable_tie()`, `get_warnings()` (the already-preprocessed
  warning is readable now), `clear_warnings()`. Family getters:
  `TLARS_Solver$get/set_kahan_refresh_interval()`, `TENET_Solver$get_alpha()`,
  `TIKNET_Solver` / `TIKOMP_Solver$get_tikhonov_k()` (a `dgCMatrix`),
  `TCDSolver_Base$get_kkt_violation()`. `reconnect` is not bound in R: an
  R6 solver is bound to its matrices at construction.
- Selectors: active bindings `tFDR`, `K`, `max_num_dummies`, `max_T_stop`,
  `stagnation_check`, `num_selected`, `num_dummies` (new C++ getter
  `TRexSelector::getNumDummies()`, also bound in Python as
  `TRexSelector.num_dummies`); DA `cor_coef`; GVS `hc_method_used`.
- ML: `MLElasticNet$n_nonconverged()`; Python `RidgeCV.interpolates()`,
  `RidgeCV.selection_at_grid_edge()`, `ElasticNetCV.selection_at_grid_edge()`,
  `ElasticNetCV.grid_truncated()`.
- Utils: `trex_get_max_threads()` / `trex_set_num_threads()` (the OpenMP
  thread count of the backend); `mmap_matrix_size` returns a double (an
  int overflowed above 2^31 - 1 elements); `mem_matrix_from()` coerces an
  integer matrix explicitly instead of through a hidden binding copy.

#### Selectors: memory-mapped Biobank and SPCA, restore(), results wrapped once

- Biobank screening and sparse PCA on a memory-mapped X: the R6 classes
  used to read the whole file into RAM through `as.matrix()`. New entry
  points `trex_biobank_screening_{1d,2d}_mmap_create` and
  `trex_spca_select_mmap` map the file directly (readwrite: the inner
  selectors normalize / center it in place and restore it); the R6 classes
  route `mmap_matrix` inputs to them. New R tests compare the mmap and
  in-memory selections exactly.
- `TRexSelector::restore()` (C++, R `$restore()`, Python `restore()`):
  gives X back before `select()` — the constructor normalizes X in place,
  and until now an unselected selector handed it back only from its
  destructor (in R and Python: the garbage collector). Idempotent, a no-op
  after `select()`; a restored selector refuses to run (`select()` throws),
  including the screening selector's own `select()`. Inherited by every R6
  selector class.
- Results wrapped once: the selector's `phi_mat`, `fdp_hat_mat` and `r_mat`
  (T_stop x p) were copied twice per access; the DA `*_array_BT` results
  (rho x T_stop x p) three times (a vector of matrices by value, the wrap,
  `simplify2array`). They are now wrapped once from the const reference,
  the DA arrays built directly as R arrays of the same shape.
  `TRexSPCA::select()` returns a const reference to the stored result, so
  the bindings wrap Z and V once (`auto res = spca.select()` still copies).
- `trex_da_to_base_ptr` (unused from R) pins the owning pointer in its
  protection slot instead of returning an unprotected non-owning view.

#### T-solvers: shared-buffer claim and the in-place writer/reader contract

- Finding (Examples repo, `demo_ts_05_tstagewise` in C++): the checkpoint
  round trip reported "Paths differ" from step 21 on. Not a solver bug: the
  demo constructed the reference solver, the checkpoint solver and the
  resumed solver with `normalize = true, intercept = true` on the SAME
  Eigen maps, so X and D were centered and scaled in place three times. A
  second in-place transformation is not idempotent (it moves the data by
  rounding, RSS difference 9e-13), and the stagewise near-ties downstream
  resolve differently. Pre-scaling once in front and running the solvers
  as readers gives bit-identical paths.
- The contract, now written on `TSolver_Base::restore()` (C++), on every
  Python constructor and on `tsolver_options()` / `TSolver_Base` (R): a
  solver constructed with `normalize = true` or `intercept = true` is a
  WRITER that transforms X and D in place, keeps the statistics, and holds
  a claim on both buffers until `restore()` or destruction, both of which
  write X and D back; a solver constructed with both false is a READER that
  runs as constructed and never touches the buffers; a solver resumed from
  a checkpoint is a reader that expects the buffers in the saved state (so
  the writer that saved it stays alive, or the pipeline pre-scales once in
  front and runs readers); y is never written (the solver centers an
  internal copy).
- The claim is executable: `tsolvers/ad/solver_utils/buffer_registry` (a
  mutex-guarded set of buffer pointers, RAII `BufferClaim`) is shared by the
  T-solvers and the T-Rex selector, which moved its private registry there.
  A second writer on a claimed buffer, T-solver or selector, throws
  `std::runtime_error` (Python `RuntimeError`, R error) with the remedy:
  pre-scale once and construct readers, `restore()` the holder, or let it
  go out of scope. A writer that finds its input already centered and
  unit-norm (every column within 1e-10) records a warning (`getWarnings()`)
  and continues.
- `restore()` is idempotent and a no-op for readers and resumed solvers
  (`holdsPreprocessed()`); the former unconditional inverse scaled and
  shifted the data a second time on a second call. It now also checks that
  the passed X and D are the buffers the solver was constructed on. The
  destructor writes back a writer that was never restored, which also
  covers a derived constructor that throws after the base preprocessed
  (`TENET` with a negative `lambda2` used to leave X normalized). The
  write-back duty travels with a moved solver.
- R binding rewritten zero-copy: the R6 T-solvers used to copy X and D into a
  private holder (`TSolverRcpp`, also for a memory-mapped X, which was read
  into RAM). They now map the R matrix (double storage) or the mmap file
  directly, pin the R objects in the XPtr protection slot, and follow the
  contract: a writer preprocesses the caller's matrix or file in place, the
  new `$restore()` writes it back and releases the claim, the GC finalizer
  writes back a writer that was never restored, `$load()` into a writer
  keeps its write-back duty (the archived statistics are the same). The
  stale copying holder in `TRexSelectorNeo_types.h` (a second definition of
  the same struct, included by `RcppExports.cpp`) is the single definition
  now.
- R `get_actions()` fixed: the C++ solver already records the 1-based
  signed lars codes, and the R6 method shifted them by one more, so
  `r_index = TRUE` returned 2-based codes and `r_index = FALSE` the 1-based
  ones it documented as 0-based. `r_index = TRUE` (default) now returns the
  codes as recorded; `r_index = FALSE` returns the 0-based magnitudes
  (`abs(a) - 1`, sign dropped: a 0-based signed code cannot distinguish the
  addition from the removal of the first variable).
- Tests: `test_tsolvers_buffer_guard.cpp` (registry primitives; second
  writer throws across families; intercept-only and normalize-only claim;
  `restore()` idempotent, readers and resumed solvers never write back;
  destruction writes back, also after a throwing derived constructor; moves
  transfer claim and duty; readers and loaded solvers never claim; the
  already-preprocessed warning; the selector/solver shared registry; the
  pre-scaled shared-buffer checkpoint round trip is bit-identical). The CD
  and Tikhonov checkpoint tests keep the writer alive across the resumed
  run. R: `test-tsolvers-inplace.R` (in-place change visible on the R matrix
  and in the mmap file, second writer errors, readers coexist, `$restore()`
  idempotent and path kept, GC write-back, `$load()` into a writer restores,
  pre-scaled readers bit-identical); `test-tsolvers.R` restores each writer
  and runs the 22-solver checkpoint round trip on pre-scaled readers with
  exact equality. Three Python tests that rebound a writer in a loop now
  release it first.

#### Linux verification (OrbStack): two platform-bound tests, broken Rd pages

- The three Ubuntu CI jobs were reproduced in an `ubuntu:24.04` container
  (GCC 13, Eigen 3.4, R 4.3.3, Python 3.12) on the committed tree, plus the
  install smoke, `R CMD check` and the cibuildwheel manylinux_2_28 x86_64
  wheels with their in-container test run. Everything passed except:
- Two C++ tests calibrated on macOS data: `std::normal_distribution` and
  `std::rand` (behind `Eigen::Random`) are implementation-defined, so the
  same seed draws different data under libc++ / macOS libc and libstdc++ /
  glibc. `EnetCvCcdTest.KeepsGlmnetLambdaRuleAndDetectsGridFloorPin` now
  asserts the wide-versus-pinned ratio with a margin both platforms clear
  (1.2; observed 1.9 on macOS, 1.32 on Linux), and
  `TSolverTikTest.SerializeLoadContinueAcrossDropRecovery_TIKNET` probes a
  small seed range until a design yields a checkpoint mid-recovery instead
  of relying on seed 12. Both pass on both platforms.
- Three R help pages were unparseable: an apostrophe in
  `\code{lambda2 * beta' K beta}` opens an R string literal in Rd (the
  TIKNET, TIKOMP and TCDTIKNET pages; now `t(beta) K beta`), a bare
  `[0.25, 4]` became a broken link and `\ldots` sat inside a `\code{}` on
  the `trex_control` page. CI never saw this because the R job installs
  with `--no-docs` and does not run `R CMD check`; the check step is added
  to the workflow.

#### Python: one subpackage per dummy mechanism, and the Tikhonov selector

- `trex_selector_neo.tsolvers.ad.{lars_based, omp_based, afs_based, cd_based}`
  with `SolverOptions`, `ScalingMode` and `GroupPenaltyVariant` under
  `tsolvers.ad`, and `trex_selector_neo.trex_selector_methods.ad` with every
  selector, mirroring `cpp/src/{tsolvers, trex_selector_methods}/ad`; the
  compiled extension carries the same level (`_core.tsolvers.ad`,
  `_core.trex_selector_methods.ad`). The package root keeps re-exporting the
  AD selectors as the everyday spelling; the `vd` / `sd` families will be
  reachable only through their subpackage. Breaking for
  `from trex_selector_neo.tsolvers.lars_based import ...` and
  `from trex_selector_neo.tsolvers import SolverOptions` (insert `.ad`); R
  needs nothing (flat namespace, group words).
- `TRexTikhonovSelector` is bound: `TRexTikhonovControlParameter` (the p x p
  penalty as a dense `tikhonov_K` array converted to sparse internally, or
  `GROUP_CENTERING`; the lambda_2 rule; the pilot knobs; cluster-aware
  dummies), the enums `TikLambda2Method`, `PilotStatistic`,
  `TikPenaltyStructure`, the records `Lambda2PilotResult` / `Lambda2PickInfo`,
  the builders `gamma_to_k` and `build_group_centering_k`, and the wrapper
  with `lambda2_used`, `lambda2_cv_diagnostics` and `lambda2_pilot`. The base
  control merges as in R (a non-Tikhonov `solver_type` warns and becomes
  TCDTIKNET). Every selector R has is now in Python.

#### GVS K-loop: per-experiment tie-break seed

- The in-memory GVS K-loop now derives a deterministic `tie_seed` per
  experiment from the selector seed, as the base ExperimentRunner does.
  Inert unless the exchangeable-tie band is on (`exch_tie_alpha > 0`); with
  it on, a seeded GVS run is reproducible instead of drawing tie-breaks from
  `random_device`.

### 2026-09-24

#### Voting ties: exact vote counts and a single-division voting grid

- The classical experiment runner accumulated Phi as K additions of 1/K. The
  rounding residue (15 x 0.05 = 0.75 + 1 ulp) decided the strict `Phi > v`
  comparison at the ties the voting lattice makes routine, because the grid
  `0.5, 0.5 + 1/K, ...` lives on the same 1/K lattice by design (R's
  `seq(0.5, 1 - eps, by = 1/K)`). The runner therefore selected a variable
  that sat exactly on a grid point; the GVS K-loop, which divides exact
  counts once, and R's `rowMeans()` did not. One such tie flips the L-loop
  stop, which is how `TRexTikhonovSelector` and `TRexGVSSelector` came to
  stop at different L on the same dummies (the divergence pinned since
  2026-08-25; the dummies and the per-experiment paths were bit-identical).
- Fix: `ExperimentRunner::aggregateResults` accumulates exact counts and
  divides by K once (Phi is the correctly rounded m/K), and
  `createVotingGrid` computes every point as one division (K + 2i)/(2K), so
  a mathematical tie yields the identical double and is excluded, as the
  definition and R do. Affects every selector on the runner (classical, DA,
  screening, SPCA, Tikhonov) only where a variable's vote share equals a
  grid point exactly; such variables are no longer selected at that
  threshold. The equivalence test beyond L = p is enabled (decision 9 in the
  plan is closed: there was no convention difference).

#### Beyond architecture: one folder per dummy mechanism (breaking, C++ only)

- `cpp/src/tsolvers/` and `cpp/src/trex_selector_methods/` each hold three peer
  folders, `ad/` (augmented dummies: the classical T-solvers and every T-Rex
  selector, with `trex_utils`), `vd/` (virtual dummies: `vd/` solvers, `trex_vd`)
  and `sd/` (sparse dummies: `sd/` solvers, `trex_sd`). The families share no
  code, and no class serves more than one mechanism; a fourth mechanism would be a
  fourth sibling. Placement rule: a file lives at the level of its lowest common
  user (`TSolver_Base`, `SolverOptions`, `solver_utils/` and all of `trex_utils/`
  are AD-only today; `VD_Base` and `SDTSolver_Base` are standalone roots).
- Namespaces follow the folders: `trex::tsolvers::ad::{linear_model::…,
  solver_utils::…, TSolver_Base, SolverOptions, ScalingMode, …}`,
  `trex::tsolvers::vd`, `trex::tsolvers::sd` (the SD solvers were declared inside
  the AD subnamespaces before; `sd_calibration` is `sd::calibration`),
  `trex::trex_selector_methods::ad::{trex_core, trex_da, trex_gvs, trex_screening,
  trex_spca, trex_tikhonov, utils::…}`, `…::vd::trex_vd`, `…::sd::trex_sd`.
- Umbrellas: `tsolvers/ad.hpp` (was `lm_tsolvers.hpp`), `tsolvers/vd.hpp`,
  `tsolvers/sd.hpp`; `tsolvers/tsolvers.hpp` still means the classical family.
  Installed header paths move with the tree (`TRexSelector/tsolvers/ad/…`,
  `TRexSelector/trex_selector_methods/ad/…`).
- The AD alias for the solver dispatch is `disp` (was `sd`, which now names the
  sparse-dummy family). Tests and validation programs mirror the layout
  (`test_tsolvers/{ad,sd}/`, `test_trex_selector_methods/{ad,vd}/`,
  `validation/{tsolvers,trex_selector_methods}/ad/`).
- R and Python: the glue includes and namespaces follow; no exported name changes.
  `bundle_for_cran.sh` excludes the `vd` and `sd` folders and umbrellas (four
  patterns, was five).

#### R namespace: one prefix per group (breaking)

- The exported symbols now follow one rule: functions and constant lists carry the
  group prefix (`trex_`, `tsolver_`, `ml_`, `mem_`), R6 classes the group word
  (`TRex*Selector`, `T*_Solver`, `ML*`). Renamed: the ML classes `ZScoreScaler`,
  `LpNormScaler`, `RidgeSolver`, `RidgeCV`, `ElasticNet`, `ElasticNetCV`,
  `SVDSolver`, `PCA` to `MLZScoreScaler`, `MLLpNormScaler`, `MLRidge`, `MLRidgeCV`,
  `MLElasticNet`, `MLElasticNetCV`, `MLSVD`, `MLPCA`; `agglomerative_cluster`,
  `cut_tree`, `LinkageMethod`, `DistanceMetric` to `ml_agglomerative_cluster`,
  `ml_cut_tree`, `ml_linkage_method`, `ml_distance_metric`; the metrics
  `compute_fdp`, `compute_tpp`, `compute_precision`, `compute_recall`,
  `compute_fdp_dense`, `compute_tpp_dense` to `trex_fdp`, `trex_tpp`,
  `trex_precision`, `trex_recall`, `trex_fdp_dense`, `trex_tpp_dense`;
  `gamma_to_k` to `trex_gamma_to_k`; `mmap_matrix` and `convert_to_memory_mapped`
  to `mem_matrix` and `mem_matrix_from` (the S3 class `mmap_matrix` and its methods
  are unchanged). `TSolver_Base`, `TCDSolver_Base` and `get_ptr` are internal now.
- Every help page lists its group ("Other T-solvers: ..."), `?TRexSelectorNeo`
  presents the four groups, and `_pkgdown.yml` indexes the reference by them.

#### The CD family, native TIENET, TIKNET and TIKOMP as R6 and Python solvers

- R: `TCDLASSO_Solver`, `TCDENET_Solver`, `TCDIENET_Solver`, `TCDTIKNET_Solver`
  (sharing the CD getters through `TCDSolver_Base`), `TIENET_Solver`,
  `TIKNET_Solver`, `TIKOMP_Solver`, all on `$new(X, D, y, [groups | tikhonov_K,]
  options)`; a Tikhonov matrix is a `Matrix::dgCMatrix`, another `Matrix` sparse
  class, or a dense matrix.
- Python: the `tsolvers.cd_based` submodule with the four CD solvers;
  `TIENET_Solver` and `TIKNET_Solver` in `lars_based`, `TIKOMP_Solver` in
  `omp_based`; a Tikhonov matrix is a dense array (converted internally, no SciPy
  dependency).

#### C++ install layout (breaking for manual include paths)

- The SDK headers now install under a namespaced root,
  `<prefix>/include/TRexSelector/{ml_methods,trex_selector_methods,tsolvers,utils}`,
  instead of directly under `<prefix>/include/`. Consumers that link the exported
  targets (`find_package(TRexSelector)`, `TRexSelector::trex_selector_methods`, ...)
  need no change; a manual `-I<prefix>/include` becomes `-I<prefix>/include/TRexSelector`.
  The pre-install prune that kept the header tree free of stale files now removes only
  that namespaced root, and only when it carries the `TRexSelector.install-stamp` written
  by a previous install: a foreign directory of the same name is refused, never deleted.
  Previously the prune deleted the generic `include/utils`, `include/ml_methods`, ...
  under the prefix, which in a shared prefix (conda, /usr/local) could belong to
  another package. `cpp/tests/install_smoke/run_install_smoke.sh` installs into a
  scratch prefix and builds a `find_package` consumer against it. Its first run
  found seven headers that the code includes but no install file set carried
  (`lambda2_cv_diagnostics.hpp`, `post_selection_criteria.hpp`,
  `trex_cluster_dummies.hpp`, `trex_cluster_hac.hpp`, three `sd_tsolvers` headers):
  an installed SDK could not compile a consumer since August. Fixed.
- New CMake option `TREX_INSTALL_CPP_SDK` (default ON): headers, static archives and
  the CMake package files. The Python build forces it OFF, so wheels no longer ship
  about 16 MB of C++ SDK into `site-packages/include` and `site-packages/lib`.

#### Hardening

- Agglomerative clustering rejects non-finite input at its entry with a message naming
  the entry (a NaN distance made the nearest-neighbour search index out of range);
  the NN chain also refuses its own sentinel. A zero-norm column keeps its documented
  convention (maximally dissimilar, distance 1).
- `cut_tree` / `cut_tree_by_height` validate a merge matrix before indexing it: row
  count, cluster ids that exist at their step, no self-merge, no id consumed twice,
  finite and (for the height cut) non-decreasing distances. Reachable through the
  public R and Python `cut_tree`.
- Overflow-checked byte counts (`utils/checked_arith.hpp`) for the memory-mapped
  matrix, the tiled distance matrix and the binary Eigen serialization: a wrapping
  `rows * cols * sizeof` throws `std::overflow_error` before any file is touched.
- Eigen's process-global thread count is set once per experiment runner from the
  calling thread and restored on exit, never written inside an OpenMP team. Before,
  the sequential runners left it changed in the host process (an R or Python session)
  and the workers raced on it. Only observable in builds without `EIGEN_USE_BLAS`
  (the R package and the Python wheel); with a BLAS, Eigen's count is inert. The
  thread-hygiene tests watch two globals: the host thread's OpenMP thread count in
  every build (the parallel runners set a per-worker count inside their team, and it
  must not reach the caller), and Eigen's count in no-BLAS builds, where it is live.

#### Scaling mode NONE: the user-managed design

- `ScalingMode::NONE` joins `L2` (default, unchanged) and `ZSCORE`. Under NONE the
  selectors take `X` and `y` exactly as supplied: no centering, no scaling, the
  stored normalization parameters are the identity (means 0, divisors 1, `y_mean` 0)
  and nothing is written back into the caller's buffer. The caller pre-centers and
  pre-scales as the model demands, or supplies an intercept-free dictionary
  (steering vectors, Vandermonde, Fourier atoms) whose columns must not be centered
  and may contain a constant atom. The dummies keep the unit-L2 standard, so a NONE
  design should be unit-L2-commensurate; one warning names how many columns have an
  L2 norm outside [0.25, 4] (z-scored or raw-amplitude data). A design pre-normalized
  with the library's own helper and run under NONE selects exactly what the raw design
  selects under L2 (bitwise, test-enforced in C++, R and Python). The solver-side
  convention and the lambda_2 tuners treat NONE as unit-L2; the correlation HAC of the
  cluster-dummy selectors then measures cosine coherence on the columns as supplied.
- The selectors' preprocessing is callable: `trex_normalize_design(X, y, scaling_mode)`
  in R and `trex_selector_neo.normalize_design(X, y, scaling_mode)` in Python return
  the normalized copies and the parameters. R exposes the mode as
  `trex_control(scaling_mode = "L2" | "ZSCORE" | "NONE")` (it had no scaling switch).
- The normalizer rejects bad input in every mode with the column named (0-based):
  a non-finite entry in `X` or `y`, and a dead column — constant under L2 / ZSCORE,
  all-zero under NONE (`std::invalid_argument`; `ValueError` in Python; an R error).
  Before, a constant column was silently pinned to scale 1 and could never be selected.
  Designs that carry constant columns must drop them now (breaking, deliberate:
  PLAN WP4 4.1 (d)). The matrix is unmodified when the constructor throws.

#### T-solver constructors take one options structure (breaking)

- Every T-solver is constructed as `Solver(X, D, y, options)`, with
  `trex::tsolvers::SolverOptions` (`tsolvers/solver_options.hpp`) carrying every other
  input: the preprocessing flags (`normalize`, `intercept`, `scaling_mode`, `verbose`),
  the tolerance, the tie-break seed and the exchangeable-tie configuration, and the
  family hyperparameters (`lambda2`, `rho_afs`, `ncgmp_variant`, `tenet_aug_use_lars`,
  `kappa_dummy`, `fold_dummies`, `group_penalty`, the `cd_*` knobs). The structural
  inputs that are data, not options, stay explicit: the group vector of the IEN solvers
  (`TIENET`, `TIENETAug`, `TCDIENET`) and the Tikhonov matrix of the sparse-K solvers
  (`TIKNET`, `TIKOMP`, `TCDTIKNET`) follow `y`. Each family reads the fields it needs;
  the defaults are the former solver-level defaults (`lambda2` defaults to 0, the lasso
  path, where it used to be a required argument). The positional constructors with
  their boolean triple and trailing enum are gone; the `NCGMPVariant` and
  `GroupPenaltyVariant` enums moved next to the structure (their old names still
  resolve). The T-Rex layer maps its `SolverHyperparameters` onto the structure in one
  place of the dispatch, which lost its per-family constructor table and its
  post-construction setters; selections are unchanged (the full suite is bitwise the
  same on the anchors).
- R: `tsolver_options()` builds the list every solver's `$new(X, D, y, options = ...)`
  takes (`TIENETAug_Solver$new(X, D, y, groups, options = ...)`); the per-solver
  arguments `normalize`, `intercept`, `verbose`, `lambda2`, `rho`, `variant` and
  `use_lars_inner` are gone (`rho` is `rho_afs`, `variant` is `ncgmp_variant`).
- Python: `trex_selector_neo.tsolvers.SolverOptions` with a keyword constructor
  (`SolverOptions(lambda2=0.5, normalize=False)`), passed as the fourth positional or
  the `options=` keyword; `GroupPenaltyVariant` exported next to it. The former
  keyword arguments of the solver constructors are gone.

#### lambda_2 tuner boundary guard (warning-only)

- The GVS and Tikhonov selectors warn once when a CV tuner returns a
  boundary-dominated lambda_2, and return the pick unchanged. Two modes are
  covered: the returned lambda_2 sits on the first or last point of its grid (the
  grid boundary, not a CV optimum: a flat or still-improving curve, typically low
  SNR), and, for the IEN tuner, the profiled lambda_1 optimum of the chosen
  lambda_2 sits on the truncation edge of the lambda_1 path. Before, both were
  visible only in `n_profile_at_edge`, which is blind to the grid-ceiling mode
  (measured in August on the demo designs: a pick at the ceiling with the counter at
  zero). The warning names the selector, the method, the mode and the remedies
  (1SE over MIN, a fixed lambda_2, a longer lambda_1 path, a finer grid). For the
  EN CCD tuner a last-point pick on a path that glmnet's deviance rules truncated
  is benign and is not flagged.
- The lambda_2 diagnostics record (C++ `lambda2_cv_diagnostics`, R
  `lambda2_cv_diagnostics`, Python `Lambda2CVDiagnostics`) carries the pick:
  `grid_index` / `grid_size`, `at_grid_edge`, `profile_at_edge` and
  `boundary_pick`. Every CV method fills the grid fields now, not only the IEN
  tuner. `ienet_cv_ccd` gained `selection_at_grid_edge()` and a per-lambda_2
  `profile_at_edge(k)`; its `diagnostics()` takes the rule (MIN or 1SE) whose pick
  the caller uses.

### 2026-07-21

#### R package refreshed: vendored core synced, exchangeable tie-break exposed

- The vendored C++ core is byte-identical with `cpp/src` again (tie-break
  additions in `tsolver_base`, the solver dispatcher, TOMP/TAFS, plus the
  cosmetic trex_da drift). `trex_control()` gains `exch_tie_alpha` (default 0)
  and `exch_tie_floor` (default 0.5), parsed into `SolverHyperparameters` by
  the Rcpp layer; man pages regenerated via roxygen. New testthat coverage:
  control defaults and tie-broken TOMP/TAFS selection runs on a correlated
  cluster design. Full testthat suite green.

#### Python package refreshed for the exchangeable tie-break API

- `SolverHyperparameters` now exposes `exch_tie_alpha` / `exch_tie_floor`, and the
  solver wrappers (TOMP/TAFS/TLARS families) gain `setExchangeableTie(alpha, floor=0.5)`;
  type stubs updated accordingly (the PRIOR_GROUPS binding changes shipped with the
  core rewrite). New tests: hyperparameter defaults/assignment, reproducible
  tie-broken paths + alpha = 0 no-op on TOMP/TAFS, data-driven rho grid
  (ascending, rho = 1 anchor) and negative-label rejection on the prior-groups
  DA path. Suite 453/453.

### 2026-07-17

#### TLARS "No valid gamma" warning classified and demoted

- Investigated the frequent `No valid gamma found, using fallback` warnings seen in the
  trex_da demos. Instrumentation on the exact demo DGPs/seeds showed **100 % of observed
  events are the benign path-exhaustion case**: on small designs
  (p + num_dummies <= n - 1, e.g. demo 04's p=25 cells) the T-loop drives the LARS path
  until every real + dummy variable is active, and `computeStepSize()` then runs with an
  empty inactive set — the terminal LS step `Cmax / A_A` is the correct step. A Gaussian
  control reproduced the events at the same rate as the heavy-tailed DGP (173 vs. 166 per
  40 trials), and the events vanish exactly at the exhaustion boundary — not heavy tails,
  not collinearity, not the demo RNG (seeds are `base + mc`, fresh mt19937 per trial).
  Demo 05 (p_total=500) emits none; its earlier attribution was console interleaving.
- `computeStepSize()` now logs the empty-inactive-set case verbose-only
  (`logInfo("Inactive set empty; taking terminal LS step.")`) and keeps a `logWarning`
  (with candidate count) only for the never-observed all-candidates-filtered case.
  Fallback step unchanged. Suite 291/291; R mirror synced. Details in
  `TRex_Research/documentation/CTRexSelector_Open_Tasks.md` §7 (CTLARS should gain the
  same classified logging when synced — it is currently silent in both cases).

#### PRIOR_GROUPS redefined: constrained sub-clustering (deflation-mismatch fix)

- `setupDA_PriorGroups()` no longer uses the user-supplied groups as deflation
  neighbourhoods. The nearest-partner penalty `delta = 2 - min|Phi_j - Phi_k|`
  (Definition 1 of the T-Rex+DA paper) assumes tight neighbourhoods in which
  the closest-Phi member is a genuine shadow; inside a large prior group nearly
  every variable has SOME member with an almost identical Phi, so delta -> 2
  uniformly and the shadow/active ranking is left invariant (the FDR-plateau /
  TPR-cap signature documented in
  `TRex_Research/documentation/Prior_Groups_Deflation_Mismatch_DA_TRex.md`).
- New semantics: the prior groups are hard merge CONSTRAINTS on the dependency
  structure. The finest common refinement of all supplied `prior_groups` levels
  partitions the variables; HAC (`hc_linkage`, correlation distance) runs
  within each constraint group (merges across boundaries are structurally
  impossible); the rho grid (length `hc_grid_length`) is subsampled from the
  pooled within-group dendrogram heights and terminated by the conservative
  rho = 1 singleton anchor (the paper's Psi = 1/2 empty-group penalty). The
  BT-style 3D calibration then operates on tight data-driven neighbourhoods
  INSIDE the known groups — the coarsest grid point recovers at most the whole
  prior group, the finest is provably conservative.
- `rho_grid_labels` removed from `TRexDAControlParameter`, the R
  `trex_da_control()`, and the Python bindings — the grid is now data-driven.
  `prior_groups` labels must be non-negative (validated). Linkage validation
  now covers PRIOR_GROUPS as well as BT.
- New guard tests in `cpp/tests/test_trex_selector_methods/test_trex_da.cpp`:
  boundary containment + rho = 1 anchor, finest-common-refinement semantics,
  full select() smoke run, negative-label rejection.

### 2026-07-15

#### Exchangeable tie-breaking for greedy solvers (DA FDR fix)

- New solver hyperparameters `exch_tie_alpha` (default 0 = off) and
  `exch_tie_floor` (default 0.5) in `SolverHyperparameters`, consumed by
  TOMP/TAFS via `TSolver_Base::applyExchangeableTieBreak()`. When alpha > 0
  and the step's top candidate is a real predictor, inactive real predictors
  with |corr| >= floor to the top candidate whose correlation gap lies within
  `alpha * ||x_j* -/+ x_k|| * ||r|| / sqrt(n)` are treated as statistically
  exchangeable, and ONE member is picked uniformly at random per step (via
  the existing tie-break RNG, so `setTieSeed()` keeps runs reproducible).
  Dummy picks are never randomized — dummy budget and early stopping are
  unchanged.
- Motivation: greedy solvers are winner-take-all — within a T-Rex random
  experiment the in-sample winner of a highly correlated cluster is a
  deterministic function of (X, y) (dummies never perturb the real-variable
  ranking), so occurrence mass concentrates on one cluster member and the
  DA deflation `delta = 2 - min|Phi_j - Phi_k|` degenerates to the identity:
  collinear shadows that beat their active in-sample survive as false
  positives (DA-AR1 FDR ~ intrinsic shadow-win rate, e.g. ~34% per cluster at
  rho = 0.9, SNR = 2, n = 300 — above tFDR = 0.2). Randomizing statistically
  indistinguishable picks restores the within-trial occurrence spread across
  exchangeable cluster members that the DA FDR control relies on (LARS-type
  solvers produce this spread through their path geometry and are unaffected).
- Validated on the demo_trex_da_01 AR(1) cells (rho = 0.9, gap = 100 and
  Random support): alpha = 0.25 restores FDR < tFDR for TOMP/TAFS while
  retaining more power than TLARS; the safe regimes (rho <= 0.7) are
  essentially unchanged. Recommended: `exch_tie_alpha = 0.25` for greedy
  solvers under trex+DA.
- The dispatcher re-applies the parameters on all three solver lifecycle
  paths (fresh, serialized warm start, in-memory warm start), mirroring
  `setTolerance`; the fields are intentionally not serialized.

### 2026-07-08

#### R <-> Python binding parity

- BREAKING (R): `TRexBiobankScreeningSelector$select()` now returns the same
  per-phenotype record shape as the Python binding. A matrix `Y` previously
  returned a list with a `$statistics` data.frame plus parallel top-level
  `$selected_indices` / `$selected_indices_screen_ordinary` /
  `$selected_indices_screen_bootstrap` lists; it now returns a plain list of
  records, one per phenotype, where `res[[i]]` has the identical shape to the
  single-phenotype (vector `Y`) return. No information is lost — every
  `$statistics` column and index list is a field of each record — and the
  summary table is a one-line projection:
  `do.call(rbind, lapply(res, function(r) as.data.frame(r[c("phenotype_index", "estimated_FDR", "method_used", "used_fallback_trex")])))`.
  This matches Python, which returns one `BiobankScreenTRexResult` (vector) or a
  list of them (matrix).
- R gains `prior_groups` + `rho_grid_labels` in `trex_da_control()` — the
  prior-groups dependency-aware path (`setupDA_PriorGroups`), which Python
  already exposed; the R glue previously could not route to it. `hc_grid_length`
  is now actually read by the DA glue (it was silently inert before).
- R gains `cv_seed` / `cv_n_folds` / `cv_n_lambda` in `trex_gvs_control()` — the
  cross-validation fold controls for the CV `lambda2_method` paths (`cv_seed`
  enables per-trial fold reproducibility), at parity with Python.
- `TRexGVSSelector` gains read-only `$groups` (1-based per-variable cluster
  labels, length `p`) and `$group_labels`, so cluster diagnostics (e.g. purity)
  can be computed R-side; Python exposes the same via `GVSSelectionResult`.

### 2026-07-07

#### Shared-buffer safety guard (new runtime check)

- `TRexSelector` (and every derived selector) normalizes the input matrix `X`
  in place via an `Eigen::Map` and restores it when `select()` finishes or on
  destruction. Two selectors mapped onto the SAME buffer therefore silently
  corrupted each other: the first to restore left the buffer un-normalized, so
  the second ran on raw data (typically selecting everything). The core now
  claims a buffer when it normalizes it and throws a clear error if a second
  live selector is constructed on the same matrix. Fix each caller by giving
  each selector its own copy of `X`, or by running selectors one at a time
  (construct, `select()`, and let each go out of scope before the next). The
  guard is in the C++ core, so C++, R, and Python all benefit; internal
  sequential orchestrators (SPCA per-PC, biobank screen->fallback) are
  unaffected because `select()` releases the claim before the next selector
  is built.

#### R <-> Python binding parity

- Elastic-net model selection is now exposed in BOTH languages: `ElasticNet`
  (glmnet-style coordinate-descent coefficient path, backed by `enet_gaussian`)
  and `ElasticNetCV` (K-fold CV lambda.min/1se selection, backed by
  `enet_cv_ccd`). `RidgeCV` (SVD-path ridge CV, `ridge_cv_svd`) is added to
  Python, which previously exposed no model-selection CV at all. Verified
  R and Python produce bit-identical elastic-net paths on shared data.
- R gains `da_method = "NN"` in `trex_da_control()` — the nearest-neighbour
  correlation-threshold DA sweep (`DAMethod::NN`), which Python already
  exposed; the R glue previously rejected it. `hc_grid_length = 0` lets the
  core auto-select `min(20, p)`.
- R gains `tenet_aug_use_lars` in `trex_gvs_control()` and
  `trex_spca_control()` — the LARS-vs-LASSO inner-path switch for `TENET_AUG`
  (Python already had it). Documentation clarifies that `en_solver`
  (`TENET` vs `TENET_AUG`) is a solver-EQUIVALENCE knob (both select the same
  variables), NOT the LARS-vs-LASSO distinction.
- `TENETAug_Solver` / `TIENETAug_Solver` are exported from the R package (were
  internal-only), with plain and mmap creators, 1-based group accessors, and
  save/load.

#### GVS / SPCA / clustering fixes

- BREAKING (R): `gvs_type = "IEN"` now derives the `TIENET_AUG` solver
  (matching the 2026-07-05 solver rework and the C++/Python wrappers); the R
  wrapper previously forced the old `TLASSO` solver.
- `TRexGVSSelector` now emits a warning when `control$solver` / `solver_type`
  is set to a value that does not match the solver derived from `gvs_type`
  (the field is silently overridden).
- `lambda_2 == 0` (degenerate no-ridge / pure T-LASSO) now emits a warning —
  the `-1` auto / `0` degenerate / `> 0` fixed sentinel convention had tripped
  in-house callers.
- Ward + `Correlation_LSH_Approx` is now permitted in both the R and Python
  clustering glues (the core dispatcher supports it, running Ward in a 64-D
  SimHash-projected Euclidean space; the glues wrongly rejected every
  non-Euclidean metric under Ward).

#### Python API

- `TRexBiobankScreeningSelector`: the control parameter `biosctrex_ctrl` is
  renamed to `bio_ctrl`, and `BiobankScreenTRexControl` now exposes its nested
  `trex_screen_ctrl` (and that control's nested `trex_ctrl`), so Python can set
  K, solver, screening method, etc., at parity with the R knobs.
- Scalers gain returning, layout-agnostic `transform()` / `fit_transform()` /
  `inverse_transform()` methods alongside the strict zero-copy
  `*_inplace` variants (which require a writeable Fortran-ordered float64
  buffer), removing the sharpest edge in the Python scaler API.

#### Core fixes

- `TENETAug_Solver` / `TIENETAug_Solver` deserialization: loading into an
  already-constructed solver (the R/Python bindings' `load()` path) now rebuilds
  the augmented maps and reconnects the inner solver unconditionally on load,
  instead of gating on a null-map check that only held for fresh-object loads —
  fixing checkpoint resume through the bindings.
- `cmake --install` now prunes the package-owned header trees before copying,
  so a renamed or deleted header can no longer poison the local install tree
  (a stale header that `#include`d a since-renamed file previously broke every
  consumer until the install was cleaned by hand).

### 2026-07-05

#### T-solver review fixes (tsolvers module)

- Full code review of `cpp/src/tsolvers/` with correctness fixes: stale
  correlations after T-Stagewise NNLS drops, missing dummy de-normalization in
  T-ENET's `getBeta`/`getBetaPath`, a data race and a scale-aware
  near-constant-column drop threshold in the preprocessing, guarded
  `getBeta(step)`/`getCp()` accessors, and a consistent
  actives/inactives/dropped partition across all solvers. `restore()` no
  longer corrupts the caller's `y` (the response is copied internally and was
  never centered in place).
- API: `executeStep()` with `T_stop = 0` now computes the full solution path as
  documented; `setTieSeed()` and `getWarnings()`/`clearWarnings()` are public
  and exposed in Python; warnings print unconditionally and are recorded
  per solver; TENETAug mirrors its inner solver's diagnostics, validates
  `lambda2`, and supports `save()`/`load()`.
- Performance: RAII guard restores Eigen's global thread count after solver
  runs; beta-path storage grows geometrically instead of per-step; dead
  Cholesky backup copies removed.
- T-ACGP now applies the conjugate correction across active-set growth steps
  via zero-padded previous directions (Blumensath & Davies 2008).
- Internals deduplicated (correlation refresh, tie scans, save/load,
  cycling-ratio tracking); tie-breaking RNG state, tolerance, and warnings are
  serialized. NOTE: the checkpoint format changed — `.bin` files saved with
  earlier builds cannot be loaded.

#### DA lifecycle refactor and NN-setup parallelisation

- `TRexDASelector` no longer overrides `select()`: the dependency-structure
  setup moved into an `onSelectBegin()` override and the DA selection /
  result widening into a `finalizeSelectionResult()` override, so DA now
  runs the inherited base pipeline like GVS does (behavior-neutral; the
  duplicated scaffolding — voting-grid setup, loop drivers, cleanup,
  X denormalization — existed twice before and could drift).
- The DA-NN neighbourhood setup is OpenMP-parallelised: each variable's
  neighbour list is built independently (race-free row-wise scan instead of
  the symmetric pair loop), with an early exit over the ascending rho grid.
- `corr_max` (GVS) documents that the default 0.5 matches the R reference
  and that the EUSIPCO 2022 paper's GWAS experiments used rho_thr = 1/3.

#### T-Rex+GVS rework: solver-encapsulated augmentation and the TIENETAug solver

- New `TIENETAug_Solver` (tsolvers, LARS family): the Informed Elastic Net
  (Machkour et al., CAMSAP 2023, Theorem 1) as a terminating forward selector.
  Builds the group-informed row-augmented system `[X; B]`, `[D; B tiled per
  dummy layer]`, `[y; 0_M]` (B(m, j) = sqrt(lambda2)/sqrt(p_m) on group-m
  columns) internally from a 0-based group assignment, applies the R-faithful
  post-augmentation column scaling, and runs an inner T-LASSO (or pure T-LARS)
  — full save/load support, registered in the polymorphic solver tests.
- `TRexGVSSelector` now hands every solver the plain (n x p) design, the plain
  (n x L) cluster-MVN dummy block, and the centered response; all augmentation
  lives in the solvers. `GVSType::IEN` routes through `TIENETAug_Solver`
  (new `SolverTypeForTRex::TIENET_AUG`, also in Python; IEN previously
  required TLASSO — breaking config change), and `en_solver = TENET_AUG`
  routes through `TENETAug_Solver`, which finally makes the
  `tenet_aug_use_lars` R-parity flag functional (it was silently ignored).
  The inline augmentation machinery (X_aug_/y_aug_ members, five build
  helpers, the variable effective row count, and the EN-aug mmap
  reinterpretation) is deleted; memory-mapped dummy files now always have n
  rows. EN-aug thereby switches from post-augmentation re-normalization to
  the demo-verified as-constructed convention that reproduces the Gram-based
  TENET path exactly (asserted by a new equivalence test).
- Review fixes: the GVS solver cache is cleared BEFORE the dummy buffers it
  views are reassigned (transient-dangling-view hardening, mirroring the
  base-class contract); EN dummy normalization goes through the shared
  data-normalizer (degenerate columns warn and keep scale 1.0 instead of
  throwing).
- New end-to-end GVS tests (grouped-recovery for EN/TENET, EN/TENETAug,
  EN/TENETAug with pure-LARS inner, IEN/TIENETAug, plus TENET==TENETAug
  selection equivalence). The IEN test pins lambda_2: the ridge-CV value
  (identical to the R reference's, which shares one CV and one * p/2
  conversion between EN and IEN) can make the sqrt(lambda2)-scaled group
  rows dominate the augmented geometry on n > p designs, leaving nothing
  selectable — a documented sensitivity of the IEN track (see
  computeLambda2()).

#### T-Rex selector review fixes (trex_selector_methods module)

- Warm starts now actually engage: `WarmStartManager` never stored solvers or
  checkpoint paths, so every T-loop step silently re-solved from scratch. The
  manager now retains solvers in-memory (co-retaining the dummy matrix for
  the PERMUTATION/DIRECT strategies, whose D_k is a per-call temporary) and
  registers serialized checkpoints for the memory-mapped path; `invalidate()`
  keeps the K slots repopulatable instead of shrinking them away. The T-loop
  now continues each solver's path incrementally, as in the R reference;
  results are unchanged (warm/cold equality is asserted by new tests). The
  checkpoint temp directory is only created for serialized mode.
- DIRECT / PERMUTATION_DIRECT dummy generation honours the user seed: it
  previously derived from the PERMUTATION-only base seed (always 0, so
  different seeds produced identical dummies) and, when unseeded, drew fresh
  `random_device` entropy on every call — changing the dummies between
  T-steps. The base seed is now resolved once per `DummyGenerator` and
  reused, making unseeded runs internally consistent and seeded runs
  reproducible across instances.
- Documentation: the intentional `Phi_prime` [0, 1] clamp (deviation from the
  R reference / paper to avoid negative FDP estimates) now documents its
  conservative/anti-conservative asymmetry; the DA-BT deflation cites
  Definition 1 of the T-Rex+DA paper (empty-group penalty Ψ = 1/2) and the
  deliberate R-parity choice of terminal-Φ deltas vs. the paper's per-step
  penalties (Eq. (10)).
- Tests: warm-start engagement and warm/cold equality (in-memory, DIRECT, and
  an end-to-end memory-mapped vs. in-memory run that asserts the T-loop
  actually continues), DIRECT seed handling; fixed a checkpoint filename race
  between the typed serialization tests under parallel ctest.
- `tloop_stagnation_stop` is now tri-state (`std::optional<bool>`, Python:
  `bool | None`). The default (auto) resolves by solver family: DISABLED for
  the equiangular LARS-path solvers (TLARS/TLASSO/TENET/TENET_AUG), which
  terminate properly on their own — matching the R reference and the T-Rex
  paper, neither of which has a stagnation stop — and ENABLED for greedy
  solvers, which otherwise run into a noise trap and select until the
  maximum number of dummies is reached. An explicit true/false overrides the
  family default.
- Cleanups from the review follow-up: the voting-grid sentinel uses machine
  epsilon instead of the solver tolerance (a large user tol with K >= 100
  could make the grid non-monotone; also matches the R reference grid); the
  AR1 window fallback for an exactly-zero correlation estimate is the
  formula's limit (±1 window) instead of the full window; the
  equi-correlation estimator excludes near-constant columns from the pair
  count; warm-start solvers are invalidated before the dummy matrices they
  view are mutated (STANDARD/HCONCAT), with the buffer-reallocation contract
  documented on the DummyGenerator; dead validation branches removed; the
  Yule-Walker vs. `arima(..., method="ML")` auto-estimation difference for
  the AR1 `cor_coef` is documented on the parameter and covered by a
  known-rho recovery test.

## TRexSelector 2.0.0

### 2026-06-11

#### Architecture

- Complete architectural rework: C++ core library extracted into a standalone
  CMake project (`TRexSelector`) with full `find_package()` support; R and
  Python bindings are now thin language layers on top.
- C++20 standard adopted throughout; Eigen3, Boost, Cereal, and OpenMP are
  managed via vcpkg or system package managers.
- Simulation and demo executables separated into a standalone sibling project
  (`TRexSelector_Simulations`); the core library no longer depends on Matplot++.

#### Selector variants

- `TRexSelector` — base FDR-controlled selector for linear regression.
- `TRexDASelector` — dependency-aware variant for correlated predictors.
- `TRexGVSSelector` — grouped variable selection.
- `TRexScreeningSelector` — fast FDR-controlled variable pre-screening for
  large *p*.
- `TRexBiobankScreeningSelector` — high-throughput genomic / biobank screening.

#### T-algorithm solvers

- LARS-based family: T-LARS, T-LASSO, T-ENET, T-Stagewise, T-Stepwise, TIENET.
- OMP-based family: T-OMP, T-GP, T-ACGP, T-MP, T-NCGMP, T-OLS.
- AFS-based family: T-AFS.
- All solvers support Cereal polymorphic serialization for warm-start workflows.

#### CMake / build system

- Four exported CMake targets:
  `TRexSelector::trex_utils`, `TRexSelector::trex_ml_methods`,
  `TRexSelector::trex_tsolvers`, `TRexSelector::trex_selector_methods`.
- `cmake --install` produces a fully relocatable install tree consumed via
  `find_package(TRexSelector REQUIRED)`.
- `TREX_PORTABLE_BUILD` option disables `-march=native` for portable wheel
  builds (PyPI / CI).
- `TREX_USE_BLAS_LAPACK` option controls BLAS/LAPACK backend (default ON;
  disabled automatically for Python wheel builds).
- CMakePresets.json provides `debug`, `release`, `release-libs`, and
  `vcpkg-release` presets.

#### R package (version 2.0.0)

- New R6 object-oriented API replaces the previous functional `trex()` API.
- C++ backend is compiled via Rcpp/RcppEigen; no longer depends on the `tlars`
  R package.
- R6 classes: `TRexSelector`, `TRexDASelector`, `TRexGVSSelector`,
  `TRexScreeningSelector`, `TRexBiobankScreeningSelector`.
- Helper functions: `compute_fdp()`, `compute_tpp()`.
- Full testthat test suite and a comprehensive vignette added.
- `R CMD check --as-cran` passes with 0 errors, 0 warnings, 1 NOTE (private
  GitHub URL, expected for new submission).
- Compiler flags set to `-O3 -funroll-loops` in Makevars.in / Makevars.win for
  improved runtime performance.

#### Python package (version 0.1.2)

- pybind11 bindings built via scikit-build-core.
- API mirrors the R6 interface; indices are 0-based.
- Classes: `TRexSelector`, `TRexDASelector`, `TRexGVSSelector`,
  `TRexScreeningSelector`, `TRexBiobankScreeningSelector`.
- 355 pytest tests covering all selectors, all solvers, ML utilities, and
  memory-boundary correctness.
- GIL released during long-running C++ computations.
- `pyproject.toml` contains full PyPI metadata (authors, classifiers, URLs).

#### Documentation

- Root `README.md` updated with C++ developer installation guide, full
  `find_package` consume instructions, and published Elsevier paper references.
- R `README.Rmd` and vignette updated with published journal references
  (Signal Processing, 2025).
- `inst/CITATION` lists both published papers.

---

## TRexSelector 0.1.0

### 2026-05-11

First installable release of the C++ TRexSelector library.

#### Package infrastructure

- `cmake --install` produces a fully relocatable install tree:
  `include/`, `lib/`, `lib/cmake/TRexSelector/`.
- Downstream projects can use `find_package(TRexSelector REQUIRED)` and link
  against `TRexSelector::trex_selector_methods`.

#### CMake changes

- All four library targets declare `FILE_SET HEADERS` for correct header
  installation and IDE indexing.
- Single authoritative `install(EXPORT TRexSelectorTargets)` in root
  `CMakeLists.txt`.
- `CMakePresets.json`: added `release-libs`, `vcpkg-release` presets.

#### CPack packaging

- `cmake --build <dir> --target package` produces binary archives.
- `cmake --build <dir> --target package_source` produces source archives.

#### vcpkg dependency manifest

- `vcpkg.json` declares eigen3, boost-headers, boost-iostreams, cereal, and
  libomp (macOS only) as managed dependencies.

#### Licensing

- GNU General Public License v3 added.
