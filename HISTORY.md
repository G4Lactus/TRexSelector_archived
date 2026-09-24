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
- The selectors reject non-finite entries and constant columns of X (constant
  columns were silently unselectable before), see `ScalingMode::NONE` below.
- Every T-solver is constructed as `(X, D, y, [groups | K,] options)` with one
  options structure; the positional constructors and the per-solver R / Python
  keyword arguments are gone, see below.
- Floating-point classification is bit-based and the Release build no longer uses
  `-ffast-math`: results are IEEE 754 conformant; runs from fast-math builds are
  not bit-reproducible against these.

### 2026-09-24

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
