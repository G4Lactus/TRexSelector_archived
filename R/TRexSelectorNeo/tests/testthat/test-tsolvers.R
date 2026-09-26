# =============================================================================
# test-tsolvers.R
# =============================================================================
#
# Unit tests for the R6 T-Solver wrappers (TSolver_Base and all subclasses).
# Mirrors the structure of the C++ test bed in
#   tests/test_tsolvers/test_tsolvers_polymorphic.cpp
#   tests/test_tsolvers/test_tsolvers_serialization.cpp
#
# Test groups
#   1. Input validation
#   2. Polymorphic execution test (all 22 factory entries)
#   3. Getter return-type and contract checks (TLARS_Solver reference)
#   4. Solver-specific getter checks (lambda, num_removals, cycling_ratio)
#   5. Warm-start / serialization (all 22 factory entries)
# =============================================================================

# ---------------------------------------------------------------------------
# Shared test fixture
# ---------------------------------------------------------------------------

n <- 100L
p <- 30L
d <- 30L
set.seed(42)
X <- matrix(rnorm(n * p), n, p)
D <- matrix(rnorm(n * d), n, d)
y <- X[, 5] * 4.0 + X[, 15] * (-2.5) + X[, 25] * 3.1 + rnorm(n)

# ---------------------------------------------------------------------------
# Solver factory (22 entries: 21 solver types, TNCGMP counted for both variants)
#
# In-place contract: the R6 solvers map X and D without copying.
# A factory call with the defaults builds a WRITER that preprocesses X and D in
# place and claims them until $restore() (or GC); a test that builds several
# writers on the fixture restores each before the next. With (FALSE, FALSE)
# it builds a READER on pre-scaled data (see the serialization group).
# ---------------------------------------------------------------------------

# 1-based contiguous group ids for the IEN solvers (10 groups of 3 predictors)
.tienet_groups <- rep(seq_len(10L), each = 3L)

# Group-mean Tikhonov matrix W = sum_m 1_m 1_m^T / p_m on the same 10 groups
# (dense; the sparse-K solvers accept a dense matrix or a dgCMatrix).
.group_mean_k_30 <- local({
  K <- matrix(0, 30, 30)
  for (g in seq_len(10L)) {
    idx <- ((g - 1L) * 3L + 1L):(g * 3L)
    K[idx, idx] <- 1 / 3
  }
  K
})

.factories <- list(
  TLARS      = function(nrm = TRUE, ic = TRUE) TLARS_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TLASSO     = function(nrm = TRUE, ic = TRUE) TLASSO_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TSTEPWISE  = function(nrm = TRUE, ic = TRUE) TSTEPWISE_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TENET      = function(nrm = TRUE, ic = TRUE) TENET_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TENET_AUG  = function(nrm = TRUE, ic = TRUE) TENETAug_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TIENET_AUG = function(nrm = TRUE, ic = TRUE) TIENETAug_Solver$new(X, D, y, groups = .tienet_groups,
                                              options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TSTAGEWISE = function(nrm = TRUE, ic = TRUE) TSTAGEWISE_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TOMP       = function(nrm = TRUE, ic = TRUE) TOMP_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TGP        = function(nrm = TRUE, ic = TRUE) TGP_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TACGP      = function(nrm = TRUE, ic = TRUE) TACGP_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TMP        = function(nrm = TRUE, ic = TRUE) TMP_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TOOLS      = function(nrm = TRUE, ic = TRUE) TOOLS_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TNCGMP_ls  = function(nrm = TRUE, ic = TRUE) TNCGMP_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic, ncgmp_variant = "line_search")),
  TNCGMP_fc  = function(nrm = TRUE, ic = TRUE) TNCGMP_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic, ncgmp_variant = "fully_corrective")),
  TAFS       = function(nrm = TRUE, ic = TRUE) TAFS_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic, rho_afs = 1.0)),
  # Native IEN, the sparse-K solvers, the CD family.
  TIENET     = function(nrm = TRUE, ic = TRUE) TIENET_Solver$new(X, D, y, groups = .tienet_groups,
                                            options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TIKNET     = function(nrm = TRUE, ic = TRUE) TIKNET_Solver$new(X, D, y, tikhonov_K = .group_mean_k_30,
                                            options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TIKOMP     = function(nrm = TRUE, ic = TRUE) TIKOMP_Solver$new(X, D, y, tikhonov_K = .group_mean_k_30,
                                            options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TCDLASSO   = function(nrm = TRUE, ic = TRUE) TCDLASSO_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic)),
  TCDENET    = function(nrm = TRUE, ic = TRUE) TCDENET_Solver$new(X, D, y, options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TCDIENET   = function(nrm = TRUE, ic = TRUE) TCDIENET_Solver$new(X, D, y, groups = .tienet_groups,
                                              options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5)),
  TCDTIKNET  = function(nrm = TRUE, ic = TRUE) TCDTIKNET_Solver$new(X, D, y, tikhonov_K = .group_mean_k_30,
                                               options = tsolver_options(normalize = nrm, intercept = ic, lambda2 = 0.5,
                                                                         fold_dummies = TRUE))
)


# =============================================================================
# Group 1: Input validation
# =============================================================================

test_that("TSolver_Base rejects non-matrix X", {
  expect_error(TLARS_Solver$new("not_a_matrix", D, y))
})

test_that("TSolver_Base rejects y length mismatch", {
  expect_error(TLARS_Solver$new(X, D, y[seq_len(50)]))
})

test_that("TSolver_Base rejects D row count mismatch", {
  expect_error(TLARS_Solver$new(X, D[seq_len(50), ], y))
})

test_that("TENET_Solver rejects negative lambda2", {
  expect_error(TENET_Solver$new(X, D, y, options = tsolver_options(lambda2 = -0.1)))
})

test_that("TNCGMP_Solver rejects unknown variant string", {
  expect_error(TNCGMP_Solver$new(X, D, y, options = tsolver_options(ncgmp_variant = "bad_variant")))
})


# =============================================================================
# Group 2: Polymorphic execution test (all 22 factory entries)
# Mirrors test_tsolvers_polymorphic.cpp :: ExecutionSmokeTest
# =============================================================================

test_that("All solver types execute without error and produce a valid RSS path", {
  for (nm in names(.factories)) {
    s <- .factories[[nm]]()
    expect_no_error(s$execute_step(T_stop = 5))

    nsteps <- s$get_num_steps()
    expect_type(nsteps, "integer")
    expect_gt(nsteps, 0L, label = nm)

    rss <- s$get_rss()
    expect_length(rss, nsteps)

    # RSS must decrease from first to last step
    expect_gt(rss[[1]], rss[[nsteps]], label = nm)

    s$restore()   # hand X and D back before the next writer
  }
})


# =============================================================================
# Group 3: Getter return-type and contract checks (TLARS_Solver as reference)
# =============================================================================

local({
  s <- TLARS_Solver$new(X, D, y)
  s$execute_step(T_stop = 8)
  nsteps <- s$get_num_steps()

  test_that("get_num_steps returns a positive integer", {
    expect_type(nsteps, "integer")
    expect_gte(nsteps, 1L)
  })

  test_that("get_beta returns a numeric vector of length p", {
    b <- s$get_beta()
    expect_type(b, "double")
    expect_length(b, p)
  })

  test_that("get_beta(step) returns a numeric vector of length p", {
    b <- s$get_beta(step = 1L)
    expect_type(b, "double")
    expect_length(b, p)
  })

  test_that("get_rss returns a positive, strictly decreasing numeric vector", {
    rss <- s$get_rss()
    expect_type(rss, "double")
    expect_length(rss, nsteps)
    expect_true(all(rss > 0), label = "all RSS > 0")
    expect_true(all(diff(rss) <= 1e-12), label = "RSS is non-increasing")
  })

  test_that("get_r2 returns a numeric vector in [0,1] that is non-decreasing", {
    r2 <- s$get_r2()
    expect_type(r2, "double")
    expect_length(r2, nsteps)
    expect_true(all(r2 >= 0 & r2 <= 1), label = "R2 in [0,1]")
    expect_true(all(diff(r2) >= -1e-12), label = "R2 non-decreasing")
  })

  test_that("get_cp returns a numeric vector of length nsteps", {
    cp <- s$get_cp()
    expect_type(cp, "double")
    expect_length(cp, nsteps)
  })

  test_that("get_dof returns a non-negative numeric vector of length nsteps", {
    dof <- s$get_dof()
    expect_type(dof, "double")
    expect_length(dof, nsteps)
    expect_true(all(dof >= 0), label = "all DoF >= 0")
  })

  test_that("get_actives(r_index=TRUE) returns 1-based indices within [1, p]", {
    a <- s$get_actives(r_index = TRUE)
    expect_type(a, "double")
    expect_true(all(a >= 1L & a <= p), label = "actives 1-based in [1,p]")
  })

  test_that("get_actives(r_index=FALSE) returns 0-based indices within [0, p-1]", {
    a <- s$get_actives(r_index = FALSE)
    expect_type(a, "double")
    expect_true(all(a >= 0L & a <= p - 1L), label = "actives 0-based in [0,p-1]")
  })

  test_that("get_inactives(r_index=TRUE) returns 1-based indices within [1, p]", {
    i <- s$get_inactives(r_index = TRUE)
    expect_type(i, "double")
    expect_true(all(i >= 1L & i <= p), label = "inactives 1-based in [1,p]")
  })

  test_that("get_actives and get_inactives partition [1, p]", {
    a <- s$get_actives(r_index = TRUE)
    i <- s$get_inactives(r_index = TRUE)
    expect_setequal(c(a, i), seq_len(p))
  })

  test_that("get_intercept returns a single numeric scalar", {
    ic <- s$get_intercept()
    expect_type(ic, "double")
    expect_length(ic, 1L)
  })

  test_that("get_actions(r_index=TRUE) returns list of 1-based integer-vector elements", {
    acts <- s$get_actions(r_index = TRUE)
    expect_type(acts, "list")
    expect_length(acts, nsteps)
    for (step_acts in acts) {
      expect_type(step_acts, "integer")
      # All action indices must be >= 1 in absolute value (1-based)
      expect_true(all(abs(step_acts) >= 1L), label = "all action |index| >= 1")
    }
  })

  test_that("get_actions(r_index=FALSE) returns list of 0-based integer-vector elements", {
    acts <- s$get_actions(r_index = FALSE)
    expect_type(acts, "list")
    expect_length(acts, nsteps)
    for (step_acts in acts) {
      expect_type(step_acts, "integer")
      # All action indices must be >= 0 in absolute value (0-based)
      expect_true(all(abs(step_acts) >= 0L), label = "all action |index| >= 0")
    }
  })

  test_that("get_actions(r_index=TRUE) is the recorded lars code, r_index=FALSE its 0-based magnitude", {
    a1 <- s$get_actions(r_index = TRUE)
    a0 <- s$get_actions(r_index = FALSE)
    expect_identical(lapply(a1, function(a) abs(a) - 1L), a0)
    # The first action of T-LARS is the addition of one variable with a
    # 1-based code in [1, p + d]
    expect_length(a1[[1]], 1L)
    expect_true(a1[[1]] >= 1L && a1[[1]] <= p + d)
  })

  s$restore()
})


# =============================================================================
# Group 4: Solver-specific getter contracts
# =============================================================================

test_that("TLARS_Solver: get_lambda returns a non-increasing numeric vector of length nsteps", {
  s <- TLARS_Solver$new(X, D, y)
  s$execute_step(T_stop = 8)
  lam <- s$get_lambda()
  expect_type(lam, "double")
  expect_length(lam, s$get_num_steps())
  expect_true(all(diff(lam) <= 1e-12), label = "lambda non-increasing")
  s$restore()
})

test_that("TLASSO_Solver: get_lambda, get_num_removals, get_cycling_ratio all correct", {
  s <- TLASSO_Solver$new(X, D, y)
  s$execute_step(T_stop = 8)
  expect_type(s$get_lambda(), "double")
  nr <- s$get_num_removals()
  expect_type(nr, "integer")
  expect_gte(nr, 0L)
  cr <- s$get_cycling_ratio()
  expect_type(cr, "double")
  expect_gte(cr, 0.0)
  expect_lte(cr, 1.0)
  s$restore()
})

test_that("TSTEPWISE_Solver: get_lambda returns a numeric vector", {
  s <- TSTEPWISE_Solver$new(X, D, y)
  s$execute_step(T_stop = 8)
  expect_type(s$get_lambda(), "double")
  s$restore()
})

test_that("TENET_Solver: get_lambda, get_num_removals, get_cycling_ratio all correct", {
  s <- TENET_Solver$new(X, D, y, options = tsolver_options(lambda2 = 0.5))
  s$execute_step(T_stop = 8)
  expect_type(s$get_lambda(), "double")
  nr <- s$get_num_removals()
  expect_type(nr, "integer")
  expect_gte(nr, 0L)
  cr <- s$get_cycling_ratio()
  expect_type(cr, "double")
  expect_gte(cr, 0.0)
  expect_lte(cr, 1.0)
  s$restore()
})

test_that("TSTAGEWISE_Solver: get_num_removals non-negative, get_cycling_ratio in [0,1]", {
  s <- TSTAGEWISE_Solver$new(X, D, y)
  s$execute_step(T_stop = 8)
  nr <- s$get_num_removals()
  expect_type(nr, "integer")
  expect_gte(nr, 0L)
  cr <- s$get_cycling_ratio()
  expect_type(cr, "double")
  expect_gte(cr, 0.0)
  expect_lte(cr, 1.0)
  s$restore()
})

test_that("TIENETAug_Solver: get_lambda2, get_groups (1-based), get_num_groups", {
  s <- TIENETAug_Solver$new(X, D, y, groups = .tienet_groups,
                            options = tsolver_options(lambda2 = 0.5))
  s$execute_step(T_stop = 8)
  expect_equal(s$get_lambda2(), 0.5)
  expect_identical(s$get_groups(), as.integer(.tienet_groups))
  expect_identical(s$get_num_groups(), 10L)
  s$restore()
})


# =============================================================================
# Group 5: Warm-start / serialization — all 22 factory entries
# Mirrors test_tsolvers_polymorphic.cpp :: SerializationEquivalence and the
# C++ buffer-guard round trip: X and D are pre-scaled ONCE (the statistics
# live in the scalers), and the reference, checkpoint and resumed solvers are
# READERS on the same matrices, so the three paths see bit-identical data.
# The former per-solver writers would each re-normalize the already
# normalized X (the T-Stagewise "Paths differ" finding).
# =============================================================================

test_that("All solver types resume from checkpoint identically to a continuous run", {
  T_partial <- 4L
  T_final   <- 8L

  sx <- MLLpNormScaler$new(); sx$fit_transform_inplace(X)
  sd <- MLLpNormScaler$new(); sd$fit_transform_inplace(D)
  on.exit({ sx$inverse_transform_inplace(X); sd$inverse_transform_inplace(D) }, add = TRUE)

  for (nm in names(.factories)) {
    # --- Reference: uninterrupted run to T_final ---
    s_ref <- .factories[[nm]](FALSE, FALSE)
    s_ref$execute_step(T_stop = T_final)

    # --- Partial run, save, reload into a reader, continue ---
    ck <- tempfile(fileext = ".bin")
    on.exit(unlink(ck), add = TRUE)

    s1 <- .factories[[nm]](FALSE, FALSE)
    s1$execute_step(T_stop = T_partial)
    s1$save(ck)
    expect_true(file.exists(ck), label = paste(nm, "checkpoint created"))

    s2 <- .factories[[nm]](FALSE, FALSE)
    s2$load(ck)
    expect_equal(s2$get_num_steps(), s1$get_num_steps(),
                 label = paste(nm, "num_steps after load"))
    s2$execute_step(T_stop = T_final)

    # Number of steps must match
    expect_equal(
      s_ref$get_num_steps(), s2$get_num_steps(),
      label = paste(nm, "num_steps after resume")
    )

    # Action path must match
    expect_identical(
      s_ref$get_actions(), s2$get_actions(),
      label = paste(nm, "action path after resume")
    )

    # Final beta: the LARS / OMP / AFS families are bit-identical; the CD
    # family re-enters its gram state after load and matches to rounding
    # (the C++ test uses the same tolerance)
    expect_equal(
      s_ref$get_beta(), s2$get_beta(),
      tolerance = 1e-9,
      label = paste(nm, "beta after resume")
    )
  }
})


# =============================================================================
# Group 6: family getters
# =============================================================================

test_that("TIENET_Solver / TCDIENET_Solver: lambda2, groups (1-based), num_groups", {
  for (cls in list(TIENET_Solver, TCDIENET_Solver)) {
    s <- cls$new(X, D, y, groups = .tienet_groups,
                 options = tsolver_options(lambda2 = 0.5))
    s$execute_step(T_stop = 6)
    expect_equal(s$get_lambda2(), 0.5)
    expect_identical(s$get_groups(), as.integer(.tienet_groups))
    expect_identical(s$get_num_groups(), 10L)
    expect_gt(s$get_num_steps(), 0L)
    s$restore()
  }
})

test_that("CD family: certified localization getters and record vectors", {
  s <- TCDENET_Solver$new(X, D, y, options = tsolver_options(lambda2 = 0.5))
  s$execute_step(T_stop = 6)
  expect_true(s$is_certified_localization())
  expect_type(s$is_certified_degraded(), "logical")
  expect_equal(s$get_lambda2(), 0.5)
  expect_length(s$get_lambda(), s$get_num_steps())   # one lambda_1 per recorded step
  expect_true(all(s$get_record_gaps() >= 0))
  expect_true(all(s$get_record_kkt_violations() >= 0))
  s$restore()

  legacy <- TCDLASSO_Solver$new(X, D, y,
                                options = tsolver_options(cd_certified_localization = FALSE))
  expect_false(legacy$is_certified_localization())
  expect_equal(legacy$get_lambda2(), 0)
  legacy$restore()
})

test_that("Sparse-K solvers accept a dense K or a dgCMatrix and report lambda2", {
  skip_if_not_installed("Matrix")
  K_sparse <- Matrix::Matrix(.group_mean_k_30, sparse = TRUE)
  # Readers on pre-scaled data: dense and sparse K then see identical X, D
  sx <- MLLpNormScaler$new(); sx$fit_transform_inplace(X)
  sd <- MLLpNormScaler$new(); sd$fit_transform_inplace(D)
  on.exit({ sx$inverse_transform_inplace(X); sd$inverse_transform_inplace(D) }, add = TRUE)
  reader <- tsolver_options(lambda2 = 0.5, fold_dummies = TRUE,
                            normalize = FALSE, intercept = FALSE)
  for (cls in list(TIKNET_Solver, TIKOMP_Solver, TCDTIKNET_Solver)) {
    dense  <- cls$new(X, D, y, tikhonov_K = .group_mean_k_30, options = reader)
    sparse <- cls$new(X, D, y, tikhonov_K = K_sparse, options = reader)
    dense$execute_step(T_stop = 6)
    sparse$execute_step(T_stop = 6)
    expect_equal(dense$get_lambda2(), 0.5)
    expect_identical(dense$get_beta(), sparse$get_beta())
  }
  tik <- TCDTIKNET_Solver$new(X, D, y, tikhonov_K = .group_mean_k_30, options = reader)
  expect_true(tik$is_folded_dummy_coupling())
  expect_error(TIKNET_Solver$new(X, D, y, tikhonov_K = -diag(30), options = reader), "PSD|positive")
})


# =============================================================================
# The base methods shared with the Python binding, and the family getters.
# =============================================================================

test_that("the T-solver base methods added in the surface catch-up work", {
  s <- TLARS_Solver$new(X, D, y, options = tsolver_options(tie_seed = 3))
  on.exit(s$restore(), add = TRUE)
  expect_true(s$is_connected())
  expect_identical(s$solver_type_to_string(), "TLARS")
  expect_identical(s$get_dummy_start_index(), p + 1L)
  s$set_tolerance(1e-12)$set_tie_seed(5)$set_exchangeable_tie(0)   # chainable, no-ops for LARS
  s$execute_step(T_stop = 6)
  nsteps <- s$get_num_steps()
  expect_length(s$get_residuals(), n)
  bp <- s$get_beta_path()                     # step 0 (all zero) included, as in Python
  expect_equal(dim(bp), c(p + d, nsteps + 1L))
  expect_true(all(bp[, 1] == 0))
  expect_equal(bp[seq_len(p), nsteps + 1L], s$get_beta(), tolerance = 1e-12)
  expect_length(s$get_cp_with_sigma(1.0), nsteps)
  expect_identical(s$get_num_actives(),
                   length(s$get_active_predictor_indices()) + length(s$get_active_dummy_indices()))
  expect_setequal(s$get_active_predictor_indices(), s$get_actives())
  expect_true(all(s$get_active_dummy_indices() >= 1L & s$get_active_dummy_indices() <= d))
  expect_type(s$get_dropped_indices(), "integer")
  expect_type(s$get_warnings(), "character")
  s$clear_warnings()
  expect_length(s$get_warnings(), 0L)
})

test_that("the family getters: Kahan interval, alpha, Tikhonov K, KKT violation", {
  s <- TLARS_Solver$new(X, D, y)
  k0 <- s$get_kahan_refresh_interval()
  expect_type(k0, "integer")
  s$set_kahan_refresh_interval(7L)
  expect_identical(s$get_kahan_refresh_interval(), 7L)
  s$restore()

  e <- TENET_Solver$new(X, D, y, options = tsolver_options(lambda2 = 0.5))
  expect_true(is.numeric(e$get_alpha()) && e$get_alpha() > 0 && e$get_alpha() <= 1)
  e$restore()

  skip_if_not_installed("Matrix")
  t <- TIKNET_Solver$new(X, D, y, tikhonov_K = .group_mean_k_30,
                         options = tsolver_options(lambda2 = 0.5))
  K <- t$get_tikhonov_k()
  expect_s4_class(K, "dgCMatrix")
  expect_equal(as.matrix(K), .group_mean_k_30, tolerance = 1e-12)
  t$restore()

  cd <- TCDLASSO_Solver$new(X, D, y)
  cd$execute_step(T_stop = 4)
  expect_true(is.numeric(cd$get_kkt_violation()) && cd$get_kkt_violation() >= 0)
  cd$restore()
})
