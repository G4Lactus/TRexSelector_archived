// ========================================================================================
// test_tsolvers_execution.cpp
// ========================================================================================
/**
 * @file test_tsolvers_execution.cpp
 *
 * @brief Unit tests verifying core execution path integrity (bounds checking, monotonicity)
 *        for TSolver_Base via the concrete proxy TLARS_Solver.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// std includes
#include <map>
#include <set>

// Eigen includes
#include <Eigen/Dense>

// tsolvers includes
#include <tsolvers/ad/linear_model/lars_based/tlars_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tlasso_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tstagewise_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tienet_solver.hpp>
#include <tsolvers/ad/linear_model/lars_based/tienet_aug_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tncgmp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tmp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tomp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tgp_solver.hpp>
#include <tsolvers/ad/linear_model/omp_based/tools_solver.hpp>

// utils includes
#include <utils/datageneration/utils_datagen.hpp>

// ========================================================================================

// Embed into test namespace
namespace trex::test::tsolvers::ad::linear_model {

using namespace trex::tsolvers::ad::linear_model::lars_based;
using namespace trex::utils::datageneration::datagen;

// ========================================================================================

class TSolverExecutionTest : public ::testing::Test {
protected:
    Eigen::MatrixXd X, D;
    Eigen::VectorXd y;

    void SetUp() override {
        // Use datagen utilities to get well-conditioned design matrices
        SyntheticData data(
            200, 50, 50,
            {0, 10, 20},
            {5.0, -3.0, 2.0},
            1.0, 42, -1, -1,
            predictor_policy::Normal(),
            dummygen::Distribution::Normal(),
            noisegen::noise_policy::Normal()
        );

        X = data.getX();
        D = data.getD();
        y = data.getY();
    }
};

// ========================================================================================

/** @brief Test to validate that the solver properly respects early T_stop halting constraints. */
TEST_F(TSolverExecutionTest, ExecuteEarlyStoppingAtTStop) {
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), D.rows(), D.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TLARS_Solver solver(X_map, D_map, y_map);

    // Goal: Stop when 4 dummy variables have entered the active set.
    int T_stop = 4;
    solver.executeStep(T_stop, true);

    // It should have executed some steps...
    EXPECT_GT(solver.getNumSteps(), 0);

    // Check dummy capacity actively constrained to exactly what we requested
    auto active_dummies = solver.getActiveDummyIndices();
    EXPECT_EQ(active_dummies.size(), T_stop);
}


/** @brief Test to ensure the residual sum of squares shrinks monotonically along the inclusion path. */
TEST_F(TSolverExecutionTest, MonotonicRSSDecrease) {
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), D.rows(), D.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TLARS_Solver solver(X_map, D_map, y_map);

    // Execute multiple steps
    solver.executeStep(10, true);

    const auto& rss_path = solver.getRSS();
    EXPECT_FALSE(rss_path.empty());

    // RSS should strictly decrease (or at least remain non-increasing) over successful
    // variable inclusions
    for (std::size_t i = 1; i < rss_path.size(); ++i) {
        // We use EXPECT_LE due to minor precision limits on edge updates
        EXPECT_LE(rss_path[i], rss_path[i-1] + 1e-12);
    }
}

/** @brief MP and GP must re-select atoms (their defining difference from
 *         OMP): on a full path, some atom is picked in more than one step,
 *         while the active set keeps each index exactly once. */
TEST_F(TSolverExecutionTest, MPAndGPReSelectAtoms) {
    using trex::tsolvers::ad::linear_model::omp_based::TMP_Solver;
    using trex::tsolvers::ad::linear_model::omp_based::TGP_Solver;

    auto max_selection_count = [this](auto make_solver) {
        Eigen::MatrixXd Xc = X, Dc = D;
        Eigen::VectorXd yc = y;
        Eigen::Map<Eigen::MatrixXd> X_map(Xc.data(), Xc.rows(), Xc.cols());
        Eigen::Map<Eigen::MatrixXd> D_map(Dc.data(), Dc.rows(), Dc.cols());
        Eigen::Map<Eigen::VectorXd> y_map(yc.data(), yc.size());
        auto solver = make_solver(X_map, D_map, y_map);
        solver.executeStep(0, false); // full path

        std::map<int, int> counts;
        for (const auto& step : solver.getActions()) {
            for (int a : step) {
                if (a > 0) { counts[a]++; } // Additions are positive 1-based indices
            }
        }
        int max_count = 0;
        for (const auto& [idx, c] : counts) { max_count = std::max(max_count, c); }

        // Active set holds each index exactly once despite repeated picks
        std::set<std::size_t> unique_actives(solver.getActives().begin(),
                                             solver.getActives().end());
        EXPECT_EQ(unique_actives.size(), solver.getActives().size());
        return max_count;
    };

    int mp_max = max_selection_count([](auto& Xm, auto& Dm, auto& ym) {
        return TMP_Solver(Xm, Dm, ym);
    });
    EXPECT_GE(mp_max, 2) << "MP never re-selected an atom on a full path";

    int gp_max = max_selection_count([](auto& Xm, auto& Dm, auto& ym) {
        return TGP_Solver(Xm, Dm, ym);
    });
    EXPECT_GE(gp_max, 2) << "GP never re-selected an atom on a full path";
}


/** @brief Pin the documented equivalence (Locatello et al. 2017, Alg. 4): for
 *         the least-squares objective (L = 1), NCGMP-LineSearch is exactly
 *         T-MP and NCGMP-FullyCorrective is exactly T-OMP. */
TEST_F(TSolverExecutionTest, NCGMPVariantsMatchMPAndOMP) {
    using trex::tsolvers::ad::linear_model::omp_based::TMP_Solver;
    using trex::tsolvers::ad::linear_model::omp_based::TOMP_Solver;
    using trex::tsolvers::ad::linear_model::omp_based::TNCGMP_Solver;
    using trex::tsolvers::ad::linear_model::omp_based::NCGMPVariant;

    auto fit = [this](auto make_solver) {
        Eigen::MatrixXd Xc = X, Dc = D;
        Eigen::VectorXd yc = y;
        Eigen::Map<Eigen::MatrixXd> X_map(Xc.data(), Xc.rows(), Xc.cols());
        Eigen::Map<Eigen::MatrixXd> D_map(Dc.data(), Dc.rows(), Dc.cols());
        Eigen::Map<Eigen::VectorXd> y_map(yc.data(), yc.size());
        auto solver = make_solver(X_map, D_map, y_map);
        solver.executeStep(4, true);
        return solver.getBetaPathSparse().dense();
    };

    Eigen::MatrixXd mp_path = fit([](auto& Xm, auto& Dm, auto& ym) {
        return TMP_Solver(Xm, Dm, ym);
    });
    Eigen::MatrixXd ncg_ls_path = fit([](auto& Xm, auto& Dm, auto& ym) {
        return TNCGMP_Solver(Xm, Dm, ym);
    });
    ASSERT_EQ(mp_path.cols(), ncg_ls_path.cols());
    EXPECT_TRUE(mp_path.isApprox(ncg_ls_path, 1e-12))
        << "max abs diff: " << (mp_path - ncg_ls_path).cwiseAbs().maxCoeff();

    Eigen::MatrixXd omp_path = fit([](auto& Xm, auto& Dm, auto& ym) {
        return TOMP_Solver(Xm, Dm, ym);
    });
    Eigen::MatrixXd ncg_fc_path = fit([](auto& Xm, auto& Dm, auto& ym) {
        return TNCGMP_Solver(Xm, Dm, ym, trex::tsolvers::ad::SolverOptions{.ncgmp_variant = NCGMPVariant::FullyCorrective});
    });
    ASSERT_EQ(omp_path.cols(), ncg_fc_path.cols());
    EXPECT_TRUE(omp_path.isApprox(ncg_fc_path, 1e-12))
        << "max abs diff: " << (omp_path - ncg_fc_path).cwiseAbs().maxCoeff();
}


/** @brief Coverage for the Fully Corrective (OMP) variant of TNCGMP, which the
 *         typed suite does not reach (its factory hardcodes LineSearch). Also
 *         asserts the actions/diagnostics alignment invariant that collinear
 *         rejections must not break. */
TEST_F(TSolverExecutionTest, TNCGMPFullyCorrectiveExecutesAndStaysAligned) {
    using trex::tsolvers::ad::linear_model::omp_based::TNCGMP_Solver;
    using trex::tsolvers::ad::linear_model::omp_based::NCGMPVariant;

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), D.rows(), D.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TNCGMP_Solver solver(X_map, D_map, y_map, trex::tsolvers::ad::SolverOptions{.ncgmp_variant = NCGMPVariant::FullyCorrective});
    solver.executeStep(4, true);

    EXPECT_GT(solver.getNumSteps(), 0u);
    EXPECT_EQ(solver.getActiveDummyIndices().size(), 4u);

    // One actions entry per recorded step, aligned with the diagnostics arrays
    EXPECT_EQ(solver.getActions().size(), solver.getNumSteps());
    EXPECT_EQ(solver.getActions().size(), solver.getRSS().size() - 1);
    EXPECT_EQ(solver.getRSS().size(), solver.getDoF().size());

    // Fully corrective refits must decrease RSS monotonically
    const auto& rss = solver.getRSS();
    for (std::size_t i = 1; i < rss.size(); ++i) {
        EXPECT_LE(rss[i], rss[i - 1] + 1e-12);
    }
}


/** @brief Regression test: executeStep() with defaults (T_stop=0, early_stop=true)
 *         must compute the full solution path as documented — historically it was
 *         a silent no-op because the loop guard treated T_stop=0 as "zero dummies
 *         allowed". T_stop=0 now means "no dummy limit". */
TEST_F(TSolverExecutionTest, DefaultExecuteStepComputesFullPath) {
    auto X2 = X; auto D2 = D; auto y2 = y;

    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::MatrixXd> D_map(D.data(), D.rows(), D.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TLARS_Solver by_default(X_map, D_map, y_map);
    by_default.executeStep(); // T_stop = 0, early_stop = true

    Eigen::Map<Eigen::MatrixXd> X_map2(X2.data(), X2.rows(), X2.cols());
    Eigen::Map<Eigen::MatrixXd> D_map2(D2.data(), D2.rows(), D2.cols());
    Eigen::Map<Eigen::VectorXd> y_map2(y2.data(), y2.size());

    TLARS_Solver explicit_full(X_map2, D_map2, y_map2);
    explicit_full.executeStep(0, false); // documented full-path form

    EXPECT_GT(by_default.getNumSteps(), 0u);
    EXPECT_EQ(by_default.getNumSteps(), explicit_full.getNumSteps());
    EXPECT_EQ(by_default.getActions(), explicit_full.getActions());
}


/** @brief Pin T-LASSO cycling behavior: zero-crossing drops must keep the
 *         lambda sequence monotone, zero the dropped coefficient at the drop
 *         step, and keep the actions/diagnostics arrays aligned. */
TEST_F(TSolverExecutionTest, LassoZeroCrossingDropsBehaveConsistently) {
    std::srand(3);

    // Strongly correlated predictors with opposing effects provoke
    // coefficient sign changes and hence LASSO drops.
    Eigen::MatrixXd Xc = Eigen::MatrixXd::Random(200, 20);
    Xc.col(1) = 0.95 * Xc.col(0) + 0.05 * Xc.col(1);
    Xc.col(2) = 0.90 * Xc.col(0) + 0.10 * Xc.col(2);
    Eigen::MatrixXd Dc = Eigen::MatrixXd::Random(200, 20);
    Eigen::VectorXd yc = 5.0 * Xc.col(0) - 4.0 * Xc.col(1) + 2.0 * Xc.col(2)
                         + 0.1 * Eigen::VectorXd::Random(200);

    Eigen::Map<Eigen::MatrixXd> X_map(Xc.data(), Xc.rows(), Xc.cols());
    Eigen::Map<Eigen::MatrixXd> D_map(Dc.data(), Dc.rows(), Dc.cols());
    Eigen::Map<Eigen::VectorXd> y_map(yc.data(), yc.size());

    TLASSO_Solver solver(X_map, D_map, y_map);
    solver.executeStep(0, false); // full path

    // The design must actually trigger drops, otherwise this test is vacuous.
    ASSERT_GT(solver.getNumRemovals(), 0u);

    // Max correlation must be non-increasing along the path
    const auto& lambda = solver.getLambda();
    for (std::size_t k = 1; k < lambda.size(); ++k) {
        EXPECT_LE(lambda[k], lambda[k - 1] + 1e-8) << "lambda jumped upward at step " << k;
    }

    // One actions entry per recorded step
    const auto& actions = solver.getActions();
    EXPECT_EQ(actions.size(), solver.getRSS().size() - 1);

    // Every drop must zero the coefficient at its step (beta column k+1).
    // Actions are 1-based signed indices; decode drops with -a - 1.
    Eigen::MatrixXd path = solver.getBetaPathSparse().dense();
    for (std::size_t k = 0; k < actions.size(); ++k) {
        for (int a : actions[k]) {
            if (a < 0) {
                EXPECT_EQ(path(static_cast<Eigen::Index>(-a) - 1,
                                static_cast<Eigen::Index>(k + 1)), 0.0)
                    << "dropped variable " << (-a - 1)
                    << " has nonzero beta at step " << k + 1;
            }
        }
    }
}


/** @brief Regression test: NNLS drops in T-Stagewise must refresh the dropped
 *         variable's correlation entry. A stale (inflated) entry re-enters
 *         findTiedMaxCorrelations and shows up as an upward jump in the lambda
 *         (max correlation) sequence, which must be non-increasing. */
TEST_F(TSolverExecutionTest, StagewiseNNLSDropKeepsLambdaMonotone) {
    std::srand(7);

    // Strongly correlated predictors with opposing effects provoke
    // direction-sign violations and hence NNLS removals.
    Eigen::MatrixXd Xc = Eigen::MatrixXd::Random(200, 20);
    Xc.col(1) = 0.95 * Xc.col(0) + 0.05 * Xc.col(1);
    Xc.col(2) = 0.90 * Xc.col(0) + 0.10 * Xc.col(2);
    Eigen::MatrixXd Dc = Eigen::MatrixXd::Random(200, 20);
    Eigen::VectorXd yc = 5.0 * Xc.col(0) - 4.0 * Xc.col(1) + 2.0 * Xc.col(2)
                         + 0.1 * Eigen::VectorXd::Random(200);

    Eigen::Map<Eigen::MatrixXd> X_map(Xc.data(), Xc.rows(), Xc.cols());
    Eigen::Map<Eigen::MatrixXd> D_map(Dc.data(), Dc.rows(), Dc.cols());
    Eigen::Map<Eigen::VectorXd> y_map(yc.data(), yc.size());

    TSTAGEWISE_Solver solver(X_map, D_map, y_map);
    solver.executeStep(0, false); // full path

    // The design must actually trigger NNLS removals, otherwise this test is vacuous.
    EXPECT_GT(solver.getNumRemovals(), 0u);

    // Max correlation must be non-increasing along the path.
    const auto& lambda = solver.getLambda();
    ASSERT_FALSE(lambda.empty());
    for (std::size_t k = 1; k < lambda.size(); ++k) {
        EXPECT_LE(lambda[k], lambda[k - 1] + 1e-8) << "lambda jumped upward at step " << k;
    }
}


/** @brief The native pathwise TIENET and the row-augmented TIENETAug solve
 *  the same IEN problem (both standardize the data blocks before the group
 *  penalty enters, and neither re-normalizes the penalty rows), so with
 *  identical inputs the two must walk the identical solution path. */
TEST_F(TSolverExecutionTest, TIENETNativeMatchesAugmentedPath) {
    // 10 contiguous groups of 5; the single dummy layer (L = p = 50) shares
    // the group of its originating variable, per the T-Rex+GVS convention.
    Eigen::VectorXi groups(X.cols());
    for (Eigen::Index j = 0; j < X.cols(); ++j) {
        groups(j) = static_cast<int>(j / 5);
    }

    for (double lambda2 : {0.5, 2.0}) {
        // Native solver normalizes in place -> give it its own copies.
        Eigen::MatrixXd Xn = X, Dn = D;
        Eigen::VectorXd yn = y;
        Eigen::Map<Eigen::MatrixXd> Xn_map(Xn.data(), Xn.rows(), Xn.cols());
        Eigen::Map<Eigen::MatrixXd> Dn_map(Dn.data(), Dn.rows(), Dn.cols());
        Eigen::Map<Eigen::VectorXd> yn_map(yn.data(), yn.size());
        TIENET_Solver native(Xn_map, Dn_map, yn_map, groups, trex::tsolvers::ad::SolverOptions{.lambda2 = lambda2});
        native.executeStep(5, true);

        Eigen::MatrixXd Xa = X, Da = D;
        Eigen::VectorXd ya = y;
        Eigen::Map<Eigen::MatrixXd> Xa_map(Xa.data(), Xa.rows(), Xa.cols());
        Eigen::Map<Eigen::MatrixXd> Da_map(Da.data(), Da.rows(), Da.cols());
        Eigen::Map<Eigen::VectorXd> ya_map(ya.data(), ya.size());
        TIENETAug_Solver aug(Xa_map, Da_map, ya_map, groups, trex::tsolvers::ad::SolverOptions{.lambda2 = lambda2});
        aug.executeStep(5, true);

        // Identical step-by-step action sequences (adds and drops).
        EXPECT_EQ(native.getActions(), aug.getActions())
            << "native and augmented IEN paths diverged at lambda2 = "
            << lambda2;
        EXPECT_EQ(native.getActives(), aug.getActives());
        EXPECT_EQ(native.getActiveDummyIndices().size(),
                  aug.getActiveDummyIndices().size());
    }
}


/** @brief lambda2 == 0 collapses the native TIENET exactly to T-LASSO. */
TEST_F(TSolverExecutionTest, TIENETLambda2ZeroReducesToTLASSO) {
    Eigen::VectorXi groups(X.cols());
    for (Eigen::Index j = 0; j < X.cols(); ++j) {
        groups(j) = static_cast<int>(j / 5);
    }

    Eigen::MatrixXd Xi = X, Di = D;
    Eigen::VectorXd yi = y;
    Eigen::Map<Eigen::MatrixXd> Xi_map(Xi.data(), Xi.rows(), Xi.cols());
    Eigen::Map<Eigen::MatrixXd> Di_map(Di.data(), Di.rows(), Di.cols());
    Eigen::Map<Eigen::VectorXd> yi_map(yi.data(), yi.size());
    TIENET_Solver ien(Xi_map, Di_map, yi_map, groups);
    ien.executeStep(5, true);

    Eigen::MatrixXd Xl = X, Dl = D;
    Eigen::VectorXd yl = y;
    Eigen::Map<Eigen::MatrixXd> Xl_map(Xl.data(), Xl.rows(), Xl.cols());
    Eigen::Map<Eigen::MatrixXd> Dl_map(Dl.data(), Dl.rows(), Dl.cols());
    Eigen::Map<Eigen::VectorXd> yl_map(yl.data(), yl.size());
    TLASSO_Solver lasso(Xl_map, Dl_map, yl_map);
    lasso.executeStep(5, true);

    EXPECT_EQ(ien.getActions(), lasso.getActions());
    auto beta_ien   = ien.getBeta(-1);
    auto beta_lasso = lasso.getBeta(-1);
    EXPECT_TRUE(beta_ien.isApprox(beta_lasso, 1e-10));
}

// ========================================================================================
/**
 * @brief T-OLS with the selection score recomputed from scratch at every step:
 *        h_j = X_A^T x_j and one triangular solve per candidate, the direct
 *        form of ||P_perp x_j||^2 = ||x_j||^2 - ||R^{-T} h_j||^2.
 */
class TOOLSFromScratch : public trex::tsolvers::ad::linear_model::omp_based::TOOLS_Solver {
public:
    using TOOLS_Solver::TOOLS_Solver;

protected:
    std::pair<double, std::vector<std::size_t>> findTiedMaxCorrelations() const override {
        double Cmax = 0.0;
        std::vector<std::size_t> tied;
        if (inactives_.empty()) return {Cmax, tied};

        const auto m = static_cast<Eigen::Index>(actives_.size());
        for (std::size_t j : inactives_) {
            const double corr = correlations_(static_cast<Eigen::Index>(j));
            double score;
            if (m == 0) {
                const double nrm = getColumn(j).norm();
                score = (nrm > eps_) ? std::abs(corr) / nrm : -1.0;
            } else {
                Eigen::VectorXd hj(m);
                for (Eigen::Index k = 0; k < m; ++k) {
                    hj(k) = getColumn(actives_[static_cast<std::size_t>(k)]).dot(getColumn(j));
                }
                const Eigen::VectorXd w = backsolveT(R_, hj);
                const double orth_sq = getColumn(j).squaredNorm() - w.squaredNorm();
                score = (orth_sq < eps_) ? -1.0 : std::abs(corr) / std::sqrt(orth_sq);
            }
            if (score > Cmax + eps_) {
                Cmax = score;
                tied.clear();
                tied.push_back(j);
            } else if (score >= Cmax - eps_ && score > 0.0) {
                tied.push_back(j);
            }
        }
        return {Cmax, tied};
    }
};

/**
 * @brief T-OLS downdates the projected candidate norms by (q^T x_j)^2 per new
 *        active position instead of re-projecting every candidate: the path
 *        (actions and coefficients, bit for bit) must match the from-scratch
 *        score on centered designs run to the n - 1 row bound.
 */
TEST(TOOLSScoring, DowndatedNormsMatchFromScratchProjection) {
    using trex::tsolvers::ad::SolverOptions;
    using trex::tsolvers::ad::linear_model::omp_based::TOOLS_Solver;

    constexpr Eigen::Index n = 60, p = 120, L = 120;
    for (int seed : {1, 2, 3, 4, 5}) {
        SyntheticData data(n, p, L, {3, 17, 41, 77, 101}, {3.0, -2.0, 2.0, -3.0, 2.5},
                           10.0, seed, -1, -1, predictor_policy::Normal(),
                           dummygen::Distribution::Normal(), noisegen::noise_policy::Normal());
        Eigen::MatrixXd Xs = data.getX();
        Eigen::MatrixXd Ds = data.getD();
        Eigen::VectorXd ys = data.getY();

        // Caller-side centering and unit-L2 scaling; the solvers are readers
        auto center_scale = [](Eigen::MatrixXd& M) {
            M.rowwise() -= M.colwise().mean();
            M.array().rowwise() /= M.colwise().norm().array();
        };
        center_scale(Xs);
        center_scale(Ds);
        ys.array() -= ys.mean();

        Eigen::Map<Eigen::MatrixXd> X_map(Xs.data(), n, p);
        Eigen::Map<Eigen::MatrixXd> D_map(Ds.data(), n, L);
        Eigen::Map<Eigen::VectorXd> y_map(ys.data(), n);

        const SolverOptions opts{.normalize = false, .intercept = false, .precentered = true,
                                 .tie_seed = seed};
        TOOLS_Solver downdated(X_map, D_map, y_map, opts);
        TOOLSFromScratch reference(X_map, D_map, y_map, opts);
        for (std::size_t T : {1, 5, 20, 120}) {  // warm starts through a T grid
            downdated.executeStep(T, true);
            reference.executeStep(T, true);
        }

        EXPECT_EQ(downdated.getActions(), reference.getActions()) << "seed " << seed;
        const auto pd = downdated.getBetaPathSparse();
        const auto pr = reference.getBetaPathSparse();
        ASSERT_EQ(pd.steps.size(), pr.steps.size()) << "seed " << seed;
        for (std::size_t s = 0; s < pd.steps.size(); ++s) {
            EXPECT_EQ(pd.steps[s].idx, pr.steps[s].idx) << "seed " << seed << ", step " << s;
            EXPECT_EQ(pd.steps[s].val, pr.steps[s].val) << "seed " << seed << ", step " << s;
        }
        // The paths run into the row bound of the centered design
        EXPECT_EQ(downdated.getNumActives(), static_cast<std::size_t>(n - 1)) << "seed " << seed;
    }
}

// ========================================================================================
} /* End of namespace trex::test::tsolvers::ad::linear_model */
