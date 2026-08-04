// ========================================================================================
// test_tsolvers_tik.cpp
// ========================================================================================
/**
 * @file test_tsolvers_tik.cpp
 *
 * @brief Unit tests for the native sparse-K Tikhonov solvers (TIKNET,
 * TIKOMP): degenerate anchors (lambda2 = 0 collapses, group-mean K ==
 * native TIENET, K = I == singleton-group TIENET), crossing-support
 * equivalence against the CCD GROUP_CENTERING reference, the augmented-data
 * OMP cross-check, the generalized-ridge refit closed form, and
 * save/load/continue round trips.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// std includes
#include <filesystem>
#include <set>
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

// tsolvers includes
#include <tsolvers/linear_model/cd_based/tcdienet_solver.hpp>
#include <tsolvers/linear_model/lars_based/tiknet_solver.hpp>
#include <tsolvers/linear_model/lars_based/tienet_solver.hpp>
#include <tsolvers/linear_model/lars_based/tlasso_solver.hpp>
#include <tsolvers/linear_model/omp_based/tikomp_solver.hpp>
#include <tsolvers/linear_model/omp_based/tomp_solver.hpp>

// utils includes
#include <utils/datageneration/utils_datagen.hpp>

// ========================================================================================

// Embed into test namespace
namespace trex::test::tsolvers::linear_model {

using namespace trex::tsolvers::linear_model::cd_based;
using namespace trex::tsolvers::linear_model::lars_based;
using namespace trex::tsolvers::linear_model::omp_based;
using namespace trex::utils::datageneration::datagen;
using trex::tsolvers::SparseBetaPath;

// ========================================================================================

namespace {

/** @brief Real-variable support at the FIRST recorded step with exactly t
 *  active dummies, for t = 1..T (the quantity the T-Rex phi calibration
 *  consumes). result[t-1] is empty when no recorded step hits count t. */
std::vector<std::set<std::size_t>> supportsAtDummyCounts(
    const SparseBetaPath& path, std::size_t p, std::size_t T,
    double eps = 1e-12) {

    std::vector<std::set<std::size_t>> result(T);
    std::vector<bool> found(T, false);

    for (const auto& st : path.steps) {
        std::size_t dummies = 0;
        for (std::size_t e = 0; e < st.idx.size(); ++e) {
            if (st.idx[e] >= p && std::abs(st.val[e]) > eps) { ++dummies; }
        }
        if (dummies >= 1 && dummies <= T && !found[dummies - 1]) {
            found[dummies - 1] = true;
            for (std::size_t e = 0; e < st.idx.size(); ++e) {
                if (st.idx[e] < p && std::abs(st.val[e]) > eps) {
                    result[dummies - 1].insert(st.idx[e]);
                }
            }
        }
    }
    return result;
}


/** @brief Group-mean IEN matrix W = sum_m 1_m 1_m^T / p_m for contiguous
 *  groups of size gsize (entries 1/gsize within each group block). */
Eigen::SparseMatrix<double> groupMeanK(Eigen::Index p, Eigen::Index gsize) {
    Eigen::SparseMatrix<double> W(p, p);
    std::vector<Eigen::Triplet<double>> trip;
    for (Eigen::Index i = 0; i < p; ++i) {
        for (Eigen::Index j = 0; j < p; ++j) {
            if (i / gsize == j / gsize) {
                trip.emplace_back(i, j, 1.0 / static_cast<double>(gsize));
            }
        }
    }
    W.setFromTriplets(trip.begin(), trip.end());
    return W;
}


/** @brief Within-group centering matrix K_GC = blockdiag(I - J/p_m) for
 *  contiguous groups of size gsize (idempotent projector). */
Eigen::SparseMatrix<double> groupCenteringK(Eigen::Index p,
                                            Eigen::Index gsize) {
    Eigen::SparseMatrix<double> K(p, p);
    std::vector<Eigen::Triplet<double>> trip;
    const double w = 1.0 / static_cast<double>(gsize);
    for (Eigen::Index i = 0; i < p; ++i) {
        for (Eigen::Index j = 0; j < p; ++j) {
            if (i / gsize == j / gsize) {
                const double delta = (i == j) ? 1.0 : 0.0;
                trip.emplace_back(i, j, delta - w);
            }
        }
    }
    K.setFromTriplets(trip.begin(), trip.end());
    return K;
}


/** @brief Identity K (p x p). */
Eigen::SparseMatrix<double> identityK(Eigen::Index p) {
    Eigen::SparseMatrix<double> K(p, p);
    K.setIdentity();
    return K;
}


/** @brief Center every column (and y) and scale the columns to unit
 *  L2-norm in place — the solvers' working scale, applied externally for
 *  the augmented cross-check (penalty rows must stay raw). */
void standardizeInPlace(Eigen::MatrixXd& X, Eigen::MatrixXd& D,
                        Eigen::VectorXd& y) {
    auto scale_cols = [](Eigen::MatrixXd& M) {
        for (Eigen::Index j = 0; j < M.cols(); ++j) {
            M.col(j).array() -= M.col(j).mean();
            const double nrm = M.col(j).norm();
            if (nrm > 0.0) { M.col(j) /= nrm; }
        }
    };
    scale_cols(X);
    scale_cols(D);
    y.array() -= y.mean();
}

} // namespace

// ========================================================================================

class TSolverTikTest : public ::testing::Test {
protected:
    Eigen::MatrixXd X, D;
    Eigen::VectorXd y;
    Eigen::VectorXi groups;   // 10 contiguous groups of 5 over p = 50

    void SetUp() override {
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

        groups.resize(X.cols());
        for (Eigen::Index j = 0; j < X.cols(); ++j) {
            groups(j) = static_cast<int>(j / 5);
        }
    }

    // Fresh maps over private copies (the solvers normalize in place).
    struct Mapped {
        Eigen::MatrixXd Xc, Dc;
        Eigen::VectorXd yc;
        Eigen::Map<Eigen::MatrixXd> Xm, Dm;
        Eigen::Map<Eigen::VectorXd> ym;
        Mapped(const Eigen::MatrixXd& X, const Eigen::MatrixXd& D,
               const Eigen::VectorXd& y)
            : Xc(X), Dc(D), yc(y),
              Xm(Xc.data(), Xc.rows(), Xc.cols()),
              Dm(Dc.data(), Dc.rows(), Dc.cols()),
              ym(yc.data(), yc.size()) {}
    };

    /** @brief Serialize -> load -> continue must replicate the single-shot
     *  run (K travels as triplets; the coupling state is restored or
     *  rebuilt). */
    template <class SolverT, class MakeFn>
    void expectRoundTripMatchesContinuous(MakeFn make, const std::string& tag,
                                          std::size_t T_partial = 3,
                                          std::size_t T_final = 6) {
        Mapped mr(X, D, y);
        SolverT ref = make(mr);
        ref.executeStep(T_final, true);

        Mapped mp(X, D, y);
        const std::string file =
            ::testing::TempDir() + "tik_ckpt_" + tag + ".bin";
        {
            SolverT part = make(mp);
            part.executeStep(T_partial, true);
            part.save(file);
        }
        ASSERT_TRUE(std::filesystem::exists(file));

        SolverT res = SolverT::load(file, mp.Xm, mp.Dm);
        res.executeStep(T_final, true);

        EXPECT_EQ(ref.getActions(), res.getActions())
            << tag << ": resumed path diverged from the continuous run.";
        EXPECT_TRUE(ref.getBeta(-1).isApprox(res.getBeta(-1), 1e-9))
            << tag << ": resumed final beta diverged from the continuous run.";

        std::filesystem::remove(file);
    }
};

// ========================================================================================
// TIKNET degenerate anchors
// ========================================================================================

/** @brief lambda2 == 0 collapses T-TIKNET exactly onto T-LASSO (every
 *  coupling term multiplies to literal zero; K is carried but inert). */
TEST_F(TSolverTikTest, TiknetLambda2ZeroReducesToTlasso) {
    Mapped mt(X, D, y);
    TIKNET_Solver tik(mt.Xm, mt.Dm, mt.ym, 0.0, groupCenteringK(X.cols(), 5),
                      true, true, false);
    tik.executeStep(5, true);

    Mapped ml(X, D, y);
    TLASSO_Solver lasso(ml.Xm, ml.Dm, ml.ym, true, true, false);
    lasso.executeStep(5, true);

    EXPECT_EQ(tik.getActions(), lasso.getActions());
    EXPECT_TRUE(tik.getBeta(-1).isApprox(lasso.getBeta(-1), 1e-9));
}


/** @brief With K = sum_m 1_m 1_m^T / p_m the LARS-TIK recursion computes
 *  the identical modified correlations, Gram couplings, and joining rates
 *  as the native LARS-IEN algorithm — TIENET becomes the special case.
 *  (Accumulation orders differ in the joining rate, so knife-edge designs
 *  could flip at ULP margins; this fixture is wide-margin and the action
 *  sequence is pinned exactly.) */
TEST_F(TSolverTikTest, TiknetGroupMeanKMatchesTienet) {
    const double lambda2 = 1.0;

    Mapped mt(X, D, y);
    TIKNET_Solver tik(mt.Xm, mt.Dm, mt.ym, lambda2, groupMeanK(X.cols(), 5),
                      true, true, false);
    tik.executeStep(5, true);

    Mapped mi(X, D, y);
    TIENET_Solver ien(mi.Xm, mi.Dm, mi.ym, lambda2, groups,
                      true, true, false);
    ien.executeStep(5, true);

    EXPECT_EQ(tik.getActions(), ien.getActions());
    EXPECT_TRUE(tik.getBeta(-1).isApprox(ien.getBeta(-1), 1e-8));
}


/** @brief K = I equals TIENET with all-singleton groups: each variable is
 *  ridged jointly with its own dummy copies through the folded sum (NOT
 *  the decoupled TENET ridge — the same collapse/coupling trade-off as on
 *  the CCD side). */
TEST_F(TSolverTikTest, TiknetIdentityKMatchesTienetSingletonGroups) {
    const double lambda2 = 0.5;

    Mapped mt(X, D, y);
    TIKNET_Solver tik(mt.Xm, mt.Dm, mt.ym, lambda2, identityK(X.cols()),
                      true, true, false);
    tik.executeStep(5, true);

    Eigen::VectorXi singleton = Eigen::VectorXi::LinSpaced(
        X.cols(), 0, static_cast<int>(X.cols()) - 1);
    Mapped mi(X, D, y);
    TIENET_Solver ien(mi.Xm, mi.Dm, mi.ym, lambda2, singleton,
                      true, true, false);
    ien.executeStep(5, true);

    EXPECT_EQ(tik.getActions(), ien.getActions());
    EXPECT_TRUE(tik.getBeta(-1).isApprox(ien.getBeta(-1), 1e-9));
}


// ========================================================================================
// TIKNET vs the CCD reference: crossing-support equivalence
// ========================================================================================

/** @brief T-TIKNET with K_GC and the CCD GROUP_CENTERING mode solve the
 *  same penalized-lasso path family from opposite algorithmic sides
 *  (piecewise-linear homotopy vs exact fixed-lambda1 minimizers); the real
 *  supports at each dummy-count crossing (what T-Rex consumes) must
 *  coincide. */
TEST_F(TSolverTikTest, TiknetKgcCrossingSupportsMatchCcdGroupCentering) {
    const std::size_t T = 5;
    const std::size_t p = static_cast<std::size_t>(X.cols());
    const double lambda2 = 2.0;

    Mapped mt(X, D, y);
    TIKNET_Solver tik(mt.Xm, mt.Dm, mt.ym, lambda2,
                      groupCenteringK(X.cols(), 5), true, true, false);
    tik.executeStep(T, true);

    Mapped mc(X, D, y);
    TCDIENET_Solver ccd(mc.Xm, mc.Dm, mc.ym, lambda2, groups,
                       true, true, false, trex::tsolvers::ScalingMode::L2,
                       GroupPenaltyVariant::CENTERING);
    ccd.setLambdaRelTol(1e-6);
    ccd.executeStep(T, true);

    const auto sup_tik = supportsAtDummyCounts(tik.getBetaPathSparse(), p, T);
    const auto sup_ccd = supportsAtDummyCounts(ccd.getBetaPathSparse(), p, T);

    for (std::size_t t = 1; t <= T; ++t) {
        EXPECT_EQ(sup_tik[t - 1], sup_ccd[t - 1])
            << "TIKNET vs CCD GROUP_CENTERING real support diverged at "
               "dummy count t = " << t;
    }
}


// ========================================================================================
// TIKOMP degenerate anchors
// ========================================================================================

/** @brief lambda2 == 0 collapses T-TIKOMP exactly onto T-OMP. */
TEST_F(TSolverTikTest, TikompLambda2ZeroReducesToTomp) {
    Mapped mt(X, D, y);
    TIKOMP_Solver tik(mt.Xm, mt.Dm, mt.ym, 0.0, groupCenteringK(X.cols(), 5),
                      true, true, false);
    tik.executeStep(5, true);

    Mapped mo(X, D, y);
    TOMP_Solver omp_ref(mo.Xm, mo.Dm, mo.ym, true, true, false);
    omp_ref.executeStep(5, true);

    EXPECT_EQ(tik.getActions(), omp_ref.getActions());
    EXPECT_TRUE(tik.getBeta(-1).isApprox(omp_ref.getBeta(-1), 1e-9));
}


/** @brief Augmented-data cross-check (Theorem-1 pattern): plain T-OMP on
 *  X_aug = [X; sqrt(l2) Gamma], y_aug = [y; 0],
 *  D_aug = [D; sqrt(l2) Gamma] (one layer) with Gamma = K_GC (idempotent
 *  projector: Gamma^T Gamma = K_GC) must walk the same greedy path as the
 *  native generalized-ridge TIKOMP. Both sides standardize the DATA blocks
 *  only, so the augmented run gets externally standardized inputs and
 *  normalize = false. */
TEST_F(TSolverTikTest, TikompNativeMatchesAugmentedTomp) {
    const double lambda2 = 2.0;
    const Eigen::Index n = X.rows();
    const Eigen::Index p = X.cols();
    const double sq = std::sqrt(lambda2);

    // Shared standardized data blocks.
    Eigen::MatrixXd Xs = X, Ds = D;
    Eigen::VectorXd ys = y;
    standardizeInPlace(Xs, Ds, ys);

    const Eigen::MatrixXd Gamma =
        Eigen::MatrixXd(groupCenteringK(p, 5));

    // Native TIKOMP on the standardized blocks.
    Eigen::MatrixXd Xn = Xs, Dn = Ds;
    Eigen::VectorXd yn = ys;
    Eigen::Map<Eigen::MatrixXd> Xn_m(Xn.data(), Xn.rows(), Xn.cols());
    Eigen::Map<Eigen::MatrixXd> Dn_m(Dn.data(), Dn.rows(), Dn.cols());
    Eigen::Map<Eigen::VectorXd> yn_m(yn.data(), yn.size());
    TIKOMP_Solver tik(Xn_m, Dn_m, yn_m, lambda2, groupCenteringK(p, 5),
                      false, false, false);
    tik.executeStep(5, true);

    // Augmented system (penalty rows raw, no re-normalization).
    Eigen::MatrixXd X_aug(n + p, p);
    X_aug.topRows(n) = Xs;
    X_aug.bottomRows(p) = sq * Gamma;
    Eigen::MatrixXd D_aug(n + p, p);   // one dummy layer: Gamma tiled once
    D_aug.topRows(n) = Ds;
    D_aug.bottomRows(p) = sq * Gamma;
    Eigen::VectorXd y_aug = Eigen::VectorXd::Zero(n + p);
    y_aug.head(n) = ys;

    Eigen::Map<Eigen::MatrixXd> Xa_m(X_aug.data(), X_aug.rows(), X_aug.cols());
    Eigen::Map<Eigen::MatrixXd> Da_m(D_aug.data(), D_aug.rows(), D_aug.cols());
    Eigen::Map<Eigen::VectorXd> ya_m(y_aug.data(), y_aug.size());
    TOMP_Solver aug(Xa_m, Da_m, ya_m, false, false, false);
    aug.executeStep(5, true);

    EXPECT_EQ(tik.getActions(), aug.getActions());
    EXPECT_TRUE(tik.getBeta(-1).isApprox(aug.getBeta(-1), 1e-8));
}


/** @brief K = I ridge-OMP closed-form sanity (no dummies): the final
 *  active-set refit must equal (X_A^T X_A + lambda2 I)^{-1} X_A^T y
 *  computed directly. */
TEST_F(TSolverTikTest, TikompIdentityKRidgeRefitClosedForm) {
    const double lambda2 = 3.0;
    const Eigen::Index p = X.cols();

    Eigen::MatrixXd Xs = X;
    Eigen::MatrixXd D_empty(X.rows(), 0);
    Eigen::VectorXd ys = y;
    standardizeInPlace(Xs, D_empty, ys);

    Eigen::Map<Eigen::MatrixXd> Xm(Xs.data(), Xs.rows(), Xs.cols());
    Eigen::Map<Eigen::MatrixXd> Dm(D_empty.data(), D_empty.rows(), 0);
    Eigen::Map<Eigen::VectorXd> ym(ys.data(), ys.size());

    TIKOMP_Solver tik(Xm, Dm, ym, lambda2, identityK(p), false, false, false);
    tik.executeStep(0, true);   // full path (no dummies to terminate on)
    ASSERT_GT(tik.getNumSteps(), 0u);

    const auto& actives = tik.getActives();
    ASSERT_FALSE(actives.empty());

    Eigen::MatrixXd XA(Xs.rows(), static_cast<Eigen::Index>(actives.size()));
    for (std::size_t k = 0; k < actives.size(); ++k) {
        XA.col(static_cast<Eigen::Index>(k)) =
            Xs.col(static_cast<Eigen::Index>(actives[k]));
    }
    const Eigen::MatrixXd G =
        XA.transpose() * XA +
        lambda2 * Eigen::MatrixXd::Identity(
            static_cast<Eigen::Index>(actives.size()),
            static_cast<Eigen::Index>(actives.size()));
    const Eigen::VectorXd beta_ref = G.ldlt().solve(XA.transpose() * ys);

    const Eigen::VectorXd beta = tik.getBeta(-1);
    for (std::size_t k = 0; k < actives.size(); ++k) {
        EXPECT_NEAR(beta(static_cast<Eigen::Index>(actives[k])),
                    beta_ref(static_cast<Eigen::Index>(k)), 1e-8)
            << "refit mismatch at active variable " << actives[k];
    }
}


// ========================================================================================
// Serialization: save -> load -> continue
// ========================================================================================

/** @brief T-TIKNET round trip: K travels as triplets; the folded coupling
 *  state (s, K s) is part of the archive. */
TEST_F(TSolverTikTest, SerializeLoadContinueMatchesContinuous_TIKNET) {
    const Eigen::Index p = X.cols();
    expectRoundTripMatchesContinuous<TIKNET_Solver>(
        [p](Mapped& m) {
            return TIKNET_Solver(m.Xm, m.Dm, m.ym, 2.0,
                                 groupCenteringK(p, 5), true, true, false);
        },
        "tiknet");
}


/** @brief T-TIKOMP round trip: K travels as triplets; the coupling state
 *  is transient and rebuilt from the refit coefficients. */
TEST_F(TSolverTikTest, SerializeLoadContinueMatchesContinuous_TIKOMP) {
    const Eigen::Index p = X.cols();
    expectRoundTripMatchesContinuous<TIKOMP_Solver>(
        [p](Mapped& m) {
            return TIKOMP_Solver(m.Xm, m.Dm, m.ym, 2.0,
                                 groupCenteringK(p, 5), true, true, false);
        },
        "tikomp");
}

// ========================================================================================
} /* End of namespace trex::test::tsolvers::linear_model */
