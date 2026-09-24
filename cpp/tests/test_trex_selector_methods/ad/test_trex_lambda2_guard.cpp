// ========================================================================================
// test_trex_lambda2_guard.cpp
// ========================================================================================
/**
 * @file test_trex_lambda2_guard.cpp
 *
 * @brief The lambda_2 tuner boundary-pick guard at the selector level
 *        (PLAN WP6, TODO 4.12; warning-only by decision of 2026-09-24).
 *
 * @details Contracts pinned here:
 *   - GVS-IEN with the IEN CCD tuner on a pure-noise response: the pick is
 *     boundary-dominated (flat CV curve, grid ceiling), the diagnostics say
 *     so (boundary_pick, at_grid_edge, grid_index / grid_size), exactly ONE
 *     warning is emitted naming the selector, the method and the mode, and
 *     the same record is carried in the GVS result.
 *   - On any design the warning is emitted if and only if boundary_pick
 *     holds (signal design, IEN and EN tracks).
 *   - The EN track records the grid pick as well (profile_at_edge stays
 *     false: no nuisance dimension there).
 *   - The Tikhonov selector's IEN tuner warns the same way.
 *   - A fixed lambda_2 never warns and leaves the guard fields at their
 *     defaults.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

// std includes
#include <cmath>
#include <map>
#include <random>
#include <string>
#include <vector>

// project trex includes
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>

// ========================================================================================

namespace trex::test::trex_selector_methods::lambda2_guard {

using namespace trex::trex_selector_methods::ad::trex_core;
using namespace trex::trex_selector_methods::ad::trex_gvs;
using namespace trex::trex_selector_methods::ad::trex_tikhonov;
using namespace trex::trex_selector_methods::ad::utils::solver_dispatch;
namespace ms = trex::ml_methods::model_selection;

namespace {

/** @brief Two equicorrelated blocks (cols 0-2, 3-5) among iid nulls; the
 *  response is the planted group contrast plus noise, or pure noise. */
struct Design {
    Eigen::MatrixXd X;
    Eigen::VectorXd y;

    Design(Eigen::Index n, Eigen::Index p, double rho, bool signal, unsigned seed) {
        std::mt19937 rng(seed);
        std::normal_distribution<double> N01(0.0, 1.0);
        X.resize(n, p);
        const double a = std::sqrt(rho), b = std::sqrt(1.0 - rho);
        for (int g = 0; g < 2; ++g) {
            Eigen::VectorXd shared(n);
            for (Eigen::Index i = 0; i < n; ++i) { shared(i) = N01(rng); }
            for (Eigen::Index j = 0; j < 3; ++j) {
                for (Eigen::Index i = 0; i < n; ++i) {
                    X(i, 3 * g + j) = a * shared(i) + b * N01(rng);
                }
            }
        }
        for (Eigen::Index j = 6; j < p; ++j) {
            for (Eigen::Index i = 0; i < n; ++i) { X(i, j) = N01(rng); }
        }
        y.resize(n);
        for (Eigen::Index i = 0; i < n; ++i) {
            double s = 0.0;
            if (signal) {
                for (Eigen::Index j = 0; j < 3; ++j) { s += X(i, j); }
                for (Eigen::Index j = 3; j < 6; ++j) { s -= X(i, j); }
            }
            y(i) = s + N01(rng);
        }
    }
};

/** @brief Group-mean (informed elastic net) matrix W = sum_m 1_m 1_m^T / p_m. */
Eigen::SparseMatrix<double> groupMeanK(const std::vector<Eigen::Index>& groups) {
    const auto p = static_cast<Eigen::Index>(groups.size());
    std::map<Eigen::Index, std::vector<Eigen::Index>> members;
    for (Eigen::Index j = 0; j < p; ++j) { members[groups[static_cast<std::size_t>(j)]].push_back(j); }
    std::vector<Eigen::Triplet<double>> trip;
    for (const auto& [g, cols] : members) {
        const double w = 1.0 / static_cast<double>(cols.size());
        for (const Eigen::Index a : cols) {
            for (const Eigen::Index b : cols) { trip.emplace_back(a, b, w); }
        }
    }
    Eigen::SparseMatrix<double> W(p, p);
    W.setFromTriplets(trip.begin(), trip.end());
    return W;
}

std::size_t countOf(const std::string& hay, const std::string& needle) {
    std::size_t n = 0;
    for (std::size_t pos = hay.find(needle); pos != std::string::npos;
         pos = hay.find(needle, pos + needle.size())) { ++n; }
    return n;
}

TRexGVSControlParameter gvsIenControl(LambdaSelectionMethod method) {
    TRexGVSControlParameter ctrl;
    ctrl.gvs_type       = GVSType::IEN;
    ctrl.lambda_2       = -1.0;
    ctrl.lambda2_method = method;
    ctrl.cv_n_lambda2   = 10;
    ctrl.cv_n_lambda1   = 20;
    ctrl.trex_ctrl.solver_type = SolverTypeForTRex::TIENET;
    ctrl.trex_ctrl.K = 3;
    ctrl.trex_ctrl.max_dummy_multiplier = 2;
    return ctrl;
}

struct GvsRun {
    ms::lambda2_cv_diagnostics diag;
    ms::lambda2_cv_diagnostics result_diag;
    double lambda2_used;
    std::string stderr_text;
};

GvsRun runGvs(Design& d, const TRexGVSControlParameter& ctrl) {
    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());
    testing::internal::CaptureStderr();
    TRexGVSSelector sel(Xm, ym, 0.2, ctrl, 42, false);
    sel.select();
    const std::string err = testing::internal::GetCapturedStderr();
    return {sel.getLambda2CVDiagnostics(), sel.getGVSResult().lambda2_cv,
            sel.getGVSResult().lambda2_used, err};
}

} // namespace

// ========================================================================================

TEST(Lambda2Guard, GvsIenPureNoiseIsABoundaryPickAndWarnsOnce) {
#ifndef CRAN_BUILD
    Design d(150, 40, 0.75, /*signal=*/false, 7);
    const GvsRun r = runGvs(d, gvsIenControl(LambdaSelectionMethod::CV_1SE_IEN_CCD));

    EXPECT_TRUE(r.diag.ran);
    EXPECT_EQ(r.diag.grid_size, 10);
    EXPECT_TRUE(r.diag.boundary_pick)
        << "grid_index " << r.diag.grid_index << ", profile_at_edge " << r.diag.profile_at_edge;
    EXPECT_TRUE(r.diag.at_grid_edge);
    EXPECT_EQ(r.diag.grid_index, r.diag.grid_size - 1) << "expected the grid ceiling";

    EXPECT_EQ(countOf(r.stderr_text, "BOUNDARY pick"), 1u) << r.stderr_text;
    EXPECT_NE(r.stderr_text.find("TRexGVSSelector::computeLambda2: CV_1SE_IEN_CCD"),
              std::string::npos) << r.stderr_text;
    EXPECT_NE(r.stderr_text.find("last point of the lambda_2 grid"), std::string::npos);
    EXPECT_NE(r.stderr_text.find("returned unchanged"), std::string::npos);

    // The pick is returned unchanged and the record rides in the result.
    EXPECT_GT(r.lambda2_used, 0.0);
    EXPECT_EQ(r.result_diag.boundary_pick, r.diag.boundary_pick);
    EXPECT_EQ(r.result_diag.grid_index, r.diag.grid_index);
#else
    SUCCEED();
#endif
}


TEST(Lambda2Guard, GvsIenSignalDesignWarnsIffBoundaryPick) {
#ifndef CRAN_BUILD
    for (auto method : {LambdaSelectionMethod::CV_1SE_IEN_CCD,
                        LambdaSelectionMethod::CV_MIN_IEN_CCD}) {
        Design d(150, 40, 0.75, /*signal=*/true, 7);
        const GvsRun r = runGvs(d, gvsIenControl(method));
        EXPECT_TRUE(r.diag.ran);
        EXPECT_GE(r.diag.grid_index, 0);
        EXPECT_LT(r.diag.grid_index, r.diag.grid_size);
        EXPECT_EQ(r.diag.boundary_pick, r.diag.at_grid_edge || r.diag.profile_at_edge);
        EXPECT_EQ(countOf(r.stderr_text, "BOUNDARY pick"), r.diag.boundary_pick ? 1u : 0u)
            << r.stderr_text;
    }
#else
    SUCCEED();
#endif
}


TEST(Lambda2Guard, GvsEnTrackRecordsTheGridPickAndWarnsIffAtEdge) {
#ifndef CRAN_BUILD
    for (auto method : {LambdaSelectionMethod::CV_1SE_CCD, LambdaSelectionMethod::CV_MIN_CCD,
                        LambdaSelectionMethod::CV_1SE_SVD, LambdaSelectionMethod::CV_MIN_SVD}) {
        Design d(150, 40, 0.75, /*signal=*/true, 7);
        TRexGVSControlParameter ctrl;
        ctrl.gvs_type       = GVSType::EN;
        ctrl.lambda_2       = -1.0;
        ctrl.lambda2_method = method;
        ctrl.cv_n_lambda    = 50;
        ctrl.trex_ctrl.solver_type = SolverTypeForTRex::TENET;
        ctrl.trex_ctrl.K = 3;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        const GvsRun r = runGvs(d, ctrl);
        EXPECT_TRUE(r.diag.ran);
        EXPECT_GT(r.diag.grid_size, 0);
        EXPECT_GE(r.diag.grid_index, 0);
        EXPECT_LT(r.diag.grid_index, r.diag.grid_size);
        EXPECT_FALSE(r.diag.profile_at_edge);   // no nuisance dimension on the EN track
        EXPECT_EQ(r.diag.boundary_pick, r.diag.at_grid_edge);
        EXPECT_EQ(countOf(r.stderr_text, "BOUNDARY pick"), r.diag.boundary_pick ? 1u : 0u)
            << r.stderr_text;
    }
#else
    SUCCEED();
#endif
}


TEST(Lambda2Guard, GvsFixedLambda2NeverWarns) {
#ifndef CRAN_BUILD
    Design d(150, 40, 0.75, /*signal=*/false, 7);
    auto ctrl = gvsIenControl(LambdaSelectionMethod::CV_1SE_IEN_CCD);
    ctrl.lambda_2 = 1.0;
    const GvsRun r = runGvs(d, ctrl);
    EXPECT_FALSE(r.diag.ran);
    EXPECT_EQ(r.diag.grid_index, -1);
    EXPECT_EQ(r.diag.grid_size, 0);
    EXPECT_FALSE(r.diag.boundary_pick);
    EXPECT_EQ(r.stderr_text.find("BOUNDARY pick"), std::string::npos) << r.stderr_text;
#else
    SUCCEED();
#endif
}


TEST(Lambda2Guard, TikhonovIenPureNoiseIsABoundaryPickAndWarnsOnce) {
#ifndef CRAN_BUILD
    const Eigen::Index n = 150, p = 40;
    Design d(n, p, 0.75, /*signal=*/false, 7);
    std::vector<Eigen::Index> groups(static_cast<std::size_t>(p));
    for (Eigen::Index j = 0; j < p; ++j) { groups[static_cast<std::size_t>(j)] = j / 5; }

    TRexTikhonovControlParameter ctrl;
    ctrl.penalty_structure = TikPenaltyStructure::SPARSE_GENERAL;
    ctrl.tikhonov_K        = groupMeanK(groups);
    ctrl.prior_groups      = groups;
    ctrl.cluster_dummies   = true;
    ctrl.lambda_2          = -1.0;
    ctrl.lambda2_method    = TikLambda2Method::CV_1SE_IEN_CCD;
    ctrl.cv_n_lambda2      = 10;
    ctrl.cv_n_lambda1      = 20;
    ctrl.trex_ctrl.solver_type = SolverTypeForTRex::TCDTIKNET;
    ctrl.trex_ctrl.K = 3;
    ctrl.trex_ctrl.max_dummy_multiplier = 2;

    Eigen::Map<Eigen::MatrixXd> Xm(d.X.data(), d.X.rows(), d.X.cols());
    Eigen::Map<Eigen::VectorXd> ym(d.y.data(), d.y.size());
    testing::internal::CaptureStderr();
    TRexTikhonovSelector sel(Xm, ym, 0.2, ctrl, 42, false);
    sel.select();
    const std::string err = testing::internal::GetCapturedStderr();
    const auto diag = sel.getLambda2CVDiagnostics();

    EXPECT_TRUE(diag.ran);
    EXPECT_EQ(diag.grid_size, 10);
    EXPECT_TRUE(diag.boundary_pick)
        << "grid_index " << diag.grid_index << ", profile_at_edge " << diag.profile_at_edge;
    EXPECT_EQ(countOf(err, "BOUNDARY pick"), 1u) << err;
    EXPECT_NE(err.find("TRexTikhonovSelector::computeLambda2: CV_1SE_IEN_CCD"),
              std::string::npos) << err;
    EXPECT_GT(sel.getLambda2Used(), 0.0);
#else
    SUCCEED();
#endif
}

} // namespace trex::test::trex_selector_methods::lambda2_guard
