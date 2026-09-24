// ========================================================================================
/**
 * @file test_trex_thread_hygiene.cpp
 *
 * @brief Eigen's process-global thread count must be exactly what it was
 *  before a selection once select() returns, on every runner path (PLAN WP4
 *  4.5, Arnau): the library sets its budget once from the calling thread,
 *  never inside an OpenMP team, and restores it. A leak here changes the
 *  threading of the host process (an R or Python session) behind its back.
 */
// ========================================================================================

// google test includes
#include <gtest/gtest.h>

// Eigen includes
#include <Eigen/Dense>

// std includes
#include <vector>

// project includes
#include <trex_selector_methods/trex_core/trex.hpp>
#include <trex_selector_methods/trex_gvs/trex_gvs.hpp>
#include <utils/datageneration/utils_datagen.hpp>
#include <utils/openmp/utils_openmp.hpp>

// ========================================================================================

namespace trex::test::trex_selector_methods::trex_core {

using namespace trex::trex_selector_methods::trex_core;
using namespace trex::utils::datageneration::datagen;
namespace sd  = trex::trex_selector_methods::utils::solver_dispatch;
namespace gvs = trex::trex_selector_methods::trex_gvs;

namespace {

/** @brief A sentinel no code path would choose on its own, restored afterwards. */
constexpr int kSentinelThreads = 3;

struct SentinelScope {
    int previous;
    SentinelScope() : previous(Eigen::nbThreads()) { Eigen::setNbThreads(kSentinelThreads); }
    ~SentinelScope() { Eigen::setNbThreads(previous); }
};

/**
 * @brief Whether Eigen's global thread count is live in this build.
 *
 * @details With EIGEN_USE_BLAS (the default whenever a BLAS is found — Apple
 *  Accelerate here) Eigen delegates its products to the BLAS and never
 *  spawns its own OpenMP teams: setNbThreads() is a no-op and nbThreads()
 *  is always 1, so a leak cannot be observed and cannot happen. The global
 *  is live in the no-BLAS configuration (TREX_USE_BLAS_LAPACK=OFF — the
 *  Python wheel's), which is where these tests bite. Skip otherwise.
 */
bool eigenThreadCountIsLive() {
    const int previous = Eigen::nbThreads();
    Eigen::setNbThreads(kSentinelThreads);
    const bool live = (Eigen::nbThreads() == kSentinelThreads);
    Eigen::setNbThreads(previous);
    return live;
}

#define TREX_SKIP_UNLESS_EIGEN_THREADS_LIVE()                                        \
    if (!eigenThreadCountIsLive()) {                                                  \
        GTEST_SKIP() << "Eigen's thread count is inert in this build (EIGEN_USE_BLAS); " \
                        "run the no-BLAS configuration (TREX_USE_BLAS_LAPACK=OFF).";  \
    }

void runCore(LLoopStrategy strategy, bool parallel) {
    std::vector<std::size_t> support = {1, 2, 3};
    std::vector<double> coefs = {5.0, -3.0, 2.0};
    SyntheticData data(60, 20, support, coefs, 1.0, 42);
    auto X = data.getX();
    auto y = data.getY();
    Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
    Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

    TRexControlParameter ctrl;
    ctrl.K = 3;
    ctrl.max_dummy_multiplier = 2;
    ctrl.lloop_strategy = strategy;
    if (parallel) {
        ctrl.parallel_rnd_experiments = true;
        ctrl.max_outer_threads = 2;
        ctrl.max_inner_threads = 1;
    }
    TRexSelector trex(X_map, y_map, 0.1, ctrl, 42, false);
    trex.select();
}

}  // namespace

// ========================================================================================

/** @brief Every sequential runner (stored, on-demand and permutation
 *  strategies) leaves Eigen's thread count as it found it. */
TEST(TRexThreadHygieneTest, SequentialRunnersRestoreEigenThreadCount) {
    TREX_SKIP_UNLESS_EIGEN_THREADS_LIVE();
    SentinelScope sentinel;
    for (const LLoopStrategy s : {LLoopStrategy::HCONCAT, LLoopStrategy::STANDARD,
                                  LLoopStrategy::SKIPL, LLoopStrategy::SEEDED,
                                  LLoopStrategy::PERMUTATION,
                                  LLoopStrategy::PERMUTATION_SEEDED}) {
        runCore(s, /*parallel=*/false);
        EXPECT_EQ(Eigen::nbThreads(), kSentinelThreads)
            << "strategy " << static_cast<int>(s) << " leaked the Eigen thread count";
    }
}


/** @brief The parallel-experiments runner sets its budget once and restores it. */
TEST(TRexThreadHygieneTest, ParallelRunnerRestoresEigenThreadCount) {
    TREX_SKIP_UNLESS_EIGEN_THREADS_LIVE();
    SentinelScope sentinel;
    for (const LLoopStrategy s : {LLoopStrategy::HCONCAT, LLoopStrategy::STANDARD}) {
        runCore(s, /*parallel=*/true);
        EXPECT_EQ(Eigen::nbThreads(), kSentinelThreads)
            << "strategy " << static_cast<int>(s) << " leaked the Eigen thread count";
    }
}


/** @brief The GVS in-memory K-loop, which runs its own parallel region. */
TEST(TRexThreadHygieneTest, GvsRunnerRestoresEigenThreadCount) {
    TREX_SKIP_UNLESS_EIGEN_THREADS_LIVE();
    SentinelScope sentinel;
    std::vector<std::size_t> support = {1, 2, 3};
    std::vector<double> coefs = {5.0, -3.0, 2.0};
    for (const bool parallel : {false, true}) {
        SyntheticData data(60, 20, support, coefs, 1.0, 42);
        auto X = data.getX();
        auto y = data.getY();
        Eigen::Map<Eigen::MatrixXd> X_map(X.data(), X.rows(), X.cols());
        Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());

        gvs::TRexGVSControlParameter ctrl;
        ctrl.gvs_type = gvs::GVSType::EN;
        ctrl.lambda_2 = 1.0;
        ctrl.trex_ctrl.K = 3;
        ctrl.trex_ctrl.max_dummy_multiplier = 2;
        if (parallel) {
            ctrl.trex_ctrl.parallel_rnd_experiments = true;
            ctrl.trex_ctrl.max_outer_threads = 2;
            ctrl.trex_ctrl.max_inner_threads = 1;
        }
        gvs::TRexGVSSelector sel(X_map, y_map, 0.1, ctrl, 42, false);
        sel.select();
        EXPECT_EQ(Eigen::nbThreads(), kSentinelThreads)
            << (parallel ? "parallel" : "serial") << " GVS run leaked the Eigen thread count";
    }
}

// ========================================================================================
} /* End of namespace trex::test::trex_selector_methods::trex_core */
