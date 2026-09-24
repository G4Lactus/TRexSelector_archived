// ===================================================================================
// tik_lambda2_power_rule.hpp
// ===================================================================================
#ifndef TREX_SELECTOR_METHODS_TREX_TIKHONOV_TIK_LAMBDA2_POWER_RULE_HPP
#define TREX_SELECTOR_METHODS_TREX_TIKHONOV_TIK_LAMBDA2_POWER_RULE_HPP
// ===================================================================================
/**
 * @file tik_lambda2_power_rule.hpp
 *
 * @brief Data-driven lambda_2 for the Tikhonov T-Rex (and for the elastic-net
 *        ridge weight of the GVS T-Rex) by the selector's own principle: on a
 *        grid of lambda_2 values, take the one at which the FDR-controlled
 *        selection is LARGEST.
 *
 * @details
 *   **Why not a fit criterion.** In the Tikhonov T-Rex the quadratic weight
 *   lambda_2 does not decide how much noise enters the selection — the dummies
 *   do that, at every lambda_2 — it decides how many true members of a
 *   structured group the path reaches before the dummies stop it. Prediction
 *   criteria (cross-validation of the ridge or of a post-selection refit, BIC
 *   and its relatives) judge a fitted estimate and are known to land two
 *   orders of magnitude above the selection-optimal range on K_center, or to
 *   fall back to the lasso support wherever the extra members carry
 *   coefficients an individual fit cannot justify (Scientific Notes,
 *   tik_lambda2_determination.md). The T-Rex calibration itself chooses its
 *   internal parameters (the stopping time and the voting level) by
 *   maximizing the number of selected variables subject to the estimated FDP;
 *   this header applies the same rule to lambda_2.
 *
 *   **The rule.** For every lambda_2 on a user grid one full selection is run
 *   with the SAME seed, so the probes differ by the penalty alone (common
 *   dummies: the base DummyGenerator's HCONCAT streams are prefix-stable in
 *   the seed). With `n_k = |selection at grid point k|`:
 *
 *     - MAX_SELECTION : k* = argmax_k n_k, ties resolved toward the SMALLEST
 *                       lambda_2 (the least reliance on the prior);
 *     - KNEE          : the smallest k with n_k >= max_k n_k - knee_tolerance
 *                       (the "minimum sufficient shrinkage" reading of the
 *                       same profile).
 *
 *   The result carries the whole profile — every probe's selection, stopping
 *   time and cost — so criteria that need the selections (a refit, an
 *   information criterion, a prediction error) can be evaluated on the same
 *   runs without repeating them.
 *
 *   **FDR status.** Every probe is an FDR-controlled selection at the target.
 *   The pick over the grid depends on y and is therefore OUTSIDE the
 *   fixed-lambda_2 guarantee; under the global null the union bound gives
 *   P(any selection) <= |grid| x P(any selection at one lambda_2), and under
 *   an alternative the exposure is the spread of the counts across the grid,
 *   which on structured designs is dominated by true positives. The rule is
 *   validated empirically (demo_14 in the Tikhonov T-Rex workspace); a
 *   sample-splitting use (pick on one half of the rows, select on the other)
 *   restores the guarantee at a power cost.
 */
// ===================================================================================

// std includes
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <utility>
#include <vector>

// Eigen includes
#include <Eigen/Core>

// T-Rex selector includes
#include <trex_selector_methods/ad/trex_gvs/trex_gvs.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/tik_lambda2_pick.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>

// ===================================================================================

namespace trex {
namespace trex_selector_methods {
namespace ad {
namespace trex_tikhonov {

// ===================================================================================
// Configuration (the rule enum and the pick live in tik_lambda2_pick.hpp)
// ===================================================================================

/** @brief The grid and the rule. */
struct Lambda2PowerRuleConfig {

    /** @brief lambda_2 candidates (>= 0, finite; sorted ascending and
     *  de-duplicated by the driver). Include 0 when "no penalty" must be
     *  an admissible answer — a prior that misfits the truth loses power
     *  monotonically and the rule then switches it off. */
    std::vector<double> grid;

    /** @brief The pick rule. Default: MAX_SELECTION. */
    Lambda2PowerRule rule = Lambda2PowerRule::MAX_SELECTION;

    /** @brief KNEE only: how many selections below the maximum still count
     *  as "on the plateau". Default: 1. */
    double knee_tolerance = 1.0;
};


// ===================================================================================
// Results
// ===================================================================================

/** @brief One grid point's selection. */
struct Lambda2Probe {

    /** @brief The lambda_2 this probe ran at. */
    double lambda2 = 0.0;

    /** @brief The selected original variables (0-based, ascending). */
    std::vector<std::size_t> selected;

    /** @brief The calibrated stopping time T_stop of the run. */
    std::size_t T_stop = 0;

    /** @brief Wall-clock seconds of the run. */
    double seconds = 0.0;
};

/** @brief The profile and the pick. */
struct Lambda2PowerRuleResult {

    /** @brief The chosen lambda_2. */
    double lambda2 = 0.0;

    /** @brief Index of the pick into `probes`. */
    std::size_t index = 0;

    /** @brief Every grid point's selection, ascending in lambda_2. */
    std::vector<Lambda2Probe> probes;

    /** @brief The pick diagnostics. */
    Lambda2PickInfo pick;

    /** @brief The pick sits on the first / last grid point (a boundary
     *  solution: widen the grid before trusting it). */
    bool at_lower_end = false;
    bool at_upper_end = false;

    /** @brief Wall-clock seconds of the whole profile. */
    double seconds = 0.0;

    /** @brief The selection at the chosen lambda_2. */
    [[nodiscard]] const std::vector<std::size_t>& selection() const { return probes[index].selected; }

    /** @brief The selection counts, ascending in lambda_2. */
    [[nodiscard]] std::vector<std::size_t> counts() const
    {
        std::vector<std::size_t> c(probes.size());
        for (std::size_t k = 0; k < probes.size(); ++k) { c[k] = probes[k].selected.size(); }
        return c;
    }
};


// ===================================================================================
// The grid
// ===================================================================================

/** @brief The validated, ascending, de-duplicated grid (tik_lambda2_pick.hpp). */
inline std::vector<double> lambda2_power_rule_grid(const Lambda2PowerRuleConfig& cfg)
{
    return validate_lambda2_grid(cfg.grid);
}


// ===================================================================================
// The driver
// ===================================================================================

/**
 * @brief Profile a selector over the grid and pick.
 *
 * @tparam Probe  Callable `Lambda2Probe(double lambda2)` running ONE full
 *                selection at the given lambda_2 (same data, same seed) and
 *                returning at least `selected` and `T_stop`; `lambda2` and
 *                `seconds` are filled in here.
 */
template <typename Probe>
Lambda2PowerRuleResult run_lambda2_power_rule(const Lambda2PowerRuleConfig& cfg, Probe&& probe)
{
    const std::vector<double> grid = lambda2_power_rule_grid(cfg);
    const auto t_all = std::chrono::steady_clock::now();

    Lambda2PowerRuleResult res;
    res.probes.reserve(grid.size());
    for (const double l2 : grid) {
        const auto t0 = std::chrono::steady_clock::now();
        Lambda2Probe pr = probe(l2);
        pr.lambda2 = l2;
        pr.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::sort(pr.selected.begin(), pr.selected.end());
        res.probes.push_back(std::move(pr));
    }

    res.pick         = pick_lambda2_index(res.counts(), cfg.rule, cfg.knee_tolerance);
    res.index        = res.pick.index;
    res.lambda2      = res.probes[res.index].lambda2;
    res.at_lower_end = (res.index == 0);
    res.at_upper_end = (res.index + 1 == res.probes.size());
    res.seconds      = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_all).count();
    return res;
}


// ===================================================================================
// Convenience overloads for the two selectors
// ===================================================================================

/**
 * @brief The rule on a TRexTikhonovSelector.
 *
 * @param X     Design (n x p); copied once, the copy is normalized in place
 *              per probe and restored (base-class contract).
 * @param y     Response (n).
 * @param tFDR  Target FDR.
 * @param ctrl  Tikhonov control; `lambda_2` is overwritten per probe, the
 *              rest (penalty structure, K, groups, dummies, base control)
 *              is shared by every probe.
 * @param seed  Selector seed, the SAME for every probe (common dummies).
 * @param cfg   Grid and rule.
 * @param verbose Selector verbosity.
 */
inline Lambda2PowerRuleResult tikhonov_lambda2_power_rule(const Eigen::MatrixXd& X,
                                                          const Eigen::VectorXd& y,
                                                          double tFDR,
                                                          TRexTikhonovControlParameter ctrl,
                                                          int seed,
                                                          const Lambda2PowerRuleConfig& cfg,
                                                          bool verbose = false)
{
    Eigen::MatrixXd Xc = X;
    Eigen::VectorXd yc = y;
    Eigen::Map<Eigen::MatrixXd> X_map(Xc.data(), Xc.rows(), Xc.cols());
    Eigen::Map<Eigen::VectorXd> y_map(yc.data(), yc.size());
    return run_lambda2_power_rule(cfg, [&](double l2) {
        ctrl.lambda_2 = l2;
        TRexTikhonovSelector selector(X_map, y_map, tFDR, ctrl, seed, verbose);
        selector.select();
        Lambda2Probe pr;
        pr.selected = selector.getSelectedIndices();
        pr.T_stop   = selector.getTStop();
        return pr;
    });
}

/**
 * @brief The rule on a TRexGVSSelector (the elastic-net / informed-elastic-net
 *        ridge weight). Same contract as the Tikhonov overload; `lambda_2 = 0`
 *        on the grid is the plain lasso T-Rex.
 */
inline Lambda2PowerRuleResult gvs_lambda2_power_rule(const Eigen::MatrixXd& X,
                                                     const Eigen::VectorXd& y,
                                                     double tFDR,
                                                     trex_gvs::TRexGVSControlParameter ctrl,
                                                     int seed,
                                                     const Lambda2PowerRuleConfig& cfg,
                                                     bool verbose = false)
{
    Eigen::MatrixXd Xc = X;
    Eigen::VectorXd yc = y;
    Eigen::Map<Eigen::MatrixXd> X_map(Xc.data(), Xc.rows(), Xc.cols());
    Eigen::Map<Eigen::VectorXd> y_map(yc.data(), yc.size());
    return run_lambda2_power_rule(cfg, [&](double l2) {
        ctrl.lambda_2 = l2;
        trex_gvs::TRexGVSSelector selector(X_map, y_map, tFDR, ctrl, seed, verbose);
        selector.select();
        Lambda2Probe pr;
        pr.selected = selector.getSelectedIndices();
        pr.T_stop   = selector.getGVSResult().T_stop;
        return pr;
    });
}

}  // namespace trex_tikhonov
}  // namespace ad
}  // namespace trex_selector_methods
}  // namespace trex

#endif  // TREX_SELECTOR_METHODS_TREX_TIKHONOV_TIK_LAMBDA2_POWER_RULE_HPP
