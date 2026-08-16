// ===================================================================================
// lambda2_cv_diagnostics.hpp
// ===================================================================================
/**
 * @file lambda2_cv_diagnostics.hpp
 *
 * @brief Health record of one lambda_2 CV tuner run, reported alongside the
 *        value the tuner returned.
 *
 * @details The selectors resolve lambda_2 internally and used to expose only
 *  the pick, so a boundary solution or a binding path truncation was
 *  invisible to the caller. This struct is what `computeLambda2()` fills in
 *  addition to the value; see `TRexGVSSelector::getLambda2CVDiagnostics()`
 *  and `TRexTikhonovSelector::getLambda2CVDiagnostics()`.
 */
// ===================================================================================
#ifndef TREX_ML_METHODS_MODEL_SELECTION_LAMBDA2_CV_DIAGNOSTICS_HPP
#define TREX_ML_METHODS_MODEL_SELECTION_LAMBDA2_CV_DIAGNOSTICS_HPP

// Eigen includes
#include <Eigen/Core>

namespace trex {
namespace ml_methods {
namespace model_selection {

/**
 * @brief What a lambda_2 CV tuner did.
 *
 * `ran == false` means lambda_2 was user-supplied and no tuner executed.
 * The structural fields are specific to the 2D IEN-geometry tuner
 * (`ienet_cv_ccd`); the other methods fill only `ran` and `seconds`.
 */
struct lambda2_cv_diagnostics {

    /** @brief A CV tuner ran (false: lambda_2 came from the caller). */
    bool ran = false;

    /** @brief Wall-clock seconds spent resolving lambda_2. */
    double seconds = 0.0;

    /** @brief Number of lambda_2 values whose profiled lambda_1 optimum
     *  landed on a truncation edge. Zero is the healthy case: the path
     *  truncation never touched the answer. Nonzero means boundary
     *  solutions — raise `dfmax_ratio` or lower `lambda1_min_ratio`. */
    Eigen::Index n_profile_at_edge = 0;

    /** @brief Support cap (glmnet `dfmax`) applied, or -1 when disabled. */
    Eigen::Index support_cap = -1;

    /** @brief p <= `support_cap`: the cap cannot bind, so the run pays the
     *  full cost of the deviance rules alone. Explains an expensive tuner
     *  call on a small-p design. */
    bool low_dimensional = false;
};

}  // namespace model_selection
}  // namespace ml_methods
}  // namespace trex

#endif  // TREX_ML_METHODS_MODEL_SELECTION_LAMBDA2_CV_DIAGNOSTICS_HPP
