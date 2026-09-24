// ==============================================================================
#ifndef RCPP_TREX_TIKHONOV_H
#define RCPP_TREX_TIKHONOV_H
// ==============================================================================

// std includes
#include <memory>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// TRex includes
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>
#include "rcpp_trex_selector.h"

// ==============================================================================

using namespace trex::trex_selector_methods::ad::trex_tikhonov;
using namespace trex::trex_selector_methods::ad::trex_core;

// ==============================================================================

/**
 * @brief Rcpp wrapper class for TRexTikhonovSelector (general-Tikhonov T-Rex).
 *
 * Maps R objects to the C++ structures for the sparse-K selector. The nested
 * base control is merged from the separately-supplied trex_control; a
 * Tikhonov-family solver_type (TCDTIKNET / TIKNET / TIKOMP — the only
 * solvers with a general-K mode) is honored, anything else becomes the
 * TCDTIKNET default (the R level warns first on an explicit mismatch).
 */
class RTRexTikhonovSelector : public RTRexSelector {
public:
    RTRexTikhonovSelector(
        Eigen::Map<Eigen::MatrixXd> X,
        Eigen::Map<Eigen::VectorXd> y,
        double tFDR,
        const TRexTikhonovControlParameter& tik_control,
        const tc::TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X);
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y);

        namespace tik_sd = trex::trex_selector_methods::ad::utils::solver_dispatch;
        TRexTikhonovControlParameter tik_ctrl = tik_control;
        tik_ctrl.trex_ctrl = trex_control;
        // The C++ constructor rejects any non-Tikhonov solver; the
        // caller-facing default (TLARS) is neutral, so it becomes TCDTIKNET.
        const auto st = tik_ctrl.trex_ctrl.solver_type;
        const bool tik_family =
            st == tik_sd::SolverTypeForTRex::TCDTIKNET ||
            st == tik_sd::SolverTypeForTRex::TIKNET ||
            st == tik_sd::SolverTypeForTRex::TIKOMP;
        if (!tik_family) {
            tik_ctrl.trex_ctrl.solver_type =
                tik_sd::SolverTypeForTRex::TCDTIKNET;
        }

        this->selector_ = std::make_unique<TRexTikhonovSelector>(
            *(this->X_map_), *(this->y_map_), tFDR, tik_ctrl, seed, verbose
        );
    }

    TRexTikhonovSelector* get() const {
        return static_cast<TRexTikhonovSelector*>(this->selector_.get());
    }
};

#endif // RCPP_TREX_TIKHONOV_H
