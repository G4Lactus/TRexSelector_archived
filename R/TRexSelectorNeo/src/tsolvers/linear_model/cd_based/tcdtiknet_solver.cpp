// ============================================================================
// tcdtiknet_solver.cpp
// ============================================================================
/**
 * @file tcdtiknet_solver.cpp
 *
 * @brief Implementation of the Terminating CCD Tikhonov Net (T-CDTIKNET)
 * solver class, extending the T-CDENET solver.
 */
// ============================================================================

// tsolvers includes
#include <tsolvers/linear_model/cd_based/tcdtiknet_solver.hpp>

// ============================================================================

// Embedded into namespace trex::tsolvers::linear_model::cd_based
namespace trex::tsolvers::linear_model::cd_based {

// ============================================================================
// Specialized penalty hooks (FOLDED dummy coupling)
// ============================================================================

void TCDTIKNET_Solver::ensureCdState() {
    TCDLASSO_Solver::ensureCdState();     // base caches + configurePenaltyDiag()

    if (!fold_dummies_) { return; }   // INDEPENDENT_RIDGE: base machinery

    // Rebuild the folded coupling state w = lambda2 * K * s from the
    // (possibly deserialized) coefficients: s_i sums variable i's
    // coefficient and those of all its dummy copies.
    const auto p_idx = static_cast<Eigen::Index>(p_original_);
    Eigen::VectorXd s = Eigen::VectorXd::Zero(p_idx);
    const std::size_t p_total = p_original_ + num_dummies_;
    for (std::size_t t = 0; t < p_total; ++t) {
        const double b = beta_cd_(static_cast<Eigen::Index>(t));
        if (b != 0.0) {
            s(static_cast<Eigen::Index>(foldedBase(t))) += b;
        }
    }
    folded_w_ = lambda2_ * (tikhonov_K_ * s);
}


void TCDTIKNET_Solver::configurePenaltyDiag() {
    if (!fold_dummies_) {
        TCDLASSO_Solver::configurePenaltyDiag();   // Tikhonov diag + dummy ridge
        return;
    }

    // FOLDED coupling: every unified column carries the diagonal of its
    // base variable, lambda2 * K_bb — dummies included (kappa_dummy_
    // is ignored; the dummy block is coupled through K, not ridged).
    const std::size_t p_total = p_original_ + num_dummies_;
    penalty_diag_.setZero(static_cast<Eigen::Index>(p_total));
    if (lambda2_ <= 0.0) { return; }

    Eigen::VectorXd kdiag(static_cast<Eigen::Index>(p_original_));
    for (std::size_t j = 0; j < p_original_; ++j) {
        kdiag(static_cast<Eigen::Index>(j)) =
            tikhonov_K_.coeff(static_cast<Eigen::Index>(j),
                              static_cast<Eigen::Index>(j));
    }
    for (std::size_t t = 0; t < p_total; ++t) {
        penalty_diag_(static_cast<Eigen::Index>(t)) =
            lambda2_ * kdiag(static_cast<Eigen::Index>(foldedBase(t)));
    }
}


void TCDTIKNET_Solver::penaltyRankUpdate(std::size_t j, double delta) {
    if (!fold_dummies_) {
        TCDLASSO_Solver::penaltyRankUpdate(j, delta);
        return;
    }

    // FOLDED coupling: s_base(j) += delta, maintained through
    // w += lambda2 * delta * K.col(base(j)) — O(nnz(K col)).
    const auto b = static_cast<Eigen::Index>(foldedBase(j));
    const double f = lambda2_ * delta;
    for (Eigen::SparseMatrix<double>::InnerIterator it(tikhonov_K_, b);
         it; ++it) {
        folded_w_(it.row()) += f * it.value();
    }
}

// ============================================================================
// Serialization
// ============================================================================

void TCDTIKNET_Solver::save(const std::string& filename) const {
    saveImpl(*this, filename);
}


TCDTIKNET_Solver TCDTIKNET_Solver::load(const std::string& filename,
                                        Eigen::Map<Eigen::MatrixXd>& X,
                                        Eigen::Map<Eigen::MatrixXd>& D) {
    return loadImpl<TCDTIKNET_Solver>(filename, X, D);
}

// ============================================================================
}  /* End of namespace trex::tsolvers::linear_model::cd_based */
