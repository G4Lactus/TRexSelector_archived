// ===================================================================================
// tikomp_solver.cpp
// ===================================================================================
/**
 * @file tikomp_solver.cpp
 *
 * @brief Implementation of the Terminating Tikhonov-K OMP (T-TIKOMP) solver
 * (generalized-ridge refits, FOLDED dummy coupling).
 */
// ===================================================================================

// std includes
#include <algorithm>
#include <cmath>
#include <stdexcept>

// tsolvers includes
#include <tsolvers/linear_model/omp_based/tikomp_solver.hpp>

// utils includes
#include <utils/openmp/utils_openmp.hpp>

// ==================================================================================

// Embedded into trex::tsolvers::linear_model::omp_based namespace
namespace trex::tsolvers::linear_model::omp_based {

// ==================================================================================
// Constructor Implementation
// ==================================================================================

TIKOMP_Solver::TIKOMP_Solver(
    Eigen::Map<Eigen::MatrixXd>& X,
    Eigen::Map<Eigen::MatrixXd>& D,
    Eigen::Map<Eigen::VectorXd>& y,
    double lambda2,
    const Eigen::SparseMatrix<double>& tikhonov_matrix,
    bool normalize,
    bool intercept,
    bool verbose,
    ScalingMode scaling_mode)
    : TOMP_Solver(X, D, y, normalize, intercept, verbose,
                  SolverTypeOMPBased::TIKOMP, scaling_mode),
      lambda2_(lambda2) {

    if (lambda2_ < 0.0) {
        throw std::invalid_argument(
            "TIKOMP_Solver: quadratic penalty lambda2 must be non-negative.");
    }

    const std::size_t p = p_original_;
    const std::size_t L = num_dummies_;
    const auto p_idx = static_cast<Eigen::Index>(p);

    // ---- Validate the penalty matrix (mirrors TIKNET / TCDTIKNET sparse) --
    if (tikhonov_matrix.rows() != p_idx || tikhonov_matrix.cols() != p_idx) {
        throw std::invalid_argument(concatMsg(
            "TIKOMP_Solver: Tikhonov matrix is ", tikhonov_matrix.rows(),
            "x", tikhonov_matrix.cols(), ", expected ", p, "x", p));
    }
    if (lambda2_ > 0.0 && L > 0 && L % p != 0) {
        throw std::invalid_argument(
            "TIKOMP_Solver: D.cols() must be a multiple of X.cols() "
            "(FOLDED layered dummy convention).");
    }
    tikhonov_K_ = tikhonov_matrix;
    tikhonov_K_.makeCompressed();

    folded_s_ = Eigen::VectorXd::Zero(p_idx);
    ks_       = Eigen::VectorXd::Zero(p_idx);

    // ---- State limits: rank(X^T X + lambda2 K-tilde) <= effective_n +
    //      rank(K); the count of structurally nonzero K columns is a cheap
    //      upper bound (the Cholesky pivot guard handles the boundary). ----
    std::size_t pcols_available = p + L - dropped_indices_.size();
    if (lambda2_ > 0.0) {
        std::size_t k_col_bound = 0;
        for (int j = 0; j < tikhonov_K_.outerSize(); ++j) {
            if (Eigen::SparseMatrix<double>::InnerIterator(tikhonov_K_, j)) {
                ++k_col_bound;
            }
        }
        max_actives_ = std::min(pcols_available, effective_n_ + k_col_bound);
    } else {
        max_actives_ = std::min(pcols_available, effective_n_);
    }
    maxSteps_ = 8 * max_actives_;

    // The TOMP base constructor already computed correlations_ = X^T y over
    // all inactive columns; at initialization s = 0, so the penalized
    // correlations coincide with the plain ones — no rescaling needed.
}


// ============================================================================
// Core Algorithm Implementation
// ============================================================================

void TIKOMP_Solver::executeStep(std::size_t T_stop, bool early_stop) {
    validateConnected();
    // Set openMP dominance for the entire algorithm
    EigenSingleThreadGuard eigen_single_thread_guard;

    // Mirrors TOMP_Solver::executeStep with the active-set bound raised to
    // max_actives_ (the penalized Gram G_A = X_A^T X_A + lambda2 K-tilde_A
    // stays positive definite beyond rank n whenever K contributes rank).
    while (currentStep_ < maxSteps_ && !inactives_.empty() &&
           actives_.size() < max_actives_ &&
           (!early_stop || T_stop == 0 || count_active_dummies_ < T_stop)) {

        // ========================================================
        // STEP 1: Max penalized correlation
        // ========================================================
        auto [Cmax, new_vars] = findTiedMaxCorrelations();
        if (new_vars.empty() || Cmax < 100 * eps_) {
            logInfo("Max |corr| approx 0 or no variables to add; exiting.");
            break;
        }
        pruneTiedDummies(new_vars, T_stop, early_stop);
        applyExchangeableTieBreak(new_vars);

        // ========================================================
        // STEP 2: Active set update (penalized Cholesky)
        // ========================================================
        actions_.emplace_back(updateActiveSet(new_vars));
        if (actives_.empty() || R_.size() == 0) {
            logWarning("Active set or Cholesky is empty; exiting.");
            break;
        }

        // ========================================================
        // STEP 3: Generalized-ridge refit and residuals
        // ========================================================
        currentStep_++;
        updateBetaPath();
        updateResiduals();

        double rss = r_.dot(r_);
        RSS_.push_back(rss);
        R2_.push_back(1.0 - rss / RSS_[0]);
        updateDummyTracking();
        DoF_.push_back(actives_.size() + (intercept_ ? 1 : 0));

        // ========================================================
        // STEP 4: Inactives and penalized correlations update
        // ========================================================
        updateInactiveSet();
        updateCorrelations();
    }
}

// ==================================================================================
// Math Helpers (Routed via getColumn)
// ==================================================================================

void TIKOMP_Solver::updateCorrelations() {
    // Rebuild the coupling state from the current refit coefficients: the
    // full OLS-style refit re-solves ALL active coefficients each step, so
    // incremental maintenance has no advantage here (and this keeps the
    // state correct after deserialization/reconnect as well).
    const auto p_idx = static_cast<Eigen::Index>(p_original_);
    folded_s_.setZero(p_idx);
    for (const std::size_t j : actives_) {
        folded_s_(static_cast<Eigen::Index>(foldedBase(j))) += runningBeta(j);
    }
    ks_ = tikhonov_K_ * folded_s_;

    const std::size_t num_inactives = inactives_.size();
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < num_inactives; ++i) {
        const std::size_t j = inactives_[i];
        correlations_(static_cast<Eigen::Index>(j)) =
            getColumn(j).dot(r_) -
            lambda2_ * ks_(static_cast<Eigen::Index>(foldedBase(j)));
    }
}


std::vector<int> TIKOMP_Solver::updateActiveSet(
    const std::vector<std::size_t>& new_vars) {

    std::vector<int> actions_this_step;
    any_dropped_ = false;

    for (auto j_new : new_vars) {
        // updateRTIK never mutates R_, so only the rank needs restoring on
        // rejection (mirror of TOMP_Solver::updateActiveSet).
        int rankR_backup = last_updateR_rank_;
        std::size_t active_size_backup = actives_.size();

        Eigen::MatrixXd newR = updateRTIK(j_new, eps_);

        if (last_updateR_rank_ == static_cast<int>(active_size_backup)) {
            last_updateR_rank_ = rankR_backup;
            dropped_indices_.push_back(j_new);
            actions_this_step.push_back(actionDrop(j_new));
            any_dropped_ = true;
            logWarning(concatMsg("Variable ", j_new, " collinear; dropped."));
        } else {
            R_ = newR;
            actives_.push_back(j_new);
            actions_this_step.push_back(actionAdd(j_new));
            num_additions_++;
            std::erase(inactives_, j_new);
            if (j_new >= dummy_start_idx_) count_active_dummies_++;
        }
    }
    return actions_this_step;
}


Eigen::MatrixXd TIKOMP_Solver::updateRTIK(std::size_t j_new, double eps) {
    std::size_t m = actives_.size();
    const auto b_new = static_cast<Eigen::Index>(foldedBase(j_new));

    // Gather K column b_new into the dense workspace (K is symmetric, so
    // the column entries indexed by row b(u) are exactly K_{b(t), b(u)}).
    if (kcol_gather_.size() != static_cast<Eigen::Index>(p_original_)) {
        kcol_gather_ = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(p_original_));
    }
    for (Eigen::SparseMatrix<double>::InnerIterator it(tikhonov_K_, b_new);
         it; ++it) {
        kcol_gather_(it.row()) = it.value();
    }

    Eigen::Ref<const Eigen::VectorXd> xnew = getColumn(j_new);
    // Penalized diagonal entry: x_j^T x_j + lambda2 * K_bb.
    double xtx = xnew.dot(xnew) + lambda2_ * kcol_gather_(b_new);

    Eigen::MatrixXd newR;
    if (R_.size() == 0 || m == 0) {
        newR.resize(1, 1);
        newR(0, 0) = std::sqrt(xtx);
        last_updateR_rank_ = 1;
    } else {
        // Cross-products with the K coupling on the active bases.
        Eigen::VectorXd cross_prod(m);
        for (std::size_t k = 0; k < m; ++k) {
            double v = xnew.dot(getColumn(actives_[k])) +
                       lambda2_ * kcol_gather_(static_cast<Eigen::Index>(
                           foldedBase(actives_[k])));
            cross_prod(static_cast<Eigen::Index>(k)) = v;
        }

        Eigen::VectorXd r = backsolveT(R_, cross_prod);
        double rpp_sq = xtx - r.dot(r);
        double rpp{};
        int new_rank = static_cast<int>(R_.rows());

        if (rpp_sq < eps * xtx) {
            rpp = eps;
        } else {
            rpp = std::sqrt(rpp_sq);
            ++new_rank;
        }

        newR = Eigen::MatrixXd::Zero(R_.rows() + 1, R_.cols() + 1);
        newR.topLeftCorner(R_.rows(), R_.cols()) = R_;
        newR.block(0, R_.cols(), R_.rows(), 1) = r;
        newR(R_.rows(), R_.cols()) = rpp;

        last_updateR_rank_ = new_rank;
    }

    // Reset the gathered column (keep the workspace all-zero between calls).
    for (Eigen::SparseMatrix<double>::InnerIterator it(tikhonov_K_, b_new);
         it; ++it) {
        kcol_gather_(it.row()) = 0.0;
    }
    return newR;
}


// ==================================================================================
// Serialization
// ==================================================================================

void TIKOMP_Solver::save(const std::string& filename) const {
    saveImpl(*this, filename);
}


TIKOMP_Solver TIKOMP_Solver::load(const std::string& filename,
                                  Eigen::Map<Eigen::MatrixXd>& X,
                                  Eigen::Map<Eigen::MatrixXd>& D) {
    return loadImpl<TIKOMP_Solver>(filename, X, D);
}

// ==================================================================================
} /* End of namespace trex::tsolvers::linear_model::omp_based */
