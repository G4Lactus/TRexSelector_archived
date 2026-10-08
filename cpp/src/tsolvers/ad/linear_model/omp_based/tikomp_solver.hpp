// ===================================================================================
// tikomp_solver.hpp
// ===================================================================================
#ifndef TSOLVERS_LINEAR_MODEL_OMP_BASED_TIKOMP_SOLVER_HPP
#define TSOLVERS_LINEAR_MODEL_OMP_BASED_TIKOMP_SOLVER_HPP
// ===================================================================================
/**
 * @file tikomp_solver.hpp
 *
 * @brief Terminating Tikhonov-K OMP (T-TIKOMP) solver: greedy orthogonal
 * matching pursuit under a general sparse-K quadratic penalty with FOLDED
 * dummy coupling, extending the T-OMP solver.
 *
 * @details
 *  Native generalized-ridge OMP for the objective family
 *
 *      F(beta) = 1/2 ||y - X beta||^2 + (lambda2/2) s^T K s,
 *      s = fold(beta),
 *
 *  with a user-supplied sparse PSD matrix K = Gamma^T Gamma (p x p) and the
 *  FOLDED dummy convention (dummy copy t couples as its base
 *  b(t) = (t - p) mod p — the T-Rex+GVS layered convention, so dummy/null
 *  exchangeability holds for arbitrary K):
 *
 *  - Selection by penalized correlation
 *    c_j = x_j^T r - lambda2 * (K s)_{b(j)} — the negative gradient of F at
 *    the current iterate (the same modified correlation as TIKNET / the
 *    TCDTIKNET FOLDED mode).
 *  - Active-set refit beta_A = (X_A^T X_A + lambda2 K-tilde_A)^{-1} X_A^T y
 *    with K-tilde_{tu} = K_{b(t), b(u)}, via the same penalized Cholesky
 *    rank-1 update as TIKNET (diagonal x_t^T x_t + lambda2 K_bb, cross
 *    terms x_t^T x_u + lambda2 K_{b(t), b(u)}).
 *  - Residual r = y - X_A beta_A; the coupling state s, v = K s is rebuilt
 *    from the refit coefficients once per step (O(nnz(K))).
 *  - Forward-only (OMP never drops variables); termination on the
 *    dummy-count contract like every t-solver.
 *
 *  Equivalence: the path coincides with plain T-OMP run on the augmented
 *  system X_aug = [X; sqrt(lambda2) Gamma], y_aug = [y; 0],
 *  D_aug = [D; sqrt(lambda2) Gamma tiled per layer] WITHOUT re-normalizing
 *  the augmented columns (both sides standardize the DATA blocks only) —
 *  the Theorem-1 pattern of TIENET_AUG, used as the cross-check anchor.
 *
 *  Degenerate cases: lambda2 == 0 reduces exactly to T-OMP. K = I ridges
 *  each variable jointly with its own dummy copies through the folded sum
 *  (ridge-OMP with folded coupling).
 */
// ===================================================================================

// std includes
#include <string>
#include <vector>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Cereal includes
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/base_class.hpp>

// Base class
#include <tsolvers/ad/linear_model/omp_based/tomp_solver.hpp>

// ===================================================================================

// Embedded into trex::tsolvers::ad::linear_model::omp_based namespace
namespace trex::tsolvers::ad::linear_model::omp_based {

// ===================================================================================

/**
 * @brief Terminating Tikhonov-K OMP (T-TIKOMP) solver with FOLDED dummy
 * coupling and generalized-ridge active-set refits.
 *
 *  References:
 *  -----------
 *   A. N. Tikhonov, "Solution of Incorrectly Formulated Problems and the Regularization Method," Soviet Math.
 *     Dokl., vol. 4, pp. 1035-1038, 1963 (English translation of Dokl. Akad. Nauk SSSR 151 (1963), no. 3,
 *     501-504).
 *   Y. Park, L. Reichel, G. Rodriguez, and X. Yu, "Parameter Determination for Tikhonov Regularization
 *     Problems in General Form," J. Comput. Appl. Math., vol. 343, pp. 12-25, 2018,
 *     doi: 10.1016/j.cam.2018.04.049.
 *   Y. C. Pati, R. Rezaiifar, and P. S. Krishnaprasad, "Orthogonal Matching Pursuit: Recursive Function
 *     Approximation with Applications to Wavelet Decomposition," in Proc. 27th Asilomar Conf. Signals, Syst.,
 *     Comput., 1993, pp. 40-44.
 *   F. Scheidt and M. Muma, "Terminating Greedy Pursuits for High-Dimensional False Discovery Rate Control,"
 *     in Proc. Eur. Signal Process. Conf. (EUSIPCO), 2026, pp. 2511-2515.
 */
class TIKOMP_Solver : public TOMP_Solver {
protected:
    // ==========================================================================
    // TIK specific state variables
    // ==========================================================================

    /** @brief Quadratic penalty weight lambda2 (>= 0). */
    double lambda2_{0.0};

    /** @brief Sparse PSD Tikhonov matrix K = Gamma^T Gamma (p x p) on the
     *  standardized scale. */
    Eigen::SparseMatrix<double> tikhonov_K_{};

    /** @brief Folded coefficient sums s in R^p (transient: rebuilt from the
     *  active refit after every step). */
    Eigen::VectorXd folded_s_{};

    /** @brief Coupling state v = K s (UNSCALED by lambda2; transient). */
    Eigen::VectorXd ks_{};

    /** @brief Dense gather workspace for one K column during the Cholesky
     *  rank-1 update (transient; lazily sized). */
    Eigen::VectorXd kcol_gather_{};

    /** @brief Base variable of unified column t (identity for reals,
     *  (t - p) mod p for dummy copies — FOLDED layered convention). */
    std::size_t foldedBase(std::size_t t) const noexcept {
        return (t < p_original_) ? t : (t - p_original_) % p_original_;
    }

public:
    // ==========================================================================
    // Constructor/Destructor
    // ==========================================================================

    /**
     * @brief Construct a new TIKOMP_Solver object.
     *
     * @param X Map to predictor matrix (n x p).
     * @param D Map to the dummy predictor matrix (n x L) with L a multiple
     *        of p for lambda2 > 0 (FOLDED layered convention).
     * @param y Map to the response vector (n x 1).
     * @param lambda2 Quadratic penalty weight (>= 0).
     * @param tikhonov_matrix Sparse PSD Tikhonov matrix K = Gamma^T Gamma
     *        (p x p) applied to the folded sums on the standardized scale.
     * @param normalize Whether to normalize the predictors.
     * @param intercept Whether to include an intercept term.
     * @param verbose Whether to enable verbose output.
     * @param scaling_mode Column-scaling convention (L2 or z-score).
     */
    TIKOMP_Solver(Eigen::Map<Eigen::MatrixXd>& X,
                  Eigen::Map<Eigen::MatrixXd>& D,
                  Eigen::Map<Eigen::VectorXd>& y,
                  const Eigen::SparseMatrix<double>& tikhonov_matrix,
                  const SolverOptions& options = SolverOptions{});

    /** @brief Default constructor for serialization. */
    TIKOMP_Solver() : TOMP_Solver(SolverTypeOMPBased::TIKOMP) {}

    /** @brief Virtual destructor. */
    ~TIKOMP_Solver() override = default;

    /** @brief Deleted copy operations. */
    TIKOMP_Solver(const TIKOMP_Solver&) = delete;

    /** @brief Deleted copy assignment operator. */
    TIKOMP_Solver& operator=(const TIKOMP_Solver&) = delete;

    /** @brief Defaulted move constructor. */
    TIKOMP_Solver(TIKOMP_Solver&&) = default;

    /** @brief Deleted move assignment operator. */
    TIKOMP_Solver& operator=(TIKOMP_Solver&&) = delete;

    // ==========================================================================
    // Main algorithm
    // ==========================================================================

    /**
     * @brief Execute the T-TIKOMP path with early stopping based on the
     * dummy-count contract. Mirrors T-OMP with the active-set bound raised
     * to max_actives_ (the penalized Gram can exceed rank n).
     *
     * @param T_stop Number of dummies for early stopping threshold
     *        (0 = full path).
     * @param early_stop If true, stop when count_active_dummies_ >= T_stop.
     */
    void executeStep(std::size_t T_stop = 0, bool early_stop = true) override;

    // ==========================================================================
    // Output & Accessors
    // ==========================================================================

    /** @brief Quadratic penalty weight used. */
    double getLambda2() const noexcept { return lambda2_; }

    /** @brief The sparse Tikhonov matrix K. */
    const Eigen::SparseMatrix<double>& getTikhonovK() const noexcept {
        return tikhonov_K_;
    }

    /** @brief TOMP_Solver::memoryBytes() plus the copy of K, the folded sums,
     *  the coupling state and the gather workspace. */
    std::size_t memoryBytes() const noexcept override;

    // ==========================================================================
    // Serialization
    // ==========================================================================

    /** @brief Friend class for Cereal serialization. */
    friend class cereal::access;

    /**
     * @brief Serialize all internal state except input matrices. K passes
     * through triplets; the coupling state s, v is transient (rebuilt from
     * the refit coefficients on the next step).
     *
     * @tparam Archive Cereal archive type.
     * @param archive Output/input archive object.
     */
    template <class Archive>
    void serialize(Archive& archive) {
        // Sparse K passes through triplets (same convention as TCDLASSO).
        std::vector<int> K_rows, K_cols;
        std::vector<double> K_vals;
        for (int j = 0; j < tikhonov_K_.outerSize(); ++j) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(tikhonov_K_, j);
                 it; ++it) {
                K_rows.push_back(static_cast<int>(it.row()));
                K_cols.push_back(j);
                K_vals.push_back(it.value());
            }
        }
        std::size_t K_dim = static_cast<std::size_t>(tikhonov_K_.rows());

        archive(
            cereal::base_class<TOMP_Solver>(this),
            CEREAL_NVP(lambda2_),
            cereal::make_nvp("K_dim", K_dim),
            cereal::make_nvp("K_rows", K_rows),
            cereal::make_nvp("K_cols", K_cols),
            cereal::make_nvp("K_vals", K_vals)
        );

        // Rebuild sparse K from the (de)serialized triplets (no-op on save).
        if (K_dim > 0) {
            std::vector<Eigen::Triplet<double>> trip;
            trip.reserve(K_vals.size());
            for (std::size_t t = 0; t < K_vals.size(); ++t) {
                trip.emplace_back(K_rows[t], K_cols[t], K_vals[t]);
            }
            tikhonov_K_.resize(static_cast<Eigen::Index>(K_dim),
                               static_cast<Eigen::Index>(K_dim));
            tikhonov_K_.setFromTriplets(trip.begin(), trip.end());
            tikhonov_K_.makeCompressed();
        }
    }

    /**
     * @brief Save the current TIKOMP solver/model state to file (binary
     * serialization).
     *
     * @param filename Output file path.
     */
    void save(const std::string& filename) const override;

    /**
     * @brief Load a TIKOMP solver state from file (binary serialization).
     *
     * @param filename Path to input file.
     * @param X Map to feature matrix.
     * @param D Map to dummy matrix.
     *
     * @return Deserialized TIKOMP_Solver object connected to X and D.
     */
    static TIKOMP_Solver load(const std::string& filename,
                              Eigen::Map<Eigen::MatrixXd>& X,
                              Eigen::Map<Eigen::MatrixXd>& D);

protected:

    // ==========================================================================
    // Internal TIK-OMP math helpers
    // ==========================================================================

    /**
     * @brief Rebuild the coupling state s (folded sums of the current refit
     * coefficients) and v = K s, then refresh the penalized correlations
     * c_j = x_j^T r - lambda2 * (K s)_{b(j)} over all inactives.
     */
    void updateCorrelations() override;

    /**
     * @brief Active set update with the penalized Cholesky rank-1 update
     * (diagonal + lambda2 K_bb, cross terms + lambda2 K_{b(t), b(u)}).
     *
     * @param new_vars Candidate indices for addition.
     *
     * @return Vector of actions for this step: positive index = accepted,
     *         negative = dropped.
     */
    std::vector<int> updateActiveSet(
        const std::vector<std::size_t>& new_vars) override;

    /**
     * @brief Penalized Cholesky rank-1 update for
     * G_A = X_A^T X_A + lambda2 K-tilde_A (needs the entering variable's
     * index for its base, hence not an override of the column-only
     * updateR()).
     *
     * @param j_new Entering unified column index.
     * @param eps   Numerical tolerance for collinearity detection.
     * @return Updated (k+1 x k+1) upper-triangular Cholesky factor.
     */
    Eigen::MatrixXd updateRTIK(std::size_t j_new, double eps);

    // updateBetaPath() and updateResiduals() are inherited unchanged: the
    // penalty enters the refit exclusively through the Cholesky factor R_
    // (normal equations G_A beta = X_A^T y), and the residual is the plain
    // observation-space r = y - X_A beta_A.
};

// =============================================================================
} /* End of namespace trex::tsolvers::ad::linear_model::omp_based */

// =============================================================================
// Cereal polymorphic registration
// =============================================================================
#include <cereal/types/polymorphic.hpp>
CEREAL_REGISTER_TYPE(trex::tsolvers::ad::linear_model::omp_based::TIKOMP_Solver)

// =============================================================================
#endif /* TSOLVERS_LINEAR_MODEL_OMP_BASED_TIKOMP_SOLVER_HPP */
