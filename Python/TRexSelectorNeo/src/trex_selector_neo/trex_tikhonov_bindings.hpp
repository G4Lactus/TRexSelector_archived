// =====================================================================================
// trex_tikhonov_bindings.hpp - Pybind11 bindings for TRexTikhonovSelector
// =====================================================================================
#ifndef TREX_PYTHON_TREX_TIKHONOV_BINDINGS_HPP
#define TREX_PYTHON_TREX_TIKHONOV_BINDINGS_HPP
// =====================================================================================
/**
 * @file trex_tikhonov_bindings.hpp
 * @brief Pybind11 bindings for TRexTikhonovSelector (general-Tikhonov T-Rex).
 *
 * @details Mirrors the R surface (rcpp_trex_tikhonov.cpp): the control struct
 *          with its nested base control merged from a separate `trex_control`
 *          argument, the lambda_2 getters (used value, CV diagnostics, pilot
 *          record) and the two K builders. The penalty matrix K crosses the
 *          boundary as a dense p x p array and is converted with
 *          sparseView() (no SciPy dependency), the same convention as the
 *          Tikhonov T-solvers.
 */
// =====================================================================================

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// trex includes
#include <ml_methods/model_selection/lambda2_cv_diagnostics.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/trex_tikhonov.hpp>
#include <trex_selector_methods/ad/trex_tikhonov/tik_lambda2_pick.hpp>
#include "trex_core_bindings.hpp"

// =====================================================================================

namespace py = pybind11;
namespace tik    = trex::trex_selector_methods::ad::trex_tikhonov;
namespace tik_sd = trex::trex_selector_methods::ad::utils::solver_dispatch;
using trex::trex_selector_methods::ad::trex_core::TRexControlParameter;

// =====================================================================================

/**
 * @brief Python wrapper class for TRexTikhonovSelector leveraging zero-copy references.
 */
class PyTRexTikhonovSelector : public PyTRexSelector {
public:
    /**
     * @brief Zero-copy constructor for the PyTRexTikhonovSelector wrapper.
     *
     * @param X Reference to the design matrix (n × p).
     * @param y Reference to the response vector (n × 1).
     * @param tFDR Target False Discovery Rate threshold.
     * @param tik_control Tikhonov-specific control (K or penalty structure, lambda_2 rule, dummies).
     * @param trex_control General configuration parameters for the T-Rex algorithm.
     * @param seed Random seed for reproducible computations.
     * @param verbose If true, prompts standard out messaging during select().
     */
    PyTRexTikhonovSelector(
        Eigen::Ref<Eigen::MatrixXd> X,
        Eigen::Ref<Eigen::VectorXd> y,
        double tFDR,
        const tik::TRexTikhonovControlParameter& tik_control,
        const TRexControlParameter& trex_control,
        int seed,
        bool verbose
    ) {
        this->X_map_ = std::make_unique<Eigen::Map<Eigen::MatrixXd>>(X.data(),
                                                                  X.rows(),
                                                                  X.cols());
        this->y_map_ = std::make_unique<Eigen::Map<Eigen::VectorXd>>(y.data(),
                                                                  y.size());

        // The C++ constructor takes a single control struct that nests the
        // base algorithmic parameters as `trex_ctrl`. Merge the separately-
        // supplied `trex_control`, then keep only a Tikhonov-family solver
        // (TCDTIKNET / TIKNET, the general-K pair): the caller-facing
        // default (TLARS) is neutral and becomes TCDTIKNET, mirroring the R
        // glue (the Python wrapper warns first on an explicit mismatch).
        tik::TRexTikhonovControlParameter tik_ctrl = tik_control;
        tik_ctrl.trex_ctrl = trex_control;
        const auto st = tik_ctrl.trex_ctrl.solver_type;
        const bool tik_family =
            st == tik_sd::SolverTypeForTRex::TCDTIKNET ||
            st == tik_sd::SolverTypeForTRex::TIKNET;
        if (!tik_family) {
            tik_ctrl.trex_ctrl.solver_type = tik_sd::SolverTypeForTRex::TCDTIKNET;
        }

        this->selector_ = std::make_unique<tik::TRexTikhonovSelector>(
            *(this->X_map_), *(this->y_map_), tFDR, tik_ctrl, seed, verbose
        );
    }

    /** @brief lambda_2 used by the last select() (solver scale); -1 before select(). */
    double getLambda2Used() const {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->getLambda2Used();
    }

    /** @brief The lambda_2 tuner health record of the last select(). */
    const trex::ml_methods::model_selection::lambda2_cv_diagnostics&
    getLambda2CVDiagnostics() const {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->getLambda2CVDiagnostics();
    }

    /** @brief The lambda_2 pilot record of the last select() (PILOT_* rules). */
    const tik::Lambda2PilotResult& getLambda2PilotResult() const {
        return static_cast<tik::TRexTikhonovSelector*>(this->selector_.get())
            ->getLambda2PilotResult();
    }
};


/**
 * @brief Binds the TRexTikhonovSelector and related configurations.
 *
 * @param m The Python module to bind against (the AD selector submodule;
 *          bind_trex_gvs must have registered Lambda2CVDiagnostics before).
 */
inline void bind_trex_tikhonov(py::module& m) {
    py::enum_<tik::TikLambda2Method>(m, "TikLambda2Method",
        "Rule for resolving lambda_2 when tik_control.lambda_2 < 0.")
        .value("CV_1SE_IEN_CCD", tik::TikLambda2Method::CV_1SE_IEN_CCD, "IEN-geometry coordinate-descent CV (lambda_1 profiled), 1-SE rule (default).")
        .value("CV_MIN_IEN_CCD", tik::TikLambda2Method::CV_MIN_IEN_CCD, "IEN-geometry coordinate-descent CV, minimum rule.")
        .value("CV_1SE_TIK_SVD", tik::TikLambda2Method::CV_1SE_TIK_SVD, "Generalized-Tikhonov SVD ridge CV, 1-SE rule.")
        .value("CV_MIN_TIK_SVD", tik::TikLambda2Method::CV_MIN_TIK_SVD, "Generalized-Tikhonov SVD ridge CV, minimum rule.")
        .value("PILOT_MAX",      tik::TikLambda2Method::PILOT_MAX,      "Power rule read off a few random experiments: the grid value with the largest selection count.")
        .value("PILOT_KNEE",     tik::TikLambda2Method::PILOT_KNEE,     "Power rule: the smallest grid value within pilot_knee_tolerance of the maximum count.")
        .export_values();

    py::enum_<tik::PilotStatistic>(m, "PilotStatistic",
        "Selection count the PILOT_* rules read per grid value.")
        .value("MAJORITY",       tik::PilotStatistic::MAJORITY,       "Variables selected in more than half of the pilot experiments.")
        .value("CALIBRATED",     tik::PilotStatistic::CALIBRATED,     "Calibrated count (soft curve).")
        .value("GROUP_MAJORITY", tik::PilotStatistic::GROUP_MAJORITY, "Groups whose selected fraction exceeds pilot_group_fraction.")
        .export_values();

    py::enum_<tik::TikPenaltyStructure>(m, "TikPenaltyStructure",
        "Where the penalty matrix comes from.")
        .value("SPARSE_GENERAL",  tik::TikPenaltyStructure::SPARSE_GENERAL,  "tik_control.tikhonov_K is the p x p PSD penalty matrix (dense array, converted to sparse).")
        .value("GROUP_CENTERING", tik::TikPenaltyStructure::GROUP_CENTERING, "K = I - J/p_m per group, built from prior_groups (or the HAC clusters).")
        .export_values();

    py::class_<tik::Lambda2PickInfo>(m, "Lambda2PickInfo",
        "How the pilot rule picked its grid value.")
        .def(py::init<>())
        .def_readonly("index",        &tik::Lambda2PickInfo::index,        "0-based index of the pick on the pilot grid.")
        .def_readonly("max_selected", &tik::Lambda2PickInfo::max_selected, "Largest selection count on the grid.")
        .def_readonly("n_ties",       &tik::Lambda2PickInfo::n_ties,       "Grid values tied at the maximum.");

    py::class_<tik::Lambda2PilotResult>(m, "Lambda2PilotResult",
        "Record of a lambda_2 pilot run (PILOT_MAX / PILOT_KNEE): the grid, the per-value counts and curves, the pick, and the pilot's size.")
        .def(py::init<>())
        .def_readonly("ran",               &tik::Lambda2PilotResult::ran,               "A pilot ran (False: lambda_2 came from the caller or a CV tuner).")
        .def_readonly("grid",              &tik::Lambda2PilotResult::grid,              "The lambda_2 grid probed.")
        .def_readonly("counts",            &tik::Lambda2PilotResult::counts,            "Selection count per grid value under pilot_statistic.")
        .def_readonly("curve_majority",    &tik::Lambda2PilotResult::curve_majority,    "Per grid value: majority count per pilot T-step.")
        .def_readonly("curve_soft",        &tik::Lambda2PilotResult::curve_soft,        "Per grid value: calibrated (soft) count per pilot T-step.")
        .def_readonly("calibrated",        &tik::Lambda2PilotResult::calibrated,        "Calibrated count per grid value.")
        .def_readonly("group_majority",    &tik::Lambda2PilotResult::group_majority,    "Group-majority count per grid value.")
        .def_readonly("groups_qualifying", &tik::Lambda2PilotResult::groups_qualifying, "Groups qualifying per grid value.")
        .def_readonly("probe_seconds",     &tik::Lambda2PilotResult::probe_seconds,     "Wall-clock seconds per grid value.")
        .def_readonly("pick",              &tik::Lambda2PilotResult::pick,              "Lambda2PickInfo of the pick.")
        .def_readonly("index",             &tik::Lambda2PilotResult::index,             "0-based index of the picked grid value.")
        .def_readonly("lambda2",           &tik::Lambda2PilotResult::lambda2,           "The picked lambda_2.")
        .def_readonly("at_lower_end",      &tik::Lambda2PilotResult::at_lower_end,      "The pick is the first grid value.")
        .def_readonly("at_upper_end",      &tik::Lambda2PilotResult::at_upper_end,      "The pick is the last grid value.")
        .def_readonly("K_pilot",           &tik::Lambda2PilotResult::K_pilot,           "Pilot experiments per grid value.")
        .def_readonly("T_pilot",           &tik::Lambda2PilotResult::T_pilot,           "Pilot T-steps.")
        .def_readonly("L_0",               &tik::Lambda2PilotResult::L_0,               "Pilot dummy multiplier.")
        .def_readonly("dummy_tag",         &tik::Lambda2PilotResult::dummy_tag,         "Dummy stream tag the pilot drew under.")
        .def_readonly("seconds",           &tik::Lambda2PilotResult::seconds,           "Wall-clock seconds of the whole pilot.");

    py::class_<tik::TRexTikhonovControlParameter>(m, "TRexTikhonovControlParameter",
        "Control parameters for TRexTikhonovSelector: the penalty (a dense p x p PSD matrix tikhonov_K, or GROUP_CENTERING built from groups), the lambda_2 rule, and the cluster-aware dummies.")
        .def(py::init<>())
        .def_readwrite("penalty_structure", &tik::TRexTikhonovControlParameter::penalty_structure, "SPARSE_GENERAL (tikhonov_K given) or GROUP_CENTERING (K from the groups).")
        .def_property("tikhonov_K",
            [](const tik::TRexTikhonovControlParameter& c) {
                return Eigen::MatrixXd(c.tikhonov_K);
            },
            [](tik::TRexTikhonovControlParameter& c, const Eigen::MatrixXd& K) {
                c.tikhonov_K = K.sparseView();
            },
            "The p x p PSD penalty matrix K as a dense array (converted to sparse internally; read back dense). Required under SPARSE_GENERAL.")
        .def_readwrite("lambda_2",        &tik::TRexTikhonovControlParameter::lambda_2,        "Tikhonov penalty (solver scale): < 0 = resolve by lambda2_method, 0 = no ridge, > 0 = fixed.")
        .def_readwrite("lambda2_method",  &tik::TRexTikhonovControlParameter::lambda2_method,  "Rule for resolving lambda_2 when lambda_2 < 0 (default CV_1SE_IEN_CCD).")
        .def_readwrite("cv_n_folds",      &tik::TRexTikhonovControlParameter::cv_n_folds,      "Folds of the CV tuners.")
        .def_readwrite("cv_n_lambda",     &tik::TRexTikhonovControlParameter::cv_n_lambda,     "SVD tuners: points on the lambda_2 grid.")
        .def_readwrite("cv_n_lambda1",    &tik::TRexTikhonovControlParameter::cv_n_lambda1,    "IEN-geometry CD tuners: points on the profiled lambda_1 path per lambda_2 value.")
        .def_readwrite("cv_dfmax_ratio",  &tik::TRexTikhonovControlParameter::cv_dfmax_ratio,  "IEN-geometry CD tuners: support cap as a fraction of the smallest fold's training size (default 0.8; <= 0 disables the cap; larger values allow larger supports at a higher tuner cost; must be finite).")
        .def_readwrite("cv_n_lambda2",    &tik::TRexTikhonovControlParameter::cv_n_lambda2,    "IEN-geometry CD tuners: points on the lambda_2 grid.")
        .def_readwrite("cv_seed",         &tik::TRexTikhonovControlParameter::cv_seed,         "Seed for the CV fold permutation (-1 = derive from the T-Rex seed).")
        .def_readwrite("pilot_grid",      &tik::TRexTikhonovControlParameter::pilot_grid,      "lambda_2 grid of the PILOT_* rules.")
        .def_readwrite("pilot_K",         &tik::TRexTikhonovControlParameter::pilot_K,         "Pilot experiments per grid value.")
        .def_readwrite("pilot_T",         &tik::TRexTikhonovControlParameter::pilot_T,         "Pilot T-steps.")
        .def_readwrite("pilot_L_multiplier", &tik::TRexTikhonovControlParameter::pilot_L_multiplier, "Pilot dummy multiplier (in [1, trex_control.max_dummy_multiplier]).")
        .def_readwrite("pilot_knee_tolerance", &tik::TRexTikhonovControlParameter::pilot_knee_tolerance, "PILOT_KNEE: count tolerance below the maximum.")
        .def_readwrite("pilot_statistic", &tik::TRexTikhonovControlParameter::pilot_statistic, "Selection count the pilot reads (MAJORITY, CALIBRATED, GROUP_MAJORITY).")
        .def_readwrite("pilot_group_fraction", &tik::TRexTikhonovControlParameter::pilot_group_fraction, "GROUP_MAJORITY: fraction of a group that must be selected.")
        .def_readwrite("pilot_shared_dummies", &tik::TRexTikhonovControlParameter::pilot_shared_dummies, "Pilot draws its dummies under the main run's first block tag.")
        .def_readwrite("fold_dummy_coupling", &tik::TRexTikhonovControlParameter::fold_dummy_coupling, "FOLDED dummy coupling: dummies take the K block of their real counterpart (default True).")
        .def_readwrite("cluster_dummies", &tik::TRexTikhonovControlParameter::cluster_dummies, "Cluster-aware MVN dummies (default True); False = i.i.d. dummies (the K = I anchor).")
        .def_readwrite("prior_groups",    &tik::TRexTikhonovControlParameter::prior_groups,    "Group assignment (0-based, length p; empty = HAC clusters from corr_max / hc_linkage).")
        .def_readwrite("corr_max",        &tik::TRexTikhonovControlParameter::corr_max,        "HAC cut: maximum within-cluster correlation when prior_groups is empty.")
        .def_readwrite("hc_linkage",      &tik::TRexTikhonovControlParameter::hc_linkage,      "HAC linkage (ml_methods.clustering.LinkageMethod).")
        .def_readwrite("trex_ctrl",       &tik::TRexTikhonovControlParameter::trex_ctrl,       "Nested base T-Rex control (overridden by the separate trex_control argument of TRexTikhonovSelector).");

    m.def("gamma_to_k",
          [](const Eigen::MatrixXd& gamma) {
              return Eigen::MatrixXd(tik::TRexTikhonovSelector::gammaToK(gamma.sparseView()));
          },
          py::arg("gamma"),
          "Build the penalty matrix K = Gamma^T Gamma from a (q x p) Tikhonov operator Gamma (dense in, dense out).");

    m.def("build_group_centering_k",
          [](const std::vector<Eigen::Index>& prior_groups) {
              return Eigen::MatrixXd(tik::TRexTikhonovSelector::buildGroupCenteringK(prior_groups));
          },
          py::arg("prior_groups"),
          "Build the GROUP_CENTERING penalty K = blockdiag(I - J/p_m) from a 0-based group assignment (dense out).");

    py::class_<PyTRexTikhonovSelector, PyTRexSelector>(m, "TRexTikhonovSelector",
        "T-Rex selector with a general Tikhonov penalty (TCDTIKNET / TIKNET) and "
        "cluster-aware dummies.\n"
        "\n"
        "References\n"
        "----------\n"
        "A. N. Tikhonov, \"Solution of Incorrectly Formulated Problems and the Regularization\n"
        "  Method,\" Soviet Math. Dokl., vol. 4, pp. 1035-1038, 1963 (English translation of Dokl.\n"
        "  Akad. Nauk SSSR 151 (1963), no. 3, 501-504).\n"
        "Y. Park, L. Reichel, G. Rodriguez, and X. Yu, \"Parameter Determination for Tikhonov\n"
        "  Regularization Problems in General Form,\" J. Comput. Appl. Math., vol. 343, pp. 12-25,\n"
        "  2018, doi: 10.1016/j.cam.2018.04.049.\n"
        "A. J. Smola and R. Kondor, \"Kernels and Regularization on Graphs,\" in Learning Theory\n"
        "  and Kernel Machines (COLT/Kernel 2003), Lecture Notes in Artificial Intelligence\n"
        "  vol. 2777, 2003, pp. 144-158.\n"
        "F. R. K. Chung, Spectral Graph Theory. CBMS Regional Conference Series in Mathematics,\n"
        "  no. 92. American Mathematical Society, 1997.\n"
        "A. E. Hoerl and R. W. Kennard, \"Ridge Regression: Biased Estimation for Nonorthogonal\n"
        "  Problems,\" Technometrics, vol. 42, no. 1, pp. 80-86, 2000 (Reprint of the 1970 article\n"
        "  in the Special 40th Anniversary Issue).\n"
        "J. Machkour, M. Muma, and D. P. Palomar, \"The Terminating-Random Experiments Selector:\n"
        "  Fast High-Dimensional Variable Selection with False Discovery Rate Control,\" Signal\n"
        "  Process., vol. 231, Art. no. 109894, 2025, doi: 10.1016/j.sigpro.2025.109894.\n")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::VectorXd>, double,
                      const tik::TRexTikhonovControlParameter&, const TRexControlParameter&, int, bool>(),
             py::arg("X"), py::arg("y"), py::arg("tFDR") = 0.1,
             py::arg("tik_control") = tik::TRexTikhonovControlParameter(),
             py::arg("trex_control") = TRexControlParameter(),
             py::arg("seed") = -1, py::arg("verbose") = true,
             py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
             "Construct the Tikhonov selector. X and y are accessed zero-copy via Eigen::Map.")
        .def("getLambda2Used", &PyTRexTikhonovSelector::getLambda2Used, "lambda_2 used by the last select() (solver scale); -1 before select().")
        .def("getLambda2CVDiagnostics", &PyTRexTikhonovSelector::getLambda2CVDiagnostics, "Return the lambda_2 tuner health record of the last select().")
        .def("getLambda2PilotResult", &PyTRexTikhonovSelector::getLambda2PilotResult, "Return the lambda_2 pilot record of the last select().");
}
// =====================================================================================
#endif /* End of TREX_PYTHON_TREX_TIKHONOV_BINDINGS_HPP */
