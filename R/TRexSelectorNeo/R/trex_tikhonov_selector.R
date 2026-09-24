#' @name TRexTikhonovSelector
#'
#' @title T-Rex Selector with a general Tikhonov penalty
#'
#' @description FDR-controlled variable selection with a Tikhonov penalty
#'   matrix \code{K = t(Gamma) \%*\% Gamma} driving the Tikhonov-net solvers
#'   (TCDTIKNET by default, or TIKNET / TIKOMP). \code{K = I} collapses to
#'   the plain elastic net; a group-mean K reproduces the GVS-IEN penalty;
#'   banded / Laplacian K encode smoothness-informed selection;
#'   \code{penalty_structure = "GROUP_CENTERING"} builds the centering
#'   penalty from the dummy clustering. lambda_2 comes from the caller, a
#'   CV tuner, or the pilot rule (see \code{\link{trex_tikhonov_control}}).
#'   Dummies are cluster-aware and FOLDED-coupled by default.
#'
#' @importFrom R6 R6Class
#'
#' @examples
#' \donttest{
#' data(Gauss_data)
#' X <- Gauss_data$X
#' y <- Gauss_data$y
#' ctrl <- trex_tikhonov_control(tikhonov_K = diag(ncol(X)), lambda_2 = 0.5,
#'                               cluster_dummies = FALSE,
#'                               fold_dummy_coupling = FALSE)
#' sel <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, verbose = FALSE,
#'                                 tik_control = ctrl)
#' sel$select()
#' sel$selected_indices
#' }
#' @family T-Rex selectors and controls
#' @export
TRexTikhonovSelector <- R6::R6Class("TRexTikhonovSelector",
  inherit = TRexSelector,
  public = list(
    #' @description Create a new TRexTikhonovSelector object.
    #'
    #' @param X Feature matrix or MemoryMappedMatrix.
    #' @param y Response vector.
    #' @param tFDR Target FDR level (default: 0.1).
    #' @param seed Random seed (default: -1).
    #' @param verbose Whether to print progress (default: TRUE).
    #' @param tik_control A Tikhonov control list from
    #'   \code{\link{trex_tikhonov_control}()} (required — it carries the
    #'   penalty matrix K).
    #' @param control A control list from \code{\link{trex_control}()}
    #'   (default: \code{trex_control()}). The \code{solver} field selects
    #'   the Tikhonov-family solver: \code{"TCDTIKNET"} (default; the
    #'   neutral \code{trex_control()} default \code{"TLARS"} maps to it),
    #'   \code{"TIKNET"} (native pathwise), or \code{"TIKOMP"}; any other
    #'   value warns and falls back to \code{"TCDTIKNET"}.
    initialize = function(X, y,
                          tFDR = 0.1,
                          seed = -1,
                          verbose = TRUE,
                          tik_control,
                          control = trex_control()) {

      if (missing(tik_control)) {
        stop("tik_control is required (build it with trex_tikhonov_control(); ",
             "it carries the penalty matrix K or the penalty structure).")
      }

      # Only the Tikhonov family has a general-K mode. "TLARS" is the
      # trex_control() default and cannot be distinguished from an explicit
      # choice, so it never warns; the C++ glue maps it to TCDTIKNET.
      tik_family <- c("TCDTIKNET", "TIKNET", "TIKOMP")
      if (!(control$solver %in% c("TLARS", tik_family))) {
        warning("control$solver = \"", control$solver,
                "\" is not a Tikhonov-family solver; TRexTikhonovSelector ",
                "accepts \"TCDTIKNET\", \"TIKNET\", or \"TIKOMP\". ",
                "Using \"TCDTIKNET\".")
        control$solver <- "TLARS"  # the C++ glue defaults it to TCDTIKNET
      }

      if (inherits(X, "MemoryMappedMatrix")) {
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_tikhonov_mmap_create(
          get_ptr(X), y, tFDR, tik_control, control, seed, verbose
        )
      } else {
        X <- .as_double_matrix(X)
        y <- .as_double_vector(y)
        private$refs <- list(X = X, y = y)
        private$ptr <- trex_tikhonov_create(
          X, y, tFDR, tik_control, control, seed, verbose
        )
      }
    },

    #' @description Run the Tikhonov T-Rex algorithm.
    select = function() {
      trex_tikhonov_select(private$ptr)
      invisible(self)
    }
  ),
  active = list(
    #' @field lambda2_used The lambda_2 actually used by the last
    #'   \code{select()} run (user-supplied or CV-resolved; solver scale).
    #'   \code{-1} before \code{select()}.
    lambda2_used = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_used(private$ptr)
    },

    #' @field lambda2_cv_diagnostics Health record of the lambda_2 tuner of
    #'   the last \code{select()}: a list with \code{ran} (a tuner or the
    #'   pilot ran), \code{seconds}, \code{n_profile_at_edge} (lambda_2
    #'   values whose profiled lambda_1 optimum sat on a truncation edge; 0
    #'   is healthy), \code{support_cap} (-1 = disabled),
    #'   \code{low_dimensional}, and the boundary-pick guard fields
    #'   \code{grid_index} / \code{grid_size}, \code{at_grid_edge},
    #'   \code{profile_at_edge} and \code{boundary_pick} (see
    #'   \code{TRexGVSSelector}; the CV methods warn once on a boundary pick
    #'   and return it unchanged; the pilot rules leave them at their
    #'   defaults).
    lambda2_cv_diagnostics = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_cv_diagnostics(private$ptr)
    },

    #' @field lambda2_pilot Record of the lambda_2 pilot rule of the last
    #'   \code{select()} (\code{ran == FALSE} when no pilot ran): the
    #'   \code{grid}, the \code{counts} the pick was read off, the
    #'   per-grid-point curves over the pilot stopping times
    #'   (\code{curve_majority}, \code{curve_soft}; G x T matrices), the
    #'   \code{calibrated} and \code{group_majority} counts, the 1-based
    #'   \code{index} of the pick, the picked \code{lambda2},
    #'   \code{max_selected}, \code{n_ties}, the boundary flags
    #'   \code{at_lower_end} / \code{at_upper_end} (widen the grid before
    #'   trusting a boundary pick), the budget actually used (\code{K_pilot},
    #'   \code{T_pilot}, \code{L_0}) and the timings.
    lambda2_pilot = function(value) {
      if (!missing(value)) {
        stop("Field is read-only.")
      }
      trex_tikhonov_get_lambda2_pilot_result(private$ptr)
    }
  )
)


#' @name trex_gamma_to_k
#'
#' @title Build a Tikhonov penalty matrix K from an operator Gamma
#'
#' @description Computes \code{K = t(Gamma) \%*\% Gamma} (pruned and
#'   compressed) for use as \code{tikhonov_K} in
#'   \code{\link{trex_tikhonov_control}}. For example, a first-difference
#'   Gamma yields the smoothness-encoding graph Laplacian of a chain.
#'
#' @param gamma Regularization operator (m x p): a base \code{matrix} or a
#'   sparse \code{Matrix::dgCMatrix}.
#'
#' @return A sparse \code{Matrix::dgCMatrix} (p x p). Requires the Matrix
#'   package (Suggests) for the sparse return class.
#'
#' @examples
#' # First-difference operator over 5 variables -> chain Laplacian
#' p <- 5
#' gamma <- diff(diag(p))
#' K <- trex_gamma_to_k(gamma)
#'
#' @family T-Rex selectors and controls
#' @export
trex_gamma_to_k <- function(gamma) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("trex_gamma_to_k() returns a sparse Matrix::dgCMatrix and needs the ",
         "Matrix package installed.")
  }
  trex_tikhonov_gamma_to_k(gamma)
}
