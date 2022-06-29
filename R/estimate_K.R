

#' Estimation of Intrinsic Mode Function (IMF) Number in Variational Mode Decomposition
#'
#' @param signal_1d a numeric vector specifying the 1-dimensional input signal
#' @param cor_thresh a numeric value specifying the minimum (positive or negative) correlation coefficient threshold where decomposition will be stopped (a value between 0.0 and 1.0)
#' @param default_vmd_params a list of parameters consisting of the (remaining) Variational Mode Decomposition default parameters (except for 'data' and 'K')
#' @param min_K a numeric value specifying the minimum value of the K (modes) parameter (from which decomposition starts)
#' @param seed a numeric value specifying the seed (for reproducibility purposes)
#' @param verbose a boolean. If TRUE then information will be printed in the console
#'
#' @return
#'
#' a numeric value specifying the optimal K parameter
#'
#' @importFrom glue glue
#' @importFrom data.table data.table
#'
#' @details
#'
#' Correlation Coefficient Method:
#' * Correlation coefficient (CC) between the mode components and the original signal will be obtained. Decomposition will be stopped when the minimum correlation coefficient is less than the given threshold, and then the value of K will be determined
#' @md
#'
#' @references
#'
#' https://doi.org/10.1155/2020/8304903
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' require(VMDecomp)
#' data(arrhythmia)
#'
#' default_vmd_params = list(alpha = 2000,
#'                           tau = 0,
#'                           DC = FALSE,
#'                           init = 1,
#'                           tol = 1e-6)
#'
#' res_k = estimate_k_modes(signal_1d = arrhythmia[['MLII']],
#'                          cor_thresh = 0.1,
#'                          default_vmd_params = default_vmd_params,
#'                          min_K = 2,
#'                          seed = 1,
#'                          verbose = TRUE)
#' res_k
#' }

estimate_k_modes = function(signal_1d,
                            cor_thresh,
                            default_vmd_params,
                            min_K = 2,
                            seed = 1,
                            verbose = FALSE) {

  if (min_K < 2) stop("The 'min_K' parameter must be at least 2!", call. = F)                        # so that in case the 'min_K' = 2 gives a correlation coefficient less than 'cor_thresh' the optimal K will be 1 (min_K - 1)
  if (cor_thresh < 0.0) stop("The 'cor_thresh' parameter must be between 0.0 and 1.0!", call. = F)
  if (verbose) t_start = proc.time()
  k_iter = min_K

  while (TRUE) {

    cp_params = default_vmd_params
    cp_params = append(cp_params, list(data = signal_1d, K = k_iter))
    if (verbose) {
      cat(glue::glue("VMD based on a K of '{k_iter}' will be computed ..."), '\n')
    }

    set.seed(seed)
    vmd_iter = do.call(VMDecomp::vmd, cp_params)
    imfs_iter = data.table::data.table(vmd_iter$u)
    imfs_cors = as.vector(apply(imfs_iter, 2, function(x) abs(stats::cor(signal_1d, x))))         # use abs() to account for negative correlation

    if (any(imfs_cors <= cor_thresh)) {
      if (verbose) {
        message(glue::glue("Optimal K parameter: '{k_iter - 1}'  Pre-specified correlation coefficient threshold: '{cor_thresh}'"))
      }
      break
    }

    k_iter = k_iter + 1
  }

  if (verbose) compute_elapsed_time(time_start = t_start)
  return(k_iter - 1)
}

