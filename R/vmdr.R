
#' Variational Mode Decomposition (1- or 2-dimensional)
#'
#' @param data either a vector or a matrix (of type numeric or integer)
#' @param alpha a numeric value specifying the balancing parameter of the data-fidelity constraint
#' @param tau a numeric value specifying the time-step of the dual ascent ( pick 0 for noise-slack )
#' @param K a numeric value specifying the number of modes to be recovered
#' @param DC a boolean. If true the first mode is put and kept at DC (0-freq)
#' @param init a numeric value. This parameter differs depending on the input 'data' parameter (1-dimensional and 2-dimensional). See the details section for more information
#' @param tol a numeric value specifying the tolerance of convergence criterion (typically this parameter is around 1e-6 for the 1-dimensional and 1e-7 for the 2-dimensional data)
#' @param verbose a boolean. If TRUE then information will be printed in the console
#'
#' @return a list object of length three which includes the
#' * 'u' (collection of decomposed modes)
#' * 'u_hat' (spectra of the modes)
#' * 'omega' (estimated mode center-frequencies)
#' objects
#' @md
#'
#' @details
#' The 'init' parameter takes the following values for,
#' * 1-dimensional data:
#'     + 0 = all omegas start at 0
#'     + 1 = all omegas start uniformly distributed
#'     + 2 = all omegas initialized randomly
#' * 2-dimensional data:
#'     + 0 = all omegas start at 0
#'     + 1 = all omegas start initialized randomly
#' @md
#'
#' @export
#'
#' @references
#'
#' https://math.montana.edu/dzosso/code/
#'
#' @examples
#'
#' require(VMDecomp)
#'
#' #..............
#' # 1-dimensional
#' #..............
#'
#' N = 250
#'
#' set.seed(1)
#' rand_unif = runif(n = N, min = 0, max = 1.0)
#'
#' f_sig1 = 6 * rand_unif
#' f_sig2 = cos(x = 8 * pi * rand_unif)
#' f_sig3 = 0.5 * cos(x = 40 * pi * rand_unif)
#'
#' f_sig = f_sig1 + f_sig2 + f_sig3
#'
#' alpha = 2000
#' tau = 0
#' K = 3
#' DC = FALSE
#' init = 1
#' tol = 1e-6
#'
#' set.seed(2)
#' res_1d = vmd(data = f_sig,
#'              alpha = alpha,
#'              tau = tau,
#'              K = K,
#'              DC = DC,
#'              init = init,
#'              tol = tol,
#'              verbose = FALSE)
#'
#' #..............
#' # 2-dimensional
#' #..............
#'
#' rows_cols = 10
#'
#' set.seed(3)
#' data = matrix(runif(rows_cols^2), rows_cols, rows_cols)
#' alpha = 5000
#' tau = 0.25
#' K = 2
#' DC = TRUE
#' init = 1
#' tol = 1e-7
#'
#' set.seed(4)
#' res_2d = vmd(data = data,
#'              alpha = alpha,
#'              tau = tau,
#'              K = K,
#'              DC = DC,
#'              init = init,
#'              tol = tol,
#'              verbose = FALSE)

vmd = function(data,
               alpha,
               tau,
               K,
               DC,
               init,
               tol,
               verbose = FALSE) {

  if (!inherits(DC, 'logical')) stop("The 'DC' parameter must be either TRUE or FALSE!", call. = F)

  if (verbose) t_start = proc.time()
  if (inherits(data, 'matrix')) {
    if (!init %in% 0:1) stop("In case of the 2-dimensional data the 'init' parameter must be either 0 or 1", call. = F)
    vmd_obj = vmd_2d(signal = data,
                     alpha = alpha,
                     tau = tau,
                     K = K,
                     DC = DC,
                     init = init,
                     tol = tol,
                     verbose = verbose)
  }
  else if (is.vector(data, mode = 'numeric') || is.vector(data, mode = 'integer')) {
    if (!init %in% 0:2) stop("In case of the 1-dimensional data the 'init' parameter must be 0, 1 or 2", call. = F)
    if (verbose) cat("1-dimensional VMD starts ...\n")
    vmd_obj = vmd_1d(signal = data,
                     alpha = alpha,
                     tau = tau,
                     K = K,
                     DC = DC,
                     init = init,
                     tol = tol,
                     verbose = verbose)
  }
  else {
    stop("The input data must be either of type numeric matrix (2-dimensional) or numeric vector (1-dimensional)!", call. = F)
  }
  if (verbose) compute_elapsed_time(time_start = t_start)
  return(vmd_obj)
}

