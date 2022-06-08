/**
 * @file VMDheader.h
 * @author Lampros Sp. Mouselimis
 * @date June 2022
 *
 * Implementation of 1-Dimensional and 2-Dimensional 'Variational Mode Decomposition'
 * See COPYRIGHTS, CITATION and papers in the 'inst' folder for more information
 *
 */


#ifndef __VMDheader__
#define __VMDheader__

#ifdef _OPENMP
#include <omp.h>
#endif


namespace vmdR {

  class VarModeDecomp {

    public:

      VarModeDecomp() { }


      /**
       * "fftshift" (input 'vector' without using the 'dim' parameter)
       * Reference: https://github.com/gnu-octave/octave/blob/default/scripts/signal/fftshift.m
       *            /inst/Octave/fftshift.m
       */

      arma::cx_vec fftshift_vec(arma::cx_vec x, bool inverse = false) {
        arma::uword xl = x.n_elem;
        double xl_div = xl / 2.0;
        arma::uword xx;
        if (inverse) {
          xx = std::floor(xl_div);
        }
        else {
          xx = std::ceil(xl_div);
        }

        arma::uvec start_shift = arma::regspace<arma::uvec>(xx, 1, xl-1);
        arma::uvec end_shift = arma::regspace<arma::uvec>(0, 1, xx-1);
        arma::uvec shift = arma::join_cols(start_shift, end_shift);
        return x(shift);
      }


      /**
       * "fftshift" (input 'matrix' without using the 'dim' parameter)
       */

      arma::cx_mat fftshift_mat(arma::cx_mat x, bool inverse = false) {
        arma::uword nd = 2;                                              // number of dimensions in matrix by default is 2 (rows, columns)
        arma::uword nr_rows = x.n_rows;
        arma::uword nr_cols = x.n_cols;
        arma::Col<arma::uword> sz(2);
        sz(0) = nr_rows;
        sz(1) = nr_cols;

        arma::uword shift_nrows, shift_ncols;

        if (inverse) {
          shift_nrows = std::floor(nr_rows / 2.0);
          shift_ncols = std::floor(nr_cols / 2.0);
        }
        else {
          shift_nrows = std::ceil(nr_rows / 2.0);
          shift_ncols = std::ceil(nr_cols / 2.0);
        }
        arma::Col<arma::uword> sz2(2);
        sz2(0) = shift_nrows;
        sz2(1) = shift_ncols;

        arma::cx_mat y(x);
        arma::uvec idx_y_rows(nr_rows);
        arma::uvec idx_y_cols(nr_cols);

        for (arma::uword i = 0; i < nd; i++) {

          arma::uvec idx_start = arma::regspace<arma::uvec>(sz2(i), 1, sz(i)-1);
          arma::uvec idx_end = arma::regspace<arma::uvec>(0, 1, sz2(i)-1);
          arma::uvec shift = arma::join_cols(idx_start, idx_end);

          if (i == 0) {                              // modify rows (as index nd(0) will correspond to the rows)
            idx_y_rows = shift;
          }
          else {                                     // modify cols (as index nd(1) will correspond to the columns)
            idx_y_cols = shift;
          }
        }
        return y(idx_y_rows, idx_y_cols);
      }


      /**
       *  Convert a boolean to integer
       */

      int bool_to_int(bool input) {
        int val = (input) ? 1 : 0;
        return val;
      }


      /**
       * switch between 'init' cases (1-dimensional)
       */

      arma::mat switch_1D(arma::uword init,
                          arma::uword N,
                          arma::uword K,
                          double fs) {

        arma::mat omega_k = arma::zeros<arma::mat>(N, K);
        switch (init) {
          case 1: {
            for (arma::uword i = 0; i < K; i++) {
            omega_k(0,i) = (0.5 / K) * i;
          }
          } break;
          case 2: {
            arma::rowvec k_rand = arma::randu<arma::rowvec>(K);
            k_rand *= (std::log(0.5) - std::log(fs));              // element-wise multiplication
            k_rand += std::log(fs);                                // element-wise addition
            k_rand = arma::exp(k_rand);
            k_rand = arma::sort(k_rand);
            omega_k.row(0) = k_rand;
          } break;
          default: {
            omega_k.row(0) = arma::zeros<arma::rowvec>(K);
          }
        }
        return omega_k;
      }


      /**
       * switch between 'init' cases (2-dimensional)
       */

      arma::cube switch_2D(arma::uword init,
                           arma::uword N,
                           arma::uword K,
                           bool DC) {

        arma::cube omega_inner = arma::zeros<arma::cube>(N, 2, K);
        arma::uword maxK;

        switch (init) {
          case 0: {                          // Case 0: spread omegas radially
            if (DC) {                        // if DC, keep first mode at 0,0
            maxK = K-1;
          }
            else {
              maxK = K;
            }
            arma::vec range = arma::regspace<arma::vec>(0, 1, maxK);
            arma::uword dc_val = bool_to_int(DC);                      // if DC is TRUE then min is 1 otherwise it's 0
            range = range + dc_val;
            arma::uword k_min = arma::as_scalar(arma::min(range));
            arma::uword k_max = arma::as_scalar(arma::max(range));
            for (arma::uword k = k_min; k < k_max; k++) {
              omega_inner(0,0,k) = 0.25 * std::cos(arma::datum::pi * k / maxK);
              omega_inner(0,1,k) = 0.25 * std::sin(arma::datum::pi * k / maxK);
            }
          } break;
          case 1: {
            for (arma::uword k = 0; k < K; k++) {
            omega_inner(0,0,k) = arma::randu() - 1.0 / 2.0;
            omega_inner(0,1,k) = arma::randu() / 2.0;
          }
            if (DC) {                           // DC component (if expected)
              omega_inner(0,0,0) = 0.0;
              omega_inner(0,1,0) = 0.0;
            }
          } break;
        }
        return omega_inner;
      }


      /**
       * drop the row of a cube and create a new matrix
       * based on the columns and cube-slices
       */

      arma::cx_mat squeeze_row(arma::cx_cube x) {

        arma::uword N_slices = x.n_slices;
        arma::cx_mat y(N_slices, x.n_cols);

        for (arma::uword s = 0; s < N_slices; s++) {
          y.row(s) = x.slice(s).row(0);
        }
        return y.t();
      }


      /**
       * similar to the 'x' matrix of the matlab 'meshgrid' function  [ populating the matrices column-wise, is faster ]
       * Armadillo uses the "column-major ordering" which means read and write by column is faster than by row, see:
       * https://stackoverflow.com/a/53498100
       */

      arma::mat meshgrid_x(int rows, int cols) {

        arma::mat out(rows, cols, arma::fill::zeros);
        arma::colvec tmp_vec(out.col(0).n_elem);

        for (int i = 0; i < cols; i++) {
          out.col(i) = tmp_vec.fill(i);
        }
        return(out);
      }


      /**
       * similar to the 'y' matrix of the matlab 'meshgrid' function  [ populating the matrices column-wise, is faster ]
       * Armadillo uses the "column-major ordering" which means read and write by column is faster than by row, see:
       *
       // https://stackoverflow.com/a/53498100
       */

      arma::mat meshgrid_y(int rows, int cols) {

        arma::mat out(rows, cols, arma::fill::zeros);

        for (int i = 0; i < cols; i++) {
          out.col(i) = arma::regspace<arma::colvec>(0, rows - 1);
        }
        return(out);
      }


      /**
       * Variational Mode Decomposition (1-dimensional)
       */

      Rcpp::List VMD_1D(arma::vec& signal,
                        double alpha,
                        double tau,
                        arma::uword K,
                        bool DC,
                        arma::uword init,
                        double tol,
                        bool verbose = false) {
        
        bool flag_isfinite = signal.is_finite();
        if (!flag_isfinite) {
          Rcpp::stop("The input data includes missing values! You have to replace these missing values before using the 1-dimensional Variational Mode Decomposition function!");
        }

        // Period and sampling frequency of input signal
        arma::uword save_T = signal.n_elem;
        double fs = 1.0 / save_T;

        // extend the signal by mirroring
        arma::uword T = save_T;
        arma::vec f_mirror(signal.n_elem * 2);
        arma::uword T_div_two = T / 2.0;

        for (arma::uword i = 0; i < T_div_two; i++) {
          f_mirror(i) = signal(T_div_two - 1 - i);
        }

        arma::uword iter = 0;
        for (arma::uword j = T_div_two; j < (3 * T_div_two); j++) {
          f_mirror(j) = signal(iter);
          iter++;
        }

        arma::uword rev_iter = T - 1;
        for (arma::uword k = (3 * T_div_two); k < 2 * T; k++) {
          f_mirror(k) = signal(rev_iter);
          rev_iter--;
        }

        arma::vec f(f_mirror);

        // Time Domain 0 to T (of mirrored signal)
        T = f.n_elem;
        arma::rowvec t = arma::regspace<arma::rowvec>(1, 1, T);
        t = t / T;

        // Spectral Domain discretization
        arma::uword len_t = t.n_elem;
        arma::rowvec freqs(len_t);
        freqs = t - 0.5 - (1.0 / T);

        // Maximum number of iterations (if not converged yet, then it won't anyway)
        arma::uword N = 500;

        // For future generalizations: individual alpha for each mode
        arma::vec rep_ones = arma::ones<arma::vec>(K);
        arma::vec Alpha = alpha * rep_ones;

        // Construct and center f_hat
        arma::cx_vec f_hat = fft(f);
        f_hat = fftshift_vec(f_hat);

        arma::cx_vec f_hat_plus(f_hat);

        T_div_two = T / 2.0;                             // overwrite 'T_div_two'

        for (arma::uword s = 0; s < T_div_two; s++) {
          f_hat_plus(s) = 0.0;
        }

        // matrix keeping track of every iterant, could be discarded for mem
        arma::uword len_freqs = freqs.n_elem;
        arma::cx_cube u_hat_plus = arma::zeros<arma::cx_cube>(N, len_freqs, K);

        // update omega_plus
        arma::mat omega_plus = switch_1D(init, N, K, fs);

        // if DC mode imposed, set its omega to 0
        if (DC) {
          omega_plus(0,0) = 0.0;
        }

        // start with empty dual variables
        arma::cx_mat lambda_hat = arma::zeros<arma::cx_mat>(N, len_freqs);

        // other inits
        double uDiff = tol + arma::datum::eps;                               // update step
        arma::uword n = 0;                                                   // loop counter
        arma::cx_rowvec sum_uk = arma::zeros<arma::cx_rowvec>(len_freqs);    // accumulator

        if (verbose) {
          Rcpp::Rcout << "--------------------------------" << std::endl;
          Rcpp::Rcout << "The 1-dimensional VMD starts ..." << std::endl;
          Rcpp::Rcout << "--------------------------------" << std::endl;
        }

        // Main loop for iterative updates
        while (uDiff > tol && n < N - 1) {                 // not converged and below iterations limit

          // update first mode accumulator
          int k = 0;
          sum_uk = u_hat_plus.slice(K-1).row(n) + sum_uk - u_hat_plus.slice(0).row(n);

          // update spectrum of first mode through Wiener filter of residuals
          arma::cx_rowvec u_hat_iter_1st = (f_hat_plus.t() - sum_uk - lambda_hat.row(n) / 2.0);      // .t() transpose the Col-vec to match the dimensions of the rowvec
          arma::rowvec u_hat_iter_2nd = freqs - arma::as_scalar(omega_plus(n,k));
          u_hat_iter_2nd = arma::pow(u_hat_iter_2nd, 2);
          u_hat_iter_2nd = u_hat_iter_2nd * arma::as_scalar(Alpha(k));
          u_hat_iter_2nd = u_hat_iter_2nd + 1.0;

          u_hat_plus.slice(k).row(n+1) = u_hat_iter_1st / u_hat_iter_2nd;

          if (!DC) {
            arma::rowvec dc_iter = arma::abs(u_hat_plus.slice(k).row(n+1).subvec(T_div_two, T-1));
            double dc_scalar_iter = arma::accu(arma::pow(dc_iter, 2));
            double cross_pr_iter = arma::dot(freqs.subvec(T_div_two, T-1), arma::pow(dc_iter, 2));
            omega_plus(n+1,k) = cross_pr_iter / dc_scalar_iter;
          }

          // update of any other mode
          for (arma::uword k = 1; k < K; k++) {

            // accumulator
            sum_uk = u_hat_plus.slice(k-1).row(n+1) + sum_uk - u_hat_plus.slice(k).row(n);

            // mode spectrum
            arma::cx_rowvec u_hat_loop_1st = (f_hat_plus.t() - (sum_uk - lambda_hat.row(n) / 2.0));
            arma::rowvec u_hat_loop_2nd = freqs - arma::as_scalar(omega_plus(n,k));
            u_hat_loop_2nd = arma::pow(u_hat_loop_2nd, 2);
            u_hat_loop_2nd = u_hat_loop_2nd * arma::as_scalar(Alpha(k));
            u_hat_loop_2nd = u_hat_loop_2nd + 1.0;
            u_hat_plus.slice(k).row(n+1) = u_hat_loop_1st / u_hat_loop_2nd;

            // center frequencies
            arma::rowvec dc_loop = arma::abs(u_hat_plus.slice(k).row(n+1).subvec(T_div_two, T-1));
            double dc_scalar_loop = arma::accu(arma::pow(dc_loop, 2));
            double cross_pr_loop = arma::dot(freqs.subvec(T_div_two, T-1), arma::pow(dc_loop, 2));
            omega_plus(n+1,k) = cross_pr_loop / dc_scalar_loop;
          }

          arma::cx_cube dual_item = u_hat_plus.row(n+1);
          dual_item = arma::sum(dual_item, 2);                     // returns an arma::cube
          arma::cx_rowvec dual_vec = dual_item.slice(0).row(0);    // convert to rowvec

          // Dual ascent
          lambda_hat.row(n+1) = lambda_hat.row(n) + tau * (dual_vec - f_hat_plus.t());

          // loop counter
          n++;

          // converged yet?
          arma::cx_double uDiff_inner = arma::datum::eps;

          for (arma::uword z = 0; z < K; z++) {
            arma::cx_rowvec rep_cx = u_hat_plus.slice(z).row(n) - u_hat_plus.slice(z).row(n-1);
            arma::cx_rowvec rep_cx_conj = arma::conj(rep_cx);
            arma::cx_double iter_cx = 1.0 / T * arma::dot(rep_cx, rep_cx_conj);
            // arma::cx_double iter_cx = 1.0 / T * arma::accu(rep_cx % rep_cx_conj);      // alternative to dot-product
            uDiff_inner += iter_cx;
          }

          uDiff = std::abs(uDiff_inner);
          if (verbose) {
            if (((n+1) % 10) == 0) {
              Rcpp::Rcout << "Iteration: " << n + 1 << "   uDiff: " << uDiff << std::endl;
            }
          }
        }

        // Post-processing and cleanup

        // discard empty space if converged early
        N = std::min(N,n);

        if (verbose) {
          Rcpp::Rcout << "-----------------------------------------" << std::endl;
          Rcpp::Rcout << "The algorithm converged in iteration: " << N + 1 << std::endl;
          Rcpp::Rcout << "-----------------------------------------" << std::endl;
        }

        arma::mat omega = omega_plus.rows(0,N);

        // Signal reconstruction
        arma::cx_cube u_hat_plus_subs = u_hat_plus.cols(T_div_two, T-1);
        u_hat_plus_subs = u_hat_plus_subs.row(N-1);
        arma::cx_mat u_hat_plus_sq = squeeze_row(u_hat_plus_subs);                // squeezed
        arma::cx_mat u_hat_plus_sq_cnj = arma::conj(u_hat_plus_sq);               // conjugated squeezed

        arma::cx_mat u_hat = arma::zeros<arma::cx_mat>(T, K);
        u_hat.rows(T_div_two, T-1) = u_hat_plus_sq;                               // overwrite lower part

        arma::uvec non_contig_idx = arma::regspace<arma::uvec>(1, 1, T_div_two);
        non_contig_idx = arma::reverse(non_contig_idx);
        u_hat.rows(non_contig_idx) = u_hat_plus_sq_cnj;                            // overwrite upper part

        u_hat.row(0) = arma::conj(u_hat.row(u_hat.n_rows - 1));                    // replace first with conj. last row

        arma::mat u = arma::zeros<arma::mat>(K, len_t);
        for (arma::uword k_iter = 0; k_iter < K; k_iter++) {
          u.row(k_iter) = arma::real(arma::ifft(fftshift_mat(u_hat.col(k_iter), true))).t();    // inverse fftshift (matrix)
        }

        // remove mirror part
        u = u.cols(T/4, (3 * T/4)-1);

        // recompute spectrum
        u_hat.clear();
        u_hat.set_size(u.n_cols, K);

        for (arma::uword k_reit = 0; k_reit < K; k_reit++) {
          u_hat.col(k_reit) = fftshift_mat(arma::fft(u.row(k_reit))).t();
        }

        return Rcpp::List::create(Rcpp::Named("u") = u.t(),                        // the collection of decomposed modes
                                  Rcpp::Named("u_hat") = u_hat,                    // spectra of the modes
                                  Rcpp::Named("omega") = omega);                   // estimated mode center-frequencies
      }


      /**
       * Variational Mode Decomposition (2-dimensional)
       * Convergence (or number of iterations) might differ from execution to execution
       */

      Rcpp::List VMD_2D(arma::mat& signal,
                        double alpha,
                        double tau,
                        arma::uword K,
                        bool DC,
                        arma::uword init,
                        double tol,
                        bool verbose = false) {
        
        bool flag_isfinite = signal.is_finite();
        if (!flag_isfinite) {
          Rcpp::stop("The input data includes missing values! You have to replace these missing values before using the 2-dimensional Variational Mode Decomposition function!");
        }

        // Resolution of image
        arma::uword Hy = signal.n_rows;
        arma::uword Hx = signal.n_cols;

        if (Hy != Hx) {
          Rcpp::stop("The rows and columns (X- and Y- dimensions) of the input image ('signal' parameter) must be the same!\n  Transposing the matrices and arrays inside the while loop will increase the computation time considerably!\n  Thus, you have to resize the input image before using the function!");
        }

        arma::mat X = meshgrid_x(Hx, Hy);
        X = (1.0 + X) / Hx;
        arma::mat Y = meshgrid_y(Hx, Hy);
        Y = (1.0 + Y) / Hy;

        // Spectral Domain discretization
        double fx = 1.0 / Hx;
        double fy = 1.0 / Hy;
        arma::mat freqs_1 = X - 0.5 - fx;
        arma::mat freqs_2 = Y - 0.5 - fy;

        // N is the maximum number of iterations
        arma::uword N = 3000;

        // For future generalizations: alpha might be individual for each mode
        arma::mat rep_ones = arma::ones<arma::mat>(K, 1);
        arma::mat Alpha = alpha * rep_ones;

        // Construct f and f_hat
        arma::cx_mat f_hat = fftshift_mat(fft2(signal));

        // Storage matrices for (Fourier) Lagrange multiplier.
        arma::cx_mat mu_hat = arma::zeros<arma::cx_mat>(Hy, Hx);

        // N iterations at most, 2 spatial coordinates, K clusters
        if (K < 1 & DC) {
          Rcpp::stop("The 'K' parameter must be greater than 1 if the 'DC' parameter is TRUE!");
        }

        arma::cube omega = switch_2D(init, N, K, DC);

        // Stopping criteria tolerances
        double uDiff = tol + arma::datum::eps;
        double omegaDiff = tol + arma::datum::eps;

        // Storage matrices for (Fourier) modes. All iterations are not recorded.
        arma::cx_cube u_hat = arma::zeros<arma::cx_cube>(Hy, Hx, K);
        arma::cx_cube u_hat_old = arma::zeros<arma::cx_cube>(Hy, Hx, K);
        arma::cx_mat sum_uk = arma::zeros<arma::cx_mat>(Hy, Hx);        // accumulator
        arma::mat HilbertMask;

        if (verbose) {
          Rcpp::Rcout << "--------------------------------" << std::endl;
          Rcpp::Rcout << "The 2-dimensional VMD starts ..." << std::endl;
          Rcpp::Rcout << "--------------------------------" << std::endl;
        }

        // first run
        arma::uword n = 0;

        // run until convergence or max number of iterations
        while ((uDiff > tol || omegaDiff > tol) && n < N - 1) {

          // first things first
          int k = 0;

          // compute the halfplane mask for the 2D "analytic signal"
          HilbertMask = (arma::sign(freqs_1 * omega(n,0,k) + freqs_2 * omega(n,1,k)) + 1);

          // update first mode accumulator
          sum_uk = u_hat.slice(K-1) + sum_uk - u_hat.slice(k);

          // update first mode's spectrum through wiener filter (on half plane)
          arma::cx_mat first_eq_term = ((f_hat - sum_uk - mu_hat / 2.0) % HilbertMask);
          arma::mat second_eq_term = (1.0 + arma::as_scalar(Alpha(k,0)) * (arma::pow(freqs_1 - omega(n,0,k), 2.0) + arma::pow(freqs_2 - omega(n,1,k), 2.0)));
          u_hat.slice(k) = first_eq_term / second_eq_term;

          // update first mode's central frequency as spectral center of gravity
          if (!DC) {
            omega(n+1,0,k) = arma::sum(arma::sum(freqs_1 % (arma::abs(arma::pow(u_hat.slice(k), 2.0))))) / arma::sum(arma::sum(arma::abs(arma::pow(u_hat.slice(k), 2.0))));
            omega(n+1,1,k) = arma::sum(arma::sum(freqs_2 % (arma::abs(arma::pow(u_hat.slice(k), 2.0))))) / arma::sum(arma::sum(arma::abs(arma::pow(u_hat.slice(k), 2.0))));

            // keep omegas on same halfplane
            if (omega(n+1,1,k) < 0) {
              omega.slice(k).row(n+1) = -omega.slice(k).row(n+1);
            }
          }

          // recover full spectrum from analytic signal
          u_hat.slice(k) = fftshift_mat(arma::fft2(arma::real(arma::ifft2(fftshift_mat(u_hat.slice(k), true)))));

          // work on other modes
          for (arma::uword k = 1; k < K; k++) {

            // recompute Hilbert mask
            HilbertMask = (arma::sign(freqs_1 * omega(n,0,k) + freqs_2 * omega(n,1,k)) + 1);

            // update accumulator
            sum_uk = u_hat.slice(k-1) + sum_uk - u_hat.slice(k);

            // update signal spectrum
            arma::cx_mat first_eq_term = ((f_hat - sum_uk - mu_hat / 2.0) % HilbertMask);
            arma::mat second_eq_term = (1.0 + arma::as_scalar(Alpha(k,0)) * (arma::pow(freqs_1 - omega(n,0,k), 2.0) + arma::pow(freqs_2 - omega(n,1,k), 2.0)));
            u_hat.slice(k) = first_eq_term / second_eq_term;

            // update signal frequencies
            omega(n+1,0,k) = arma::sum(arma::sum(freqs_1 % (arma::abs(arma::pow(u_hat.slice(k), 2.0))))) / arma::sum(arma::sum(arma::abs(arma::pow(u_hat.slice(k), 2.0))));
            omega(n+1,1,k) = arma::sum(arma::sum(freqs_2 % (arma::abs(arma::pow(u_hat.slice(k), 2.0))))) / arma::sum(arma::sum(arma::abs(arma::pow(u_hat.slice(k), 2.0))));

            // keep omegas on same halfplane
            if (omega(n+1,1,k) < 0) {
              omega.slice(k).row(n+1) = -omega.slice(k).row(n+1);
            }

            // recover full spectrum from analytic signal
            u_hat.slice(k) = fftshift_mat(arma::fft2(arma::real(arma::ifft2(fftshift_mat(u_hat.slice(k), true)))));
          }

          // Gradient ascent for augmented Lagrangian
          arma::cx_mat u_hat_mat = arma::sum(u_hat,2);
          mu_hat += tau * (u_hat_mat - f_hat);

          // increment iteration counter
          n++;

          // convergence?
          arma::cx_double uDiff_inner = arma::datum::eps;
          omegaDiff = arma::datum::eps;

          for (arma::uword z = 0; z < K; z++) {
            omegaDiff += arma::accu(arma::sum(arma::abs(arma::pow(omega.row(n) - omega.row(n-1), 2.0))));
            uDiff_inner += arma::accu(arma::sum(1.0 / (Hx * Hy) * (u_hat.slice(k) - u_hat_old.slice(k)) % arma::conj((u_hat.slice(k) - u_hat_old.slice(k)))));
          }

          uDiff = std::abs(uDiff_inner);
          u_hat_old = u_hat;

          if (verbose) {
            if (((n+1) % 10) == 0) {
              Rcpp::Rcout << "Iteration: " << n + 1 << "  uDiff: " << uDiff << "  omegaDiff: " << omegaDiff << std::endl;
            }
          }
        }

        // discard empty space if converged early
        N = std::min(N,n);

        if (verbose) {
          Rcpp::Rcout << "-----------------------------------------" << std::endl;
          Rcpp::Rcout << "The algorithm converged in iteration: " << N + 1 << std::endl;
          Rcpp::Rcout << "-----------------------------------------" << std::endl;
        }

        // Signal Reconstruction
        arma::cube u = arma::zeros<arma::cube>(Hy,Hx,K);
        for (arma::uword w = 0; w < K; w++) {
          u.slice(w) = arma::real(arma::ifft2(fftshift_mat(u_hat.slice(w), true)));
        }

        // Should the omega-history be returned, or just the final results?
        // omega = omega(n,:,:);
        omega = omega.rows(0,n);

        return Rcpp::List::create(Rcpp::Named("u") = u,                        // the collection of decomposed modes
                                  Rcpp::Named("u_hat") = u_hat,                // spectra of the modes
                                  Rcpp::Named("omega") = omega);               // estimated mode center-frequencies
      }


      ~VarModeDecomp() { }
  };

}

#endif
