#' Stress-Test the Hyperbolic Distribution
#'
#' This function implements a stress-test in the HYP by stretching one the parameters
#' (mean, sigma or skewness) by 20 percent in both directions (up and down). Then, analyses
#' the impact that the shaked parameter causes in all the other parameters.
#'
#' @param .invariant An univariate data structure.
#' @param .what A \code{character} indicating which parameter should be stressed.
#' @param .simulations A \code{numeric} scalar. The number of simulations to be
#' executed for each level of correlation.
#' @param .symmetric A \code{character} flag. Should the distributions be symmetric? Default to FALSE.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' #
stress_test_univariate_hyp_distribution <- function(.invariant, .what = c("mu", "sigma", "skew"), .simulations = 300, .symmetric = FALSE) {
  UseMethod("stress_test_univariate_hyp_distribution", .invariant)
}

#' @rdname stress_test_univariate_hyp_distribution
#' @export
stress_test_univariate_hyp_distribution.default <- function(.invariant, .what = c("mu", "sigma", "skew"), .simulations = 300, .symmetric = FALSE) {
  stop(".invariant must be a tibble, xts or a matrix.", call. = FALSE)
}

#' @rdname stress_test_univariate_hyp_distribution
#' @export
stress_test_univariate_hyp_distribution.tbl <- function(.invariant, .what = c("mu", "sigma", "skew"), .simulations = 300, .symmetric = FALSE) {
  if (any(purrr::map_lgl(.invariant, lubridate::is.Date))) {
    .invariant <- .invariant |>
      timetk::tk_xts(silent = TRUE)
    .invariant <- matrix(.invariant, nrow = nrow(.invariant), ncol = ncol(.invariant))
  } else {
    .invariant <- as.matrix(.invariant)
  }
  stress_test_univariate_hyp_distribution_(.invariant = .invariant, .what = .what, .simulations = .simulations, .symmetric = .symmetric)
}

#' @rdname stress_test_univariate_hyp_distribution
#' @export
stress_test_univariate_hyp_distribution.xts <- function(.invariant, .what = c("mu", "sigma", "skew"), .simulations = 300, .symmetric = FALSE) {
  stress_test_univariate_hyp_distribution_(.invariant = as.matrix(.invariant), .what = .what, .simulations = .simulations, .symmetric = .symmetric)
}

#' @rdname stress_test_univariate_hyp_distribution
#' @export
stress_test_univariate_hyp_distribution.matrix <- function(.invariant, .what = c("mu", "sigma", "skew"), .simulations = 300, .symmetric = FALSE) {
  stress_test_univariate_hyp_distribution_(.invariant = .invariant, .what = .what, .simulations = .simulations, .symmetric = .symmetric)
}

#' @rdname stress_test_univariate_hyp_distribution
#' @export
stress_test_univariate_hyp_distribution.numeric <- function(.invariant, .what = c("mu", "sigma", "skew"), .simulations = 300, .symmetric = FALSE) {
  stress_test_univariate_hyp_distribution_(.invariant = .invariant, .what = .what, .simulations = .simulations, .symmetric = .symmetric)
}

#' @keywords internal
stress_test_univariate_hyp_distribution_ <- function(.invariant, .what = c("mu", "sigma", "skew"), .simulations = 300, .symmetric = FALSE) {

  .what <- match.arg(.what, choices = c("mu", "sigma", "skew"))

  assertthat::assert_that(assertthat::is.string(.what))
  assertthat::assert_that(assertthat::is.number(.simulations))
  assertthat::assert_that(assertthat::is.flag(.symmetric))

  n <- nrow(.invariant)

  hyp_fit <- ghyp::fit.hypuv(data = .invariant, symmetric = .symmetric, silent = TRUE)
  params <- ghyp::coefficients(hyp_fit)
  mu     <- params$mu
  sigma  <- params$sigma
  skew   <- params$gamma

  # Stress-Test
  if (.what == "mu") {
    Mus <- seq(mu * 0.8, mu * 1.2, length.out = 5)
  } else if (.what == "sigma") {
    Sigmas <- seq(sigma * 0.8, sigma * 1.2, length.out = 5)
  } else if (.what == "skew") {
    Skews <- seq(skew * 0.8, skew * 1.2, length.out = 5)
  } else {
    stop("Parameter currently not implemented. Valid params are: 'mu', 'sigma' and 'skew'.", call. = FALSE)
  }


  stress_loss_mu  <- matrix(0, nrow = .simulations, ncol = 5)
  stress_inef_mu  <- vector("numeric", 5)
  stress_bias_mu  <- vector("numeric", 5)
  stress_error_mu <- vector("numeric", 5)

  stress_loss_sigma  <- matrix(0, nrow = .simulations, ncol = 5)
  stress_inef_sigma  <- vector("numeric", 5)
  stress_bias_sigma  <- vector("numeric", 5)
  stress_error_sigma <- vector("numeric", 5)

  stress_loss_skew  <- matrix(0, nrow = .simulations, ncol = 5)
  stress_inef_skew  <- vector("numeric", 5)
  stress_bias_skew  <- vector("numeric", 5)
  stress_error_skew <- vector("numeric", 5)

  # each cycle represents a different stress-test scenario
  for (i in 1:5) {

    CyclesToGo <- 5 - i + 1
    cat("Cycles To Go:", CyclesToGo, "\n")

    if (.what == "mu") {
      mu_ <- Mus[i]
    } else if (.what == "sigma") {
      sigma_ <- Sigmas[i]
    } else if (.what == "skew") {
      skews_ <- Skews[i]
    }

    Mu_hats    <- vector("numeric", length = .simulations)
    Sigma_hats <- vector("numeric", length = .simulations)
    Skew_hats  <- vector("numeric", length = .simulations)

    # each cycle represents a simulation under a given stress-test scenario
    for (j in 1:.simulations) {

      if (.what == "mu") {
        W <- ghyp::rgig(n = n, lambda = hyp_fit@lambda, chi = hyp_fit@chi, psi = hyp_fit@psi)
        X <- mu_ + W * hyp_fit@gamma + sqrt(W) * as.vector(hyp_fit@sigma) * stats::rnorm(n)
        fit_simul <- ghyp::fit.hypuv(X, silent = TRUE)
      } else if (.what == "sigma") {
        W <- ghyp::rgig(n = n, lambda = hyp_fit@lambda, chi = hyp_fit@chi, psi = hyp_fit@psi)
        X <- hyp_fit@mu + W * hyp_fit@gamma + sqrt(W) * as.vector(sigma_) * stats::rnorm(n)
        fit_simul <- ghyp::fit.hypuv(X, silent = TRUE)
      } else {
        W <- ghyp::rgig(n = n, lambda = hyp_fit@lambda, chi = hyp_fit@chi, psi = hyp_fit@psi)
        X <- hyp_fit@mu + W * skews_ + sqrt(W) * as.vector(hyp_fit@sigma) * stats::rnorm(n)
        fit_simul <- ghyp::fit.hypuv(X, silent = TRUE)
      }

      coeffs <- ghyp::coefficients(fit_simul)

      Mu_hats[j]    <- coeffs$mu
      Sigma_hats[j] <- coeffs$sigma
      Skew_hats[j]  <- coeffs$gamma

    }

    # Mu Analysis
    stress_loss_mu[ , i] <- (Mu_hats - mu) ^ 2
    stress_inef_mu[i]   <- (stats::sd(Mu_hats) * ((.simulations - 1) / .simulations)) ^ 2
    stress_bias_mu[i]   <- (mean(Mu_hats) - mu) ^ 2
    stress_error_mu[i]  <- mean(stress_loss_mu[ , i])

    # Sigma Analysis
    stress_loss_sigma[ , i] <- (Sigma_hats - sigma) ^ 2
    stress_inef_sigma[i]   <- (stats::sd(Sigma_hats) * ((.simulations - 1) / .simulations)) ^ 2
    stress_bias_sigma[i]   <- (mean(Sigma_hats) - sigma) ^ 2
    stress_error_sigma[i]  <- mean(stress_loss_sigma[ , i])

    # Skew Analysis
    stress_loss_skew[ , i] <- (Skew_hats - skew) ^ 2
    stress_inef_skew[i]    <- (stats::sd(Skew_hats) * ((.simulations - 1) / .simulations)) ^ 2
    stress_bias_skew[i]    <- (mean(Skew_hats) - skew) ^ 2
    stress_error_skew[i]   <- mean(stress_loss_skew[ , i])

  }

  # exclude outliers
  exclusion <- stress_loss_mu > stats::quantile(stress_loss_mu, probs = 0.95)
  stress_loss_mu <- stress_loss_mu[apply(exclusion, 1, sum) == 0, ]

  exclusion <- stress_loss_sigma > stats::quantile(stress_loss_sigma, probs = 0.95)
  stress_loss_sigma <- stress_loss_sigma[apply(exclusion, 1, sum) == 0, ]

  exclusion <- stress_loss_skew > stats::quantile(stress_loss_skew, probs = 0.95)
  stress_loss_skew <- stress_loss_skew[apply(exclusion, 1, sum) == 0, ]


  data_plots <- list(mu    = list(stress_loss_mu, stress_bias_mu, stress_inef_mu),
                     sigma = list(stress_loss_sigma, stress_bias_sigma, stress_inef_sigma),
                     skew  = list(stress_loss_skew, stress_bias_skew, stress_inef_skew))


  if (.what == "mu") {
    x_axis_param <- Mus
  } else if (.what == "sigma") {
    x_axis_param <- Sigmas
  } else {
    x_axis_param <- Skews
  }

  plots <- vector("list", 3)
  names <- c("Location", "Dispersion", "Skewness")
  for (i in seq_along(plots)) {

    p1 <- data_plots[[i]][[1]] |>
      as.data.frame() |>
      tibble::as_tibble() |>
      `colnames<-`(factor(round(x_axis_param, 5))) |>
      tidyr::pivot_longer(cols = dplyr::everything()) |>
      dplyr::mutate_if(is.character, as.numeric) |>
      dplyr::mutate(name = scales::percent(.data$name)) |>
      dplyr::mutate_if(is.character, as.factor) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$value * 100, y = .data$name)) +
      ggridges::geom_density_ridges(scale = 1, stat = "binline", bins = 100, fill = "#03333e") +
      ggplot2::scale_x_continuous(labels = scales::percent_format()) +
      ggplot2::coord_flip() +
      ggplot2::labs(subtitle = paste0(names[i], ": error distribution"),
                    y        = NULL,
                    x        = NULL)


    p2 <- tibble::tibble(mu           = round(x_axis_param, 5),
                         Bias         = as.vector(data_plots[[i]][[2]]),
                         Inefficiency = as.vector(data_plots[[i]][[3]])) |>
      tidyr::pivot_longer(cols = -mu) |>
      dplyr::mutate(mu = scales::percent(mu)) |>
      dplyr::mutate_if(is.character, as.factor) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$mu, y = .data$value * 100, fill = .data$name)) +
      ggplot2::geom_col() +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::labs(subtitle = paste0(names[i], ": average error"),
                    y        = NULL,
                    x        = NULL,
                    fill     = "Decomposition") +
      ggplot2::scale_fill_manual(values = c("#03333e", "#9F9573")) +
      ggplot2::theme(legend.position = "bottom")

    plots[[i]] <- patchwork::wrap_plots(p1 / p2) +
      patchwork::plot_annotation(
        title    = "Estimation Error: Correlation StressTesting",
        subtitle = "Effects on the parameters of the Hyperbolic Distribution"
      )
    print(plots[[i]])

  }

}
