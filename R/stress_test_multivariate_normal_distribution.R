#' Stress-Testing on Correlations
#'
#' This function stresses the first two moments of a normal multivariate distribution
#' in a variety of correlation environments.
#'
#' @param .num_assets A \code{numeric} scalar. It's the dimension of the of simulated
#' market.
#' @param .sample_size A \code{numeric} scalar. It's the number of realizations of
#' each series.
#' @param .simulations A \code{numeric} scalar. The number of simulations to be
#' executed for each level of correlation.
#'
#' @importFrom rlang .data
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' #
stress_test_multivariate_normal_distribution <- function(.num_assets = 5, .sample_size = 200, .simulations = 500) {

  assertthat::assert_that(assertthat::is.number(.num_assets))
  assertthat::assert_that(assertthat::is.number(.sample_size))
  assertthat::assert_that(assertthat::is.number(.simulations))

  Mu  <- rep(0, .num_assets)
  Sig <- rep(1, .num_assets)

  Min_Theta <- 0
  Max_Theta <- 0.9
  Steps     <- 10

  Step   <- (Max_Theta - Min_Theta) / (Steps - 1)
  Thetas <- seq(from = Min_Theta, to = Max_Theta, by = Step)

  Stress_Loss_Mu      <- matrix(NA_real_, nrow = .simulations, ncol = Steps)
  Stress_Inef2_Mu     <- matrix(NA_real_, nrow = 1, ncol = Steps)
  Stress_Bias2_Mu     <- matrix(NA_real_, nrow = 1, ncol = Steps)
  Stress_Error2_Mu    <- matrix(NA_real_, nrow = 1, ncol = Steps)
  Stress_Loss_Sigma   <- matrix(NA_real_, nrow = .simulations, ncol = Steps)
  Stress_Inef2_Sigma  <- matrix(NA_real_, nrow = 1, ncol = Steps)
  Stress_Bias2_Sigma  <- matrix(NA_real_, nrow = 1, ncol = Steps)
  Stress_Error2_Sigma <- matrix(NA_real_, nrow = 1, ncol = Steps)

  Mu_hats    <- matrix(NA_real_, nrow = .simulations, ncol = .num_assets)
  Sigma_hats <- matrix(NA_real_, nrow = .simulations, ncol = .num_assets * .num_assets)
  l <- matrix(1, .simulations, 1) # used for matrix multiplication

  # each i represents a different stress-test scenario
  for (i in 1:Steps) {

    CyclesToGo <- Steps - i + 1
    cat("Cycles To Go:", CyclesToGo, "\n")

    Theta <- Thetas[i]
    C <- (1 - Theta) * diag(.num_assets) + Theta * matrix(1, .num_assets, .num_assets)
    Sigma <- diag(Sig) %*% C %*% diag(Sig)

    # specification <- ghyp::student.t(mu = Mu, sigma = Sigma)
    # true_data <- ghyp::rghyp(n = .sample_size, object = specification)

    # each j represents a simulation under a giver stress-test scenario
    for (j in 1:.simulations) {

      X <- ghyp::rghyp(
        n      = .sample_size,
        object = ghyp::gauss(mu = Mu, sigma = Sigma)
      )
      X_hat <- ghyp::fit.gaussmv(X)

      Mu_hats[j, ]    <- X_hat@mu
      Sigma_hats[j, ] <- as.vector(X_hat@sigma)

    }

    # Mu Analysis
    Stress_Loss_Mu[ , i]   <- rowSums((Mu_hats - l %*% Mu) ^ 2)
    Stress_Inef2_Mu[, i]   <- apply(Mu_hats, 2, stats::sd) %*% apply(Mu_hats, 2, stats::sd)
    Stress_Bias2_Mu[ , i]  <- sum((colMeans(Mu_hats) - Mu) ^ 2)
    Stress_Error2_Mu[ , i] <- mean(Stress_Loss_Mu[ , i])

    # Sigma Analysis
    Stress_Loss_Sigma[ , i]   <- rowSums((Sigma_hats - l %*% as.vector(Sigma)) ^ 2)
    Stress_Inef2_Sigma[ , i]  <- apply(Sigma_hats, 2, stats::sd) %*% apply(Sigma_hats, 2, stats::sd)
    Stress_Bias2_Sigma[ , i]  <- sum((colMeans(Sigma_hats) - as.vector(Sigma)) ^ 2)
    Stress_Error2_Sigma[ , i] <- mean(Stress_Loss_Sigma[ , i])

  }

  # exclude outliers
  exclusion <- Stress_Loss_Mu > stats::quantile(Stress_Loss_Mu, probs = 0.95)
  Stress_Loss_Mu <- Stress_Loss_Mu[apply(exclusion, 1, sum) == 0, ]

  exclusion <- Stress_Loss_Sigma > stats::quantile(Stress_Loss_Sigma, probs = 0.95)
  Stress_Loss_Sigma <- Stress_Loss_Sigma[apply(exclusion, 1, sum) == 0, ]

  data_plots <- list(mu    = list(Stress_Loss_Mu, Stress_Bias2_Mu, Stress_Inef2_Mu),
                     sigma = list(Stress_Loss_Sigma, Stress_Bias2_Sigma, Stress_Inef2_Sigma))

  # Start plotting
  plots <- vector("list", 2)
  names <- c("Location", "Dispersion")
  for (i in seq_along(plots)) {

    p1 <- data_plots[[i]][[1]] |>
      as.data.frame() |>
      tibble::as_tibble() |>
      `colnames<-`(factor(Thetas)) |>
      tidyr::pivot_longer(cols = dplyr::everything()) |>
      dplyr::mutate_if(is.character, as.numeric) |>
      dplyr::mutate(name = scales::percent(.data$name)) |>
      dplyr::mutate_if(is.character, as.factor) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$value, y = .data$name)) +
      ggridges::geom_density_ridges(scale = 1, stat = "binline", bins = 100, fill = "#03333e") +
      ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
      ggplot2::coord_flip() +
      ggplot2::labs(subtitle = paste0(names[i], ": error distribution"),
                    y        = "Correlation",
                    x        = NULL)

    p2 <- tibble::tibble(Correlations = Thetas,
                         Bias         = as.vector(data_plots[[i]][[2]]),
                         Inefficiency = as.vector(data_plots[[i]][[3]])) |>
      tidyr::pivot_longer(cols = -.data$Correlations) |>
      dplyr::mutate(Correlations = scales::percent(.data$Correlations)) |>
      dplyr::mutate_if(is.character, as.factor) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$Correlations, y = .data$value, fill = .data$name)) +
      ggplot2::geom_col() +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
      ggplot2::labs(subtitle = paste0(names[i], ": average error"),
                    y        = NULL,
                    x        = "Correlation",
                    fill     = "Decomposition") +
      ggplot2::scale_fill_manual(values = c("#03333e", "#9F9573")) +
      ggplot2::theme(legend.position = "bottom")

    plots[[i]] <- patchwork::wrap_plots(p1 / p2) +
      patchwork::plot_annotation(
        title    = "Estimation Error: Correlation StressTesting",
        subtitle = "Effects on the parameters of the Normal Distribution"
      )
    print(plots[[i]])

  }

}
