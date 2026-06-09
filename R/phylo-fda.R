#' Phylogenetic Flexible Discriminant Analysis (phylo.FDA)
#'
#' A phylogenetically-adjusted variant of `mda::fda()` inspired by Motani & Schmitz (2011, Evolution),
#' intended to generate discriminant predictions for fossil (or otherwise non-extant) taxa
#' while accounting for expected covariance among taxa under Brownian motion.
#'
#' @details
#' This implementation follows the `mda::fda` workflow (coding groups as an indicator matrix,
#' fitting regressions to Helmert-contrast-coded responses, then extracting discriminant
#' variates via an eigen/SVD step), but applies a phylogenetic whitening transform based on
#' the tree variance-covariance matrix `W = ape::vcv.phylo(tree)`.
#'
#' Two training functions are provided:
#' * `phylo.fda()` transforms both predictors `X` *and* the indicator matrix `Y` using the
#'   whitening transform.
#' * `phylo.fda_conventional()` transforms predictors only, leaving `Y` untransformed.
#'
#' @importFrom graphics abline axis layout legend par plot.new plot.window points rect text
#' @importFrom stats fitted predict runif
#'
#' @name phylo_fda
#' @keywords multivariate
NULL

# -----------------------------------------------------------------------------
# Internal helpers (not exported)
# -----------------------------------------------------------------------------

#' @noRd
.check_namespace_available <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
  }
}

#' @noRd
.check_inputs_and_align <- function(data, grp, tree) {
  # Basic checks
  if (is.null(tree)) stop("'tree' must be provided", call. = FALSE)
  if (!inherits(tree, "phylo")) stop("'tree' must be a 'phylo' object", call. = FALSE)
  
  X <- as.matrix(data)
  if (nrow(X) < 2L) stop("'data' must have at least 2 rows (taxa)", call. = FALSE)
  
  if (length(grp) != nrow(X)) {
    stop("'grp' must have the same length as nrow(data)", call. = FALSE)
  }
  
  tips <- tree$tip.label
  rn <- rownames(X)
  
  # Helper: can we map grp by names?
  grp_has_names <- !is.null(names(grp)) && all(tips %in% names(grp))
  
  if (length(tips) != nrow(X)) {
    # If labels allow alignment, do it. Otherwise stop.
    if (!is.null(rn) && all(tips %in% rn)) {
      X <- X[tips, , drop = FALSE]
      if (grp_has_names) {
        g <- as.factor(grp[tips])
      } else {
        g <- as.factor(grp)[match(tips, rn)]
      }
    } else if (grp_has_names) {
      ord <- match(tips, names(grp))
      X <- X[ord, , drop = FALSE]
      g <- as.factor(grp[ord])
    } else {
      stop(
        "Tree tips (length ", length(tips), ") do not match the number of rows in data (", nrow(X), "). ",
        "Provide rownames(data) that match tree$tip.label, or a named 'grp'.",
        call. = FALSE
      )
    }
  } else {
    # Same n: ensure order is correct if we can.
    if (!is.null(rn) && setequal(rn, tips) && !identical(rn, tips)) {
      X <- X[tips, , drop = FALSE]
    }
    
    # Build g in an alignment-safe way
    if (grp_has_names) {
      g <- as.factor(grp[tips])
    } else {
      g <- as.factor(grp)
      if (!is.null(rn) && setequal(rn, tips) && !identical(rn, tips)) {
        g <- g[match(tips, rn)]
      }
    }
  }
  
  list(X = X, g = g, tree = tree)
}

#' Apply Pagel's lambda directly to a phylogenetic variance–covariance matrix.
#'
#' This implements the standard lambda transform at the covariance level:
#' off-diagonal elements are multiplied by lambda, while the diagonal is unchanged.
#'
#' @noRd
.apply_lambda_to_vcv <- function(W, lambda) {
  if (is.null(lambda) || isTRUE(all.equal(lambda, 1))) return(W)

  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda)) {
    stop("lambda must be a single numeric value.")
  }

  Wl <- W
  off <- row(Wl) != col(Wl)
  Wl[off] <- lambda * Wl[off]
  Wl
}

.rescale_tree_lambda <- function(tree, val, treetrans = NULL) {
  # Deprecated: lambda is now applied directly to the phylogenetic VCV matrix.
  # This function is retained only for backward compatibility.
  warning("'treetrans'/.rescale_tree_lambda() is deprecated and ignored; lambda is applied to the VCV matrix.",
          call. = FALSE)
  tree
}

#' @noRd
.make_phylo_transform <- function(W) {
  # Whitening/transform matrix N derived from a phylogenetic covariance matrix W.
  if (is.null(W) || !is.matrix(W)) stop("'W' must be a covariance matrix", call. = FALSE)
  invW <- solve(W)

  # W should be symmetric; force symmetric eigen for numeric stability.
  eig <- eigen(invW, symmetric = TRUE)
  if (any(eig$values <= 0)) {
    stop("Non-positive eigenvalues encountered in inv(VCV). Tree/VCV may be singular.", call. = FALSE)
  }

  # N = V * sqrt(D) * V^T  (since symmetric)
  N <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

  # Determinant normalization to avoid arbitrary scaling.
  # Using det(N) can underflow/overflow for moderate n, producing 0/Inf and then Inf/NaN in downstream QR.
  # Since det(N) = prod(sqrt(eig$values)), we can compute the geometric mean scale safely on the log scale.
  log_divnum <- 0.5 * mean(log(eig$values))  # (1/n) * log(det(N))
  divnum <- exp(log_divnum)
  if (!is.finite(divnum) || divnum <= 0) {
    # Fallback: skip normalization if numerical issues arise (scaling is arbitrary anyway).
    return(N)
  }
  N / divnum
}

#' @noRd
.normalize_prior <- function(prior, cnames) {
  # Normalize and align class priors.
  #
  # Accepts either:
  #  - an unnamed numeric vector of length J (assumed to be in the order of cnames), or
  #  - a named numeric vector whose names include all cnames (order arbitrary; extras allowed).
  #
  # Returns a length-J numeric vector in cnames order that sums to 1.
  if (is.null(prior)) return(NULL)

  # Back-compat: if cnames is given as a single number, treat it as J.
  if (length(cnames) == 1 && is.numeric(cnames)) {
    cnames <- seq_len(as.integer(cnames))
  }
  cnames <- as.character(cnames)
  J <- length(cnames)

  if (!is.numeric(prior)) stop("'prior' must be a numeric vector", call. = FALSE)

  if (!is.null(names(prior))) {
    if (!all(cnames %in% names(prior))) {
      miss <- setdiff(cnames, names(prior))
      stop(
        "Named 'prior' must include all group names. Missing: ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    prior <- prior[cnames]
  } else {
    if (length(prior) != J) stop("'prior' has wrong length", call. = FALSE)
    names(prior) <- cnames
  }

  if (any(!is.finite(prior))) stop("'prior' must be finite", call. = FALSE)
  if (any(prior <= 0)) stop("'prior' must be strictly positive", call. = FALSE)

  s <- sum(prior)
  if (!is.finite(s) || s <= 0) stop("'prior' must sum to a positive finite value", call. = FALSE)
  prior / s
}

#' @noRd
.contr_fda <- function(p, contrast.default = stats::contr.helmert(length(p))) {
  # Internal analogue of mda's internal contr.fda, with safe argument handling.
  if (is.null(p) || length(p) < 2L) stop("'p' must be a probability vector of length >= 2", call. = FALSE)
  d <- dim(contrast.default)
  sqp <- sqrt(p / sum(p))
  x <- cbind(1, contrast.default) * outer(sqp, rep(1, d[2] + 1))
  qx <- base::qr(x)
  J <- qx$rank
  base::qr.qy(qx, diag(d[1])[, seq(2, J)]) / outer(sqp, rep(1, J - 1))
}

#' @noRd
.polyreg_modified <- function(x, y, w, ...) {
  y <- as.matrix(y)
  x <- as.matrix(x)

  if (iswt <- !missing(w)) {
    if (any(w <= 0)) stop("only positive weights", call. = FALSE)
    w <- sqrt(w)
    y <- y * w
    x <- x * w
  }

  qrx <- base::qr(x)
  coef <- as.matrix(base::qr.coef(qrx, y))
  fitted_vals <- base::qr.fitted(qrx, y)

  if ((df <- qrx$rank) < ncol(x)) coef[qrx$pivot, ] <- coef
  if (iswt) fitted_vals <- fitted_vals / w

  structure(
    list(fitted.values = fitted_vals, coef = coef, df = df),
    class = "polyreg.modified"
  )
}

#' @export
fitted.polyreg.modified <- function(object, ...) object$fitted.values

#' @export
predict.polyreg.modified <- function(object, newdata, ...) {
  if (missing(newdata)) {
    z <- fitted(object)
    if (is.null(z)) stop("need to supply newdata", call. = FALSE)
    return(z)
  }
  as.matrix(newdata) %*% object$coef
}

#' @noRd
.make_dummy_response <- function(g) {
  # n x J dummy matrix (one-hot) for group membership
  stats::model.matrix(~ g - 1L)
}

#' @noRd
.make_N_from_invW <- function(invW) {
  # N such that N %*% t(N) = invW (symmetric sqrt of invW)
  ev <- eigen(invW, symmetric = TRUE)

  # Guard tiny negative eigenvalues from numerical noise
  vals <- pmax(ev$values, 0)
  N <- ev$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(ev$vectors)

  # Optional scale normalization on log scale (avoid det under/overflow)
  pos <- vals > 0
  if (any(pos)) {
    divnum <- exp(0.5 * mean(log(vals[pos])))
    if (is.finite(divnum) && divnum > 0) N <- N / divnum
  }
  N
}


#' @noRd
.phylo_rss <- function(data, grp, tree, lambda = 1, intercept = TRUE) {
  aligned <- .check_inputs_and_align(data, grp, tree)
  X <- aligned$X
  g <- aligned$g
  tree <- aligned$tree
  
  if (intercept) X <- cbind(Intercept = rep(1, nrow(X)), X)
  
  # Phylogenetic covariance with lambda applied
  W <- ape::vcv.phylo(tree)
  W <- .apply_lambda_to_vcv(W, lambda)
  invW <- solve(W)
  
  # Dummy response (one-hot groups)
  y <- .make_dummy_response(g)
  
  # GLS fit: B = (X' invW X)^{-1} X' invW y
  Xt_iW <- t(X) %*% invW
  bhat <- solve(Xt_iW %*% X) %*% (Xt_iW %*% y)
  yhat <- X %*% bhat
  
  E <- y - yhat
  RSS_mat <- t(E) %*% invW %*% E
  list(RSS = sum(diag(RSS_mat)))
}

#' @noRd
.phylo_rss_pred <- function(data, grp, tree, test_idx, lambda = 1, intercept = TRUE) {
  # Prediction-style RSS for training set only (drop test rows, fit in transformed space)
  aligned <- .check_inputs_and_align(data, grp, tree)
  X <- aligned$X
  g <- aligned$g
  tree <- aligned$tree
  
  n <- nrow(X)
  if (any(test_idx < 1 | test_idx > n)) stop("'test_idx' out of bounds", call. = FALSE)
  test_idx <- unique(as.integer(test_idx))
  train_idx <- setdiff(seq_len(n), test_idx)
  
  if (length(train_idx) < 2L) stop("Too few training taxa after dropping test set", call. = FALSE)
  
  if (intercept) X <- cbind(Intercept = rep(1, nrow(X)), X)
  
  # Apply lambda to covariance, build N from invW, transform X and y
  W <- ape::vcv.phylo(tree)
  W <- .apply_lambda_to_vcv(W, lambda)
  invW <- solve(W)
  N <- .make_N_from_invW(invW)
  
  y <- .make_dummy_response(g)
  
  Xs <- N %*% X
  Ys <- N %*% y
  
  # Train only on training rows
  Xtr <- Xs[train_idx, , drop = FALSE]
  Ytr <- Ys[train_idx, , drop = FALSE]
  
  # OLS in transformed space
  bhat <- solve(t(Xtr) %*% Xtr) %*% (t(Xtr) %*% Ytr)
  Etr <- Ytr - Xtr %*% bhat
  
  list(RSS = sum(Etr^2))
}

#' Optimize lambda by residual sum of squares
#'
#' Fits phylogenetic FDA across a grid of Pagel's lambda values and returns
#' the lambda value that minimizes residual sum of squares. If multiple lambda
#' values tie, the largest lambda among the minima is returned.
#'
#' @param measurements Numeric matrix or data frame of predictor measurements,
#'   with taxa in rows.
#' @param grps Factor or vector giving group membership for each row of
#'   `measurements`.
#' @param mytree Phylogeny of class `"phylo"`.
#' @param sstep Step size for the lambda grid.
#' @param srange Numeric length-2 vector giving the lambda search range.
#' @param intercept Logical; if `TRUE`, include an intercept in the model matrix.
#'
#' @return A list with `optlambda`, the selected lambda value, and `rslt`, a
#'   data frame containing lambda values and residual sum of squares.
#'
#' @export
optLambda <- function(measurements,
                      grps,
                      mytree,
                      sstep = 0.01,
                      srange = c(0, 1),
                      intercept = TRUE) {
  lambdas <- seq(min(srange), max(srange), by = sstep)
  rss <- numeric(length(lambdas))
  
  for (i in seq_along(lambdas)) {
    rss[i] <- .phylo_rss(measurements, grps, mytree, lambda = lambdas[i], intercept = intercept)$RSS
  }
  
  # If tie, follow Motani/Schmitz convention: choose the maximum lambda among minima
  min_rss <- min(rss)
  opt <- max(lambdas[rss == min_rss])
  
  rslt <- data.frame(Lambda = lambdas, RSS = rss)
  list(optlambda = opt, rslt = rslt)
}

#' Optimize lambda by predictive residual sum of squares
#'
#' Fits the phylogenetic FDA whitening step across a grid of Pagel's lambda
#' values and evaluates residual sum of squares on a specified test subset.
#' If multiple lambda values tie, the largest lambda among the minima is returned.
#'
#' @param measurements Numeric matrix or data frame of predictor measurements,
#'   with taxa in rows.
#' @param grps Factor or vector giving group membership for each row of
#'   `measurements`.
#' @param mytree Phylogeny of class `"phylo"`.
#' @param test_idx Integer or logical index specifying the test rows.
#' @param sstep Step size for the lambda grid.
#' @param srange Numeric length-2 vector giving the lambda search range.
#' @param intercept Logical; if `TRUE`, include an intercept in the model matrix.
#'
#' @return A list with `optlambda`, the selected lambda value, and `rslt`, a
#'   data frame containing lambda values and predictive residual sum of squares.
#'
#' @export
optLambda.pred <- function(measurements,
                           grps,
                           mytree,
                           test_idx,
                           sstep = 0.01,
                           srange = c(0, 1),
                           intercept = TRUE) {
  lambdas <- seq(min(srange), max(srange), by = sstep)
  rss <- numeric(length(lambdas))
  
  for (i in seq_along(lambdas)) {
    rss[i] <- .phylo_rss_pred(measurements, grps, mytree, test_idx = test_idx,
                              lambda = lambdas[i], intercept = intercept)$RSS
  }
  
  min_rss <- min(rss)
  opt <- max(lambdas[rss == min_rss])
  
  rslt <- data.frame(Lambda = lambdas, RSS = rss)
  list(optlambda = opt, rslt = rslt)
}


# -----------------------------------------------------------------------------
# S3 methods
# -----------------------------------------------------------------------------

#' @export
print.phylo.fda <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nDimension:", format(x$dimension), "\n")
  cat("\nPercent Between-Group Variance Explained:\n")
  print(round(x$percent.explained, 2))

  if (!is.null(x$confusion)) {
    n <- as.integer(sum(x$confusion))
    err <- format(round(attr(x$confusion, "error"), 5))
    cat("\nTraining Misclassification Error:", err, "( N =", n, ")\n")
  }
  invisible(x)
}

#' @export
predict.phylo.fda <- function(object,
                              newdata,
                              type = c("class", "variates", "posterior", "hierarchical", "distances"),
                              prior = NULL,
                              dimension = NULL,
                              ...) {
  type <- match.arg(type)
  
  means <- object$means
  Jk <- dim(means)
  J <- Jk[1]
  k <- Jk[2]
  
  # Default dimension: trained dimension
  if (is.null(dimension)) dimension <- k
  if (type == "hierarchical") {
    dimension.set <- dimension[dimension <= k]
    if (!length(dimension.set)) dimension.set <- k
    dimension <- max(dimension.set)
  } else {
    dimension <- min(max(dimension), k)
  }
  # Build design and get fitted values from regression
  if (missing(newdata)) {
    y <- predict(object$fit)
  } else {
    # Accept newdata as data.frame/list (build model matrix) OR as a numeric matrix.
    if (inherits(newdata, "data.frame") || is.list(newdata)) {
      Terms <- stats::delete.response(object$terms)
      attr(Terms, "intercept") <- 0

      # If training expected an explicit Intercept column, ensure it exists in newdata
      if (isTRUE(object$intercept) && !"Intercept" %in% names(newdata)) {
        newdata$Intercept <- 1
      }

      newdata <- stats::model.matrix(Terms, newdata)
    } else {
      newdata <- as.matrix(newdata)

      # If the fit expects an Intercept column but newdata doesn't have it, add it.
      if (isTRUE(object$intercept)) {
        coef_rows <- nrow(object$fit$coef)
        if (ncol(newdata) == coef_rows - 1L) {
          newdata <- cbind(Intercept = 1, newdata)
        }
      }

      # If we can, reorder columns to match the fitted coefficient rows.
      coef_names <- rownames(object$fit$coef)
      if (!is.null(coef_names) && !is.null(colnames(newdata)) &&
          setequal(coef_names, colnames(newdata)) &&
          !identical(coef_names, colnames(newdata))) {
        newdata <- newdata[, coef_names, drop = FALSE]
      }

      # Final guard: dimensions must match.
      if (ncol(newdata) != nrow(object$fit$coef)) {
        stop(
          "newdata has ", ncol(newdata), " column(s) but the fitted regression expects ",
          nrow(object$fit$coef), ". Ensure newdata contains the same trait columns used in training, ",
          "and that an Intercept column is included only if it was used during fitting.",
          call. = FALSE
        )
      }
    }

    y <- predict(object$fit, newdata)
  }
  
  # Project into discriminant space
  y <- y %*% object$theta.mod[, seq(dimension), drop = FALSE]
  
  lambda <- object$values
  alpha <- sqrt(lambda[seq(dimension)])
  sqima <- sqrt(1 - lambda[seq(dimension)])
  newZ <- scale(y, FALSE, sqima * alpha)
  

  # Preserve training order for alignment-safe outputs
  if (missing(newdata) && !is.null(object$train_order) && length(object$train_order) == nrow(newZ)) {
    rownames(newZ) <- object$train_order
  }
  
  if (is.null(prior)) {
    prior <- object$prior
  } else {
    prior <- .normalize_prior(prior, dimnames(object$means)[[1]])
  }
  
  means <- means[, seq(dimension), drop = FALSE]
  
  dist <- function(x, mean, m = ncol(mean)) (scale(x, mean, FALSE)^2) %*% rep(1, m)
  
  switch(type,
         variates = return(newZ),
         
         class = {
           n <- nrow(newZ)
           prior_log <- 2 * log(prior)
           mindist <- dist(newZ, means[1, ], dimension) - prior_log[1]
           pclass <- rep(1, n)
           
           for (i in seq(2, J)) {
             ndist <- dist(newZ, means[i, ], dimension) - prior_log[i]
             l <- ndist < mindist
             pclass[l] <- i
             mindist[l] <- ndist[l]
           }
           
           out <- factor(pclass, levels = seq(J), labels = dimnames(means)[[1]])
           
           # Name predictions if we can
           if (!is.null(rownames(newZ)) && length(rownames(newZ)) == length(out)) {
             names(out) <- rownames(newZ)
           }
           
           return(out)
         },
         
         posterior = {
           pclass <- matrix(0, nrow(newZ), J)
           for (i in seq(J)) {
             pclass[, i] <- exp(-0.5 * dist(newZ, means[i, ], dimension)) * prior[i]
           }
           dimnames(pclass) <- list(dimnames(newZ)[[1]], dimnames(means)[[1]])
           return(pclass / drop(pclass %*% rep(1, J)))
         },
         
         hierarchical = {
           prior_log <- 2 * log(prior)
           Pclass <- vector("list", length(dimension.set))
           names(Pclass) <- paste("D", dimension.set, sep = "")
           
           for (ad in seq(along = dimension.set)) {
             d <- dimension.set[ad]
             dd <- seq(d)
             mindist <- dist(newZ[, dd, drop = FALSE], means[1, dd, drop = FALSE], d) - prior_log[1]
             pclass <- rep(1, nrow(newZ))
             
             for (i in seq(2, J)) {
               ndist <- dist(newZ[, dd, drop = FALSE], means[i, dd, drop = FALSE], d) - prior_log[i]
               l <- ndist < mindist
               pclass[l] <- i
               mindist[l] <- ndist[l]
             }
             
             Pclass[[ad]] <- factor(pclass, levels = seq(J), labels = dimnames(means)[[1]])
           }
           
           rns <- dimnames(newZ)[[1]]
           if (is.null(rns)) rns <- paste(seq(nrow(newZ)))
           return(structure(Pclass, class = "data.frame", row.names = rns, dimensions = dimension.set))
         },
         
         distances = {
           dclass <- matrix(0, nrow(newZ), J)
           for (i in seq(J)) dclass[, i] <- dist(newZ, means[i, ], dimension)
           dimnames(dclass) <- list(dimnames(newZ)[[1]], dimnames(means)[[1]])
           return(dclass)
         }
  )
}

# -----------------------------------------------------------------------------
# Public training functions
# -----------------------------------------------------------------------------

#' Fit phylogenetic flexible discriminant analysis (transform X and Y)
#'
#' @param data Numeric matrix/data.frame of predictors (rows = taxa).
#' @param grp Factor (or vector coercible to factor) giving group membership per row.
#' @param tree A `phylo` object whose tip labels correspond to taxa.
#' @param lambda Pagel's lambda value used for tree rescaling.
#' @param treetrans Function used to rescale the tree (default `geiger::rescale`).
#' @param dimension Number of discriminant dimensions to retain (default J-1).
#' @param eps Small value for numerical truncation.
#' @param keep.fitted Keep fitted values (default heuristic `n*dimension < 1000`).
#' @param method Regression backend (default internal `.polyreg_modified`).
#' @param intercept Add intercept column to predictors.
#' @param eqprior Use equal priors.
#' @param prior Optional prior vector (will be normalized).
#'
#' @return An object of class `phylo.fda`.
#' @export
phylo.fda <- function(data, grp, tree,
                      lambda = 1,
                      treetrans = NULL,
                      dimension = NULL,
                      eps = .Machine$double.eps,
                      keep.fitted = NULL,
                      method = .polyreg_modified,
                      intercept = TRUE,
                      eqprior = FALSE,
                      prior = NULL) {

  .check_namespace_available("ape")
  .check_namespace_available("mda")

  this.call <- match.call()
  aligned <- .check_inputs_and_align(data, grp, tree)
  X <- aligned$X
  g <- aligned$g
  tree <- aligned$tree

  if (intercept) X <- cbind(Intercept = rep(1, nrow(X)), X)
  if (!is.null(treetrans) && !missing(treetrans)) {
    warning("'treetrans' is deprecated and ignored; lambda is applied to the VCV matrix.", call. = FALSE)
  }

  W <- ape::vcv.phylo(tree)
  W <- .apply_lambda_to_vcv(W, lambda)

  ng <- nlevels(g)
  cnames <- levels(g)
  J <- length(cnames)
  N <- .make_phylo_transform(W)
  Xstar <- N %*% X
  rownames(Xstar) <- rownames(X)
  n <- nrow(Xstar)

  # Indicator matrix
  Y <- matrix(0, n, ng)
  for (i in seq_len(n)) Y[i, g[i]] <- 1

  # Motani & Schmitz-style: transform Y as well
  Ystar <- N %*% Y

  # Group proportions
  weights <- rep(1, n)
  dp <- tapply(weights, as.numeric(g), sum) / n

  theta <- .contr_fda(dp, stats::contr.helmert(J))
  Theta <- Ystar %*% theta

  fit <- method(Xstar, Theta, weights)

  rss <- t(Theta - fitted(fit)) %*% (Theta - fitted(fit))
  ssm <- t(Theta) %*% fitted(fit) / n

  ed <- svd(ssm, nu = 0)
  thetan <- ed$v
  lambda_vals <- ed$d
  lambda_vals[lambda_vals > 1 - eps] <- 1 - eps

  discr.eigen <- lambda_vals / (1 - lambda_vals)
  pe <- (100 * cumsum(discr.eigen)) / sum(discr.eigen)

  if (is.null(dimension)) dimension <- J - 1
  dimension <- min(dimension, sum(lambda_vals > eps))

  if (dimension == 0) {
    warning("degenerate problem; no discrimination")
    return(structure(list(dimension = 0, fit = fit, call = this.call), class = "phylo.fda"))
  }

  thetan <- thetan[, seq_len(dimension), drop = FALSE]
  pe <- pe[seq_len(dimension)]

  alpha <- sqrt(lambda_vals[seq_len(dimension)])
  sqima <- sqrt(1 - lambda_vals[seq_len(dimension)])

  vnames <- paste0("v", seq_len(dimension))
  means <- scale(theta %*% thetan, FALSE, sqima / alpha)
  dimnames(means) <- list(cnames, vnames)

  # Prior
  pr <- colSums(Ystar) / sum(colSums(Ystar))
  if (eqprior) pr <- rep(1 / ng, ng)
  if (!is.null(prior)) pr <- .normalize_prior(prior, cnames)

  # Terms for data-frame newdata
  dset <- as.data.frame(X)
  dset$grp <- g
  Terms <- stats::terms(grp ~ ., data = dset)

  if (is.null(keep.fitted)) keep.fitted <- (n * dimension < 1000)

  obj <- structure(
    list(
      percent.explained = pe,
      values = lambda_vals,
      means = means,
      theta.mod = thetan,
      dimension = dimension,
      prior = pr,
      theta = theta,
      fit = fit,
      call = this.call,
      terms = Terms,
      N = N,
      lambda = lambda,
      train_order = rownames(Xstar),
      intercept = intercept,
      eqprior = eqprior,
      rss = sum(diag(rss))
    ),
    class = "phylo.fda"
  )

  obj$confusion <- mda::confusion(predict(obj), g)
  if (!keep.fitted) obj$fit$fitted.values <- NULL

  obj
}

#' Fit phylogenetic FDA (transform X only; conventional Y)
#'
#' @inheritParams phylo.fda
#' @return An object of class `phylo.fda`.
#' @export
phylo.fda_conventional <- function(data, grp, tree,
                                   lambda = 1,
                                   treetrans = NULL,
                                   dimension = NULL,
                                   eps = .Machine$double.eps,
                                   keep.fitted = NULL,
                                   method = .polyreg_modified,
                                   intercept = TRUE,
                                   eqprior = FALSE,
                                   prior = NULL) {

  .check_namespace_available("ape")
  .check_namespace_available("mda")

  this.call <- match.call()
  aligned <- .check_inputs_and_align(data, grp, tree)
  X <- aligned$X
  g <- aligned$g
  tree <- aligned$tree

  if (intercept) X <- cbind(Intercept = rep(1, nrow(X)), X)
  if (!is.null(treetrans) && !missing(treetrans)) {
    warning("'treetrans' is deprecated and ignored; lambda is applied to the VCV matrix.", call. = FALSE)
  }

  W <- ape::vcv.phylo(tree)
  W <- .apply_lambda_to_vcv(W, lambda)

  ng <- nlevels(g)
  cnames <- levels(g)
  J <- length(cnames)
  N <- .make_phylo_transform(W)
  Xstar <- N %*% X
  rownames(Xstar) <- rownames(X)
  n <- nrow(Xstar)

  Y <- matrix(0, n, ng)
  for (i in seq_len(n)) Y[i, g[i]] <- 1

  weights <- rep(1, n)
  dp <- tapply(weights, as.numeric(g), sum) / n

  theta <- .contr_fda(dp, stats::contr.helmert(J))
  Theta <- Y %*% theta

  fit <- method(Xstar, Theta, weights)

  rss <- t(Theta - fitted(fit)) %*% (Theta - fitted(fit))
  ssm <- t(Theta) %*% fitted(fit) / n

  ed <- svd(ssm, nu = 0)
  thetan <- ed$v
  lambda_vals <- ed$d
  lambda_vals[lambda_vals > 1 - eps] <- 1 - eps

  discr.eigen <- lambda_vals / (1 - lambda_vals)
  pe <- (100 * cumsum(discr.eigen)) / sum(discr.eigen)

  if (is.null(dimension)) dimension <- J - 1
  dimension <- min(dimension, sum(lambda_vals > eps))

  if (dimension == 0) {
    warning("degenerate problem; no discrimination")
    return(structure(list(dimension = 0, fit = fit, call = this.call), class = "phylo.fda"))
  }

  thetan <- thetan[, seq_len(dimension), drop = FALSE]
  pe <- pe[seq_len(dimension)]

  alpha <- sqrt(lambda_vals[seq_len(dimension)])
  sqima <- sqrt(1 - lambda_vals[seq_len(dimension)])

  vnames <- paste0("v", seq_len(dimension))
  means <- scale(theta %*% thetan, FALSE, sqima / alpha)
  dimnames(means) <- list(cnames, vnames)

  pr <- colSums(Y) / sum(colSums(Y))
  if (eqprior) pr <- rep(1 / ng, ng)
  if (!is.null(prior)) pr <- .normalize_prior(prior, cnames)

  dset <- as.data.frame(X)
  dset$grp <- g
  Terms <- stats::terms(grp ~ ., data = dset)

  if (is.null(keep.fitted)) keep.fitted <- (n * dimension < 1000)

  obj <- structure(
    list(
      percent.explained = pe,
      values = lambda_vals,
      means = means,
      theta.mod = thetan,
      dimension = dimension,
      prior = pr,
      theta = theta,
      fit = fit,
      call = this.call,
      terms = Terms,
      N = N,
      lambda = lambda,
      train_order = rownames(Xstar),
      intercept = intercept,
      eqprior = eqprior,
      rss = sum(diag(rss))
    ),
    class = "phylo.fda"
  )

  obj$confusion <- mda::confusion(predict(obj), g)
  if (!keep.fitted) obj$fit$fitted.values <- NULL

  obj
}

#' Prepare a dataset for phylogenetic FDA
#'
#' Aligns data to tree, optionally resolves polytomies, applies simple trait transforms,
#' and splits training vs test taxa based on an "unknown" label (and NA).
#'
#' @param tree A phylo object.
#' @param dat  A data.frame with taxa, group labels, and trait columns.
#' @param taxon_col Name of the taxon column in dat.
#' @param group_col Name of the grouping column in dat.
#' @param trait_cols Character vector of trait column names.
#' @param unknown_label Label in group_col identifying test taxa (default "unknown").
#' @param log10_traits Logical; log10-transform trait matrix.
#' @param signif_digits Integer or NULL; apply signif() to traits if not NULL.
#' @param force_binary Logical; if TRUE, resolve polytomies via multi2di(random=TRUE).
#'
#' @return A list with aligned tree/data and training/test splits.
#' @export
prep_phylo_fda_dataset <- function(tree,
                                   dat,
                                   taxon_col = "taxon",
                                   group_col = "groups",
                                   trait_cols,
                                   unknown_label = "unknown",
                                   log10_traits = TRUE,
                                   signif_digits = 4,
                                   force_binary = TRUE) {
  
  # ---- checks ----
  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop("'tree' must be a 'phylo' object", call. = FALSE)
  }
  dat <- as.data.frame(dat)
  
  if (!taxon_col %in% names(dat)) stop("Missing taxon_col in dat", call. = FALSE)
  if (!group_col %in% names(dat)) stop("Missing group_col in dat", call. = FALSE)
  if (missing(trait_cols) || !length(trait_cols)) stop("You must supply trait_cols", call. = FALSE)
  if (!all(trait_cols %in% names(dat))) stop("Some trait_cols are missing in dat", call. = FALSE)
  
  # ---- tree: enforce binary if requested ----
  if (isTRUE(force_binary) && !ape::is.binary(tree)) {
    tree <- ape::multi2di(tree, random = TRUE)
  }
  
  # ---- set rownames to taxa ----
  rownames(dat) <- dat[[taxon_col]]
  
  tips <- tree$tip.label
  rn   <- rownames(dat)
  
  # ---- prune to intersection ----
  keep <- intersect(tips, rn)
  if (length(keep) < 2L) {
    stop("After aligning tree and data, fewer than 2 taxa remain.", call. = FALSE)
  }
  
  # drop tips not in data
  drop_from_tree <- setdiff(tips, keep)
  if (length(drop_from_tree)) tree <- ape::drop.tip(tree, drop_from_tree)
  
  # drop rows not in tree
  drop_from_data <- setdiff(rn, keep)
  if (length(drop_from_data)) dat <- dat[setdiff(rn, drop_from_data), , drop = FALSE]
  
  # reorder data to tree tip order
  dat <- dat[tree$tip.label, , drop = FALSE]
  
  # ---- extract group + traits ----
  grp_all <- dat[[group_col]]
  names(grp_all) <- rownames(dat)
  
  X_all <- as.matrix(dat[, trait_cols, drop = FALSE])
  
  # ---- transform traits (optional) ----
  if (isTRUE(log10_traits)) X_all <- log10(X_all)
  if (!is.null(signif_digits)) X_all <- signif(X_all, signif_digits)
  
  # ---- split training vs test ----
  is_test <- (grp_all == unknown_label) | is.na(grp_all)
  
  taxa_test  <- names(grp_all)[is_test]
  taxa_train <- names(grp_all)[!is_test]
  
  X_train    <- X_all[!is_test, , drop = FALSE]
  grp_train  <- factor(grp_all[!is_test])
  
  X_test     <- X_all[is_test, , drop = FALSE]
  grp_test   <- grp_all[is_test]
  
  tree_train <- tree
  if (length(taxa_test)) tree_train <- ape::drop.tip(tree_train, taxa_test)
  
  list(
    tree_full  = tree,
    tree_train = tree_train,
    X_all      = X_all,
    grp_all    = grp_all,
    taxa_train = taxa_train,
    taxa_test  = taxa_test,
    X_train    = X_train,
    grp_train  = grp_train,
    X_test     = X_test,
    grp_test   = grp_test
  )
}



#' Plot class predictions across lambda values
#'
#' Convenience wrapper to (i) fit phylogenetic FDA across a grid of lambda values,
#' (ii) predict classes and posterior probabilities for a set of test taxa, and
#' (iii) visualize results as a lambda-by-taxon matrix.
#'
#' Colors encode predicted class; opacity encodes max posterior probability.
#'
#' @param X_train Training trait matrix (rows = taxa).
#' @param grp_train Training grouping factor.
#' @param tree_train Training phylogeny (must match training taxa).
#' @param X_test Test trait matrix (rows = taxa).
#' @param lambda_grid Numeric vector of lambda values.
#' @param lambda_star Optional numeric scalar; draw a vertical dashed line at the
#'   nearest lambda grid value.
#' @param prior Optional prior vector (unnamed or named). If named, names should
#'   match levels(grp_train). If unnamed, it is assumed to be in the order of
#'   levels(grp_train). Will be normalized to sum to 1.
#' @param taxa Optional character vector of taxa (subset of rownames(X_test)) to
#'   plot, in the exact order desired.
#' @param test_subset Optional integer indices into rownames(X_test) to plot.
#' @param taxa_idx Optional deprecated alias for `test_subset`; retained for backward compatibility.
#' @param class_cols Optional named vector mapping class level -> color.
#' @param alpha_floor Minimum opacity between 0 and 1.
#' @param alpha_ceil Maximum opacity between 0 and 1.
#' @param sep_every Draw horizontal white separators every N taxa (default 1).
#'   Set NULL to disable.
#' @param legend_panel Logical; if TRUE, draw a dedicated legend panel.
#' @param eqprior Passed to phylo.fda (default FALSE).
#' @param ... Passed to phylo.fda.
#'
#' @return Invisibly returns a list with pred_mat, conf_mat, lambda_grid, and
#'   taxa plotted.
#' @export
plot_lambda_matrix <- function(X_train,
                               grp_train,
                               tree_train,
                               X_test,
                               lambda_grid = seq(0, 1, by = 0.05),
                               lambda_star = NULL,
                               prior = NULL,
                               taxa = NULL,
                               test_subset = NULL,
                               taxa_idx = NULL,
                               class_cols = NULL,
                               alpha_floor = 0.10,
                               alpha_ceil = 1.00,
                               sep_every = 1,
                               legend_panel = TRUE,
                               eqprior = FALSE,
                               ...) {

  if (missing(X_train) || missing(grp_train) || missing(tree_train) || missing(X_test)) {
    stop("You must supply X_train, grp_train, tree_train, and X_test.", call. = FALSE)
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(oldpar); graphics::layout(1) }, add = TRUE)


  X_train <- as.matrix(X_train)
  X_test  <- as.matrix(X_test)

  if (!inherits(tree_train, "phylo")) stop("'tree_train' must be a 'phylo' object", call. = FALSE)
  if (nrow(X_train) < 2L) stop("X_train must have at least 2 rows", call. = FALSE)
  if (nrow(X_test) < 1L) stop("X_test must have at least 1 row", call. = FALSE)
  if (ncol(X_train) != ncol(X_test)) stop("X_train and X_test must have the same number of columns", call. = FALSE)

  grp_train <- as.factor(grp_train)

  # ---- prior handling (vector form allowed) ----
  prior <- .prior_from_vector_or_null(prior, grp_train)

  # ---- compute prediction matrices (full test set, then subset) ----
  mats <- .lambda_prediction_mats(
    X_train = X_train,
    grp_train = grp_train,
    tree_train = tree_train,
    X_test = X_test,
    lambda_grid = lambda_grid,
    prior = prior,
    eqprior = eqprior,
    ...
  )

  pred_mat <- mats$pred_mat
  conf_mat <- mats$conf_mat

  sub <- .subset_taxa_mats(pred_mat, conf_mat, taxa = taxa, test_subset = test_subset, taxa_idx = taxa_idx)
  pred_plot <- sub$pred
  conf_plot <- sub$conf

  taxa_plot <- rownames(pred_plot)
  n_taxa <- nrow(pred_plot)
  n_lam  <- ncol(pred_plot)

  class_cols <- .default_class_cols(class_cols, grp_train)

  if (!is.numeric(alpha_floor) || alpha_floor < 0 || alpha_floor > 1) stop("alpha_floor must be in [0,1]", call. = FALSE)
  if (!is.numeric(alpha_ceil)  || alpha_ceil  < 0 || alpha_ceil  > 1) stop("alpha_ceil must be in [0,1]", call. = FALSE)
  if (alpha_floor > alpha_ceil) stop("alpha_floor must be <= alpha_ceil", call. = FALSE)

  alpha <- alpha_floor + (alpha_ceil - alpha_floor) * conf_plot
  alpha[alpha < alpha_floor] <- alpha_floor
  alpha[alpha > alpha_ceil]  <- alpha_ceil

  j_star <- NULL
  if (!is.null(lambda_star) && is.finite(lambda_star)) {
    j_star <- which.min(abs(lambda_grid - lambda_star))
  }

  if (isTRUE(legend_panel)) {
    layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  } else {
    layout(1)
  }

  # Panel 1: matrix
  par(mar = c(6, 14, 3, if (isTRUE(legend_panel)) 1 else 2))
  plot(NA,
       xlim = c(0.5, n_lam + 0.5),
       ylim = c(0.5, n_taxa + 0.5),
       xlab = "", ylab = "", axes = FALSE,
       main = "Predicted class across \u03BB")

  for (i in seq_len(n_taxa)) {
    y <- n_taxa - i + 1
    for (j in seq_len(n_lam)) {
      cls <- pred_plot[i, j]
      col0 <- class_cols[cls]
      colA <- .apply_alpha(col0, alpha[i, j])
      rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5, col = colA, border = NA)
    }
  }

  if (!is.null(sep_every) && sep_every > 0) {
    y_lines <- seq(0.5, n_taxa + 0.5, by = sep_every)
    abline(h = y_lines, col = "white", lwd = 1.5)
  }

  axis(1, at = seq_len(n_lam), labels = lambda_grid, las = 2, cex.axis = 0.8)
  axis(2, at = seq_len(n_taxa), labels = rev(taxa_plot), las = 2, cex.axis = 0.55)

  if (!is.null(j_star)) {
    abline(v = j_star, lty = 2, lwd = 2)
    #mtext(paste0("\u03BB* = ", signif(lambda_star, 3)), side = 3, line = 0.2, adj = 1)
  }

  # Panel 2: legends
  if (isTRUE(legend_panel)) {
    par(mar = c(6, 0, 3, 2))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))

    legend("topleft",
           legend = names(class_cols),
           fill   = unname(class_cols),
           bty = "n",
           cex = 0.95)

    conf_levels <- c(0.4, 0.7, 0.95)
    conf_alpha <- alpha_floor + (alpha_ceil - alpha_floor) * conf_levels

    y0 <- 0.55
    text(0, y0, "Max posterior", adj = c(0, 1), cex = 0.9)

    swatch_col <- unname(class_cols)[min(2L, length(class_cols))]
    for (k in seq_along(conf_levels)) {
      y <- y0 - 0.08 - k * 0.12
      rect(0.02, y, 0.22, y + 0.08,
           col = .apply_alpha(swatch_col, conf_alpha[k]),
           border = NA)
      text(0.27, y + 0.04, sprintf("%.2f", conf_levels[k]), adj = c(0, 0.5), cex = 0.85)
    }
  }

  layout(1)

  invisible(list(
    pred_mat = pred_mat,
    conf_mat = conf_mat,
    lambda_grid = lambda_grid,
    taxa_plotted = taxa_plot
  ))
}


#' Plot discriminant variate scores for training and test taxa
#'
#' Produces a simple 2D scatter plot of discriminant variate scores (DV1 vs DV2
#' by default) for the training set and an optional set of test taxa.
#'
#' Training points are drawn first (optionally faded), then test points are drawn
#' on top for visibility.
#'
#' @param fit A fitted object of class `"phylo.fda"`.
#' @param grp_train Factor or vector giving the training group membership.
#' @param X_test Optional numeric matrix or data frame of test observations to
#'   project into discriminant space.
#' @param train_color Character string specifying how training points are
#'   colored: `"true"` for observed classes or `"predicted"` for fitted classes.
#' @param dims Length-2 integer vector indicating which variates to plot.
#' @param plot_type Character string specifying the plotting mode: `"auto"`,
#'   `"2d"`, or `"jitter"`.
#' @param test_subset Optional subset of test taxa: integer indices (recommended)
#'   or a character vector of row names.
#' @param class_cols Optional named vector mapping group level -> color.
#' @param prior Optional prior class probabilities passed to `predict.phylo.fda()`.
#' @param use_posterior_alpha Logical; if `TRUE`, scale test-point opacity by
#'   posterior probability.
#' @param alpha_floor Minimum opacity for posterior-scaled test points.
#' @param alpha_ceil Maximum opacity for posterior-scaled test points.
#' @param alpha_train Opacity for training points.
#' @param cex_train Point size for training points.
#' @param cex_test Point size for test points.
#' @param pch_train Plotting character for training points.
#' @param pch_test Plotting character for test points.
#' @param main Plot title.
#' @param ... Passed to base `plot()`.
#'
#' @return Invisibly returns a list with training/test variate matrices and
#'   class vectors used for coloring.
#' @export
plot_phylo_fda_scores <- function(fit,
                                  grp_train,
                                  X_test = NULL,
                                  train_color = c("true", "predicted"),
                                  dims = c(1, 2),
                                  plot_type = c("auto", "2d", "jitter"),
                                  test_subset = NULL,
                                  class_cols = NULL,
                                  prior = NULL,
                                  use_posterior_alpha = FALSE,
                                  alpha_train = 0.30,
                                  alpha_floor = 0.15,
                                  alpha_ceil  = 1.00,
                                  cex_train = 1.1,
                                  cex_test = 1.3,
                                  pch_train = 16,
                                  pch_test = 23,
                                  main = "phylo.fda discriminant variates",
                                  ...) {
  
  train_color <- match.arg(train_color)
  plot_type <- match.arg(plot_type)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  if (length(dims) != 2L) stop("dims must be a length-2 integer vector", call. = FALSE)
  
  grp_train <- as.factor(grp_train)
  class_cols <- .default_class_cols(class_cols, grp_train)
  
  # normalize prior if supplied (used only for posterior-based alpha)
  if (!is.null(prior)) prior <- .prior_from_vector_or_null(prior, grp_train)
  
  Z_train <- predict(fit, type = "variates")
  if (!is.matrix(Z_train)) Z_train <- as.matrix(Z_train)
  nv <- ncol(Z_train)
  
  # training classes (true or predicted)
  g_true <- droplevels(as.factor(grp_train))
  if (train_color == "true") {
    cls_train <- g_true
  } else {
    cls_train <- droplevels(as.factor(predict(fit, type = "class")))
    cls_train <- factor(cls_train, levels = levels(g_true))
  }
  
  # optional test
  Z_test <- NULL
  cls_test <- NULL
  X_test_used <- NULL
  if (!is.null(X_test)) {
    X_test_used <- as.matrix(X_test)
    
    if (!is.null(test_subset)) {
      if (is.numeric(test_subset)) {
        X_test_used <- X_test_used[test_subset, , drop = FALSE]
      } else {
        keep <- as.character(test_subset)
        keep <- keep[keep %in% rownames(X_test_used)]
        X_test_used <- X_test_used[keep, , drop = FALSE]
      }
    }
    
    Z_test <- predict(fit, newdata = X_test_used, type = "variates")
    if (!is.matrix(Z_test)) Z_test <- as.matrix(Z_test)
    cls_test <- as.factor(predict(fit, newdata = X_test_used, type = "class", prior = prior))
  }
  
  # decide plot mode
  do_jitter <- FALSE
  if (plot_type == "jitter") {
    do_jitter <- TRUE
  } else if (plot_type == "auto") {
    if (nv < 2L || dims[1] == dims[2]) do_jitter <- TRUE
  }
  
  # alpha for training (either constant, or derived from posterior confidence)
  if (isTRUE(use_posterior_alpha)) {
    post_train <- predict(fit, type = "posterior", prior = prior)
    conf_train <- apply(post_train, 1, max)
    a_train <- alpha_floor + (alpha_ceil - alpha_floor) * conf_train
    a_train[a_train < alpha_floor] <- alpha_floor
    a_train[a_train > alpha_ceil]  <- alpha_ceil
  } else {
    a_train <- rep(alpha_train, nrow(Z_train))
  }
  cols_train <- .apply_alpha(class_cols[as.character(cls_train)], a_train)
  
  # alpha for test points (posterior confidence)
  bg_cols <- NULL
  if (!is.null(Z_test)) {
    if (isTRUE(use_posterior_alpha)) {
      post_test <- predict(fit, newdata = X_test_used, type = "posterior", prior = prior)
      conf_test <- apply(post_test, 1, max)
      a_test <- alpha_floor + (alpha_ceil - alpha_floor) * conf_test
      a_test[a_test < alpha_floor] <- alpha_floor
      a_test[a_test > alpha_ceil]  <- alpha_ceil
      bg_cols <- .apply_alpha(class_cols[as.character(cls_test)], a_test)
    } else {
      bg_cols <- class_cols[as.character(cls_test)]
    }
  }
  
  if (isTRUE(do_jitter)) {
    
    par(mar = c(5, 8, 4, 2))
    
    dv <- dims[1]
    if (dv > nv) stop("Requested dims exceed available variates in fit", call. = FALSE)
    
    # x values
    x_train <- drop(Z_train[, dv])
    x_all <- x_train
    if (!is.null(Z_test)) x_all <- c(x_all, drop(Z_test[, dv]))
    
    pad <- 0.04
    xlim <- range(x_all, finite = TRUE)
    if (diff(xlim) == 0) xlim <- xlim + c(-1, 1)
    xlim <- xlim + c(-1, 1) * diff(xlim) * pad
    
    # --- TRUE group factor (force clean levels)
    g_true <- factor(grp_train)
    g_levels <- unique(as.character(g_true))  # preserve order of appearance
    g_true <- factor(g_true, levels = g_levels)
    
    # y bands
    jit <- 0.18
    y_base  <- as.numeric(g_true)
    y_train <- y_base + runif(length(y_base), -jit, jit)
    
    # test band
    y_test <- NULL
    if (!is.null(Z_test)) {
      y_test_center <- length(g_levels) + 1
      y_test <- y_test_center + runif(nrow(Z_test), -jit, jit)
    }
    
    # fixed equal spacing
    ymax <- length(g_levels) + if (!is.null(Z_test)) 1 else 0
    ylim <- c(0.5, ymax + 0.5)
    
    plot(x_train, y_train,
         xlim = xlim, ylim = ylim,
         yaxt = "n", ylab = "",
         xlab = paste0("DV", dv, " score"),
         pch = pch_train, cex = cex_train,
         col = cols_train,
         main = main,
         ...)
    
    # axis labels
    at   <- seq_len(length(g_levels))
    labs <- g_levels
    
    if (!is.null(Z_test)) {
      at   <- c(at, y_test_center)
      labs <- c(labs, "test")
    }
    print(g_levels)
    print(labs)
    axis(2, at = at, labels = labs, las = 2)
    
    if (!is.null(Z_test)) {
      points(drop(Z_test[, dv]), y_test,
             pch = pch_test, cex = cex_test, lwd = 1.4,
             bg = bg_cols,
             col = "black")
    }
    
    legend("topright", bty = "n",
           legend = names(class_cols),
           pch = pch_test, pt.bg = unname(class_cols), col = "black")
    
    return(invisible(list(
      Z_train = Z_train,
      Z_test  = Z_test,
      cls_train = cls_train,
      cls_test  = cls_test,
      class_cols = class_cols
    )))
  }
  
  # 2D scatter plot
  if (max(dims) > nv) stop("Requested dims exceed available variates in fit", call. = FALSE)
  
  # axis limits from both sets (padding)
  x_all <- Z_train[, dims[1]]
  y_all <- Z_train[, dims[2]]
  if (!is.null(Z_test)) {
    x_all <- c(x_all, Z_test[, dims[1]])
    y_all <- c(y_all, Z_test[, dims[2]])
  }
  pad <- 0.04
  xlim <- range(x_all, finite = TRUE)
  ylim <- range(y_all, finite = TRUE)
  xlim <- xlim + c(-1, 1) * diff(xlim) * pad
  ylim <- ylim + c(-1, 1) * diff(ylim) * pad
  
  plot(Z_train[, dims[1]], Z_train[, dims[2]],
       xlim = xlim, ylim = ylim,
       col = cols_train,
       pch = pch_train, cex = cex_train,
       xlab = paste0("DV", dims[1], " score"),
       ylab = paste0("DV", dims[2], " score"),
       main = main,
       ...)
  
  if (!is.null(Z_test)) {
    points(Z_test[, dims[1]], Z_test[, dims[2]],
           pch = pch_test, cex = cex_test, lwd = 1.4,
           bg = bg_cols,
           col = "black")
  }
  
  legend("topright", bty = "n",
         legend = names(class_cols),
         pch = pch_test, pt.bg = unname(class_cols), col = "black")
  
  invisible(list(Z_train = Z_train,
                 Z_test = Z_test,
                 cls_train = cls_train,
                 cls_test = cls_test,
                 class_cols = class_cols))
}


#' @noRd
#' @noRd
.subset_taxa_mats <- function(pred_mat, conf_mat, taxa = NULL, test_subset = NULL, taxa_idx = NULL) {
  if (!is.matrix(pred_mat) || !is.matrix(conf_mat)) stop("pred_mat and conf_mat must be matrices", call. = FALSE)
  if (!identical(dim(pred_mat), dim(conf_mat))) stop("pred_mat and conf_mat must have identical dimensions", call. = FALSE)

  # Backward compatible alias: taxa_idx
  if (!is.null(test_subset) && !is.null(taxa_idx)) {
    stop("Provide only one of 'test_subset' or 'taxa_idx'.", call. = FALSE)
  }
  if (is.null(test_subset) && !is.null(taxa_idx)) test_subset <- taxa_idx

  if (!is.null(taxa)) {
    taxa <- taxa[taxa %in% rownames(pred_mat)]
    return(list(pred = pred_mat[taxa, , drop = FALSE],
                conf = conf_mat[taxa, , drop = FALSE]))
  }

  if (!is.null(test_subset)) {
    return(list(pred = pred_mat[test_subset, , drop = FALSE],
                conf = conf_mat[test_subset, , drop = FALSE]))
  }

  list(pred = pred_mat, conf = conf_mat)
}


#' @noRd
.lambda_prediction_mats <- function(X_train,
                                   grp_train,
                                   tree_train,
                                   X_test,
                                   lambda_grid,
                                   prior,
                                   eqprior = FALSE,
                                   ...) {

  lambda_grid <- as.numeric(lambda_grid)
  if (!length(lambda_grid)) stop("lambda_grid must have length >= 1", call. = FALSE)
  if (any(!is.finite(lambda_grid))) stop("lambda_grid contains non-finite values", call. = FALSE)

  taxa_test <- rownames(X_test)
  if (is.null(taxa_test)) taxa_test <- paste0("taxon_", seq_len(nrow(X_test)))

  pred_mat <- matrix(NA_character_, nrow = length(taxa_test), ncol = length(lambda_grid),
                     dimnames = list(taxa_test, lambda_grid))
  conf_mat <- matrix(NA_real_, nrow = length(taxa_test), ncol = length(lambda_grid),
                     dimnames = list(taxa_test, lambda_grid))

  for (j in seq_along(lambda_grid)) {
    lam <- lambda_grid[j]

    fit_lam <- phylo.fda(
      data = X_train,
      grp  = grp_train,
      tree = tree_train,
      lambda = lam,
      eqprior = eqprior,
      ...
    )

    post <- predict(fit_lam, newdata = X_test, type = "posterior", prior = prior)
    conf <- apply(post, 1, max)

    cls <- predict(fit_lam, newdata = X_test, type = "class", prior = prior)

    pred_mat[, j] <- as.character(cls[taxa_test])
    conf_mat[, j] <- conf[taxa_test]
  }

  list(pred_mat = pred_mat, conf_mat = conf_mat)
}


#' @noRd
#.apply_alpha: robust alpha application without grDevices::adjustcolor quirks
#' @noRd
.apply_alpha <- function(cols, alpha) {
  cols <- as.character(cols)
  if (length(alpha) == 1L) alpha <- rep(alpha, length(cols))
  if (length(alpha) != length(cols)) {
    stop("'alpha' must have length 1 or the same length as 'cols'", call. = FALSE)
  }
  rgb <- grDevices::col2rgb(cols) / 255
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ], alpha = alpha)
}

.default_class_cols <- function(class_cols = NULL, grp_train) {
  lev <- levels(grp_train)
  if (is.null(class_cols)) {
    if (length(lev) == 3 && all(c("nocturnal", "cathemeral", "diurnal") %in% lev)) {
      class_cols <- c(nocturnal = "#1A1A1A", cathemeral = "#6BAED6", diurnal = "#F0E442")
    } else {
      cols <- grDevices::hcl.colors(length(lev), palette = "Okabe-Ito")
      class_cols <- stats::setNames(cols, lev)
    }
  }

  if (is.null(names(class_cols))) {
    if (length(class_cols) != length(lev)) {
      stop("class_cols must be named, or have length equal to nlevels(grp_train)", call. = FALSE)
    }
    names(class_cols) <- lev
  }

  missing_cols <- setdiff(lev, names(class_cols))
  if (length(missing_cols)) {
    stop("class_cols is missing colors for: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  class_cols[lev]
}


#' @noRd
.prior_from_vector_or_null <- function(prior, grp_train) {
  J <- nlevels(grp_train)
  lev <- levels(grp_train)

  if (is.null(prior)) return(rep(1 / J, J))

  if (!is.numeric(prior)) stop("prior must be numeric", call. = FALSE)
  if (length(prior) != J) stop("prior must have length equal to nlevels(grp_train)", call. = FALSE)

  if (is.null(names(prior))) {
    names(prior) <- lev
  }

  prior <- prior[lev]
  prior <- prior / sum(prior)
  prior
}
