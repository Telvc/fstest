# function to compute the selective type I error

test_clusters_approx_tensor <- function(X, k1, k2, iso=TRUE, sig=NULL, SigInv=NULL, ndraws=2000, cl_fun, cl=NULL) {
  if(!is.matrix(X)) stop("X should be a matrix")
  n <- nrow(X)
  q <- ncol(X)

  if(is.null(cl)) cl <- cl_fun(X)
  K <- length(unique(cl))

  if(!is_integer_between_a_b(K, 2, n)) stop("number of clusters (K) should be between 2 and n")
  if(!is_integer_between_a_b(k1, 1, K) | !is_integer_between_a_b(k2, 1, K)) stop(paste("cluster indices should be between 1 and K", sep=""))
  if((iso != TRUE) & (iso != FALSE)) stop("iso should be TRUE or FALSE")

  n1 <- sum(cl == k1)
  n2 <- sum(cl == k2)
  squared_norm_nu <- 1/n1 + 1/n2
  diff_means <- colMeans(X[cl == k1, , drop=FALSE]) - colMeans(X[cl == k2, , drop=F])

  prop_k2 <- n2/(n1+n2)


  if(iso) {
    if(is.null(sig)) {
      sig <- sqrt(sum(scale(X, scale=FALSE)^2)/(n*q - q))
    }

    scale_factor <- squared_norm_nu*sig^2
    # compute test statistic
    stat <- norm_vec(diff_means)
  } else {
    if(is.null(SigInv)) {
      Sig <- stats::cov(scale(X, scale=FALSE))
      SigInv <- solve(Sig)
    }

    scale_factor <- squared_norm_nu

    # compute test statistic
    stat <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
  }

  scale_factor <- sqrt(scale_factor)
  phi <- stats::rnorm(ndraws)*scale_factor + stat


  k1_constant <- prop_k2*diff_means/stat
  k2_constant <- (prop_k2 - 1)*diff_means/stat
  orig_k1 <- t(X[cl == k1, ])
  orig_k2 <- t(X[cl == k2, ])

  Xphi <- X

  log_survives <- unlist(future.apply::future_lapply(X = 1:ndraws, FUN = function(j) {
    if(phi[j] < 0) return(NA)

    # Compute perturbed data set
    Xphi <- X
    Xphi[cl == k1, ] <- t(orig_k1 + (phi[j] - stat)*k1_constant)
    Xphi[cl == k2, ] <- t(orig_k2 + (phi[j] - stat)*k2_constant)

    # Recluster the perturbed data set
    cl_Xphi <- cl_fun(Xphi)
    if(preserve_cl(cl, cl_Xphi, k1, k2)) {
      log_survives <- -(phi[j]/scale_factor)^2/2 + (q-1)*log(phi[j]/scale_factor) - (q/2 - 1)*log(2) - lgamma(q/2) - log(scale_factor) -
        stats::dnorm(phi[j], mean=stat, sd=scale_factor, log=TRUE)
      return(log_survives)
    }

    return(NA)

  }, future.seed=TRUE))

  # Trim down to only survives
  phi <- phi[!is.na(log_survives)]
  log_survives <- log_survives[!is.na(log_survives)]

  survives <- length(log_survives)

  # Return nothing if nothing survives
  if(survives == 0) {
    warning("Oops - we didn't generate any samples that preserved the clusters! Try re-running with a larger value of ndraws.")
    return(list(stat=stat, pval=NA, stderr=NA, clusters=cl))
  }

  #  Approximate p-values
  log_survives_shift <- log_survives - max(log_survives)
  props <- exp(log_survives_shift)/sum(exp(log_survives_shift))
  pval <- sum(props[phi >= stat])

  var_pval <- (1 - pval)^2*sum(props[phi >= stat]^2) + pval^2*sum(props[phi < stat]^2)

  return(list(stat=stat, pval=pval, stderr=sqrt(var_pval), clusters=cl, n_survive = survives, phi = phi))
}

# ----- general purpose helper functions -----

norm_vec <- function(x) {
  sqrt(sum(x^2))
}

is_integer_between_a_b <- function(x, a, b) {
  (x>= min(c(a, b))) && (x %% 1 == 0) && (x <= max(c(a, b)))
}

same_cl <- function(cl1, cl2, K) {
  tab <- table(cl1, cl2)
  sum(tab != 0) == K
}

preserve_cl <- function(cl, cl_phi, k1, k2) {
  tab <- table(cl, cl_phi)

  k1_in <- (sum(tab[k1, ] != 0) == 1) & (sum(tab[, k1] != 0) == 1)
  k2_in <- (sum(tab[k2, ] != 0) == 1) & (sum(tab[, k2] != 0) == 1)

  k1_in & k2_in
}
