#' @title Nonparametric partial coherence matrix estimation
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
nonpar.pcn <- function(object, x, eta = NULL, symmetric = c('min', 'max', 'avg', 'none'), 
                       pcn.cv.args = list(n.folds = 1, path.length = 10),
                       n.cores = min(parallel::detectCores() - 1, 3)){
  
  xx <- x - object$mean.x
  p <- dim(x)[1]
  
  symmetric <- match.arg(symmetric, c('min', 'max', 'avg', 'none'))
  GG <- Re(object$spec$Sigma_i[,, 1])

  if(is.null(eta)){
    dcv <- direct.cv(object, xx, target = 'spec', path.length = pcn.cv.args$path.length, n.folds = pcn.cv.args$n.folds,
                     q = object$q, kern.bandwidth.const = object$kern.bandwidth.const, n.cores = n.cores)
    eta <- dcv$eta
  } 
  DD <- direct.inv.est(GG, eta = eta, symmetric = symmetric, n.cores = n.cores)$DD
  out <- list(Omega = DD, eta = eta)
  attr(out, 'class') <- 'fnets.pcn'
  
  return(out)
  
}

#' @title Parametric partial coherence matrix estimation
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
param.pcn <- function(object, x, eta = NULL, symmetric = c('min', 'max', 'avg', 'none'), 
                      pcn.cv.args = list(n.folds = 1, path.length = 10),
                      n.cores = min(parallel::detectCores() - 1, 3)){
  
  xx <- x - object$mean.x
  p <- dim(x)[1]
  
  symmetric <- match.arg(symmetric, c('min', 'max', 'avg', 'none'))
  GG <- object$idio.var$Gamma
  A <- t(object$idio.var$beta)
  d <- dim(A)[2]/p
  
  A1 <- diag(1, p)
  for(ll in 1:d) A1 <- A1 - A[, (ll - 1) * p + 1:p]
  
  if(is.null(eta)){
    dcv <- direct.cv(object, xx, target = 'acv', 
                     path.length = pcn.cv.args$path.length, n.folds = pcn.cv.args$n.folds,
                     q = object$q, kern.bandwidth.const = object$kern.bandwidth.const, n.cores = n.cores)
    eta <- dcv$eta
  } 
  Delta <- direct.inv.est(GG, eta = eta, symmetric = symmetric, n.cores = n.cores)$DD
  Omega <- 2 * pi * t(A1) %*% Delta %*% A1
  
  out <- list(Delta = Delta, Omega = Omega, eta = eta)
  attr(out, 'class') <- 'fnets.pcn'
  return(out)
  
}

#' @title Cross-validation for the constrained l1-minimisation problem for inverse matrix estimation
#' @export
direct.cv <- function(object, xx, target = c('spec', 'acv'), symmetric = c('min', 'max', 'avg', 'none'),
                      path.length = 10, n.folds = 1, q = 0, kern.bandwidth.const = 4, n.cores = min(parallel::detectCores() - 1, 3)){
  
  n <- ncol(xx)
  p <- nrow(xx)
  
  target <- match.arg(target, c('spec', 'acv'))
  if(target == 'spec'){
    GG <- Re(object$spec$Sigma_i[,, 1])
    eta.max <- max(abs(GG)) 
    eta.path <- round(exp(seq(log(eta.max), log(eta.max * .1), length.out = path.length)), digits = 10)
  }
  if(target == 'acv'){
    A <- t(object$idio.var$beta)
    d <- dim(A)[2]/p
    GG <- object$idio.var$Gamma
    eta.max <- max(abs(GG)) 
    eta.path <- round(exp(seq(log(eta.max), log(eta.max * .01), length.out = path.length)), digits = 10)
  }
  
  
  cv.err <- rep(0, length = path.length)
  ind.list <- split(1:n, ceiling(n.folds*(1:n)/n)) 
  for(fold in 1:n.folds){ 
    train.ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
    train.x <- xx[, ind.list[[fold]][train.ind]]
    test.x  <- xx[, ind.list[[fold]][- train.ind]] 
    if(target == 'spec'){
      train.GG <- Re(dyn.pca(train.x, q = q, kern.bandwidth.const = kern.bandwidth.const)$spec$Sigma_i[,, 1])
      test.GG <- Re(dyn.pca(test.x, q = q, kern.bandwidth.const = kern.bandwidth.const)$spec$Sigma_i[,, 1])
    }
    if(target == 'acv'){
      train.G0 <- dyn.pca(train.x, q = q, kern.bandwidth.const = kern.bandwidth.const)$acv$Gamma_i
      test.G0 <- dyn.pca(test.x, q = q, kern.bandwidth.const = kern.bandwidth.const)$acv$Gamma_i
      train.GG <- train.G0[,, 1]
      test.GG <- test.G0[,, 1]
      for(ll in 1:d){
        train.GG <- train.GG - A[, (ll - 1) * p + 1:p] %*% train.G0[,, ll + 1]
        test.GG <- test.GG - A[, (ll - 1) * p + 1:p] %*% test.G0[,, ll + 1]
      }
    }
    # sv <- svd(test.GG)
    # test.GG0 <- sv$u %*% diag((sv$d - 1e-10) * (sv$d > 1e-10) + 1e-10) %*% t(sv$u)
    for(ii in 1:path.length){
      DD <- direct.inv.est(train.GG, eta = eta.path[ii], symmetric = symmetric, n.cores = n.cores)$DD
      # sv <- svd(DD)
      # DD <- sv$u %*% diag((sv$d - 1e-10) * (sv$d > 1e-10) + 1e-10) %*% t(sv$u)
      DG <- DD %*% test.GG
      sv <- svd(DG, nu = 0, nv = 0)
      cv.err[ii] <- cv.err[ii] + sum(sv$d) - sum(log(sv$d)) - p # sum(diag(DD %*% test.GG)) - log(det(DD)) # sum(diag(DD %*% test.GG)) - log(det(DD %*% test.GG)) - p
    }
  }
  
  cv.err 
  eta.min <- eta.path[which.min(cv.err)]

  plot(eta.path, cv.err, type = 'b', col = 1, pch = 1, log = 'x', xlab = 'eta (log scale)', ylab = 'CV error')
  abline(v = eta.min)
  rp <- recordPlot()
  
  out <- list(eta = eta.min, cv.error = cv.err, eta.path = eta.path)
  return(out) 
  
}

#' @internal
direct.inv.est <- function(GG, eta = NULL, symmetric = c('min', 'max',  'avg', 'none'), n.cores = min(parallel::detectCores() - 1, 3)){
  
  p <- dim(GG)[1]
  f.obj <- rep(1, 2 * p)
  f.con <- rbind(-GG, GG)
  f.con <- cbind(f.con,-f.con)
  f.dir <- rep('<=', 2 * p)
  
  cl <- parallel::makePSOCKcluster(n.cores)
  doParallel::registerDoParallel(cl)
  
  DD <- foreach::foreach(ii = 1:p, .combine = 'cbind', .multicombine = TRUE, .export = c('lp')) %dopar% {
    ee <- rep(0, p)
    ee[ii] <- 1
    b1 <- rep(eta, p) - ee
    b2 <- rep(eta, p) + ee
    f.rhs <- c(b1, b2)
    lpout <- lp('min', f.obj, f.con, f.dir, f.rhs)
    lpout$solution[1:p] - lpout$solution[-(1:p)]
  }
  parallel::stopCluster(cl)
  
  DD <- make.symmetric(DD, symmetric)
  
  out <- list(DD = DD, eta = eta, symmetric = symmetric)
  return(out)
  
}

#' @internal
make.symmetric <- function(DD, symmetric){
  symmetric <- match.arg(symmetric, c('min', 'max', 'avg', 'none'))
  if(symmetric != 'none'){
    if(symmetric == 'min'){
      DD[abs(DD) > t(abs(DD))] <- 0
      DD <- DD + t(DD)
      diag(DD) <- diag(DD)/2
    }
    if(symmetric == 'max'){
      DD[abs(DD) < t(abs(DD))] <- 0
      DD <- DD + t(DD)
      diag(DD) <- diag(DD)/2
    }
    if(symmetric == 'avg') DD <- (DD + t(DD))/2
  }
  DD
}


#' @title Plot fnets.pcn object
#' @method plot fnets.pcn
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
plot.fnets.pcn <- function(object, names = NULL, groups = NULL, threshold = 0, size = NULL, ...){
  O <- object$Omega
  p <- dim(O)[1]
  
  mark.groups <- list()
  perm <- c()
  if(!is.null(groups)){
    for (ii in unique(groups)) {
      permii <- which(groups == ii)
      mark.groups[[ii]] <- permii
      perm <- c(perm, permii)
    }
    # perm <- Matrix::invPerm(perm)
  } else perm <- 1:p
  # Granger causal network
  
  # Long-Run Partial Correlation network
  O <- O[perm,perm]
  omega <- igraph::graph_from_adjacency_matrix(O, mode = "undirected", weighted=TRUE, diag = FALSE) 
  if(!is.null(size)) {degreeO <- igraph::degree(omega, mode = size, normalized = F)
  degreeO <- degreeO/max(degreeO)*15 }else degreeO <- rep(15, p)
  if(!is.null(names)) V(omega)$name <- names[perm]
  l_omega <- igraph::layout_in_circle(omega) 

  # Contemporaneous
  if(!is.null(object$Delta)){
   contemp <- igraph::graph_from_adjacency_matrix(object$Delta, mode = "undirected", weighted=TRUE, diag = FALSE)
   if(!is.null(size)) {degreeC <- igraph::degree(contemp, mode = size, normalized = F)
   degreeC <- degreeC/max(degreeC)*15 }else degreeC <- rep(15, p)
   if(!is.null(names)) V(contemp)$name <- names[perm]
   l_contemp <- igraph::layout_in_circle(contemp) 

   par(mfrow = c(1,2))
   plot.igraph(contemp, main = "Contemporaneous", vertex.label = V(contemp)$name, layout = l_contemp, vertex.label.font = 2,
               vertex.label.color = "black", edge.color = "gray40", edge.arrow.size = 0.1, vertex.size = degreeC,
               vertex.shape ="circle", vertex.color = groups[perm], vertex.label.cex = 0.8)
  } else par(mfrow = c(1,1))
  if(!is.null(names)) V(omega)$name <- names 
  plot.igraph(omega, main = "Long-Run Partial Correlation", vertex.label = V(omega)$name, layout = l_omega, vertex.label.font = 2,
              vertex.label.color = "black", edge.color = "gray40", edge.arrow.size = 0.1, vertex.size = degreeO,
              vertex.shape ="circle", vertex.color = groups[perm], vertex.label.cex = 0.8)
}
