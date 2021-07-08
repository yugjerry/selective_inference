require(palmerpenguins)
require(ggplot2)
options(ggplot2.discrete.colour=list(RColorBrewer::brewer.pal(6, "Dark2")))

dat <- penguins[complete.cases(penguins), ]
dat <- dat[dat$sex == "female", c(1, 3, 5)]

ggplot(dat) + geom_point(aes(x=flipper_length_mm, y = bill_length_mm, 
                             shape=as.factor(species)), size = 3, fill="grey", colour="black") + 
  scale_shape_manual(name="Species", values=c(21, 24, 22)) + 
  ylab("Bill length (mm)") + xlab("Flipper length (mm)") + coord_fixed() + 
  theme_bw(base_size=22) + ggtitle("Penguins") + theme(legend.position="right")

km_cluster <- function(X) { 
  km <- kmeans(X, 3, nstart=50)
  return(km$cluster)
}

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

sortE <- function(E) {
  E.sorted <- lapply(1:nrow(E), function(i){
    temp <- as.numeric(E[i, ])
    if (temp[1] <= 0 & temp[2] <= 0) {
      return(sort(-temp))
    }
    if (temp[1] >= 0 & temp[2] >= 0) {
      return(sort(temp))
    }
    # we know temp[1] < 0, temp[2] > 0 OR temp[1] > 0, temp[2] < 0
    temp <- abs(temp)
    return(rbind(c(0, temp[1]), c(0, temp[2])))
  })
  E.sorted <- do.call(rbind, E.sorted)
  # in order to use the approximation, we translate Inf to a large number
  return(finiteE(E.sorted))
}

finiteE <- function(E) {
  ind.inf <- which(E == Inf)
  if (length(ind.inf) == 0) return(E)
  # we know there are some infinite entries
  E.max <- max(E[-ind.inf])
  E[which(E == Inf)] <- max(10000, E.max * 2)
  return(E)
}

isSameIntervals <- function(int1, int2) {
  
  # first make int1, int2 to the default order
  int1 <- intervals::reduce(int1)
  int2 <- intervals::reduce(int2)
  
  if (nrow(int1) != nrow(int2)) return(FALSE)
  
  # int1 and int2 has the same number of intervals
  
  if (sum(int1 != int2) > 0) return(FALSE)
  
  # int1 and int2 has the same elements
  return(TRUE)
}

TNRatioApprox <- function(E1, E2, scale = NULL) {
  
  if (is.null(scale)) {
    temp <- (c(E1, E2))^2
    scale.grid <- stats::quantile(temp, probs = seq(0, 1, 0.2))
    
    for(scale in scale.grid) {
      temp <- TNRatioApprox(E1, E2, scale = scale)
      if (!is.na(temp)) {
        return(temp)
      }
      # if temp is NaN, proceed to the next loop
    }
    
    # if all scale.grid does not work, then return NaN
    return(NaN)
  }
  num1 <- magicfun(E1[, 1]) * exp(-(E1[, 1]^2 - scale)/2)
  num2 <- magicfun(E1[, 2]) * exp(-(E1[, 2]^2 - scale)/2)
  denom1 <- magicfun(E2[, 1]) * exp(-(E2[, 1]^2 - scale)/2)
  denom2 <- magicfun(E2[, 2]) * exp(-(E2[, 2]^2 - scale)/2)
  res <- sum(num1-num2)/sum(denom1-denom2)
  return(res)
}

TNProbEachInt <- function(lo, up) {
  if (up == Inf) {
    return(stats::pnorm(lo, 0, 1, lower.tail = FALSE))
  }
  # we know up < Inf, want P(lo <= X <= up).
  # we want small - small (big - big will mask small numbers),
  
  try1 <- stats::pnorm(lo, 0, 1, lower.tail = FALSE) - stats::pnorm(up, 0, 1, lower.tail = FALSE)
  if (try1 != 0) return(try1)
  
  try2 <- stats::pnorm(up, 0, 1, lower.tail = TRUE) - stats::pnorm(lo, 0, 1, lower.tail = TRUE)
  return(try2)
  
}

TNProb <- function(E) {
  # sum cdf over each disjoint interval of E
  res <- sum(sapply(1:nrow(E), function(v) {
    return(TNProbEachInt(E[v, 1], E[v, 2]))
  }))
  return(res)
}

TNSurv <- function(q, mean, sd, E, approx = FALSE) {
  # check if truncation is empty (i.e. truncated to the empty set)
  if (nrow(E) == 0) {
    stop("truncation set is empty")
  }
  
  # check if truncation is the whole real line
  if (isSameIntervals(E, intervals::Intervals(c(-Inf, Inf)))) {
    return(stats::pnorm(q, mean, sd, lower.tail = FALSE))
  }
  
  # E is not empty and is not the whole real line,
  # i.e. 0 < P(X in E) < 1
  
  # we want P(X > q | X in E) = P(X >= q AND X in E) / P(X in E)
  # {X >= q} = {Z >= (q-mean)/sd}
  # {X in E} = {Z in (E-mean)/sd}
  # Z ~ N(0, 1)
  q <- (q-mean)/sd
  E <- (E-mean)/sd
  mean <- 0
  sd <- 1
  q2 <- q*q
  region <- suppressWarnings(intervals::interval_intersection(E, intervals::Intervals(c(q, Inf))))
  # check if the result is 0 or 1
  if(nrow(region) == 0) return(0)
  if (isSameIntervals(E, region)) return(1)
  
  # transform region and E so that intervals have positive endpoints
  region <- sortE(region)
  E <- sortE(E)
  
  # we want P(Z in region) / P(Z in E)
  # try approximate calculation
  if (approx) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }
  
  # try exact calculation
  denom <- TNProb(E)
  num <- TNProb(region)
  
  if (denom < 1e-100 || num < 1e-100) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }
  
  # we know denom and num are both reasonably > 0
  
  res <- num / denom
  # force the result to lie in [0, 1]
  return(max(0, min(1, res)))
}

magicfun = function(z){
  z2 <- z*z
  z3 <- z*z*z
  temp <- (z2 + 5.575192695 * z + 12.77436324) /
    (sqrt(2*pi) * z3 + 14.38718147*z2 + 31.53531977*z + 2*12.77436324)
  return(temp)
}

result <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

delta<-function(a,b,c){
  b^2-4*a*c
}

TChisqRatioApprox <- function(df, E1, E2) {
  
  # the transform that makes x into a N(0, 1) r.v. such that
  # P(X >= x) = P(Z >= Chisq2N(x)), X ~ chisq(df), Z ~ N(0, 1)
  # this function can take either scaler, vector or matrix
  Chisq2N <- function(x, df, tol = 1e-6) {
    
    if (is.numeric(x) && length(x) == 1) {
      if (x <= tol) { # x <= 0
        return(-Inf)
      }
      if (x == Inf) {
        return(Inf)
      }
      # we know x > 0 and x is finite
      x <- (x/df)^(1/6) - (1/2) * (x/df)^(1/3) + (1/3) * (x/df)^(1/2)
      mu <- 5/6 - 1/(9*df) - 7/(648*df^2) + 25/(2187*df^3)
      sig <- sqrt(1/(18*df) + 1/(162*df^2) - 37/(11664*df^3))
      return((x-mu)/sig)
    }
    
    if (is.vector(x)) {
      return(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df))
    }
    
    if (is.matrix(x)) {
      return(structure(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df), dim = dim(x)))
    }
    
    return(intervals::Intervals())
    
  }
  
  
  E1 <- Chisq2N(E1, df)
  E1 <- sortE(E1) # notice that Chisq2N can be negative
  E2 <- Chisq2N(E2, df)
  E2 <- sortE(E2)
  
  # now we want P(Z in E1) / P(Z in E2), Z ~ N(0, 1)
  return(TNRatioApprox(E1, E2))
}

cond_const <- function(X, cl, k1, k2){
  k=length(unique(cl))
  tol=1e-10
  centr = NULL
  for(j in 1:k){
    centr[j] = mean(X[cl==j,])
  }
  root_mat = NULL
  a1 = prop_k2*diff_means/stat
  a2 = (prop_k2 - 1)*diff_means/stat
  c1 = prop_k2*diff_means
  c2 = (prop_k2 - 1)*diff_means
  # type 1, 2
  for(j in setdiff((1:k),c(k1,k2))){
    for(i in (cl==k1)){
      l1 = norm_vec(a1)^2
      l2 = sum(a1*(X[i,]-c1-centr[j]))
      l3 = norm_vec(X[i,]-c1-centr[j])^2 - norm_vec(X[i,]-centr[k1])^2
      res = result(l1,l2,l3)
      if(length(res)==1){
        root_mat = rbind(root_mat, rep(res,2))
      }
      if(length(res)==2){
        root_mat = rbind(root_mat, res)
      }
    }
    for(i in (cl==k2)){
      l1 = norm_vec(a2)^2
      l2 = sum(a2*(X[i,]-c2-centr[j]))
      l3 = norm_vec(X[i,]-c2-centr[j])^2 - norm_vec(X[i,]-centr[k2])^2
      res = result(l1,l2,l3)
      if(length(res)==1){
        root_mat = rbind(root_mat, rep(res,2))
      }
      if(length(res)==2){
        root_mat = rbind(root_mat, res)
      }
    }
  }
  # type 3
  for(i in (cl==k1)){
    l1 = norm_vec(a1-a2)^2
    l2 = sum((a1-a2)*(X[i,]-c1-centr[k2]))
    l3 = norm_vec(X[i,]-c1-centr[k2])^2 - norm_vec(X[i,]-centr[k1])^2
    res = result(l1,l2,l3)
    if(length(res)==1){
      root_mat = rbind(root_mat, rep(res,2))
    }
    if(length(res)==2){
      root_mat = rbind(root_mat, res)
    }
  }
  # type 4
  for(i in (cl==k2)){
    l1 = norm_vec(a2-a1)^2
    l2 = sum((a2-a1)*(X[i,]-c2-centr[k1]))
    l3 = norm_vec(X[i,]-c2-centr[k1])^2 - norm_vec(X[i,]-centr[k2])^2
    res = result(l1,l2,l3)
    if(length(res)==1){
      root_mat = rbind(root_mat, rep(res,2))
    }
    if(length(res)==2){
      root_mat = rbind(root_mat, res)
    }
  }
  root1 = min(root_mat[,1])
  root2 = max(root_mat[,2])
  
  if(root2 < -tol) {
    return(c(0, Inf))
  }
  
  if(root1 > tol) {
    return(c(0, root1, root2, Inf))
  }
  
  return(c(root2, Inf))
}

test_clusters_approx <- function(X, k1, k2, iso=TRUE, sig=NULL, SigInv=NULL, ndraws=2000, cl_fun, cl=NULL) {
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
  log_survives <- rep(NA, ndraws)
  phi <- stats::rnorm(ndraws)*scale_factor + stat
  
  
  k1_constant <- prop_k2*diff_means/stat
  k2_constant <- (prop_k2 - 1)*diff_means/stat
  orig_k1 <- t(X[cl == k1, ])
  orig_k2 <- t(X[cl == k2, ])
  
  Xphi <- X
  
  for(j in 1:ndraws) {
    if(phi[j] < 0) next
    
    # Compute perturbed data set
    Xphi <- X
    Xphi[cl == k1, ] <- t(orig_k1 + (phi[j] - stat)*k1_constant)
    Xphi[cl == k2, ] <- t(orig_k2 + (phi[j] - stat)*k2_constant)
    
    # Recluster the perturbed data set
    cl_Xphi <- cl_fun(Xphi)
    if(preserve_cl(cl, cl_Xphi, k1, k2)) {
      log_survives[j] <- -(phi[j]/scale_factor)^2/2 + (q-1)*log(phi[j]/scale_factor) - (q/2 - 1)*log(2) - log(gamma(q/2)) - log(scale_factor) -
        stats::dnorm(phi[j], mean=stat, sd=scale_factor, log=TRUE)
    }
  }
  
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
  
  return(list(stat=stat, pval=pval, stderr=sqrt(var_pval), clusters=cl))
}

test_exact <- function(X, K, k1, k2, iso=TRUE, sig=NULL, SigInv=NULL) {
  # error checking 
  if(!is.matrix(X)) stop("X should be a matrix")
  
  n <- nrow(X)
  q <- ncol(X)
  
  # hierarchical clustering with squared Euclidean distance and specified linkage
  hcl_at_K <- kmeans(X, K, nstart=50)
  hcl_at_K <- hcl_at_K$cluster
  
  n1 <- sum(hcl_at_K == k1)
  n2 <- sum(hcl_at_K == k2)
  squared_norm_nu <- 1/n1 + 1/n2
  diff_means <- colMeans(X[hcl_at_K == k1, , drop=FALSE]) - colMeans(X[hcl_at_K == k2, , drop=FALSE])
  
  if(iso) {
    if(is.null(sig)) {
      sig <- sqrt(sum(scale(X, scale=FALSE)^2)/(n*q - q))
    }
    
    # compute test statistic
    stat <- norm_vec(diff_means)
    
    # compute truncation set
    S = cond_const(X, hcl_at_K, k1, k2)
    S = intervals::Intervals(S)
    
    # set distribution of phi
    scale_factor <- squared_norm_nu*sig^2
    
  } else {
    if(is.null(SigInv)) {
      Sig <- stats::cov(scale(X, scale=F))
      SigInv <- solve(Sig)
    }
    
    # compute test statistic
    stat <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
    
    # compute truncation set
    S = cond_const(X, hcl_at_K, k1, k2)
    S = intervals::Intervals(S)
    
    # set distribution of phi
    scale_factor <- squared_norm_nu
  }
  
  # compute p-value using truncated chi-squared distribution
  gestat <- intervals::Intervals(c(stat^2/scale_factor, Inf))
  denom <- S^2/scale_factor
  numer <- suppressWarnings(intervals::interval_intersection(gestat, denom))
  pval <- TChisqRatioApprox(q, numer, denom)
  
  return(list(stat=stat, pval=pval, trunc=S))
}



X <- as.matrix(dat[, -c(1)]) # remove species and convert to matrix
set.seed(123) 
cl <- km_cluster(X)
ggplot(dat) + geom_point(aes(x=flipper_length_mm, y = bill_length_mm, 
                             shape=as.factor(species), fill=as.factor(cl)), size = 3, colour="black") + 
  scale_fill_discrete(name="Clusters", guide=guide_legend(ncol=2, override.aes=list(shape=21))) + 
  scale_shape_manual(name="Species", values=c(21, 24, 22), guide=guide_legend(override.aes=list(fill="black"))) +
  ylab("Bill length (mm)") + xlab("Flipper length (mm)") + coord_fixed() + 
  theme_bw(base_size=22) + ggtitle("Penguins") + theme(legend.position="right") 

# test_clusters_approx(X, k1=1, k2=3, cl_fun=km_cluster, cl=cl, ndraws=10000)
test_exact(X, K=3, k1=1, k2=3, iso=TRUE, sig=NULL, SigInv=NULL)

