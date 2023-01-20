vector_norm = function(x) {
  return(sqrt(sum(x^2)))
}

pdist2 = function(x,y) {
  dist.vec = matrix(nrow = nrow(x), ncol = ncol(x))
  for (i in c(1:nrow(x))) {
    dist.vec[i] = sqrt((x[i] - y)^2)
  }
  return(dist.vec)
}

wthresh = function(x,sorh,t){
  # WTHRESH Perform soft or hard thresholding.

  if (sorh == 's') {
    tmp = (abs(x) - t)
    tmp = (tmp+abs(tmp))/2
    y = sign(x) * tmp
  } else if (sorh == 'h') {
    y = x * (abs(x) > t)
  } else {
    print('Invalid argument value.')
  }
  return(y)
}

coordlasso = function(A, delta, lambda, tol, maxiter) {
  # Find solution

  # minimize 1/2 b'Ab - delta'b + m ||b||_1
  # hat b_l = S_m(delta_l - sum_{j != l} a_{lj} b_j)/a_{ll}
  # initial value

  p = nrow(A)

  b_ini = solve(A)%*%delta
  b_ini = b_ini/vector_norm(b_ini)

  bnew = b_ini

  converged = FALSE
  step = 0
  dA = as.matrix(diag(A))

  # obj_new = 1/2 *bnew'*A*bnew - delta'*bnew + lambda * sum(abs(bnew));

  while ((step < maxiter) & !converged) {
    bold = bnew
    step = step + 1

    for (j in c(1:p)) {
      x = delta[j] - A[j,]%*%bold + dA[j]*bold[j]

      x = wthresh(x,'s',lambda)
      x = x/dA[j]

      bold[j] = x
    }
    conv_criterion = vector_norm(bnew - bold)/vector_norm(bold)
    converged = conv_criterion < tol

    bnew = bold

    if (vector_norm(bnew) == 0) break

  }

  b = bnew

  if (vector_norm(b) != 0) {
    b = b / vector_norm(b)
  }

  return(b)
}

SDelta_binary = function(xtr0cv,xtr1cv) {
  # return covariance matrix and mean functions

  ntr0cv = nrow(xtr0cv)
  ntr1cv = nrow(xtr1cv)
  p = ncol(xtr1cv)

  S = (ntr0cv - 1)*cov(xtr0cv) + (ntr1cv - 1)*cov(xtr1cv)
  S = S/(ntr0cv + ntr1cv - 2)
  mu0 = as.matrix(colMeans(xtr0cv))
  mu1 = as.matrix(colMeans(xtr1cv))

  return(list('S' = S, 'mu0' = mu0, 'mu1' = mu1))
}

SDelta_multi = function(xtr) {
  # Compute pooled S

  nclass = length(xtr)

  p = ncol(xtr[[1]])

  ntr = c()

  S = matrix(0,nrow = p, ncol= p)

  for (n in seq(nclass)){

    S = S + (nrow(xtr[[n]]) - 1)*cov(xtr[[n]])

    ntr = c(ntr, nrow(xtr[[n]]))

  }

  S = S/(sum(ntr) - nclass)

  # Compute X_bar

  mu_list = vector(mode = "list", length = nclass)

  for (n in seq(nclass)){
    mu_list[[n]] = as.matrix(colMeans(xtr[[n]]))
  }

  # Compute cov(t(A)) (pxp matrix)

  A = do.call(cbind,mu_list)
  A_cov = cov(t(A))

  # Eigenvalue Decomposition and pick first nclass-1 eigenvectors as deltas

  eigen_A_cov = eigen(A_cov)

  delta_list = vector(mode = "list", length = nclass - 1)

  for (n in seq(nclass - 1)){
    delta_list[[n]] = eigen_A_cov$vectors[,n]
  }

  return(list('S' = S, 'delta_list' = delta_list))
}

is_constant = function(xb){

  const = c()

  for (n in seq(ncol(xb))){
    const = c(const, var(xb[,n]))
  }

  if (sum(const) == 0){
    return(NULL)
  } else{
    return(xb[,which(const != 0), drop = FALSE])
  }
}
