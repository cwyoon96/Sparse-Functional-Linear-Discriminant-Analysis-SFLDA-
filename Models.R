source('common_func.R')
library(pracma)
library(MASS)

SFLDA = function(data, y, tau, lambda){
  # train SFLDA model with given data

  # data should be a list of n x p matrices (list of length 1 for univariate SFLDA)

  # y should be a vector of classes where order matches with the data

  # Check whether data and y is in right form

  if (!is.list(data)){
    stop("data should be given as list of matrices")
  }

  if (length(data) == 0){
    stop('list is empty')
  }

  if (!is.vector(y)){
    stop("y should be given as a vector")
  }

  # divide each case

  y = as.character(y)

  nclass = length(unique(y))

  if (nclass == 1){
    stop('number of class equals 1')
  } else if(nclass == 2){
    if (length(data) == 1){
      # binary univariate case

      result = binary_univariate_SFLDA(data,y, tau, lambda)
    }
    else{
      # binary multivariate case

      result = binary_multivariate_SFLDA(data,y, tau, lambda)
    }
  } else{
    if (length(data) == 1){
      # multiclass univariate case

      result = multiclass_univariate_SFLDA(data,y, tau, lambda)
    } else{
      # multiclass multivariate case

      result = multiclass_multivariate_SFLDA(data,y, tau, lambda)
    }
  }

  return(result)

}

binary_univariate_SFLDA = function(data, y, tau, lambda){
  data = data[[1]]

  y_class = unique(y)

  xtr0 = data[y == y_class[1],]

  xtr1 = data[y == y_class[2],]

  SD = SDelta_binary(xtr0, xtr1)

  S = SD$S
  mu0 = SD$mu0
  mu1 = SD$mu1

  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  delta = (mu0 - mu1)
  delta = delta / vector_norm(delta)
  S = S / norm(S, type = 'F')

  b = coordlasso((1-tau)*S + tau*DD, delta, lambda, 1e-8, 500)

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = 'xb1'

  if (is.null(is_constant(xb))){

    md = table(y)

    constant = TRUE

  } else {
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }

  return(list('type' = 'bu','beta' = b, 'y_class' = y_class, 'md'=md, 'constant' = constant))
}

binary_multivariate_SFLDA = function(data, y, tau, lambda){
  n_mat = length(data)

  p_list = c()

  for (i in seq(n_mat)){
    p_list = c(p_list, ncol(data[[i]]))
  }

  data = do.call(cbind, data)

  y_class = unique(y)

  xtr0 = data[y == y_class[1],]

  xtr1 = data[y == y_class[2],]

  SD = SDelta_binary(xtr0, xtr1)

  S = SD$S
  mu0 = SD$mu0
  mu1 = SD$mu1

  p = length(mu0)
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  s = 0

  for (j in p_list[-n_mat]){
    D[s +j,] = rep(0,p)

    s = s + j
  }

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  delta = (mu0 - mu1)
  delta = delta / vector_norm(delta)
  S = S / norm(S, type = 'F')

  b = coordlasso((1-tau)*S + tau*DD, delta, lambda, 1e-8, 500)

  beta = vector('list', length = n_mat)

  s = 1

  for (j in seq(n_mat)){

    beta[[j]] = b[seq(s,s + p_list[j]-1),,drop = FALSE]

    s = s + p_list[j]

  }

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = 'xb1'

  if (is.null(is_constant(xb))){

    md = table(y)

    constant = TRUE

  } else {
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }


  return(list('type' = 'bm', 'beta' = beta,'y_class' = y_class,  'md'=md, 'constant' = constant))

}

multiclass_univariate_SFLDA = function(data, y, tau, lambda){

  data = data[[1]]

  y_class = unique(y)

  nclass = length(y_class)

  xtr = vector(mode = "list", length = nclass)

  for (n in seq(nclass)){
    xtr[[n]] = data[y == y_class[n],]
  }

  SD = SDelta_multi(xtr)

  S = SD$S

  delta_list = SD$delta_list

  p = length(delta_list[[1]])
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  S = S / norm(S, type = 'F')

  beta_list =  vector(mode = "list", length = nclass - 1)

  for (n in seq(nclass - 1)){
    beta_list[[n]] = coordlasso((1-tau)*S + tau*DD, delta_list[[n]], lambda, 1e-8, 500)
  }

  # make each beta orthogonal

  b = do.call(cbind, beta_list)

  if (pracma::Rank(b) == nclass - 1){
    b = gramSchmidt(b)$Q
  }

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = paste('xb',as.character(seq(nclass - 1)),sep = '')

  if (is.null(is_constant(xb))){

    md = table(y)

    constant = TRUE

  } else {
    xb = is_constant(xb)
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }

  return(list('type' = 'mu','beta' = b, 'y_class' = y_class, 'md' = md, 'constant' =  constant))

}

multiclass_multivariate_SFLDA = function(data, y, tau, lambda){

  n_mat = length(data)

  p_list = c()

  for (i in seq(n_mat)){
    p_list = c(p_list, ncol(data[[i]]))
  }

  data = do.call(cbind, data)

  y_class = unique(y)

  nclass = length(y_class)

  xtr = vector(mode = "list", length = nclass)

  for (n in seq(nclass)){
    xtr[[n]] = data[y == y_class[n],]
  }

  SD = SDelta_multi(xtr)

  S = SD$S

  delta_list = SD$delta_list

  p = length(delta_list[[1]])
  D = cbind(diag(p-1), rep(0,p-1)) + cbind(rep(0,p-1), - diag(p-1))

  s = 0

  for (j in p_list[-n_mat]){
    D[s +j,] = rep(0,p)

    s = s + j
  }

  DD = t(D)%*%D
  DD = DD / norm(DD, type = 'F')

  S = S / norm(S, type = 'F')

  beta_list =  vector(mode = "list", length = nclass - 1)

  for (n in seq(nclass - 1)){
    beta_list[[n]] = coordlasso((1-tau)*S + tau*DD, delta_list[[n]], lambda, 1e-8, 500)
  }

  # make each beta orthogonal

  b = do.call(cbind, beta_list)

  if (pracma::Rank(b) == nclass - 1){
    b = gramSchmidt(b)$Q
  }

  xb = data %*% b

  xb = data.frame(xb)
  colnames(xb) = paste('xb',as.character(seq(nclass - 1)),sep = '')

  if (is.null(is_constant(xb))){

    md = table(y)
    constant = TRUE

  } else {
    xb = is_constant(xb)
    xb$class = as.factor(y)

    md = lda(formula = class ~ ., data = xb)

    constant = FALSE

  }

  beta = vector('list', length = n_mat)

  s = 1

  for (j in seq(n_mat)){

    beta[[j]] = b[seq(s,s + p_list[j]-1),]

    s = s + p_list[j]

  }

  return(list('type' = 'mm','beta' = beta,'y_class' = y_class,'md'= md, 'constant' = constant))

}


predict_class = function(model,x_test){
  # predict class of given test set

  # model: SFLDA model object
  # x_test: test data in list of matrix

  type = model$type

  if (type == 'bu'){

    x_test = x_test[[1]]

    y_class = model$y_class

    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)

    } else{


      beta = model$beta

      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = 'xb1'

      pred_class = predict(model$md, pX)$class
    }




  }

  if (type == 'bm'){

    x_test = do.call(cbind, x_test)

    y_class = model$y_class


    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)
    }else {

      beta = model$beta

      beta = do.call(rbind, beta)


      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = 'xb1'

      pred_class = predict(model$md, pX)$class

    }


  }

  if (type == 'mu'){

    x_test = x_test[[1]]

    y_class = model$y_class

    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)

    } else{


      beta = model$beta

      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = paste('xb',as.character(seq(ncol(pX))),sep = '')

      pred_class = predict(model$md, pX)$class
    }

  }

  if (type == 'mm'){

    x_test = do.call(cbind,x_test)

    y_class = model$y_class

    if (model$constant){

      pred_class = factor(rep(names(model$md)[which.max(model$md)], nrow(x_test)),level = y_class)
    }else {

      beta = model$beta

      beta = do.call(rbind, beta)


      pX = x_test%*%beta

      pX = data.frame(pX)
      colnames(pX) = paste('xb',as.character(seq(ncol(pX))),sep = '')

      pred_class = predict(model$md, pX)$class

    }


  }

  return(pred_class)


}



