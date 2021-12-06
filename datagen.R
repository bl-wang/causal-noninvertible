
############ generate a pair of linear/nonlinear data: (x, y), y ~ x
# n: number of data points

pairdatagen = function(n = NULL, x = NULL, error = T, mx = NULL, sdx = NULL, noninv = T, noninv_method = "random", inv_method = "random",
                             se = 1, lo = 0.3, hi = 4, tau.p = NULL){
  
  scale01 = function(x){
    (x-min(x))/diff(range(x))
  }
  
  if(is.null(mx)) mx = 0 # mean of x
  if(is.null(sdx)) sdx = 1 # sd of x
  if(is.null(x)) x = rnorm(n, mx, sdx) 
  
  ## random error
  if(error){
    e = rnorm(n, 0, se)
  }else{
    e = rep(0, n)
  }
  
  
  
  if(noninv) {
    
    if(noninv_method == "random"){
      noninv_method = sample(c("piecewise", "nl1", "nl2", "nl3", "nl4"), 1) 
    }
    
    if(noninv_method == "piecewise") {
      sign = sample(c(-1,1), 1)
      
      tau.t = if(is.null(tau.p)){
        mean(sample(x, 10))
      }else{
        quantile(x, probs = tau.p)
      } 
      x1 = as.matrix(x[which(x < tau.t)])
      x2 = as.matrix(x[which(x >= tau.t)])
      
      b1 =  max(runif(10, lo, hi)) * sign
      b2 = max(runif(10, lo, hi) ) * (-sign)
      a1 = runif(1, lo, hi)
      a2 = a1 + tau.t * (b1 - b2)
      
      y = e
      y[which(x < tau.t)] = y[which(x < tau.t)] + a1 + b1 * x1
      y[which(x >= tau.t)] = y[which(x >= tau.t)] + a2 + b2 * x2
      
    }
    
    
    if(noninv_method == "nl1") y = sample(c(-1,1), 1)*max(runif(2, lo, hi))*x^2 + sample(c(-1,1), 1)*min(runif(2, lo, hi))*x + e
    
    if(noninv_method == "nl2") y = sample(c(-1,1), 1)*cos(max(2, 2*runif(1, lo, hi))*x) + e
    
    if(noninv_method == "nl3") y = sample(c(-1,1), 1)*runif(1, lo, hi)*x^3 + sample(c(-1,1), 1)*3.5*runif(1, lo, hi)*x^2 + e
    
    if(noninv_method == "nl4") y = sample(c(-1,1), 1)*tanh(x) + sample(c(-1,1), 1)*cos(2.5*x) + sample(c(-1,1), 1)*x^2 + e
    
    
  }else{
    
    sign = sample(c(-1,1), 1)
    
    
    if(inv_method == "random"){
      inv_method = sample(c("piecewise", "sl", "exp", "cubic", "log", "sinh", "tanh", "sigmoid"), 1) 
    }
    
    if(inv_method == "piecewise") {
      tau.t = mean(sample(x, 10))
      x1 = as.matrix(x[which(x < tau.t)])
      x2 = as.matrix(x[which(x >= tau.t)])
      
      b12 = sample(tail(abs(runif(15, lo, hi))), 2)
      b1 =  b12[1] * sign
      b2 = b12[2] * sign
      a1 = runif(1, lo, hi)
      a2 = a1 + tau.t * (b1 - b2)
      
      y = e
      y[which(x < tau.t)] = y[which(x < tau.t)] + a1 + b1 * x1
      y[which(x >= tau.t)] = y[which(x >= tau.t)] + a2 + b2 * x2
    }else if(inv_method == "sl"){
      a = runif(1, lo, hi)
      b = max(runif(10, lo, hi))
      
      y = a + b*x + e
    }else if(inv_method == "exp"){
      y = exp(x) + e
    }else if(inv_method == "cubic"){
      y = x^3 + e
    }else if(inv_method == "log"){
      x = scale01(x)+0.01
      y = log(x) + e
    }else if(inv_method == "sinh"){
      y = sinh(x) + e
    }else if (inv_method == "tanh"){
      y = tanh(x) + e
    }else if(inv_method == "sigmoid"){
      y = 1/(1+exp(-x))
    }
    
  }
  
  
  dat = NULL
  dat$x = x #scale(x)
  dat$y = y #scale(y)
  return(dat)
}


# plot(pairdatagen(n = 500, se = 0.5, noninv = T), cex = 0.5)
# plot(pairdatagen(n = 500, se = 0.5, noninv = F), cex = 0.5)






################ generate data for node Y=f(Xinv)+f(Xninv)


nodedatagen = function(n = NULL, Xinv = matrix(numeric(0),0,0), Xninv = matrix(numeric(0),0,0), 
                     nodeN.inv = NULL, nodeN.ninv = NULL, se = 1, noninv_method = "random", inv_method = "random"){
  
  
  Xinv = as.matrix(Xinv); Xninv = as.matrix(Xninv)
  nodeN.inv = ncol(Xinv)
  nodeN.ninv = ncol(Xninv)
  
  if(nodeN.inv > 0) {
    Yinv = sapply(1:nodeN.inv, function(i) pairdatagen(n, x = Xinv[, i], error = F, noninv = F, noninv_method = noninv_method, inv_method = inv_method)$y)
  }else {
    Yinv = matrix(0, n, 1)
  }
  
  if(nodeN.ninv > 0) {
    Yninv = sapply(1:nodeN.ninv, function(i) pairdatagen(n, x = Xninv[, i], error = F, noninv = T, noninv_method = noninv_method, inv_method = inv_method)$y)
  }else{
    Yninv = matrix(0, n, 1)
  }
  
  y = rowSums(Yinv) + rowSums(Yninv) + rnorm(n, 0, se)
  
  dat = NULL
  dat$Xinv = Xinv
  dat$Xninv = Xninv
  dat$y = scale(y)
  return(dat)
}



DAGdatagen = function(n = NULL, dag.bn = NULL, nodeN = NULL, labels = NULL, dagsparsityProb = NULL, ninvProb= 1, 
                      edgeInd.ninv = NULL, se = 1, noninv_method = "random", inv_method = "random") {
  
  # n - data size
  
  library(igraph)
  library(bnlearn)
  
  
  if(is.null(labels)) labels = paste("X", 1:nodeN, sep = "") 
  
  if(is.null(dag.bn)){
    if(is.null(dagsparsityProb)) dagsparsityProb = 2/(nodeN-1)
    dag.bn = bnlearn::random.graph(as.character(1:nodeN), prob = dagsparsityProb)
  }else{
    nodes(dag.bn) = as.character(1:nodeN)
  }
  
  edgeInd = apply(dag.bn$arcs, 2, as.numeric)
  
  edgeN = nrow(edgeInd)
  edgeN.ninv = round(edgeN * ninvProb) # non-invertible edge percentage, default 1

  #non-invertible edges
  if(is.null(edgeInd.ninv)){
    edgeInd.ninv = matrix(edgeInd[sample(1:edgeN, edgeN.ninv),], ncol= 2, byrow = F)
  }else{
    edgeInd.ninv = matrix(edgeInd.ninv, ncol = 2, byrow = F)
  }
  
  dat = scale(sapply(1:nodeN, function(i)  rnorm(n, 0, 1)))
  
  edge_order = as.numeric(topo_sort(as.igraph(dag.bn)))
  
  
  for(i in edge_order) {
    if(i %in% edgeInd[, 2]) {
      if(i %in% edgeInd.ninv[, 2]) {
        pa.ninv = edgeInd.ninv[edgeInd.ninv[, 2] %in% i, 1]
        pa.inv = edgeInd[edgeInd[, 2] %in% i, 1][!edgeInd[edgeInd[, 2] %in% i, 1]%in%pa.ninv]
        dat[, i] = nodedatagen(n = n, Xinv = dat[, pa.inv], Xninv = dat[, pa.ninv], se = se, 
                               noninv_method = noninv_method, inv_method = inv_method)$y
      }else {
        pa.inv = edgeInd[edgeInd[, 2] %in% i, 1]
        dat[, i] = nodedatagen(n = n, Xinv = dat[, pa.inv], se = se, 
                               noninv_method = noninv_method, inv_method = inv_method)$y
      }
    }
  }
  
  # ground truth. cpdag with fixed non-invertible edge 
  
  colnames(dat) = nodes(dag.bn) = labels
  
  wl = matrix(labels[edgeInd.ninv], ncol = 2)
  
  dag.bn$learning$whitelist = wl
  cpdag.bn = cpdag(dag.bn, wlbl = TRUE)
  cpdagAmat = amat(cpdag.bn)

  cpdagAmat.ninv = matrix(0, nodeN, nodeN)
  colnames(cpdagAmat.ninv) = rownames(cpdagAmat.ninv) = labels
  cpdagAmat.ninv[edgeInd.ninv] = 1 
  
  
  return(list(DAGdata = scale(dat), dag.bn = dag.bn, cpdag.bn = cpdag.bn, 
              cpdagAmat = cpdagAmat, cpdagAmat.ninv = cpdagAmat.ninv))
}



