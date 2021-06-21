packages = function(x){
  #x = as.character(match.call()[[2]])
  if(!require(x,character.only=TRUE)){
    install.packages(pkgs=x, quiet = TRUE)
    require(x,character.only=TRUE)
  }else{
    require(x,character.only=TRUE)
  }
}

pkgBiocLite = function(x){
  #x = as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("BiocInstaller", suppressUpdates = TRUE)
    biocLite(x, suppressUpdates = TRUE)
    require(x,character.only=TRUE)
  }else{
    require(x,character.only=TRUE)
  }
}

pkgBiocLite("graph")
pkgBiocLite("RBGL")
packages("pcalg")
packages("sparsebn")
packages("igraph")
packages("bnlearn")
packages("rstudioapi")

# source("./FunctionsV2/pwldatagenerator.R")
#source("./FunctionsV2/polydatagenerator.R")
# path = "~/OneDrive - University of California/Research/Bitbucket_Repos/causal-network"


##########################################
### require function pairdatafenerator()
##########################################

mydatagen = function(n = NULL, Xl = matrix(numeric(0),0,0), Xnl = matrix(numeric(0),0,0), 
                     Num.lin = NULL, Num.nl = NULL, se = 1, method = "nl1", lin = "pw"){ #, smoothdata = smoothdata
  
  # if(smoothdata){
  #   pairdatagenerator = polydatagenerator
  # }else{
  #   pairdatagenerator = pwldatagenerator
  # }
  
  Xl = as.matrix(Xl); Xnl = as.matrix(Xnl)
  Num.lin = ncol(Xl)
  Num.nl = ncol(Xnl)
  
  if(Num.lin > 0) {
    Yl = sapply(1:Num.lin, function(i) pairdatagenerator(n, x = Xl[, i], error = F, nl = F, method = method, lin = lin)$y)
  }else {
    Yl = matrix(0, n, 1)
  }
  
  if(Num.nl > 0) {
    Ynl = sapply(1:Num.nl, function(i) pairdatagenerator(n, x = Xnl[, i], error = F, nl = T, method = method, lin = lin)$y)
  }else{
    Ynl = matrix(0, n, 1)
  }
  
  y = rowSums(Yl) + rowSums(Ynl) + rnorm(n, 0, se)
  
  dat = NULL
  dat$x.lin = Xl
  dat$x.nl = Xnl
  dat$y = scale(y)
  return(dat)
}


DAGdatagen = function(n = NULL, nnode = NULL, nedge = NULL, nedge.nl= NULL, labels = NULL,
                      prob = NULL, seed.dag = NULL, seed.dat = NULL, se = 1, moral = TRUE,
                      dagEdgeList = NULL, edgeInd = NULL, NLedgeInd = NULL, method = "nl1", lin = "pw") {

  
  if(is.null(labels)) labels = paste("X", 1:nnode, sep = "")
  
  if(!is.null(seed.dag)){set.seed(seed.dag)}
  
  if(is.null(dagEdgeList)){
    dagEdgeList = sparsebnUtils::random.graph(nnode, nedge)
  }
  
  # names(dagEdgeList) = as.character(1:nnode)
  # dagAmat = as(dagEdgeList, "matrix")
  # dag.bn = to_bn(as(dagAmat, "graphNEL"))
  # 
  # edgeInd = apply(dag.bn$arcs, 2, as.numeric)
  edgeInd = matrix( as.numeric( unlist( strsplit( unlist( sapply( 1:nnode, function(i){
    if(sum(dagEdgeList[[i]])!=0) paste(dagEdgeList[[i]], i, sep = "-") } )), "-"))), ncol = 2, byrow = T )
  # edgeInd = matrix( labels[as.numeric( unlist( strsplit( unlist( sapply( 1:nnode, function(i){
  #   if(sum(dagEdgeList[[i]])!=0) paste(dagEdgeList[[i]], i, sep = "-") } )), "-")))], ncol = 2, byrow = T )
  
  #non linear edges
  if(is.null(NLedgeInd)){
    NLedgeInd = matrix(edgeInd[sample(1:nrow(edgeInd), nedge.nl),], ncol= 2, byrow = F)
  }else{
    NLedgeInd = matrix(NLedgeInd, ncol = 2, byrow = F)
  }
  
  
  if(!is.null(seed.dat)){set.seed(seed.dat)}
  dat = scale(sapply(1:nnode, function(i)  rnorm(n, 0, 1)))
  
  edge_order = as.numeric(topo_sort(to_igraph(dagEdgeList)))
  
  
  for(i in edge_order) {
    if(i %in% edgeInd[, 2]) {
      if(i %in% NLedgeInd[, 2]) {
        Pind.nl = NLedgeInd[NLedgeInd[, 2] %in% i, 1]
        Pind.lin = edgeInd[edgeInd[, 2] %in% i, 1][!edgeInd[edgeInd[, 2] %in% i, 1]%in%Pind.nl]
        dat[, i] = mydatagen(n = n, Xl = dat[, Pind.lin], Xnl = dat[, Pind.nl], se = se, method = method, lin = lin)$y
      }else {
        Pind = edgeInd[edgeInd[, 2] %in% i, 1]
        dat[, i] = (mydatagen(n = n, Xl = dat[, Pind], se = se, method = method, lin = lin)$y)
      }
    }
  }
  
  # ground truth. cpdag with nonlinear edge fixed
  
  names(dagEdgeList) = labels
  dagAmat = as(dagEdgeList, "matrix")
  dag.bn = to_bn(dagEdgeList)
  
  wl = matrix(labels[NLedgeInd], ncol = 2)
    #apply(NLedgeInd, 2, as.character)
  dag.bn$learning$whitelist = wl
  cpdag.bn = cpdag(dag.bn, moral = moral, wlbl = TRUE)
  cpdagAmat = amat(cpdag.bn)
  cpdagEdge = edgeList(cpdagAmat)
  
  cpdagAmat.nl = matrix(0, nnode, nnode)
  colnames(cpdagAmat.nl) = rownames(cpdagAmat.nl) = labels
  cpdagAmat.nl[NLedgeInd] = 1 
  cpdagEdge.nl = edgeList(cpdagAmat.nl)

  colnames(dat) = names(cpdagEdge) = names(cpdagEdge.nl) = labels
  
  
  return(list(DAGdata = scale(dat), dag.bn = dag.bn, cpdag.bn = cpdag.bn, 
              cpdagEdgeList = cpdagEdge, cpdagEdgeList.nl = cpdagEdge.nl, 
              cpdagAmat = cpdagAmat, cpdagAmat.nl = cpdagAmat.nl))
}



# mydatagen = function(n = NULL, Xl = NULL, Xnl = NULL, Num.lin = NULL, Num.nl = NULL, se = 1, 
#                       a1=NULL, a2=NULL, b1=NULL, b2=NULL, al = NULL, bl = NULL, ta = 0.25, tb = 0.75, ptau = NULL){
#   if(is.null(Xl) | length(Xl) == 0) Xl = matrix(0, n, 1) #replicate(Num.lin, rnorm(n, 0, 1))
#   if(is.null(Xnl) | length(Xnl) == 0) Xnl = matrix(0, n, 1)#replicate(Num.nl, rnorm(n, 0, 1))
#   Xl = as.matrix(Xl); Xnl = as.matrix(Xnl)
#   if(is.null(Num.lin)) Num.lin = ncol(Xl)
#   if(is.null(Num.nl)) Num.nl = ncol(Xnl)
#       
#   if(is.null(ptau)) ptau = round(runif(1, ta, tb), 2)
#   tau.t = round(apply(Xnl, 2, function(x) runif(1, quantile(x, 0.25), quantile(x, 0.75))), 2)
#   #tau.t = round(apply(Xnl, 2, function(x) quantile(x, probs = ptau)), 3)
#   
#   Xnl1 = as.matrix(sapply(1:Num.nl, function(i) Xnl[,i][which(Xnl[,i] < tau.t[i])]))
#   Xnl2 = as.matrix(sapply(1:Num.nl, function(i) Xnl[,i][which(Xnl[,i] >= tau.t[i])]))
#   
#   a1 = runif(Num.nl, -1, 1)
#   # b1.sign = sample(c(-1, 1), Num.nl, replace = TRUE)
#   # b1 = runif(Num.nl, 0.7, 2) * b1.sign
#   # b2 = - runif(Num.nl, 0.7, 2) * b1.sign
#   
#   b1.sign = sample(c(-1, 1), Num.nl, replace = TRUE)
#   b2.sign = - b1.sign
#   b = runif(2*Num.nl, 0.8, 5)
#   b1.ind = sample(1:(2*Num.nl), Num.nl)
#   b1 = b[b1.ind] * b1.sign 
#   b2 = b[-b1.ind] * b2.sign
#   
#   a2 = sapply(1:Num.nl, function(i) a1[i] + tau.t[i] * (b1[i] -  b2[i]))
#   
#   bl = c(runif(1, -1, 1), runif(Num.lin, 0.5, 2) * sample(c(-1, 1), Num.lin, replace = TRUE))
#   
#   ## random error
#   e = rnorm(n, 0, se)
# 
#   y = cbind(1, Xl) %*% bl + e
#   ## child node
#   for(i in 1:Num.nl) {
#     if(!identical(Xnl1[i][[1]], numeric(0))) {
#       y[which(Xnl[, i] < tau.t[i])] = y[which(Xnl[, i] < tau.t[i])] + a1[i] + b1[i] * Xnl1[[i]]
#     } 
#     if(!identical(Xnl2[i][[1]], numeric(0))) {
#       y[which(Xnl[, i] >= tau.t[i])] = y[which(Xnl[, i] >= tau.t[i])] + a2[i] + b2[i] * Xnl2[[i]]
#     }
#   }
# 
#   dat = NULL
#   dat$x.lin = Xl
#   dat$x.nl = Xnl
#   dat$y = y
#   return(dat)
# }