### source("mybootstrap.R")


packages = function(x){
  #x = as.character(match.call()[[2]])
  if(!require(x,character.only=TRUE)){
    install.packages(pkgs=x, quiet = TRUE)
    require(x,character.only=TRUE)
  }else{
    require(x,character.only=TRUE)
  }
}

pkgBioConductor = function(x){
  #x = as.character(match.call()[[2]])
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!require(x,character.only=TRUE)){
    BiocManager::install("RBGL")
    require(x,character.only=TRUE)
  }else{
    require(x,character.only=TRUE)
  }
}

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocInstaller", suppressUpdates = TRUE)
pkgBioConductor("graph")
pkgBioConductor("RBGL")
packages("pcalg")
packages("sparsebn")
packages("igraph")
packages("bnlearn")




### get model output given a breakpoint value
model.out = function(xp, xc, tau, stat){
  
  ### xp: the parent node
  ### xc the child node
  ### break point est
  
  
  if(min(xp) >= tau) {  ### if tau is lower boundary point
    
    
    xp.1 = xp[which(xp < tau)]; xp.2 = xp
    xc.1 = xc[which(xp < tau)]; xc.2 = xc
    
    mod.1 = NA; coef.1 = c(NA, NA); resids.1 = numeric(0); r.1 = 0
    
    mod.2 = lm(xc ~ xp)
    coef.2 = coef(mod.2) 
    resids.2 = resid(mod.2)
    r.2 = cor(xp.2, xc.2)
    if(stat=="spearman") r.2 = cor(xp.2, xc.2, method = "spearman") ### alternative stat: spearman
    
  }else if(max(xp) <= tau) { ### if tau is upper boundary point
    
    xp.1 = xp; xp.2 = xp[which(xp > tau)]
    xc.1 = xc; xc.2 = xc[which(xp > tau)]
    
    mod.1 = lm(xc ~ xp) 
    coef.1 = coef(mod.1)
    resids.1 = resid(mod.1)
    r.1 = cor(xp.1, xc.1)
    if(stat=="spearman") r.1 = cor(xp.1, xc.1, method = "spearman");  ### alternative stat: spearman
    
    mod.2 = NA; coef.2 = c(NA, NA); resids.2 = numeric(0); r.2 = 0
    
    
  }else{
    
    ### separate data to two parts by breakpoint tau
    xp.1 = xp[which(xp < tau)]; xp.2 = xp[which(xp >= tau)]
    xc.1 = xc[which(xp < tau)]; xc.2 = xc[which(xp >= tau)]
    
    ### run regression for both parts and otain coefficients, rss, r
    mod.1 = lm(xc.1 ~ xp.1)
    coef.1 = coef(mod.1)
    resids.1 = resid(mod.1)
    r.1 = cor(xp.1, xc.1)
    if(stat=="spearman") r.1 = cor(xp.1, xc.1, method = "spearman") ### alternative stat: spearman
    
    mod.2 = lm(xc.2 ~ xp.2)
    coef.2 = coef(mod.2)
    resids.2 = resid(mod.2)
    r.2 = cor(xp.2, xc.2)
    if(stat=="spearman") r.2 = cor(xp.2, xc.2, method = "spearman") ### alternative stat: spearman
    
  }
  
  n1 = length(xp.1)
  n2 = length(xp.2)
  
  resids = c(resids.1, resids.2)
  ### total sum of residuals
  Trss = sum(resids^2)
  
  if(is.na(r.1)){r.1 = 0}; if(is.na(r.2)){r.2 = 0};
  coef.1[is.na(coef.1)] = coef.2[is.na(coef.2)] = 0
  
  
  ###  averaged Rsquare statistic
  RsqAve = (n1 * r.1^2 + n2 * r.2^2) / (n1 + n2)
  
  param = NULL
  param$xp.1 = xp.1
  param$xp.2 = xp.2
  param$xc.1 = xc.1
  param$xc.2 = xc.2
  param$n1 = n1
  param$n2 = n2
  param$r.1 = r.1
  param$r.2 = r.2
  param$coef1 = coef.1 #coef1
  param$coef2 = coef.2 #coef1
  param$resids = resids
  param$Trss = Trss
  param$RsqAve = RsqAve
  
  return(param)
}


### function to find preferred direction between x->y and y->x
findDirec = function(x, y, ntau = 100, tau.x = NULL, tau.y = NULL, stat = "pearson"){
  
  ### ntau: number of break points candidates to compare
  ### tau.x: arbitrary break points for x
  ### tau.y: arbitrary break points for y
  ### stat: type of r2
  
  ### if breakpoint candidate not provided, choose among a set of quantiles of the variable
  if(is.null(tau.x)) tau.x = sapply(list(x), quantile, probs = round((0:ntau)/ntau, 3))
  if(is.null(tau.y)) tau.y = sapply(list(y), quantile, probs = round((0:ntau)/ntau, 3))
  
  ### get model output for each breakpoint candidate
  m.x = sapply(tau.x, function(t){model.out(xp = x, xc = y, tau = t, stat = stat)}) ### with x being the parent node
  m.y = sapply(tau.y, function(t){model.out(xp = y, xc = x, tau = t, stat = stat)}) ### with y being the parent node
  
  
  ### choose the optimal model
  
  ### with x being the parent node
  
  ### optimal model with minimum rss
  ind.x = which.min(unlist(m.x["Trss", ])) 
  ### estimated model parameters
  params.x = m.x[c("r.1", "r.2", "n1", "n2", "coef1", "coef2", "resids", "RsqAve"), ind.x]
  
  
  ### with y being the parent node
  
  ### optimal model with minimum rss
  ind.y = which.min(unlist(m.y["Trss", ]))
  ### estimated model parameters
  params.y = m.y[c("r.1", "r.2", "n1", "n2", "coef1", "coef2", "resids", "RsqAve"), ind.y]
  
  names(params.x$coef1) = names(params.y$coef1) =c("a1", "b1")
  names(params.x$coef2) = names(params.y$coef2) =c("a2", "b2")
  
  
  out = NULL
  
  ### Rsquared average for two directions
  optimRave.x = params.x$RsqAve ### x being the parent node
  optimRave.y = params.y$RsqAve ### y being the parent node
  
  
  ### test statistic: the ratio of the two RsqAve
  testStat = optimRave.x / optimRave.y
  
  if( testStat >= 1 ){
    
    ### prefered direction: x -> y 
    params.p = params.x; params.c = params.y
    
    out$parent = 1 ### first variable is the preferred parent
    out$testStat = testStat
    out$optimRave.p = optimRave.x; out$optimRave.c = optimRave.y
    out$tau.p = tau.x[ind.x]; out$tau.c = tau.y[ind.y]
    # out$rp.1 = r1.x; out$rp.2 = r2.x; out$rc.1 = r1.y; out$rc.2 = r2.y; 
  } else {
    
    ### prefered direction: y -> x 
    params.p = params.y; params.c = params.x
    
    out$parent = 2 ### second variable is the preferred parent
    out$testStat =  1/testStat   ### test statistic greater than 1
    out$optimRave.p = optimRave.y; out$optimRave.c = optimRave.x
    out$tau.p = tau.y[ind.y]; out$tau.c = tau.x[ind.x]
  }
  
  out$resids = params.p$resids
  out$r.p = unlist(params.p[c("r.1", "r.2")])
  out$n.p = unlist(params.p[c("n1", "n2")])        
  out$coef.p = c(params.p$coef1, params.p$coef2)
  out$r.c = unlist(params.c[c("r.1", "r.2")]) 
  out$n.c = unlist(params.c[c("n1", "n2")])
  out$coef.c = c(params.c$coef1, params.c$coef2)
  
  return(out)
}

### bootstrap function
bootstrapfunc = function(x, y, ntau=100, tau.x=NULL, tau.y=NULL, func=findDirec){
  n = length(x)
  s = sample(n, replace = T)
  res = func(x[s], y[s], ntau, tau.x, tau.y)
  return(res)
}


mybootstrap = function(x, y, k=100, ntau = 100, xlabel = "1", ylabel = "2", tau.x = NULL, tau.y = NULL, stat = "pearson",  
                       tau.boot.step = NULL, bootstrap.parameter = FALSE, reest_bp = TRUE, approx = TRUE, testonly = TRUE,
                       PrintProcess = FALSE) {
  func = findDirec
  
  ### x, y: data
  ### k: number of bootstrap samples
  ### tau.x, tau.y: arbitrary breakpoints
  
  n = length(x)
  ### find the prefered model
  object = func(x = x, y = y, ntau = ntau, tau.x = tau.x, tau.y = tau.y, stat = stat)
  
  if(object$parent == 1){
    xp = x; xc = y
    tx = object$tau.p; ty = object$tau.c
    parent.label = xlabel; child.label = ylabel
  }  else {
    xp = y; xc = x
    tx = object$tau.c; ty = object$tau.p
    parent.label = ylabel; child.label = xlabel
  }
  
  tau.p = object$tau.p
  tau.c = object$tau.c
  
  if(min(xp) >= tau.p) {
    xp.1 = xp[which(xp < tau.p)]; xp.2 = xp
    xc.1 = xc[which(xp < tau.p)]; xc.2 = xc
  }else if(max(xp) <= tau.p) {
    xp.1 = xp; xp.2 = xp[which(xp > tau.p)]
    xc.1 = xc; xc.2 = xc[which(xp > tau.p)]
  }else{
    xp.1 = xp[which(xp < tau.p)]; xp.2 = xp[which(xp >= tau.p)]
    xc.1 = xc[which(xp < tau.p)]; xc.2 = xc[which(xp >= tau.p)]
  }
  
  coefs = object$coef.p
  xc.1.hat = coefs[1] + coefs[2] * xp.1
  xc.2.hat = coefs[3] + coefs[4] * xp.2
  
  a = xc.1.hat[which.min(xp.1)]
  b = xc.1.hat[which.max(xp.1)]
  c = xc.2.hat[which.min(xp.2)]
  d = xc.2.hat[which.max(xp.2)]
  
  ### construct null data
  if( max(a, b) <= min(c, d) || min(a, b) >= max(c, d) ) {
    dist.c = 0
  }else if(abs(max(a, b) - min(c, d)) <= abs(min(a, b) - max(c, d))) {
    dist.c = max(a, b) - min(c, d)
  }else {
    dist.c = min(a, b) - max(c, d)
  }
  
  
  xc[which(xp < tau.p)] = xc.1 - dist.c
  
  
  if(reest_bp){
    ### reestimate breakpoint
    obj.null = func(x = xp, y = xc, ntau = ntau, stat = stat)
  }else{
    # get null model at the estimated breakpoint 202005
    obj.null = func(x = xp, y = xc, ntau = ntau, tau.x = tau.p, tau.y = tau.c, stat = stat)
  }
  
  
  if(obj.null$parent == 1) {
    tx.null = obj.null$tau.p; ty.null = obj.null$tau.c;
    x.null = xp; y.null = xc
  } else {
    tx.null = obj.null$tau.c; ty.null = obj.null$tau.p; 
    x.null = xp; y.null = xc
  }
  
  dat.null = cbind(x.null, y.null)
  
  ntau.boot = 2 * tau.boot.step + 1
  
  ### if there are too few observations, set approx to False
  if(min(c(obj.null$n.p, obj.null$n.c))<=3) approx = F
  
  res.null = NULL
  
  if(approx){ ### if using approximation
    
    n1.p = obj.null$n.p[1]; n2.p = obj.null$n.p[2]
    m1.p = atanh(obj.null$r.p[1]); m2.p = atanh(obj.null$r.p[2])
    s1.p = 1/sqrt(n1.p-3); s2.p = 1/sqrt(n2.p-3)
    
    n1.c = obj.null$n.c[1]; n2.c = obj.null$n.c[2]
    m1.c = atanh(obj.null$r.c[1]); m2.c = atanh(obj.null$r.c[2])
    s1.c = 1/sqrt(n1.c-3); s2.c = 1/sqrt(n2.c-3)
    
    phi1.p = rnorm(n = 100, mean = m1.p, sd = s1.p)
    phi2.p = rnorm(n = 100, mean = m2.p, sd = s2.p)
    
    phi1.c = rnorm(n = 100, mean = m1.c, sd = s1.c)
    phi2.c = rnorm(n = 100, mean = m2.c, sd = s2.c)
    
    # rsq1.x = tanh(phi1.x)^2
    # rsq2.x = tanh(phi2.x)^2
    # rsq1.y = tanh(phi1.y)^2
    # rsq2.y = tanh(phi2.y)^2
    # 
    # 
    # R2.x = ((r$n.x[1] * rsq1.x + r$n.x[2] * rsq2.x) / (r$n.x[1] + r$n.x[2]))
    # R2.y = ((r$n.y[1] * rsq1.y + r$n.y[2] * rsq2.y) / (r$n.y[1] + r$n.y[2]))
    
    testStat = object$testStat
    testStat.boot = (n1.p * tanh(phi1.p)^2 + n2.p * tanh(phi2.p)^2) / (n1.c * tanh(phi1.c)^2 + n2.c * tanh(phi2.c)^2)
    testStat.boot[testStat.boot<1] = 1/testStat.boot[testStat.boot<1]
    
    # testStat.boot = (n1.x * r1.x^2 + n2.x * r2.x^2) / (n1.y * r1.y^2 + n2.y * r2.y^2)
    # testStat.boot[testStat.boot<1] = 1/testStat.boot[testStat.boot<1]
    
    p.r2 = ifelse(testStat == 1, 1, sum(testStat.boot >= testStat) / k)
    
    
    ## boostrap original data for coefficients significance
    
    
    dat = cbind(x, y)
    colnames(dat) = c(xlabel, ylabel)
    
    outList = NULL
    
    outList$testStat = testStat
    outList$tau = tau.p
    outList$p.r2 = p.r2
    outList$parent = parent.label
    outList$child = child.label
    outList$obj = object
    outList$obj.null = obj.null
    outList$testStat.boot = testStat.boot
    outList$dat.null = dat.null
    outList$dat = dat
    
  }else{
    
    if(testonly){
      
      # for(i in 1:k){
      #   s = sample(n, replace = T);
      #   # func.s = suppressWarnings(func(x = x.null[s], y = y.null[s], tau.x = tx.null, tau.y = ty.null,
      #   #                                ntau = ntau, stat = stat))
      #   func.s = RRfunc(x = x.null[s], y = y.null[s], tau.x = tx.null, tau.y = ty.null)
      #   res.null = cbind(res.null, func.s);
      #   if(PrintProcess == TRUE) print(i)
      # }
      
      res.null = replicate(k, bootstrapfunc(x = x.null, y = y.null, tau.x = tx.null, tau.y = ty.null, func = func))
      
      ## r2 average pvalue
      testStat = object$testStat
      testStat.boot = as.vector(unlist(res.null["testStat",]))
      # hist(testStat.boot)
      #parent.pct = sum(as.vector(unlist(res.null["parent",])) == object$parent) / k
      
      p.r2 = ifelse(testStat == 1, 1, sum(testStat.boot >= testStat) / k)
      
      ## boostrap original data for coefficients significance
      
      
      dat = cbind(x, y)
      colnames(dat) = c(xlabel, ylabel)
      # # p value of statistic under null
      # if(p.r2 > 0.5) p.r2 = (1-p.r2)
      
      
      outList = NULL
      
      outList$testStat = testStat
      outList$tau = tau.p
      outList$p.r2 = p.r2
      outList$parent = parent.label
      outList$child = child.label
      outList$obj = object
      outList$obj.null = obj.null
      outList$testStat.boot = testStat.boot
      outList$dat.null = dat.null
      outList$dat = dat
      
    }else {
      if(is.null(tau.boot.step)){
        
        #ptm <- proc.time()
        # for(i in 1:k){
        #   s = sample(1:n, n, replace = TRUE);
        #   func.s = suppressWarnings(func(x = x.null[s], y = y.null[s], #tau.x = tx.null, tau.y = ty.null,
        #                                  ntau = ntau, stat = stat))
        #   res.null = cbind(res.null, func.s);
        #   if(PrintProcess == TRUE) print(i)
        # }
        res.null = replicate(k, bootstrapfunc(x = x.null, y = y.null, ntau = 100, tau.x = NULL, tau.y = NULL, func = func))
        
        #proc.time() - ptm
      }else{
        
        tau.null.x = quantile( x.null, seq(max(0, round(sum(x.null <= tx.null)/n - tau.boot.step / ntau, 3)),
                                           min(1, round(sum(x.null <= tx.null)/n + tau.boot.step / ntau, 3)),
                                           length.out = ntau.boot) );
        tau.null.y = quantile( y.null, seq(max(0, round(sum(y.null <= ty.null)/n - tau.boot.step / ntau, 3)),
                                           min(1, round(sum(y.null <= ty.null)/n + tau.boot.step / ntau, 3)),
                                           length.out = ntau.boot) );
        
        for(i in 1:k){
          s = sample(1:n, n, replace = TRUE);
          func.s = suppressWarnings(func(x = x.null[s], y = y.null[s], ntau = ntau,
                                         tau.x = tau.null.x, tau.y = tau.null.y, stat = stat));
          
          
          if(func.s$parent == 1) {
            null.p = x.null; null.c = y.null
            tau.null.p = tau.null.x; tau.null.c = tau.null.y; 
          } else {
            null.p = y.null; null.c = x.null
            tau.null.p = tau.null.y; tau.null.c = tau.null.x;
          }
          
          it = 0
          while( (func.s$tau.p != min(null.p)) & (func.s$tau.p != max(null.p)) &
                 (func.s$tau.p == min(tau.null.p) | func.s$tau.p == max(tau.null.p)) ){
            tau.null.p = quantile(null.p, seq(max(0, round(sum(null.p <= func.s$tau.p)/n - tau.boot.step / ntau, 3)),
                                              min(1, round(sum(null.p <= func.s$tau.p)/n + tau.boot.step / ntau, 3)), 
                                              length.out = ntau.boot))
            tau.null.c = quantile(null.c, seq(max(0, round(sum(null.c <= func.s$tau.c)/n - tau.boot.step / ntau, 3)),
                                              min(1, round(sum(null.c <= func.s$tau.c)/n + tau.boot.step / ntau, 3)), 
                                              length.out = ntau.boot))
            func.s = suppressWarnings(func(x = null.p[s], y = null.c[s], ntau = ntau,
                                           tau.x = tau.null.p, tau.y = tau.null.c, stat = stat))
            it = it + 1
            if(it > 10) break
          }
          
          res.null = cbind(res.null, func.s)
          if(PrintProcess == TRUE) print(i)
        }
      }
      
      
      ## r2 average pvalue
      testStat = object$testStat
      testStat.boot = as.vector(unlist(res.null["testStat",]))
      parent.pct = sum(as.vector(unlist(res.null["parent",])) == object$parent) / k
      
      p.r2 = ifelse(testStat == 1, 1, sum(testStat.boot >= testStat) / k)
      
      
      ## boostrap data for parameters
      if (bootstrap.parameter == TRUE) {
        
        res.orig = NULL
        if(is.null(tau.boot.step)){
          
          for(i in 1:k){
            s = sample(1:n, n, replace = TRUE);
            func.s = suppressWarnings(func(x = x[s], y = y[s], ntau = ntau, stat = stat));
            res.orig = cbind(res.orig, func.s);
            if(PrintProcess == TRUE) print(i)
          }
          
        }else{
          
          tau.orig.x = quantile( x, seq(max(0, round(sum(x <= tx)/n - tau.boot.step / ntau, 3)),
                                        min(1, round(sum(x <= tx)/n + tau.boot.step / ntau, 3)), 
                                        length.out = ntau.boot) );
          tau.orig.y = quantile( y.null, seq(max(0, round(sum(y <= ty)/n - tau.boot.step / ntau, 3)),
                                             min(1, round(sum(y <= ty)/n + tau.boot.step / ntau, 3)), 
                                             length.out = ntau.boot) );
          
          for(i in 1:k){
            s = sample(1:n, n, replace = TRUE);
            func.s = suppressWarnings(func(x = x[s], y = y[s], ntau = ntau, stat = stat,
                                           tau.x = tau.orig.x, tau.y = tau.orig.y));
            
            if (func.s$parent == 1) {
              orig.p = x; orig.c = y
              tau.orig.p = tau.orig.x; tau.orig.c = tau.orig.y; 
            } else {
              orig.p = y; orig.c = x
              tau.orig.p = tau.orig.y; tau.orig.c = tau.orig.x;
            }
            
            it = 0
            while( (func.s$tau.p != min(orig.p)) & (func.s$tau.p != max(orig.p)) &
                   (func.s$tau.p == min(tau.orig.p) | func.s$tau.p == max(tau.orig.p)) ){
              tau.orig.p = quantile(orig.p, seq(max(0, round(sum(orig.p <= func.s$tau.p)/n - tau.boot.step / ntau, 3)),
                                                min(1, round(sum(orig.p <= func.s$tau.p)/n + tau.boot.step / ntau, 3)), 
                                                length.out = ntau.boot))
              tau.orig.c = quantile(orig.c, seq(max(0, round(sum(orig.c <= func.s$tau.c)/n - tau.boot.step / ntau, 3)),
                                                min(1, round(sum(orig.c <= func.s$tau.c)/n + tau.boot.step / ntau, 3)), 
                                                length.out = ntau.boot))
              func.s = suppressWarnings(func(x = orig.p[s], y = orig.c[s], ntau = ntau, 
                                             tau.x = tau.orig.p, tau.y = tau.orig.c, stat = stat))
              it = it + 1
              if(it > 10) break
            }
            
            res.orig = cbind(res.orig, func.s)
            if(PrintProcess == TRUE) print(i)
          }
        }
        
        coef.boot = matrix(unlist(res.orig["coef.p",]), ncol = 4, byrow = TRUE)
        colnames(coef.boot) = c("a1", "b1", "a2", "b2")
        
        p.a1 = sum( coef.boot[, 1] >= 0 ) / k
        p.a1 = ifelse( p.a1 <= 0.5, 2 * p.a1, 2 * (1 - p.a1) )
        
        p.b1 = sum( coef.boot[, 2] >= 0 ) / k
        p.b1 = ifelse( p.b1 <= 0.5, 2 * p.b1, 2 * (1 - p.b1) )
        
        p.a2 = sum( coef.boot[, 3] >= 0 ) / k
        p.a2 = ifelse( p.a2 <= 0.5, 2 * p.a2, 2 * (1 - p.a2) )
        
        p.b2 = sum( coef.boot[, 4] >= 0 ) / k
        p.b2 = ifelse( p.b2 <= 0.5, 2 * p.b2, 2 * (1 - p.b2) )
        
        p.b12 = sum( (coef.boot[, 4] - coef.boot[, 2]) >= 0 ) / k
        p.b12 = ifelse( p.b12 <= 0.5, 2 * p.b12, 2 * (1 - p.b12) )
      }
      
      
      
      dat = cbind(x, y)
      colnames(dat) = c(xlabel, ylabel)
      
      
      outList = NULL
      
      outList$testStat = testStat
      outList$tau = tau.p
      outList$p.r2 = p.r2
      outList$parent = parent.label
      outList$child = child.label
      outList$parent.pct = parent.pct
      
      if (bootstrap.parameter == TRUE) {
        outList$coef.boot = coef.boot
        outList$p.a1 = p.a1
        outList$p.b1 = p.b1
        outList$p.a2 = p.a2
        outList$p.b2 = p.b2
        outList$p.b12 = p.b12
      }
      
      outList$obj = object
      outList$obj.null = obj.null
      outList$testStat.boot = testStat.boot
      outList$dat.null = dat.null
      outList$dat = dat
      
    }
  }
  
  return(outList)
}





updatefunc = function(dat, labels, AdjMat){  #, EdgeList
  #print("Updating Graph")
  
  ##########################################
  # added 20200628
  wl.ind = which(AdjMat - t(AdjMat) == 1, arr.ind = T)
  from = labels[wl.ind[,1]]
  to = labels[wl.ind[,2]]
  whitelist = cbind(from, to)
  ##########################################
  
  #EdgeList.nl = edgeList(matrix(0, ncol(dat), ncol(dat)))
  AdjMat.nl = matrix(0, ncol(dat), ncol(dat))
  
  #names(EdgeList) = names(EdgeList.nl) = labels
  colnames(AdjMat) = rownames(AdjMat) = colnames(AdjMat.nl) = rownames(AdjMat.nl) = labels
  # convert to adjacency matrix
  #AdjMat =  G.adjmat 
  #AdjMat = skeleton.adjmat#G.adjmat
  
  # which edges are undirected
  undirectedInd = matrix(which((AdjMat == t(AdjMat) & (AdjMat != 0)), 
                               arr.ind = TRUE), ncol = 2)
  undirectedPairs = matrix(undirectedInd[undirectedInd[, 1] < undirectedInd[, 2], ], ncol = 2)
  
  dat.resid = dat
  # regress each node on its parents to get residuals
  for(i in 1:ncol(dat)) {
    
    # parents of node 
    parent_node = which(AdjMat[, i] - AdjMat[i, ] == 1)
    
    # residual of node i if has parents
    dat.resid[, i] = if(sum(parent_node)) {
      #print(i)
      mod = paste(labels[i], "~", paste(labels[parent_node], collapse = " + "))
      resid(lm(mod, data = as.data.frame(dat)))
    }else {
      dat[, i]
    }
  }
  
  AdjMat.update = AdjMat
  # EdgeList = EdgeList
  
  # undirectedPairs = undirectedPairs
  NpInd = NcInd =  NULL
  #dat.update = vector("list", ncol(AdjMat.update))
  # PairRes = NULL
  undirectedRes = vector("list", nrow(undirectedPairs))
  
  pairInd.update = 1:nrow(undirectedPairs)
  firstround = TRUE
  s=1
  while( s <= nrow(undirectedPairs) ) {
    
    for(k in pairInd.update) {
      # node i and node j
      i =  undirectedPairs[k, 1]; j = undirectedPairs[k, 2]
      
      Vi = dat.resid[, i]; Vj = dat.resid[, j]
      
      parent_i = which(AdjMat.update[, i] - AdjMat.update[i, ] == 1)
      parent_j = which(AdjMat.update[, j] - AdjMat.update[j, ] == 1)
      
      
      # determine reversible edge
      if( sum(NpInd %in% parent_i) | sum(NpInd %in% parent_j) | firstround ) {
        undirectedRes[[k]] = bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
                                       k = iter, ntau = ntau, tau.boot.step = tau.boot.step, 
                                       stat = stat, approx = approx, testonly = testonly)
      }
    }
    
    #colnames(PairRes) = pairnames
    
    # Sort by pvalue ascending, testStat descending
    #edge.update = order(unlist(PairRes[c("p.r2"), ]), -unlist(PairRes[c("testStat"), ]))[1]
    edge.update = order(sapply(undirectedRes, function(x) x$p.r2), 
                        -sapply(undirectedRes, function(x) x$testStat))[[1]]
    
    #undirectedRes[[2]]$p.r2
    if( undirectedRes[[edge.update]]$p.r2 <= p_value ) {
      # edge parent & child
      Np = which(labels %in% undirectedRes[[edge.update]]$parent)
      Nc = which(labels %in% undirectedRes[[edge.update]]$child)
      adjmat.valid.check = AdjMat.update
      # remove the edge from child to parent
      adjmat.valid.check[Nc, Np] = 0
      # remove undirected edges to check if there is directed cycle
      adjmat.valid.check[which(adjmat.valid.check == t(adjmat.valid.check))] = 0
      valid = isValidGraph(t(adjmat.valid.check), type = "dag")
      
      if( valid ) {
        #########################################################
        # Add conditional independence test
        #########################################################
        if(conIndTest) {
          cutoff = qnorm(1-conIndAlpha/2)
          
          # parent set of the child node
          Sind = which(AdjMat.update[, Nc]==1 & AdjMat.update[, Nc]!= t(AdjMat.update)[ , Nc])
          conIndDat = dat[, c(Np, Nc, Sind)] # test on original variable
          conIndDat.resid = dat.resid[, c(Np, Nc, Sind)] # estimates come from residuals
          if(length(Sind)!=0) {
            S = 1:length(Sind) +2
          }else {
            S = NULL
          }
          
          test.tau = undirectedRes[[edge.update]]$tau
          conIndDat1 = conIndDat[which(conIndDat.resid[, 1] <= test.tau), ]
          conIndDat2 = conIndDat[which(conIndDat.resid[, 1] > test.tau), ]
          
          test1 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat1, method=test.method), n=nrow(conIndDat1), cutoff=cutoff)
          test2 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat2, method=test.method), n=nrow(conIndDat2), cutoff=cutoff)
          test = test1 + test2
        }else{
          test = FALSE
        }
        
        
        if(!test) {
          ## if conditional test rejected for both parts
          
          NpInd = cbind(NpInd, Np)
          NcInd = cbind(NcInd, Nc)
          AdjMat.update[Nc, Np] = 0
          # update child node to residuals
          #dat.update[[Np]][[Nc]] = undirectedRes[[edge.update]]$obj$resids
          dat.resid[, Nc] = undirectedRes[[edge.update]]$obj$resids
          
          #EdgeList.nl[[Nc]] = c(EdgeList.nl[[Nc]], Np)
          AdjMat.nl[Np, Nc] = 1
          
          # undirectedRes.update = undirectedRes[[-edge.update]]
          pairInd.update = pairInd.update[-which(pairInd.update == edge.update)]
          
          ### Update after everystep            
          whitelist = rbind(whitelist, matrix(labels[which(AdjMat.nl==1, arr.ind = TRUE)], ncol = 2))
          G.bn = to_bn(edgeList(AdjMat.update))
          G.bn$learning$whitelist = whitelist
          cpdag.bn = cpdag(G.bn, moral = moral, wlbl = TRUE)
          AdjMat.update = amat(cpdag.bn)
          
        }
      }
    }
    
    
    firstround = FALSE
    s=s+1
    
    if( undirectedRes[[edge.update]]$p.r2 <= p_value ) {
      undirectedRes[[edge.update]]$p.r2 = Inf
    } else break
    # if(is.na(sum(undirectedPairs))) break
    # print(s)
    
  }
  
  
  # update Graph edgeList
  #EdgeList.update = edgeList(AdjMat.update)
  
  whitelist = rbind(whitelist, matrix(labels[which(AdjMat.nl==1, arr.ind = TRUE)], ncol = 2))
  G.bn = to_bn(edgeList(AdjMat.update))
  G.bn$learning$whitelist = whitelist
  cpdag.bn = cpdag(G.bn, moral = moral, wlbl = TRUE)
  AdjMat.update = amat(cpdag.bn)
  # EdgeList.update = edgeList(AdjMat.update)
  
  
  
  AdjMat.update.all = AdjMat.update
  AdjMat.nl.all = AdjMat.nl
  #EdgeList.nl.all = EdgeList.nl
  
  noEdgeInd = which(AdjMat.update.all + t(AdjMat.update.all) == 0, arr.ind = TRUE)
  noEdgePairs =  matrix(noEdgeInd[noEdgeInd[, 1] < noEdgeInd[, 2], ], ncol = 2)
  
  if( CheckAllPairs & length(noEdgePairs)>0 ) {
    
    # dat_resid = dat
    res_noEdgePairs = list()
    for(k in 1:nrow(noEdgePairs)) {
      # node i and node j
      i =  noEdgePairs[k, 1]; j = noEdgePairs[k, 2]
      
      Vi = dat.resid[, i]; Vj = dat.resid[, j]
      
      # determine reversible edge
      
      res_noEdgePairs[[k]] = bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
                                       k = iter, ntau = ntau, tau.boot.step = tau.boot.step, 
                                       stat = stat, testonly = testonly)
      # PairRes = cbind(PairRes, bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
      #                          k = iter, ntau = ntau, tau.boot.step = tau.boot.step))
      
    }
    
    p_noEdgePairs = sapply(res_noEdgePairs, function(x) x$p.r2)
    # fdr_noEdgePairs = round(sapply(1:length(p_noEdgePairs), FUN = function(i) {
    #   length(p_noEdgePairs) * p_noEdgePairs[i] / length(which(p_noEdgePairs <= p_noEdgePairs[i]))}), 3)
    
    
    for(k in which(p_noEdgePairs <= p_value)){
      Np = which(labels %in% res_noEdgePairs[[k]]$parent)
      Nc = which(labels %in% res_noEdgePairs[[k]]$child)
      
      adjmat.valid.check = AdjMat.update.all
      adjmat.valid.check[Np, Nc] = 1
      adjmat.valid.check[which(adjmat.valid.check == t(adjmat.valid.check))] = 0
      valid = isValidGraph(t(adjmat.valid.check), type = "dag")
      
      if( valid ) {
        #########################################################
        # Add conditional independence test
        #########################################################
        if(conIndTest) {
          cutoff = qnorm(1-conIndAlpha/2)
          
          # parent set of the child node
          Sind = which(AdjMat.update.all[, Nc]==1 & AdjMat.update.all[, Nc]!= t(AdjMat.update.all)[ , Nc])
          conIndDat = dat[, c(Np, Nc, Sind)]
          conIndDat.resid = dat.resid[, c(Np, Nc, Sind)]
          
          if(length(Sind)!=0) {
            S = 1:length(Sind) +2
          }else {
            S = NULL
          }
          
          test.tau = res_noEdgePairs[[k]]$tau
          conIndDat1 = conIndDat[which(conIndDat.resid[, 1] <= test.tau), ]
          conIndDat2 = conIndDat[which(conIndDat.resid[, 1] > test.tau), ]
          
          test1 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat1, method=test.method), n=nrow(conIndDat1), cutoff=cutoff)
          test2 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat2, method=test.method), n=nrow(conIndDat2), cutoff=cutoff)
          test = test1 + test2
        }else{
          test = FALSE
        }
        
        
        if(!test) {
          ## if conditional test rejected for both parts
          AdjMat.update.all[Np, Nc] = 1
          AdjMat.nl.all[Np, Nc] = 1
          
          
          whitelist.all = rbind(whitelist, matrix(labels[which(AdjMat.nl.all==1, arr.ind = TRUE)], ncol = 2))
          G.bn.all = to_bn(to_graphNEL(edgeList(AdjMat.update.all)))
          G.bn.all$learning$whitelist = whitelist.all
          cpdag.bn.all = cpdag(G.bn.all, moral = moral, wlbl = TRUE)
          AdjMat.update.all = amat(cpdag.bn.all)
        }
      }
    }
    
    
  } else {
    AdjMat.update.all = AdjMat.nl.all = matrix(0, ncol(dat), ncol(dat))
    #EdgeList.update.all = EdgeList.nl.all = edgeList(AdjMat.update.all)
  }
  
  
  #EdgeList.update.all = edgeList(AdjMat.update.all)
  
  
  whitelist.all = rbind(whitelist, matrix(labels[which(AdjMat.nl.all==1, arr.ind = TRUE)], ncol = 2))
  G.bn.all = to_bn(to_graphNEL(edgeList(AdjMat.update.all)))
  G.bn.all$learning$whitelist = whitelist.all
  cpdag.bn.all = cpdag(G.bn.all, moral = moral, wlbl = TRUE)
  
  AdjMat.update.all = amat(cpdag.bn.all)
  #EdgeList.update.all = edgeList(AdjMat.update.all)
  #EdgeList.nl.all = edgeList(AdjMat.nl.all)
  
  # out = list()
  # out$G_update = EdgeList
  # out$adjacency.matrix = AdjMat.update
  # return(out)
  # return(list(AdjMat.update = AdjMat.update, 
  #             EdgeList.update = EdgeList.update))
  return(list(cpdagBN = cpdag.bn, cpdagAddBN = cpdag.bn.all,
              AdjMatrix=AdjMat.update, #EdgeList=EdgeList.update, 
              AdjMatrixNL = AdjMat.nl, #EdgeListNL = EdgeList.nl, 
              AdjMatrixAdd = AdjMat.update.all, #EdgeListAdd = EdgeList.update.all,
              AdjMatrixAddNL = AdjMat.nl.all #EdgeListAddNL = EdgeList.nl.all
  ))
}





###################################################################################################
DAGupdate = function(dat = NULL, labels = NULL, sbn_alpha = 0.2, pc_alpha = 0.01, skel_alpha = NULL,
                     bootstrap = mybootstrap, p_value = 0.01, stat = "pearson", iter = 100, 
                     ntau = 100, tau.boot.step = 3, pc.maxdegree = Inf, V1 = TRUE, V2 = TRUE,
                     UsePC = FALSE, UseSparsebn = FALSE, UseSkeleton = FALSE, UseEmptyG = TRUE, 
                     matchPC = FALSE, CheckAllPairs = TRUE, approx = TRUE, testonly = TRUE, 
                     conIndAlpha = 0.01, conIndTest = TRUE, test.method = "pearson", seed = 12345,
                     lambdas = NULL, sbn.select = NULL, moral = FALSE, printProgress = FALSE) {
  
  set.seed(seed)
  
  # P = ncol(dat)
  # varible labels
  labels = if (is.null(labels)){
    paste("X", 1:ncol(dat), sep = "")
  } else {
    labels
  }
  colnames(dat) = labels
  
  func = function(dat, labels, AdjMat){  #, EdgeList
    # print("Updating Graph")
    amat.fdr.in = amat.fdr.out = matrix(NA, ncol(dat), ncol(dat))
    
    # wl.ind = which(AdjMat - t(AdjMat) == 1, arr.ind = T)
    # from = labels[wl.ind[,1]]
    # to = labels[wl.ind[,2]]
    # whitelist = cbind(from, to)
    whitelist = NULL
    
    # EdgeList.nl = edgeList(matrix(0, ncol(dat), ncol(dat)))
    AdjMat.nl = matrix(0, ncol(dat), ncol(dat))
    
    # names(EdgeList) = names(EdgeList.nl) = labels
    colnames(AdjMat) = rownames(AdjMat) = colnames(AdjMat.nl) = rownames(AdjMat.nl) =
      colnames(amat.fdr.in) = rownames(amat.fdr.in) = colnames(amat.fdr.out) = rownames(amat.fdr.out)= labels
    # convert to adjacency matrix
    # AdjMat =  G.adjmat 
    # AdjMat = skeleton.adjmat#G.adjmat
    
    # which edges are undirected
    undirectedInd = matrix(which((AdjMat == t(AdjMat) & (AdjMat != 0)), 
                                 arr.ind = TRUE), ncol = 2)
    undirectedPairs = matrix(undirectedInd[undirectedInd[, 1] < undirectedInd[, 2], ], ncol = 2)
    
    dat.resid = dat
    # regress each node on its parents to get residuals
    for(i in 1:ncol(dat)) {
      
      # parents of node 
      parent_node = which(AdjMat[, i] - AdjMat[i, ] == 1)
      
      # residual of node i if has parents
      dat.resid[, i] = if(sum(parent_node)) {
        #print(i)
        mod = paste(labels[i], "~", paste(labels[parent_node], collapse = " + "))
        resid(lm(mod, data = as.data.frame(dat)))
      }else {
        dat[, i]
      }
    }
    
    AdjMat.update = AdjMat
    # EdgeList = EdgeList
    
    # undirectedPairs = undirectedPairs
    NpInd = NcInd =  NULL
    #dat.update = vector("list", ncol(AdjMat.update))
    # PairRes = NULL
    undirectedRes = vector("list", nrow(undirectedPairs))
    
    pairInd.update = 1:nrow(undirectedPairs)
    firstround = TRUE
    s=1
    while( s <= nrow(undirectedPairs) ) {
      
      for(k in pairInd.update) {
        # node i and node j
        i =  undirectedPairs[k, 1]; j = undirectedPairs[k, 2]
        
        Vi = dat.resid[, i]; Vj = dat.resid[, j]
        
        parent_i = which(AdjMat.update[, i] - AdjMat.update[i, ] == 1)
        parent_j = which(AdjMat.update[, j] - AdjMat.update[j, ] == 1)
        
        
        # determine reversible edge
        if( sum(NpInd %in% parent_i) | sum(NpInd %in% parent_j) | firstround ) {
          undirectedRes[[k]] = bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
                                         k = iter, ntau = ntau, tau.boot.step = tau.boot.step, 
                                         stat = stat, approx = approx, testonly = testonly)
        }
      }
      
      # colnames(PairRes) = pairnames
      
      # Sort by pvalue ascending, testStat descending
      # edge.update = order(unlist(PairRes[c("p.r2"), ]), -unlist(PairRes[c("testStat"), ]))[1]
      pvalue = sapply(undirectedRes, function(x) x$p.r2)
      if(s == 1) p_all = pvalue
      pvalue = pvalue[pvalue<Inf]
      Npairs = length(pvalue)
      fdr_in = round(sapply(1:Npairs, FUN = function(i) Npairs * pvalue[i] / length(which(pvalue <= pvalue[i]))), 3)
      
      
      edge.update = order(sapply(undirectedRes, function(x) x$p.r2), 
                          -sapply(undirectedRes, function(x) x$testStat))[[1]]
      
      # undirectedRes[[2]]$p.r2
      if( undirectedRes[[edge.update]]$p.r2 <= p_value ) {
        # edge parent & child
        Np = which(labels %in% undirectedRes[[edge.update]]$parent)
        Nc = which(labels %in% undirectedRes[[edge.update]]$child)
        adjmat.valid.check = AdjMat.update
        # remove the edge from child to parent
        adjmat.valid.check[Nc, Np] = 0
        # remove undirected edges to check if there is directed cycle
        adjmat.valid.check[which(adjmat.valid.check == t(adjmat.valid.check))] = 0
        valid = isValidGraph(t(adjmat.valid.check), type = "dag")
        
        if( valid ) {
          #########################################################
          # Add conditional independence test
          #########################################################
          if(conIndTest) {
            cutoff = qnorm(1-conIndAlpha/2)
            
            # parent set of the child node
            Sind = which(AdjMat.update[, Nc]==1 & AdjMat.update[, Nc]!= t(AdjMat.update)[ , Nc])
            conIndDat = dat[, c(Np, Nc, Sind)] # test on original variable
            conIndDat.resid = dat.resid[, c(Np, Nc, Sind)] # estimates come from residuals
            if(length(Sind)!=0) {
              S = 1:length(Sind) +2
            }else {
              S = NULL
            }
            
            test.tau = undirectedRes[[edge.update]]$tau
            conIndDat1 = conIndDat[which(conIndDat.resid[, 1] <= test.tau), ]
            conIndDat2 = conIndDat[which(conIndDat.resid[, 1] > test.tau), ]
            
            test1 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat1, method=test.method), n=nrow(conIndDat1), cutoff=cutoff)
            test2 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat2, method=test.method), n=nrow(conIndDat2), cutoff=cutoff)
            test = test1 + test2
          }else{
            test = FALSE
          }
          
          
          if(!test) {
            ## if conditional test rejected for both parts
            
            NpInd = cbind(NpInd, Np)
            NcInd = cbind(NcInd, Nc)
            AdjMat.update[Nc, Np] = 0
            # update child node to residuals
            # dat.update[[Np]][[Nc]] = undirectedRes[[edge.update]]$obj$resids
            dat.resid[, Nc] = undirectedRes[[edge.update]]$obj$resids
            
            # EdgeList.nl[[Nc]] = c(EdgeList.nl[[Nc]], Np)
            AdjMat.nl[Np, Nc] = 1
            
            # undirectedRes.update = undirectedRes[[-edge.update]]
            pairInd.update = pairInd.update[-which(pairInd.update == edge.update)]
            
            ### Update after everystep            
            whitelist = rbind(whitelist, matrix(labels[which(AdjMat.nl==1, arr.ind = TRUE)], ncol = 2))
            G.bn = to_bn(edgeList(AdjMat.update))
            G.bn$learning$whitelist = whitelist
            cpdag.bn = cpdag(G.bn, moral = moral, wlbl = TRUE)
            AdjMat.update = amat(cpdag.bn)
            
            amat.fdr.in[Np, Nc] = fdr_in[edge.update]
            # print(fdr_in[edge.update])
          }
        }
      }
      
      
      firstround = FALSE
      s=s+1
      
      if( undirectedRes[[edge.update]]$p.r2 <= p_value ) {
        undirectedRes[[edge.update]]$p.r2 = Inf
      } else break
      # if(is.na(sum(undirectedPairs))) break
      # print(s)
      
    }
    
    
    # update Graph edgeList
    #EdgeList.update = edgeList(AdjMat.update)
    
    whitelist = rbind(whitelist, matrix(labels[which(AdjMat.nl==1, arr.ind = TRUE)], ncol = 2))
    G.bn = to_bn(edgeList(AdjMat.update))
    G.bn$learning$whitelist = whitelist
    cpdag.bn = cpdag(G.bn, moral = moral, wlbl = TRUE)
    AdjMat.update = amat(cpdag.bn)
    # EdgeList.update = edgeList(AdjMat.update)
    
    
    
    AdjMat.update.all = AdjMat.update
    AdjMat.nl.all = AdjMat.nl
    #EdgeList.nl.all = EdgeList.nl
    
    noEdgeInd = which(AdjMat.update.all + t(AdjMat.update.all) == 0, arr.ind = TRUE)
    noEdgePairs =  matrix(noEdgeInd[noEdgeInd[, 1] < noEdgeInd[, 2], ], ncol = 2)
    
    if( CheckAllPairs & length(noEdgePairs)>0 ) {
      
      # dat_resid = dat
      res_noEdgePairs = list()
      for(k in 1:nrow(noEdgePairs)) {
        # node i and node j
        i =  noEdgePairs[k, 1]; j = noEdgePairs[k, 2]
        
        Vi = dat.resid[, i]; Vj = dat.resid[, j]
        
        # determine reversible edge
        
        res_noEdgePairs[[k]] = bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
                                         k = iter, ntau = ntau, tau.boot.step = tau.boot.step, 
                                         stat = stat, testonly = testonly)
        # PairRes = cbind(PairRes, bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
        #                          k = iter, ntau = ntau, tau.boot.step = tau.boot.step))
        
      }
      
      p_noEdgePairs = sapply(res_noEdgePairs, function(x) x$p.r2)
      p_all = c(p_all, p_noEdgePairs)
      # fdr_noEdgePairs = round(sapply(1:length(p_noEdgePairs), FUN = function(i) {
      #   length(p_noEdgePairs) * p_noEdgePairs[i] / length(which(p_noEdgePairs <= p_noEdgePairs[i]))}), 3)
      pvalue = p_noEdgePairs
      Npairs = length(p_noEdgePairs)
      fdr_out = round(sapply(1:Npairs, FUN = function(i) Npairs * pvalue[i] / length(which(pvalue <= pvalue[i]))), 3)
      
      
      for( k in which(p_noEdgePairs <= p_value) ){
        Np = which(labels %in% res_noEdgePairs[[k]]$parent)
        Nc = which(labels %in% res_noEdgePairs[[k]]$child)
        
        adjmat.valid.check = AdjMat.update.all
        adjmat.valid.check[Np, Nc] = 1
        adjmat.valid.check[which(adjmat.valid.check == t(adjmat.valid.check))] = 0
        valid = isValidGraph(t(adjmat.valid.check), type = "dag")
        
        if( valid ) {
          #########################################################
          # Add conditional independence test
          #########################################################
          if(conIndTest) {
            cutoff = qnorm(1-conIndAlpha/2)
            
            # parent set of the child node
            Sind = which(AdjMat.update.all[, Nc]==1 & AdjMat.update.all[, Nc]!= t(AdjMat.update.all)[ , Nc])
            conIndDat = dat[, c(Np, Nc, Sind)]
            conIndDat.resid = dat.resid[, c(Np, Nc, Sind)]
            
            if(length(Sind)!=0) {
              S = 1:length(Sind) +2
            }else {
              S = NULL
            }
            
            test.tau = res_noEdgePairs[[k]]$tau
            conIndDat1 = conIndDat[which(conIndDat.resid[, 1] <= test.tau), ]
            conIndDat2 = conIndDat[which(conIndDat.resid[, 1] > test.tau), ]
            
            test1 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat1, method=test.method), n=nrow(conIndDat1), cutoff=cutoff)
            test2 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat2, method=test.method), n=nrow(conIndDat2), cutoff=cutoff)
            test = test1 + test2
          }else{
            test = FALSE
          }
          
          
          if(!test) {
            ## if conditional test rejected for both parts
            AdjMat.update.all[Np, Nc] = 1
            AdjMat.nl.all[Np, Nc] = 1
            
            
            whitelist.all = rbind(whitelist, matrix(labels[which(AdjMat.nl.all==1, arr.ind = TRUE)], ncol = 2))
            G.bn.all = to_bn(to_graphNEL(edgeList(AdjMat.update.all)))
            G.bn.all$learning$whitelist = whitelist.all
            cpdag.bn.all = cpdag(G.bn.all, moral = moral, wlbl = TRUE)
            AdjMat.update.all = amat(cpdag.bn.all)
            
            amat.fdr.out[Np, Nc] = fdr_out[k]
          }
        }
      }
      
      
    } else {
      AdjMat.update.all = AdjMat.nl.all = matrix(0, ncol(dat), ncol(dat))
      #EdgeList.update.all = EdgeList.nl.all = edgeList(AdjMat.update.all)
    }
    
    
    #EdgeList.update.all = edgeList(AdjMat.update.all)
    
    
    whitelist.all = rbind(whitelist, matrix(labels[which(AdjMat.nl.all==1, arr.ind = TRUE)], ncol = 2))
    G.bn.all = to_bn(to_graphNEL(edgeList(AdjMat.update.all)))
    G.bn.all$learning$whitelist = whitelist.all
    cpdag.bn.all = cpdag(G.bn.all, moral = moral, wlbl = TRUE)
    
    AdjMat.update.all = amat(cpdag.bn.all)
    #EdgeList.update.all = edgeList(AdjMat.update.all)
    #EdgeList.nl.all = edgeList(AdjMat.nl.all)
    
    # out = list()
    # out$G_update = EdgeList
    # out$adjacency.matrix = AdjMat.update
    # return(out)
    # return(list(AdjMat.update = AdjMat.update, 
    #             EdgeList.update = EdgeList.update))
    return(list(cpdagBN = cpdag.bn, cpdagAddBN = cpdag.bn.all,
                AdjMatrix=AdjMat.update, #EdgeList=EdgeList.update, 
                AdjMatrixNL = AdjMat.nl, #EdgeListNL = EdgeList.nl, 
                AdjMatrixAdd = AdjMat.update.all, #EdgeListAdd = EdgeList.update.all,
                AdjMatrixAddNL = AdjMat.nl.all, #EdgeListAddNL = EdgeList.nl.all
                amat.fdr.in = amat.fdr.in,
                amat.fdr.out = amat.fdr.out,
                pvalue_all = p_all
    ))
  }
  
  #out = res1 = res2 = res_sbn = res_pc = res_skeleton = nl_update= NULL
  out = NULL
  
  
  if(V1) {
    
    
    if(UseEmptyG) {
      if(printProgress) cat("NL only...")      
      G.amat = matrix(0, ncol(dat), ncol(dat))
      #G.edgelist = edgeList(G.amat)
      nl_update = func(dat = dat, labels = labels, AdjMat = G.amat)
      #plot(res_nobase$EdgeListAdd, layout = test.layout)
      
      out$bn.nl = nl_update$cpdagAddBN
      out$amat.nl = nl_update$AdjMatrixAdd
      # out$fdr.nl = nl_update$amat.fdr
    }
    
    
    
    # Use pcalg graph
    if(UsePC) {
      # pc algorithm Graph
      #print("updating PC")
      if(printProgress) cat("PC + NL...")
      ### retry until find an extendable DAG or randomly generate one
      pc.learn = pc(suffStat = list(C = cor(dat), n = nrow(dat)),
                    indepTest = gaussCItest, m.max = pc.maxdegree,## indep.test: partial correlations
                    alpha = pc_alpha, labels = labels, u2pd = "retry")
      G.learn = as.bn(pc.learn)
      G.amat = amat(G.learn)
      # G.edgelist = edgeList(t(as(pc.learn, "matrix")))
      # G.learn = to_bn(G.edgelist)
      
      #G.learn = pc.stable(as.data.frame(dat), alpha = pc_alpha)
      cpdag.learn = cpdag(G.learn, moral = moral)
      cpdag.amat = amat(cpdag.learn)
      cpdag.edgelist = edgeList(cpdag.amat)
      #isValidGraph(t(cpdag.amat), "cpdag")
      
      skeleton.learn= bnlearn::skeleton(cpdag.learn)
      skeleton.amat = amat(skeleton.learn)
      #skeleton.edgelist = edgeList(skeleton.amat)
      
      out$bn.pc = G.learn
      out$bn.pc.cpdag = cpdag.learn
      out$amat.pc = G.amat
      out$amat.pc.cpdag = cpdag.amat
      
      G_update = func(dat = dat, labels = labels, 
                      AdjMat = cpdag.amat)
      
      out$bn.pc_nl = G_update$cpdagBN
      out$bn.pc_nl_all = G_update$cpdagAddBN
      out$amat.pc_nl = G_update$AdjMatrix
      out$amat.pc_nl_all = G_update$AdjMatrixAdd
      out$NLamat.pc_nl = G_update$AdjMatrixNL
      out$NLamat.pc_nl_all = G_update$AdjMatrixAddNL
      
      
      if(UseSkeleton) {
        skeleton_update = func(dat = dat, labels = labels, AdjMat = skeleton.amat)
        
        out$bn.pcSkl_nl = skeleton_update$cpdagBN
        out$bn.pcSkl_nl_all = skeleton_update$cpdagAddBN
        out$amat.pcSkl_nl = skeleton_update$AdjMatrix
        out$amat.pcSkl_nl_all = skeleton_update$AdjMatrixAdd
        out$NLamat.pcSkl_nl = skeleton_update$AdjMatrixNL
        out$NLamat.pcSKl_nl_all = skeleton_update$AdjMatrixAddNL
      }
      
      # out$fdr.pc_nl_all = G_update$amat.fdr
      
      
      # pc.learn = pc(suffStat = list(C = cor(dat), n = nrow(dat)),
      #               indepTest = gaussCItest, ## indep.test: partial correlations
      #               alpha = pc_alpha, labels = labels, u2pd = "retry")
      # G.adjmat = t(as(pc.learn, "matrix"))
      # if(!isValidGraph(t(G.adjmat), type = "cpdag")) {
      #   G.adjmat = as(dag2cpdag(as(G.adjmat, "graphNEL")), "matrix")
      # }
      # G.edgelist = edgeList(G.adjmat)
      
      
      
      # pc algorithm skeleton
      # if(is.null(skel_alpha)) skel_alpha = pc_alpha
      # skeleton.learn = pcalg::skeleton(suffStat = list(C = cor(dat), n = nrow(dat)),
      #                           indepTest = gaussCItest,#method = "original", ## indep.test: partial correlations
      #                           alpha=skel_alpha, labels = labels)
      # skeleton.adjmat = t(as(skeleton.learn, "matrix"))
      # skeleton.edgelist = edgeList(skeleton.adjmat)
      
    }
    
    
    # Use sparsebn graph
    if(UseSparsebn) {
      #detach(package:igraph)
      if(printProgress) cat("SBN + NL...")      
      
      packages("graph")
      dat.sbn = suppressMessages(sparsebnData(x = dat, type = "continuous"))
      
      if(matchPC) {
        sbn.learn = estimate.dag(data = dat.sbn, edge.threshold = sum(as(pc.learn, "matrix")))
      } else if(is.null(lambdas)){
        sbn.learn = estimate.dag(data = dat.sbn)
      }else{
        sbn.learn = estimate.dag(data = dat.sbn, lambdas = lambdas)
      }
      
      if(is.null(sbn_alpha)) sbn_alpha = 0.1
      if(is.null(sbn.select)){
        G.learn = sbn.learn[[select.parameter(x = sbn.learn, data = dat.sbn, alpha = sbn_alpha)]]$edges
      }else{
        G.learn = sbn.learn[[sbn.select]]$edges
      }
      
      G.learn = to_bn(to_graphNEL(G.learn))
      G.amat = amat(G.learn)
      cpdag.learn = cpdag(G.learn, moral = moral)
      cpdag.amat = amat(cpdag.learn)
      #cpdag.edgelist = edgeList(cpdag.amat)
      
      skeleton.learn= bnlearn::skeleton(G.learn)
      skeleton.amat = amat(skeleton.learn)
      #skeleton.edgelist = edgeList(skeleton.amat)
      out$sbn.learn = sbn.learn
      out$bn.sbn = G.learn
      out$bn.sbn.cpdag = cpdag.learn
      out$amat.sbn = G.amat
      out$amat.sbn.cpdag = cpdag.amat
      # res_sbn$G_est = list(G.learn = G.learn, cpdag.learn = cpdag.learn, EdgeList = cpdag.edgelist, AdjMatrix = cpdag.amat, 
      #                      skeleton.EdgeList = skeleton.edgelist, skeleton.AdjMatrix = skeleton.amat)
      G_update = func(dat = dat, labels = labels, AdjMat = cpdag.amat)
      
      out$bn.sbn_nl = G_update$cpdagBN
      out$bn.sbn_nl_all = G_update$cpdagAddBN
      out$amat.sbn_nl = G_update$AdjMatrix
      out$amat.sbn_nl_all = G_update$AdjMatrixAdd
      out$NLamat.sbn_nl = G_update$AdjMatrixNL
      out$NLamat.sbn_nl_all = G_update$AdjMatrixAddNL
      
      out$pvalue_all = G_update$pvalue_all
      # out$fdr.sbn_nl_all = G_update$amat.fdr.in
      
      if(UseSkeleton) {
        skeleton_update = func(dat = dat, labels = labels, AdjMat = skeleton.amat)
        
        out$bn.sbnSkl_nl = skeleton_update$cpdagBN
        out$bn.sbnSkl_nl_all = skeleton_update$cpdagAddBN
        out$amat.sbnSkl_nl = skeleton_update$AdjMatrix
        out$amat.sbnSkl_nl_all = skeleton_update$AdjMatrixAdd
        out$NLamat.sbnSkl_nl = skeleton_update$AdjMatrixNL
        out$NLamat.sbnSkl_nl_all = skeleton_update$AdjMatrixAddNL
      }
      
      
      
      # G.amat = as(G.edgelist, "matrix")
      # # convert sparsebn graph to cpDAG using pcalg dag2cpDAG
      # G.amat = as(dag2cpdag(as(G.amat, "graphNEL")), "matrix")
      # G.edgelist = edgeList(G.amat)
      # # sparsebn skeleton
      # 
      
      # skeleton.amat = G.amat
      # for(i in 1:nrow(G.amat)) {
      #   for(j in 1:ncol(G.amat))
      #     if(G.amat[i, j] == 1) {
      #       skeleton.amat[j, i] = 1
      #     }
      # }
      # skeleton.edgelist = edgeList(skeleton.amat)
      # 
      # res_sbn$G_est = list(G.learn = sbn.learn, EdgeList = G.edgelist, AdjMatrix = G.amat, 
      #                      skeleton.EdgeList = skeleton.edgelist, skeleton.AdjMatrix = skeleton.amat)
      # res_sbn$G_update = func(dat = dat, labels = labels, 
      #                         AdjMat = G.amat, EdgeList = G.edgelist)
      # if(UseSkeleton) {
      #   res_sbn$skeleton_update = func(dat = dat, labels = labels, 
      #                                  AdjMat = skeleton.amat, EdgeList = skeleton.edgelist)
      # }
    }
    
    # res1$res_pc = res_pc
    # res1$res_sbn = res_sbn
    # res1$res_emptyG = res_emptyG
    #out$res = res
  }
  
  
  
  
  if(V2) {
    
    # estimate nonlinear edges from empty graph
    #print("updating empty graph")
    
    if(!V1){
      G.amat = matrix(0, ncol(dat), ncol(dat))
      #edgelist = edgeList(G.amat)
      nl_update = func(dat = dat, labels = labels, AdjMat = G.amat)
      out$bn.nl = nl_update$cpdagAddBN
      out$amat.nl = nl_update$AdjMatrixAdd
      # out$fdr.nl = nl_update$amat.fdr
    }
    
    amat.nl = out$amat.nl
    #G.edgeList.nl = edgeList(amat.nl)
    
    
    # apply PC algorithm
    if(UsePC) {
      
      if(printProgress) cat("NL + PC...")      
      
      # add nonlinear to fixed edge set
      whitelist = matrix(labels[which(amat.nl==1, arr.ind = TRUE)], ncol = 2)
      #colnames(whitelist) = c("from", "to")
      G.learn = pc.stable(as.data.frame(dat), alpha = pc_alpha, 
                          whitelist = whitelist, max.sx = pc.maxdegree)
      cpdag.learn = cpdag(G.learn, moral = moral, wlbl = TRUE)
      # fixedEdges = amat.nl
      # fixedEdges[t(fixedEdges)==1] = 1
      # 
      
      
      # pc.learn = pc(suffStat = list(C = cor(dat), n = nrow(dat)),
      #               indepTest = gaussCItest, ## indep.test: partial correlations
      #               alpha = pc_alpha, labels = labels, u2pd = "retry", fixedEdges = fixedEdges)
      
      cpdag.amat = amat(cpdag.learn)
      #cpdag.edgelist = edgeList(cpdag.amat)
      
      out$bn.nl_pc = cpdag.learn
      #out$bn.pc_nl_all = G_update$cpdagAddBN
      out$amat.nl_pc = cpdag.amat
      
      # out$fdr.nl_pc = G_update$amat.fdr
      # res_pc$G_update = list(cpdagBN = cpdag.learn, cpdagAddBN = cpdag.learn,
      #                        EdgeList = cpdag.edgelist, AdjMatrix = cpdag.amat, 
      #                        EdgeListNL = G.edgeList.nl, AdjMatrixNL = amat.nl,
      #                        EdgeListAdd = cpdag.edgelist, AdjMatrixAdd = cpdag.amat, 
      #                        EdgeListAddNL = G.edgeList.nl, AdjMatrixAddNL = amat.nl)
    }
    
    # apply sbn algorithm
    if(UseSparsebn) {
      #detach(package:igraph)
      if(printProgress) cat("NL + SBN...Done! \n")      
      
      packages("graph")
      
      # add nonlinear to whitelist
      dat.sbn = suppressMessages(sparsebnData(x = dat, type = "continuous"))
      whitelist = matrix(labels[which(amat.nl==1, arr.ind = TRUE)], ncol = 2)
      if(length(whitelist)==0) whitelist = NULL
      
      if(matchPC) {
        sbn.learn = estimate.dag(data = dat.sbn, edge.threshold = sum(as(pc.learn, "matrix")), whitelist = whitelist)
      } else {
        sbn.learn = estimate.dag(data = dat.sbn, whitelist = whitelist)
      }
      
      if(is.null(sbn_alpha)) sbn_alpha = 0.1
      
      ind = select.parameter(x = sbn.learn, data = dat.sbn, alpha = sbn_alpha)
      if(is.finite(ind)) {
        G.edgelist = sbn.learn[[ind]]$edges
      }else{
        G.edgelist = sbn.learn[[1]]$edges
      }
      
      #G.adjmat = as(G.edgelist, "matrix")
      G.learn = to_bn(to_graphNEL(G.edgelist))
      G.learn$learning$whitelist = whitelist
      cpdag.learn = cpdag(G.learn, moral = moral, wlbl = TRUE)
      cpdag.amat = amat(cpdag.learn)
      #cpdag.edgelist = edgeList(cpdag.amat)
      
      out$bn.nl_sbn = cpdag.learn
      #out$bn.pc_nl_all = G_update$cpdagAddBN
      out$amat.nl_sbn = cpdag.amat
      
      # out$fdr.nl_sbn = G_update$amat.fdr
      #G.adjmat = G.adjmat
      #G.edgelist = edgeList(G.adjmat)
      
      # res_sbn$G_update = list(cpdagBN = cpdag.learn, cpdagAddBN = cpdag.learn,
      #                         EdgeList = cpdag.edgelist, AdjMatrix = cpdag.amat, 
      #                         EdgeListNL = G.edgeList.nl, AdjMatrixNL = amat.nl,
      #                         EdgeListAdd = cpdag.edgelist, AdjMatrixAdd = cpdag.amat, 
      #                         EdgeListAddNL = G.edgeList.nl, AdjMatrixAddNL = amat.nl)
    }
    
    #res2 = list()
    # res2$res_pc = res_pc
    # res2$res_sbn = res_sbn
    # res2$nl_update = nl_update
    # out$res2 = res2
  }
  
  return(out)
  
}








