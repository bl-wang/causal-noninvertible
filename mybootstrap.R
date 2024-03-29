


### get model output given a breakpoint value
model.out = function(xp, xc, tau){
  
  ### xp: the parent node
  ### xc the child node
  ### break point est
  
  
  if(min(xp) >= tau) {  ### if tau is lower boundary point
    
    
    xp.1 = xp[which(xp < tau)]; xp.2 = xp
    xc.1 = xc[which(xp < tau)]; xc.2 = xc
    
    mod.1 = NA; coef.1 = c(NA, NA); resids.1 = r.1 = rs.1 = 0
    
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
    
    mod.2 = NA; coef.2 = c(NA, NA); resids.2 = r.2 = rs.2 = 0
    
    
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
    
    resids = c(resids.1, resids.2)
  }
  
  n1 = length(xp.1)
  n2 = length(xp.2)
  
  ### total sum of residuals
  Trss = sum(resids.1^2) + sum(resids.2^2)
  
  if(is.na(r.1)){r.1 = 0}; if(is.na(r.2)){r.2 = 0};
  if(is.na(rs.1)){rs.1 = 0}; if(is.na(rs.2)){rs.2 = 0};
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
  m.x = sapply(tau.x, function(t){model.out(xp = x, xc = y, tau = t)}) ### with x being the parent node
  m.y = sapply(tau.y, function(t){model.out(xp = y, xc = x, tau = t)}) ### with y being the parent node
    
    
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



#save("mybootstrap", file="~/Bitbucket_Repos/chipseqproject/Functions/mybootstrap.Rdata")


