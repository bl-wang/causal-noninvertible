### added different linear patterns 11/14/2020

pairdatagenerator = function(n = NULL, x = NULL, error = T, mx = NULL, sdx = NULL, nl = T, method = "nl1", lin = "pw",
                             se = 1, lo = 0.3, hi = 4, a = NULL, tau.p = NULL, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  # if(is.null(xl) | length(xl) == 0) xl = matrix(0, n, 1) #replicate(Num.lin, rnorm(n, 0, 1))
  # if(is.null(xnl) | length(xnl) == 0) xnl = matrix(0, n, 1)#replicate(Num.nl, rnorm(n, 0, 1))
  # xl = as.matrix(xl); xnl = as.matrix(xnl)
  # if(is.null(Num.lin)) Num.lin = ncol(xl)
  # if(is.null(Num.nl)) Num.nl = ncol(xnl)
  scale01 = function(x){
    (x-min(x))/diff(range(x))
  }
  
  if(is.null(mx)) mx = 0
  if(is.null(sdx)) sdx = 1
  if(is.null(x)) x = rnorm(n, mx, sdx) #runif(n, -1, 1)#
  ## random error
  if(error){
    e = rnorm(n, 0, se)
  }else{
    e = rep(0, n)
  }
  
  if(is.null(a)) {
    ab = runif(2, lo, hi)
    a = runif(1, lo, hi)
    # a = max(ab); b = ab[which.min(abs(ab))]
    #a = ab[which.max(abs(ab))]*1.5; b = ab[which.min(abs(ab))]
  }
  
  sign1 = sample(c(-1,1), 1)
  sign2 = sample(c(-1,1), 1)
  sign3 = sample(c(-1,1), 1)
  
  if(nl) {
    
    if(method == "pw") {
      sign = sample(c(-1,1), 1)
      #y = sign*a*x + e
      # y = sign*(a + b*x) + e
      
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
    
    if(method == "nl1") y = sign1*max(ab)*x^2 + sign2*min(ab)*x + e #y = sign1*max(ab)*x^2 + e #
    
    if(method == "nl2") y = sign1*cos(max(2, 2*a)*x) + e
    
    if(method == "nl3") y = sign1*a*x^3 + sign2*3.5*a*x^2 + e
    
    if(method == "nl4") y = sign1*tanh(x) + sign2*cos(2.5*x) + sign3*x^2 + e
    #sin(-pi/2*x) + tanh(x) + e
    
    #plot(x, y, cex = 0.1)
    #a = ab[which.max(abs(ab))]*1.5; #b = ab[which.min(abs(ab))]
    ## child node
    #y = a*x^2 + b*x + e
    
  }else{
    #a = ab[which.min(abs(ab))]; b = ab[which.max(abs(ab))]
    #y = a + b*x + e
    sign = sample(c(-1,1), 1)
    
    
    if(lin == "random"){
      lin = sample(c("pw", "sl", "exp", "cubic", "log", "sinh", "tanh", "sigmoid"), 1) 
    }
    
    if(lin == "pw") {
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
    }else if(lin == "sl"){
      a = runif(1, lo, hi)
      b = max(runif(10, lo, hi))
      
      y = a + b*x + e
    }else if(lin == "exp"){
      y = exp(x) + e
    }else if(lin == "cubic"){
      y = x^3 + e
    }else if(lin == "log"){
      x = scale01(x)+0.01
      y = log(x) + e
    }else if(lin == "sinh"){
      y = sinh(x) + e
    }else if (lin == "tanh"){
      y = tanh(x) + e
    }else if(lin == "sigmoid"){
      y = 1/(1+exp(-x))
    }
      

    
  

    #y = sign*a*x + e
    # 
    
    
  }
  
  
  dat = NULL
  dat$x = x #scale(x)
  dat$y = y #scale(y)
  return(dat)
}

# set.seed(201905)
# plot(polydatagenerator(n, se=0.5), cex = 0.5)
# ttp = NULL
# for(i in 1:100){
#   tt = polydatagenerator(n, se=0.5)
#   ttt = mybootstrap(tt$x, tt$y)
#   ttp = c(ttp, ttt$p.r2)
# }
# 
# ttp
