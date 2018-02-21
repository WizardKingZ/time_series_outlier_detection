library(MTS)

##
## getAdjustedSeries takes xt: time series, a.hat: estimated residual value
## p: AR(p), q: MA(q), k: dimension of time series, outlier: identified outlier
## io: inovation outlier omega value, ad: additive outlier omega value,
## lshift: level shift omega value, tc: temporary change omega value,
## delta: delta value applied to temporary change, delta in (0,1)
## list.phi.hat: estimated phi, list.theta.hat: estimated theta
## size: number of observation
##
getAdjustedSeries <- function(xt, a.hat, p, q, k, outlier, io, ad, lshift, tc, delta, list.phi.hat, list.theta.hat, size) { 
  type = outlier[1, 2]
  time = outlier[1, 1]
  ## type 1: IO, type 2: AO, type 3: LS, type 4: TC
  if (type == 1) {
      ## ls is MA representation for time series model
      ## to adjust type 1 outlier we subtract Psi impact since the inception of IO outlier
    ls = getARrep(list.theta.hat, list.phi.hat, q, p, size, k)
    for (i in time:n) { 
      if (i == time) {
        xt[i, ] = xt[i, ] - io[time, ]
      } else {
        impact = ls[[i-time]]%*%t(t(io[time, ]))
        xt[i, ] = xt[i, ] + t(impact)
      }
    }
  } else if (type == 2 ) {
      ## to adjust type 2 outlier we subtract omega at the inception of AO outlier
    xt[time, ] = xt[time, ] - ad[time, ]
  } else if (type == 3 ) {
      ## ls is 1/(I - B)
      ## to adjust type 3 outlier we subtract ls[time]*omega since the inception of LS outlier
    ls = getARrep(NULL, list(diag(k)), 0, 1, n, k)
    for (i in time:n) { 
      if (i == time) {
        xt[i, ] = xt[i, ] - lshift[time, ]
      } else {
        impact = ls[[i-time]]%*%t(t(lshift[time, ]))
        xt[i, ] = xt[i, ] + t(impact)
      }
    }
  } else if (type == 4) {
      ## ls is 1/(I - delta*B)
      ## to adjust type 4 outlier we subtract delta*ls[time]*omega since the inception of TC outlier
    ls = getARrep(NULL, list(delta*diag(k)), 0, 1, n, k)
    for (i in time:n) { 
      if (i == time) {
        xt[i, ] = xt[i, ] - tc[time, ]
      } else {
        impact = ls[[i-time]]%*%t(t(tc[time, ]))
        xt[i, ] = xt[i, ] + t(impact)
      }
    }
  }
  return(xt)
}

## mat2list convert an estimated parameter from matrix form to a list form
mat2list <- function(mat, k, p) { 
  if (p <= 0) {
    return(NULL)
  } else {
    l <- list()
    for (i in 1:p) { 
      l[[i]] =mat[, (k*(i-1)+1):(k*i)]
    }
    return(l)
  } 
}

## getARrep is to find the polynomial for pi (= phi^-1 * theta)
## pi = I - sum_i=1^infty pi_i B^i where pi_is are the returning list
getARrep <- function(phi, theta, p, q, size, k) {
  pi <- list()
  for (i in 1:size) { 
    pi.i = matrix(0, k, k)
    phi.i = matrix(0, k, k)
    if (i <= p) {
      phi.i = phi[[i]]
    }
    for (j in 1:i) { 
      theta.j = matrix(0, k, k)
      if (j <= q) { 
        theta.j = theta[[j]]
      }
      if (i == j) { 
        pi.i = pi.i - theta.j
      } else {
        pi.i = pi.i + theta.j%*%pi[[i-j]]
      }
    }
    pi[[i]] = phi.i + pi.i
  }
  return(pi)
}

getOmega <- function(a.hat, pi, invSigma, size, k, p) { 
  omega = matrix(0, size, k)
  result = list()
  c = list()
  n = size
  for (h in 1:n) {
    left = matrix(0, k, k)
    right = matrix(0, k, 1) 
    for (ind in 1:(n-h+1)) { 
      i = ind - 1
      pi.i = -diag(k) 
      if (i > 0 ) { 
        pi.i = pi[[i]]
      }
      left = left + t(pi.i)%*%invSigma%*%pi.i
      
      right = right + t(pi.i)%*%invSigma%*%t(t(a.hat[h+i, ]))
    }
    omega.h = -solve(left)%*%right;
    omega[h, ] = t(omega.h) 
    c[[h]] = solve(left)
  }
  result$Omega = omega
  result$Cov = c
  return(result)
}

## getOutlierResidual calculates the omega value for each type of outlier
## a.hat : estimated residual
## type: outlier type
## delta : delta value for TC
##
## *** Note: this function is specifically tailored for condition estimation of multiple time series ***
## *** MTS returns n-p by k vector for estimated residual value. The first p term is conditioned to be zero for ML estimation procedure ***
##
getOutlierResidual <- function(a.hat, pi, sigma, size, k, type, delta, p, flag) {
    ## type {1 : innovation outlier, 2 : additive outlier, 3 : level shift, 4 : temporary}
    # flag false means conditional;
    result = list()
    omega = matrix(0, size, k)
    invSigma = solve(sigma)
    if (!flag) {
        omega[(p+1):size, ] = a.hat
    } else {
        omega = a.hat
    }
    #print(a.hat)
    if (type == 1) {
        result$Omega = omega
        c = list()
        cinv = list()
        for (i in 1:size) {
            c[[i]] = sigma
            cinv[[i]] = invSigma
        }
        result$Cov = c
        result$invCov = cinv
    } else if (type == 2) {
        npi = pi
        result = getOmega(omega, npi, invSigma, size, k, p)
    } else if (type == 3) {
        npi = getARrep(pi, list(diag(k)), size, 1, size, k)
        result = getOmega(omega, npi, invSigma, size, k, p)
    } else if (type == 4) {
        npi = getARrep(pi, list(delta*diag(k)), size, 1, size, k)
        result = getOmega(omega, npi, invSigma, size, k, p)
    }
    return(result)
}

getC.i.h <- function(omega, sigma, d, k) {
    max = 0
    for (j in 1:d) {
        c = abs(omega[j])/sqrt(sigma[j, j])
        if (c >= max) { 
            max = c 
        }  
    }
    return(max)
}

## getStatistics obtain the t-values for J and C hypothesis testing
getStatistics<-function(io, ad, ls, tc, size, k) { 
  Jmax = matrix(0, size, 6)
  Cmax = matrix(0, size, 6)
  for (i in 1:size) { 
    Jio = t(io$Omega[i, ])%*%solve(io$Cov[[i]])%*%t(t(io$Omega[i, ]))
    Jad = t(ad$Omega[i, ])%*%solve(ad$Cov[[i]])%*%t(t(ad$Omega[i, ]))
    Jls = t(ls$Omega[i, ])%*%solve(ls$Cov[[i]])%*%t(t(ls$Omega[i, ]))
    Jtc = t(tc$Omega[i, ])%*%solve(tc$Cov[[i]])%*%t(t(tc$Omega[i, ]))
    Jmax[i, 1] = max(c(Jio, Jad, Jls, Jtc))
    Jmax[i, 2] = which.max(c(Jio, Jad, Jls, Jtc))
    Jmax[i, 3] = Jio
    Jmax[i, 4] = Jad
    Jmax[i, 5] = Jls
    Jmax[i, 6] = Jtc
    
    
    Cio = getC.i.h(io$Omega[i, ], io$Cov[[i]], k, k)
    Cad = getC.i.h(ad$Omega[i, ], ad$Cov[[i]], k, k)
    Cls = getC.i.h(ls$Omega[i, ], ls$Cov[[i]], k, k)
    Ctc = getC.i.h(tc$Omega[i, ], tc$Cov[[i]], k, k)
    Cmax[i, 1] = max(c(Cio, Cad, Cls, Ctc))
    Cmax[i, 2] = which.max(c(Cio, Cad, Cls, Ctc))
    Cmax[i, 3] = Cio
    Cmax[i, 4] = Cad
    Cmax[i, 5] = Cls
    Cmax[i, 6] = Ctc
  }
  return(data.frame(Jmax=Jmax, Cmax=Cmax))
}

outlierDetect<-function(data, p, q, k, n, delta, critical.j, critical.c) {
  result = list()
  est = VARMA(data, p = p, q = q, include.mean = T)
  a.hat = est$residuals
  sigma.hat = est$Sigma
  phi.hat = est$Phi
  list.phi.hat = mat2list(phi.hat, k, p)
  theta.hat = est$Theta
  list.theta.hat = mat2list(theta.hat, k, q)
  list.ar.rep = getARrep(list.phi.hat, list.theta.hat, p, q, n, k)
  #list.cov.est.ao = getCovEst(list.ar.rep, sigma.hat, n, k)
  
  
  omega.hat.innovation = getOutlierResidual(a.hat, list.ar.rep, sigma.hat, n, k, 1, 0, p, FALSE)
  omega.hat.additive = getOutlierResidual(a.hat, list.ar.rep, sigma.hat, n, k, 2, 0, p, FALSE)
  omega.hat.levelshift = getOutlierResidual(a.hat, list.ar.rep, sigma.hat, n, k, 3, 0, p, FALSE)
  omega.hat.temporary = getOutlierResidual(a.hat, list.ar.rep, sigma.hat, n, k, 4, delta, p, FALSE)
  
  stats = getStatistics(omega.hat.innovation, omega.hat.additive, omega.hat.levelshift, omega.hat.temporary, n, k)
  JmaxTable = matrix(0, 4, 2)
  CmaxTable = matrix(0, 4, 2)
  JmaxTableTemp = matrix(0, 4, 2)
  CmaxTableTemp = matrix(0, 4, 2)
  
  J = stats[, 1:6]
  C = stats[, 7:12]
  for (i in 1:4) { 
    jm = max(J[, i+2])
    jm.index = which.max(J[, i+2])
    cm = max(C[, i+2])
    cm.index = which.max(C[, i+2])
    JmaxTable[i, 1] = jm
    JmaxTable[i, 2] = jm.index
    if (jm > critical.j[i]) { 
      JmaxTableTemp[i, 1] = jm
      JmaxTableTemp[i, 2] = jm.index
    }
    CmaxTable[i, 1] = cm
    CmaxTable[i, 2] = cm.index
    if (cm > critical.c[i]) { 
      CmaxTableTemp[i, 1] = cm
      CmaxTableTemp[i, 2] = cm.index
    }
  }
  
  outlier = matrix(0, 1, 2)
  m = which.max(JmaxTableTemp[, 1])
  mc = which.max(CmaxTableTemp[, 1])
  if (max(JmaxTable[, 1]) > critical.j[m]) { 
    outlier[1, 1] = JmaxTableTemp[m, 2]
    outlier[1, 2] = m
  } else if (max(CmaxTable[, 1]) > critical.c[mc]) {
    outlier[1, 1] = CmaxTableTemp[m, 2]
    outlier[1, 2] = mc 
  }
  
  xt = est$data
  if (outlier[1, 1] != 0) {
    xt = getAdjustedSeries(xt, a.hat, p, q, k, outlier, omega.hat.inovation$Omega, omega.hat.additive$Omega, omega.hat.levelshift$Omega, omega.hat.temporary$Omega, delta, list.phi.hat, list.theta.hat, n)
  }
  result$Jmax = JmaxTable
  result$Cmax = CmaxTable
  result$Outlier = outlier
  result$xt = xt
  return(result) 
}


##################################################################################################
######### Example 1: gas furnance ################################################################
##################################################################################################

##
##---- critical value for J and C should be provided as a form of table. Here since the p-value table haven't tested, we only compare with critical test statistics value
##---- the critical value presented here is 2.5% critical benchmark for AR(6) and ith element of both vectors corresponds to ith type of outlier
##
critical.j <-c(17.29, 17.98, 11.42, 16.73) 
critical.c <-c(3.90, 4.17, 3.19, 3.79)


p <- 6
q <- 0
k <- 2
 
## type of outliers
## innovation outlier
## additive outlier
## level shift
## temporary

##Loading data
## Please place the data file under the same directory of this file.
gas.furnace <- read.csv("~/Dropbox/Research/Johnny Li/gas-furnace.csv")

rate<-gas.furnace$InputGasRate
co2 <- gas.furnace$CO2

n <- length(co2)
## delta value for temporary change is normally set to be 0.7
delta <- 0.7
xt = gas.furnace

res = list()
it = 1
while (TRUE) {
  out = outlierDetect(xt, p, q, k, n, delta, critical.j, critical.c)
  
  if (out$Outlier[1, 1] == 0) {
    break;
  } else {
  
    l = list()
    l$Outlier = out$Outlier
    l$Jmax = out$Jmax
    l$Cmax = out$Cmax
  
    xt = out$xt
    print(l)
    res[[it]] = l
    it = it + 1
  }
}
##
## report returns a summary of outlier and its statistics for each individual iteration
## each interation takes two rows
## first row contains Jmax, Cmax value for each type of outlier and the position of the identified outlier
## first 4 columns are for Jmax from type 1 to 4 and column 5 to 8 are for Cmax from type 1 to 4
## the 9th column is the type of the identified outlier and the 10th is the time index
## the second row contains time indexes for each Jmax and Cmax
##
report = matrix(0, length(res)*2, 10)

for (i in 1:length(res)) {
  item = res[[i]]
  index =2*i - 1
  report[index, 1:4] = item$Jmax[, 1]
  report[index+1, 1:4] = item$Jmax[, 2] 
  report[index, 5:8] = item$Cmax[, 1]
  report[index+1, 5:8] = item$Cmax[, 2]
  report[index, 9:10] =item$Outlier
}
