# Functions used for MCMC for SIR model
# Jan 2019
# Sarah Lukens
# ---------------------------------------

require('optimx')
require(deSolve)

############################################################
# RIGHT HAND SIDE SIR                                      #
############################################################
# Sets up the RHS of the ODE, y = f(t; p)
# where S=y1, I=y2, R=y3, and p are the parameters
# for SIR, p=c(b, g) 
# where b is the transmission and g is recovery rate
SIR <- function (t, y, p) {
  {
    S <-y[1]
    I <- y[2]
    R <- y[3]
  }
  with(as.list(p), {
    dS.dt <- -b * I * S
    dI.dt <- b * I * S - g * I
    dR.dt <- g * I
    return(list(c(dS.dt, dI.dt, dR.dt)))
  })
}

############################################################
# function that just computes cost:                        #
############################################################
# least squares
costfcn <- function(params0, data){
  # friendly name for cases at time t:
  t <- data[,1]
  cases <- data[,2]
  
  b <- params0[1]
  g <- params0[2]
  S0 <- 762
  I0 <- 1
  R0 <- 0
  out <- as.data.frame(ode(y=c(S=S0,I=I0,R=R0),times=t,SIR,parms=c(b,g),hmax=1/120))
  sse<-sum((out$I-cases)^2)
  return(sse)
}

############################################################
# Function that makes a newpoint, but keeps within bounds  #
############################################################
newpoint <- function(params, scale, U, L){
  
  paramsnew = params
  for (i in 1:length(params)){
    paramsnew[i] = params[i]*exp( scale*rnorm(1,mean=0,sd=1) )
    if ( paramsnew[i] > U[i] ){
      paramsnew[i] = ( L[i]/U[i] )*paramsnew[i]
    }
    if ( paramsnew[i] < L[i] ){
      paramsnew[i] = ( U[i]/L[i] )*paramsnew[i]
    }
  }
  
  
  return(paramsnew)
  
}