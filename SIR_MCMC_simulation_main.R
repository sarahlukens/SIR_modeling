# MCMC for SIR model
# January 2019

# This brute force MCMC takes some time because it evaluates an ODE function
# for each simulation.
# --------------------------------------
source('SIRmodel_MCMC_input_functions.R')

data = read.csv("ReportedCases_by_week.csv", header = TRUE, sep = ",",as.is = TRUE)

############################################################
# Set main parameters:
############################################################

# initial guess:
p0 = c(b=0.0025, g=0.47);
# b= 0.002567967,g=0.473856076

# Parameters for MCMC simulation:
blen = 500;
scale = 0.001; # tweak for acceptance rate, smaller for higher acceptance rate
heatbeta = 1; #0.025  # tweak to control the width of the variation
# Set bounds:
U <- c(0.1, 1)
L <- c(0.0001, 0.1)

# --------------------------------------------------------------
# Initialize for the run:
ball = rep(0, blen)
gall = rep(0, blen)
Energy = rep(0, blen)
ball[1]<- p0[1]
gall[1]<- p0[2]
Energy[1] <- costfcn(p0, data)

accept = 0;
reject = 0;
# --------------------------------------------------------------
# Next, we will program a Metropolis--Hastings scheme to sample
# from a distribution proportional to the target

# todo: create a function here
# --------------------------------------------------------------
for(i in 2:blen){
  
  currentb = ball[i-1]
  currentg = gall[i-1]
  
  proposedx = newpoint( c( b=currentb, g=currentg), scale, U, L )
  A = costfcn( proposedx, data)
  
  # Accept and reject
  deltE = A - Energy[i-1]
  h = min(1, exp( -heatbeta*deltE ) )
  
  if (runif(1) < h){   # accept
    accept = accept + 1
    ball[i] = proposedx[1]
    gall[i] = proposedx[2]
    Energy[i] = A
  } else {
    reject = reject + 1
    ball[i] = currentb
    gall[i] = currentg
    Energy[i] = Energy[i-1]
  }
  
}
print('acceptance rate')
print(accept/blen)
print('rejection rate')
print(reject/blen)

par(mfrow=c(1,1))
plot(Energy, type="l")

############################################################
#                      Graph results                       #
###########################################################
burnin = min(500, blen-100);
par(mfrow=c(2,2))

hist(ball[burnin:blen], breaks=40, col=4)
plot(ball[burnin:blen], type='l')
lines(seq(1, length(ball)), rep( 0.002567967, length(ball)), col=4, lty=2)
hist(gall[burnin:blen], breaks=40, col=4)
plot(gall[burnin:blen], type='l')
lines(seq(1, length(gall)), rep( 0.473856076, length(gall)), col=4, lty=2)

###########################################################
# Best fit:
par(mfrow=c(1,1))
bopt <- ball[which(Energy == min(Energy))][1]
gopt <- gall[which(Energy == min(Energy))][1]

#params <- c(b = 0.002567967,g=0.473856076 )
days <- seq(1,max(data$Weeks),by=0.05)  # set the time steps for running
N <- 763
I <- 1
R <- 0
S <- N - I - R
SIR.out <- data.frame(ode(c(S,I,R), days, SIR, c(b=bopt, g=gopt)))
plot(days, SIR.out[,3], type = "l", lty = 2, col = 'green', ylim=c(0,350) )
lines(data$Weeks,data$Cases, type='b', col='red')

