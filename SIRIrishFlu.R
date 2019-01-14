# Uses fit SIR model to boarding school influenza data using the optimx package
# and try to contain the epidemic via vaccination

#install.packages(
#  "https://cran.r-project.org/bin/macosx/el-capitan/contrib/3.4/optimx_2013.8.7.tgz", 
#  repos = NULL, type = "source"
#)
#https://CRAN.R-project.org/package=optimx

data = read.csv("ReportedCases_by_week.csv", header = TRUE, sep = ",",as.is = TRUE)
require('optimx')
require(deSolve)

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



# Function runs model then computes sum of squares of error:
sse.sir <- function(params0,data){
t <- data[,1]
cases <- data[,2]
b <- params0[1]
g <- params0[2]
S0 <- 762
I0 <- 1
R0 <- 0
out <- as.data.frame(ode(y=c(S=S0,I=I0,R=R0),times=t,SIR,parms=c(b,g),hmax=1/120))
sse<-sum((out$I-cases)^2)
}

# Initial guess:
params <- c(b = 0.002567967,g=0.473856076 )
days <- seq(1,max(data$Weeks),by=0.05)  # set the time steps for running
N <- 763
I <- 1
R <- 0
S <- N - I - R
SIR.out <- data.frame(ode(c(S,I,R), days, SIR, params))
plot(days, SIR.out[,3], type = "l", lty = 2, col = 'green', ylim=c(0,350) )
lines(data$Weeks,data$Cases, type='b', col='red')

# Estimate the reproductive rate using the fitted model parameters:
R_not=params[1]*N/params[2]
cat('Reproductive number R_0 = ', as.numeric(R_not) )


