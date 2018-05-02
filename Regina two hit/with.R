
#  1. Setup statements 
rm(list = ls())
#Set your working directory
setwd('~/Dropbox/R') #ONLY WORKS ON MAC
#Check that you set your working directory correctly
getwd()
require('deSolve')

# 2. Input data, set parameter values, and/or set state_variable conditions.
#state_variable population size

S0 = 1-.0001
I0 = .0001
R0 = 0
state_vars = c(SS = S0, II = I0, RR = R0)

#Generate a  series of times at which you want the ODE solver to output population sizes
tseq <- seq(0, 5000, by = 1)
#View time sequence
tseq

#Generate a vector of parameter values. This syntax gives the appropriate name to each entry in the vector, so you don't have to remember which parameter comes first, second and third later.
pars <- c(beta = 0.5, gamma = 0.4, mu = 0.0007, vac = .005)
#View variable
pars; names(pars) = NULL



# 3. FUNCTION to pass to the ODE solver

#  use to passequations to the ODE solver as a function. 
#  function MUST input the following:
# 1. a single vector of state_variable conditions, one for each state variable
# 2. a single vector of time steps, here called tseq
# 3. a single vector of parameters, here called pars
# function calculates and outputs a derivative for each state variable


#pars will be a generic term for parameters
SIR_system <- function(tseq, state_vars, pars){
  SS = state_vars[1]
  II = state_vars[2]
  RR = state_vars[3]
  
  beta = pars[1]
  gamma = pars[2]
  mu = pars[3]
  vac = pars[4]
  
  dS_dt = mu-beta*state_vars[1]*state_vars[2]-mu*state_vars[1]-vac*state_vars[1]
  dI_dt = beta*state_vars[1]*state_vars[2] - gamma*state_vars[2]-mu*state_vars[2]
  dR_dt = gamma*state_vars[2]-mu*state_vars[3]+vac*state_vars[1]
  #Hint - with one state variable, we call state_vars[1] to get the first and only entry, which represents NN
  list(c(dS = dS_dt, dI = dI_dt, dR = dR_dt))
}



# 4. Call the ODE solver to output results

#Store the results to a variable called "output"
output <- lsoda(c(S0, I0, R0), tseq, SIR_system, pars)

# Look at the format of the results.
head(output) #The first column is time, "1" is your first state variable (S), "2" is I, and "3" is R

#Plot the output
par(mfrow = c(2,2)) # Set up a plot with two rows and two columns
plot(output[,1], output[,2], col = "blue", type = "l", xlab="time", ylab="N",  main = "S") #Plot column 1 (time), against column 2 (S)
plot(output[,1], output[,3], col = 'red', main = 'I', type = 'l') #Plot column 1 (time), against column 3 (I)
plot(output[,1], output[,4], main = 'R', type = 'l') #Plot column 1 (time), against column 4 (R)

