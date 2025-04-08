## Lab 10: age-structured model

########################################################################################
## PART 1: Two AGE GROUP MODEL
########################################################################################

library(deSolve)

# sample code for the 2-age group model
SIR2ageGrs <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSC = nu - SC * (betaCC * IC + betaCA * IA) - muC * SC - lC * SC
    dIC = SC * (betaCC * IC + betaCA * IA) - gamma * IC - muC * IC - lC * IC
    dSA = lC * SC - SA * (betaAC * IC + betaAA * IA) - muA * SA
    dIA = lC * IC + SA * (betaAC * IC + betaAA * IA) - gamma * IA - muA * IA
    
    list(c(dSC, dIC, dSA, dIA))
  })
}



## parameters/initial conditions 
## note all time units are in year here
betaCC = 1100; # transmission rate among children per year
betaCA = betaAC = 110;  # transmission rate between children and adults per year
betaAA = 220; # transmission rate among adults per year
gamma = 1 / (14/365); # recovery rate: 1/ infectious period, infectious period = 2 weeks  

lC = 0.066667; # rate of leaving the children group: 1/15 per year

muC = 0; # mortality rate of children
muA = 0.0166667; # mortality rate of adults:  1/(75-15);  mean life expectancy=75yr; subtract the 15 yrs in youth with 0 death rate

# fraction of population in the two age groups
nC = 0.2;  # % children = 15/75; 
nA = 0.8; # % adults = 1-nC;
nu = 1/75; # birth rate, per year

# initial susceptibility 
SC0 = .1; # initial susceptibility in children 10%
SA0 =.1; # initial susceptibility in adults 10%
# initial infectious:
IC0 = .0001; # initial infection rate in children 
IA0 = .0001; # initial infection rate in adults
# initial recovered:
RC0 = nC - IC0 - SC0; 
RA0 = nA - IA0 - SA0; 

######  CODE ON YOUR OWN  ###### 
## to run the model, you will have to put all pieces (parameters, initial conditions, etc.) together. 
## [Q1] Run the model using parameters/initial conditions listed in the next slide for 100 years.  
## Plot the fraction of susceptibles (Si/ni) for the two groups. what do you find, which group has a higher susceptibility? Why? (1pt)

parametersSIR2=(c(betaCC=betaCC, betaCA=betaCA,
                  betaAA=betaAA, betaAC=betaAC, gamma=gamma))

state=(c(SC=SC0, IC=IC0,
         SA=SA0, IA=IA0))

times=seq(0,100,by=1)

simSIR2=ode(y=state,times=times,func=SIR2ageGrs,parms=parametersSIR2)

SC_frac=simSIR2[,'SC'] / nC * 100  # % of susceptible children
SA_frac=simSIR2[,'SA'] / nA * 100  # % of susceptible adults

plot(simSIR2[,'time'], SC_frac, type = "l", col = "blue", lwd = 2,
     xlab = "Time (years)", ylab = "Percent Susceptible (%)",
     main = "Percent Susceptible Over Time", ylim = c(0, max(SC_frac, SA_frac)))
lines(simSIR2[,'time'], SA_frac, col = "red", lwd = 2)
legend("topright", legend = c("Children", "Adults"), col = c("blue", "red"), lwd = 2)


## SAMPLE CODE TO NORMALIZE THE PROPORTIONS RELATIVE TO THE SPECIFIC GROUP (VS THE ENTIRE POPULATION)
if(F){
  ## proportion susceptible in the two groups, i.e. Si/ni
  fSC=sim[,'SC']/nC; # note: it is normalized by population size in each group
  fSA=sim[,'SA']/nA;
  
  ## proportion infectious in the two groups, i.e. Ii/ni
  fIC=sim[,'IC']/nC;
  fIA=sim[,'IA']/nA;
  
  
  ## plot # people infected
  totI=rowSums(sim[,c('IC','IA')]);  # total infectious
  totS=rowSums(sim[,c('SC','SA')]);  # total susecptible
  
}



# [Q2] What is R0 for the entire population? What would the average age of infection be, given this R0 value and a life span of 75 yr? (0.5pt)
# [Hint: use the same method we used for risk structured model, ignore death rate and aging rate]

beta= matrix(c(1100,110,110,220),2,2) # the beta matrix
nC = 0.2;  # % children = 15/75; 
nA = 0.8; # % adults = 1-nC;
nu = 1/75; # birth rate, per year

n=c(nC,nA, nu)      # n is the vector storing th proportion in each group
n.matrix=diag(n,2,2)  # matrix related to the population size in each group
# to see it:
View(n.matrix)

gamma=1 / (14/365); b = 1;
R.matrix=n.matrix %*% beta / gamma
# to see the output of the eigen function:
eigen(R.matrix)
## To find R0
R0=eigen(R.matrix)$values[1]

L=75
infection_age = L/R0


# [Q3] What is the force of infection for each of the two groups, when the disease / your simulation reaches equilibrium? (0.5pt)[Hint: What is the force of infection?]

# Extract the last row of the simulation (t = 100)
IC_eq <- tail(simSIR2[,'IC'], 1)
IA_eq <- tail(simSIR2[,'IA'], 1)

# Force of infection for children and adults
lambda_C <- betaCC * (IC_eq / nC) + betaCA * (IA_eq / nA)
lambda_A <- betaAC * (IC_eq / nC) + betaAA * (IA_eq / nA)

lambda_C
lambda_A



# [Q4] What is the force of infection for each of the two groups, given the parameters and initial state variables at the beginning of the simulation? (0.5pt)[Hint: What is the force of infection?]

# Initial infected and population proportions
IC0 <- 0.0001
IA0 <- 0.0001
nC <- 0.2
nA <- 0.8

# Transmission rates
betaCC <- 1100
betaCA <- 110
betaAC <- 110
betaAA <- 220

# Force of infection at time 0
lambda_C0 <- betaCC * (IC0 / nC) + betaCA * (IA0 / nA)
lambda_A0 <- betaAC * (IC0 / nC) + betaAA * (IA0 / nA)

lambda_C0
lambda_A0



########################################################################################
## PART 2: VACCINATION
########################################################################################
# [Q5] Code a model with 2 age groups and vaccination at birth using the equations in the slide 12
# (Submit the finished code in your lab report.)  


# Modified 2-age group SIR model with vaccination
SIR2ageVacc <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Force of infection for each group
    lambdaC <- betaCC * IC + betaCA * IA
    lambdaA <- betaAC * IC + betaAA * IA
    
    # Differential equations
    dSC = nu * (1 - p) - SC * lambdaC - muC * SC - lC * SC
    dIC = SC * lambdaC - gamma * IC - muC * IC - lC * IC
    dRC = nu * p + gamma * IC - muC * RC - lC * RC
    
    dSA = lC * SC - SA * lambdaA - muA * SA
    dIA = lC * IC + SA * lambdaA - gamma * IA - muA * IA
    dRA = lC * RC + gamma * IA - muA * RA
    
    return(list(c(dSC, dIC, dRC, dSA, dIA, dRA)))
  })
}

# Parameters (same as before)
betaCC = 1100
betaCA = betaAC = 110
betaAA = 220
gamma = 1 / (14 / 365)
lC = 1 / 15
muC = 0
muA = 1 / (75 - 15)
nu = 1 / 75
p = 0.5  # 50% of newborns vaccinated

nC = 0.2
nA = 0.8

# Initial conditions
SC0 = 0.1
IC0 = 0.0001
RC0 = nC - SC0 - IC0

SA0 = 0.1
IA0 = 0.0001
RA0 = nA - SA0 - IA0

state = c(SC = SC0, IC = IC0, RC = RC0,
          SA = SA0, IA = IA0, RA = RA0)

parameters = c(betaCC = betaCC, betaCA = betaCA, betaAA = betaAA, betaAC = betaAC,
               gamma = gamma, lC = lC, muC = muC, muA = muA, nu = nu, p = p)

times = seq(0, 100, by = 1)

# Run simulation
simVax = ode(y = state, times = times, func = SIR2ageVacc, parms = parameters)

# Fraction recovered (Ri / ni)
RC_frac = simVax[, "RC"] / nC * 100
RA_frac = simVax[, "RA"] / nA * 100

# Plotting
plot(simVax[, "time"], RC_frac, type = "l", col = "pink", lwd = 2,
     xlab = "Time (years)", ylab = "Percent Recovered (%)",
     main = "Recovered Population Over Time (p = 0.5)", ylim = c(0, max(RC_frac, RA_frac)))
lines(simVax[, "time"], RA_frac, col = "turquoise", lwd = 2)
legend("bottomright", legend = c("Children", "Adults"), col = c("pink", "turquoise"), lwd = 2)



# [Q6] Then test your model with p=0.5 and the same initial conditions/parameters as in Part 1. 
# Plot the fraction of recovered/immune (Ri/ni) for the two age groups.   (0.5 pt)






####################
## TEST DIFFERENT VACCINATION RATES
## use the R Shiny App: ShinyApp_Vac.R

R0s=15
p = 1 - (1 / R0s)
