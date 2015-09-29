# A simple SIR (susceptible, infected, recovered) model
# Stephen Beckett (stephen.beckett@biology.gatech.edu)
# 28th September 2015


# The SIR model (Kermack and McKendrick, 1927) represents the spread of disease through a fixed population in terms of those susceptible to the infection, those which are infected and those which have recovered from the disease.

# Kermack, W. O. and McKendrick, A. G. "A Contribution to the Mathematical Theory of Epidemics." Proc. Roy. Soc. Lond. A 115, 700-721, 1927.




#### Setup ####
#Simulation parameters
h <- 0.01  #step size
Slength <- 500 #number of steps to run simulation over

#State Variables
S <- 0.9 # proportion of susceptibles
I <- 0.1 # proportion of infected
R <- 0 # proportion of recovered

#Model parameters
Beta <- 0.5 # infection rate
Gamma <- 0.2 # recovery rate



#### Simulation  ####
#Run the model
for(time in 1:Slength) {
	#Find time derivatives
	dSdt = -Beta*S*I 
	dIdt = Beta*S*I - Gamma*I
	dRdt = Gamma*I

	#Update state variables using forward-Euler integration
	S = S + h*dSdt
	I = I + h*dIdt
	R = R + h*dRdt

}



