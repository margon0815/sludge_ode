library(deSolve)

# conceptual model for dynamic system of activated sludge reactor and secondary 
# clarifier with variable sludge height based on Krebs et al., 2000

# Markus Ahnert, Dresden University of Technology, Institute for Urban Water Management
# markus.ahnert@tu-dresden.de

# define ODE function
sludge_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    tlag <- t - 1
    if (tlag < 0)
      ylag <- state
    else 
      ylag <- lagvalue(tlag)  # returns a vector
    
    # ODEs
    dhs   <- Q*(1+R)*X_BB/(0.5*A_NB*X_B) - 0.7*Q*R/(0.5*A_NB)
    dtE   <- 0.5*A_NB/(0.7*Q*R)*dhs
    dX_BB <- 0.7*Q/V_BB*R*X_B - Q/V_BB*(1+R)*X_BB
    
    # Variante 1: Nutzung des vorherigen Zeitschrittes
    #dX_B  <- 1000/DSVI*(24*(tE+dtE))^(1/3) - 1000/DSVI*(24*ylag[4])^(1/3) # Umrechnung von Einheit d zu h über 24 h/d
    
    # Variante 2: Nutzung von X_B - das ist der state und der müsste ja auch dem vorherigen Zeitschritt entsprechen?
    dX_B  <- 1000/DSVI*(24*(tE+dtE))^(1/3) - X_B # Umrechnung von Einheit d zu h über 24 h/d
    
    # offensichtlich sind beide Varianten falsch, da sich jeweils kein statischer Zustand 
    # einstellt und das Ergebnis unplausibel ist
    
    # return derivatives
    return(list(c(dhs,dtE,dX_BB,dX_B)))
  })
}

# set up parameters in a list
params = c(Q = 200, R=1, V_BB=100, A_NB=50, DSVI=100 )

state_0 = c(hs=0, tE=0, X_BB=3000, X_B=5000) # initial state

t_start = 0 # start time
t_end = 10 # end time
t_points = seq(t_start, t_end, by = 0.01) # time points for solution

# solve ODE
out <- dede(y = state_0, times = t_points, func = sludge_ode, parms = params)

plot(out)