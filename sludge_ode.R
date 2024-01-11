library(deSolve)

# conceptual model for dynamic system of activated sludge reactor and secondary 
# clarifier with variable sludge height based on Krebs et al., 2000

# Markus Ahnert, Dresden University of Technology, Institute for Urban Water Management
# markus.ahnert@tu-dresden.de

# define ODE function
growth_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    Q <- influentfun(t)
    
    # empir. Funktion für X_B aus A 131
    tE   <- 0.5*A_NB*hs/(0.7*Q*R)                               # Eindickzeit ist mittlere Verweilzeit des Schlammes im Schlammbett des NKB
    X_B  <- 1000/DSVI*(24*tE)^(1/3)                             # Bodenschlammkonz. im Nachklärbecken, Umrechnung von Einheit d zu h über 24 h/d
    
    # ODEs
    dhs   <- Q*(1+R)*X_BB/(0.5*A_NB*X_B) - 0.7*Q*R/(0.5*A_NB)   # Schlammspiegelhöhe
    dX_BB <- 0.7*Q/V_BB*R*X_B - Q/V_BB*(1+R)*X_BB               # Schlammkonz. im Belebungsbecken
    
    # Hilfsberechnungen
    Ms <- 0.5*A_NB*hs*X_B + V_BB*X_BB                           # Schlamm im System = Schlamm im NKB + BB 
    MsNKB <- 0.5*A_NB*hs*X_B
    MsBB <- V_BB*X_BB

    # return derivatives
    return(list(c(dhs,dX_BB),Q=Q, X_B=X_B,X_R=0.7*X_B, Ms=Ms,MsBB=MsBB,MsNKB=MsNKB))
  })
}

# set up parameters in a list
params = c(Q = 1000, R=1, V_BB=1000, A_NB=250, DSVI=100 )

state_0 = c(hs=0.1, X_BB=3) # initial state

t_start = 0 # start time
t_end = 5 # end time
t_points = seq(t_start, t_end, by = 0.01) # time points for solution

# dynamischer Zulauf soll für Schlammverlagerung vom BB in NKB sorgen

fQ <- 2

influentdata <- data.frame(time = c(0, 1, 2, 100) ,
                           value = c(params["Q"], fQ*params["Q"], params["Q"], params["Q"]))
influentfun <- approxfun(influentdata, method="constant")

out <- ode(y = state_0, times = t_points, func = growth_ode, parms = params)

plot(out)