# conceptual model for dynamic system of activated sludge reactor and secondary 
# clarifier with variable sludge height based on Krebs et al., 2000

# Markus Ahnert, Dresden University of Technology, Institute for Urban Water Management
# markus.ahnert@tu-dresden.de

# define ODE function - in that special case a differential algebraic system
sludge_ode <- function(t, state, params,influentfcn) {
  with(as.list(c(state, params)), {
    
    Q <- influentfcn(t)
    K <- (1000*24^(1/3))/DSVI
    
    X_R <-0.7*(Ms_NKB/((1/K)^3*0.7*Q*R))^(1/4)
    
    # ODEs
    # sludge conc. in activated sludge reactor (direct calculation of concentration due to assumed constant volume of AS reactor)
    dX_BB <- Q/V_BB*R*X_R - Q/V_BB*(1+R)*X_BB 
    # sludge mass in secondary clarifier (variable volume of sludge bed as reactor volume)
    dMs_NKB <- Q*(1+R)*X_BB-Q*R*X_R             
    
    # additional calculations
    # height of sludge bed
    hs <- 0.7*Ms_NKB/(0.5*A_NKB*X_R)
    # sludge mass in AS reactor
    MsBB <- V_BB*X_BB
    # total sludge mass as control variable
    Ms <- MsBB+Ms_NKB                           
    
    # return derivatives and additional results
    return(list(c(dX_BB,dMs_NKB),Ms_BB=MsBB,X_R=X_R,hs=hs,Qin=Q))
  })
}