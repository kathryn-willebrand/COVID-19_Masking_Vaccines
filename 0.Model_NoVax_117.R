#Model Structure - Constant Masking
library(deSolve)
setwd("C:/Users/Justin/Box/Masks_Vaccines/")
source("0.Parameters.r")

###### ODEs ###### 

### Base Model - No Vaccination 
seir_noVax <- function(times,init_novax,parms_novax){
  
  with(as.list(c(init_novax,parms_novax)), {
    
    dS <- -beta_117*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*mask_prop))
    dE <-  beta_117*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*mask_prop))- l*E
    dI <- l*E*p_symp - h*I - p*I
    dA <- l*E*p_asymp - p*A
    dH <- h*I - g*H - cfr*H
    dR <- p*(I+A) + g*H 
    dM <- cfr*H 
    
    s.flow.out <- -beta_117*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*mask_prop))
    e.flow.in <- beta_117*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*mask_prop))
    e.flow.out <- - l*E
    i.flow.in <- l*E*p_symp
    i.flow.out <- - h*I - p*I
    a.flow.in <- l*E*p_asymp
    a.flow.out <- - p*A
    h.flow.in <- h*I 
    h.flow.out <- - g*H - cfr*H
    ra.flow <- p*A
    ri.flow <- p*I
    rh.flow <- g*H
    mh.flow <- cfr*H
    
    return(list(c(dS, dE, dI, dA, dH, dR, dM, 
                  s.flow.out, e.flow.in, e.flow.out, i.flow.in, i.flow.out, a.flow.in, a.flow.out, h.flow.in, h.flow.out, ra.flow, ri.flow, rh.flow, mh.flow)))
  })	
  
}

t0 = 0 #Start
tn =300 #End 
times <- seq(t0, tn, by = 1) 

mask_prop <- 0.75

# Parameters 
parms_novax <- cbind(beta_117, presymp_infect, asymp_infect, p_symp, p_asymp, mask_eff, mask_prop, n, l, h, p, g, cfr)

# initial starting values
init_novax <- c(S = S,                  # Susceptible;
                E = E,                  # Incubating;
                I = I,                  # Symptomatically Infectious; 
                A = A,                  # Asymptomatically Infectious;
                H = H,                  # Hospitalized;
                R = R,                  # Recovered;
                M = M,                  # Deceased
                
                s.flow.out = 0, 
                e.flow.in = 0, 
                e.flow.out = 0,  
                i.flow.in = 0, 
                i.flow.out = 0,  
                a.flow.in = 0, 
                a.flow.out = 0, 
                h.flow.in = 0, 
                h.flow.out = 0, 
                ra.flow = 0, 
                ri.flow = 0, 
                rh.flow = 0,
                mh.flow = 0)

model_noVax_117 <- as.data.frame(ode(init_novax, times, seir_noVax, parms_novax)) 

#### Validation
# S -> E; S -> V
#round((model_noVax$s.flow.out)+(model_noVax$e.flow.in))

# E -> A; E -> I
#round(model_noVax$e.flow.out+(model_noVax$a.flow.in)+(model_noVax$i.flow.in))

# I -> H; I -> R
#round(model_noVax$i.flow.out + model_noVax$h.flow.in + model_noVax$ri.flow)

# A -> R
#round(model_noVax$a.flow.out + model_noVax$ra.flow)

# H -> R; H -> M
#round(model_noVax$h.flow.out+model_noVax$rh.flow+model_noVax$mh.flow)


# Total Flow
#round(rowSums(model_noVax[,c("S", "E", "I", "A", "H", "R", "M")])-n)

# Visualize 
#plot(times, model_noVax$S, type = 'l')
plot(times, model_noVax_117$I, type = 'l')
plot(times, model_noVax_117$A, type = 'l')
#plot(times, model_noVax$H, type = 'l')
#plot(times, model_noVax$R, type = 'l')
#plot(times, model_noVax$M, type = 'l')
