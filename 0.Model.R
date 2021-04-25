#Model Structure - Constant Masking
library(deSolve)
setwd("C:/Users/Justin/Box/Masks_Vaccines/")
source("0.Parameters.R")

###### ODEs ###### 
### Varying Masking Proportion by vaccinated proportion

seir_maskFun <- function(times,init_maskFun, parms_maskFun){
  
  with(as.list(c(init_maskFun,parms_maskFun)), {
    
    dS <- -beta*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M)))) - (vax_fun(times)/(n-M))*S 
    dE <-  beta*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M)))) - l*E
    dI <- l*E*p_symp - h*I - p*I
    dA <- l*E*p_asymp - p*A
    dH <- h*I - g*H - cfr*H
    dV <-  (vax_fun(times)/(n-M))*S - vax_eff*beta*(V/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M)))) 
    dEV <- (vax_eff*beta*(V/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M)))) ) - l*EV
    dIV <- l*EV*p_symp - pv*IV - hv*IV
    dAV <- l*EV*p_asymp - p*AV
    dHV <- hv*IV - gv*HV - cfr_v*HV
    dR <- p*(I+A+AV) + pv*IV + g*H + gv*HV
    dM <- cfr*H + cfr_v*HV
    
    s.flow.out <- -beta*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M))))  - (vax_fun(times)/(n-M))*S 
    e.flow.in <- beta*(S/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M)))) 
    e.flow.out <- - l*E
    i.flow.in <- l*E*p_symp
    i.flow.out <- - h*I - p*I
    a.flow.in <- l*E*p_asymp
    a.flow.out <- - p*A
    h.flow.in <- h*I 
    h.flow.out <- - g*H - cfr*H
    v.flow.in <- (vax_fun(times)/(n-M))*S
    v.flow.out <- -vax_eff*beta*(V/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M)))) 
    ev.flow.in <- vax_eff*beta*(V/n)*((E*presymp_infect) + (A*asymp_infect) + I + (EV*presymp_infect) + (AV*asymp_infect) + IV) * (1-(mask_eff*(((1-(((1-m_u)*mask_fun(times))/(n-M)))*base_mask*((n-M)-mask_fun(times))+(1-(((1-m_v)*mask_fun(times))/(n-M)))*base_mask*(mask_fun(times)))/(n-M)))) 
    ev.flow.out <- - l*EV
    iv.flow.in <- l*EV*p_symp 
    iv.flow.out <- - pv*IV - hv*IV
    av.flow.in <- l*EV*p_asymp
    av.flow.out <- - p*AV
    hv.flow.in <- hv*IV 
    hv.flow.out <- - gv*HV - cfr_v*HV
    ra.flow <- p*A
    ri.flow <- p*I
    rh.flow <- g*H
    rav.flow <- p*AV
    riv.flow <- pv*IV
    rhv.flow <- gv*HV
    mh.flow <- cfr*H
    mhv.flow <- cfr_v*HV
    
    return(list(c(dS, dE, dI, dA, dH, dV, dEV, dIV, dAV, dHV, dR, dM, 
                  s.flow.out, e.flow.in, e.flow.out, i.flow.in, i.flow.out, a.flow.in, a.flow.out, h.flow.in, h.flow.out, 
                  v.flow.in, v.flow.out, ev.flow.in, ev.flow.out, iv.flow.in, iv.flow.out, av.flow.in, av.flow.out, hv.flow.in, hv.flow.out, 
                  ra.flow, ri.flow, rh.flow, rav.flow, riv.flow, rhv.flow, mh.flow, mhv.flow)))
  })	
  
}

t0 = 0 #Start
tn =300 #End 
times <- seq(t0, tn, by = 1) 


### Vaccination Rate
    vax <- c(vax_first74, rep(1200000, 132), rep(0,95))
    vax_fun <- approxfun(vax, rule=2)
    

### Masking Parameters 
    base_mask <- 0.75
    m_v <- 1
    m_u <- 1
    mask_fun <- approxfun(cumsum(vax), rule = 2)
    unvax_mask <- 0.75*(1-(mask_fun(times)/(n-M)))

# Parameters 
parms_maskFun <- cbind(m_u, m_v, base_mask, base_mask, beta, presymp_infect, asymp_infect, p_symp, p_asymp, vax_fun, mask_eff, mask_fun, n, l, h, p, g,pv, hv, gv, cfr, cfr_v)

# initial starting values
init_maskFun <- c(S = S,                  # Susceptible;
                  E = E,                  # Incubating;
                  I = I,                  # Symptomatically Infectious; 
                  A = A,                  # Asymptomatically Infectious;
                  H = H,                  # Hospitalized;
                  V = V,                  # Vaccinated;
                  EV = EV,                # Incubating - Vaccinated;
                  IV = IV,                # Symptomatically Infectious - Vaccinated;
                  AV = AV,                # Asymptomatically Infectious - Vaccinated;
                  HV = HV,                # Hospitalized;
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
                  v.flow.in = 0, 
                  v.flow.out = 0, 
                  ev.flow.in = 0, 
                  ev.flow.out = 0, 
                  iv.flow.in = 0, 
                  iv.flow.out = 0, 
                  av.flow.in = 0,  
                  av.flow.out = 0, 
                  hv.flow.in = 0, 
                  hv.flow.out = 0, 
                  ra.flow = 0, 
                  ri.flow = 0, 
                  rh.flow = 0,
                  rav.flow = 0,
                  riv.flow = 0,
                  rhv.flow = 0,
                  mh.flow = 0,
                  mhv.flow = 0)

model_maskFun <- as.data.frame(ode(init_maskFun, times, seir_maskFun, parms_maskFun))

#### Validation
# S -> E; S -> V
#round((model_maskFun$s.flow.out)+(model_maskFun$e.flow.in)+(model_maskFun$v.flow.in))

# E -> A; E -> I
#round(model_maskFun$e.flow.out+(model_maskFun$a.flow.in)+(model_maskFun$i.flow.in))

# I -> H; I -> R
#round(model_maskFun$i.flow.out + model_maskFun$h.flow.in + model_maskFun$ri.flow)

# A -> R
#round(model_maskFun$a.flow.out + model_maskFun$ra.flow)

# H -> R; H -> M
#round(model_maskFun$h.flow.out+model_maskFun$rh.flow+model_maskFun$mh.flow)

# V -> EV
#round(model_maskFun$v.flow.out+model_maskFun$ev.flow.in)

# EV -> AV; EV -> IV
#round(model_maskFun$ev.flow.out+model_maskFun$iv.flow.in+model_maskFun$av.flow.in)

# IV -> R; IV -> HV
#round(model_maskFun$iv.flow.out+model_maskFun$hv.flow.in+model_maskFun$riv.flow)

# AV -> R
#round(model_maskFun$av.flow.out+model_maskFun$rav.flow)

# HV -> R; HV -> M
#round(model_maskFun$hv.flow.out+model_maskFun$rhv.flow+model_maskFun$mhv.flow)

# Total Flow
#round(rowSums(model_maskFun[,c("S", "E", "I", "A", "H", "V", "EV", "IV", "AV", "HV", "R", "M")])-n)

# Visualize 
#plot(times, model_maskFun$S, type = 'l')
#plot(times, model_maskFun$I, type = 'l')
#plot(times, model_maskFun$IV, type = 'l')
#plot(times, model_maskFun$H, type = 'l')
#plot(times, model_maskFun$HV, type = 'l')
#plot(times, model_maskFun$R, type = 'l')

