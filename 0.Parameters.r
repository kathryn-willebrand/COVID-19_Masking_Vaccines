###### DATA ######
### Last Data Point: 3/22/21
# https://ourworldindata.org/covid-vaccinations
vax_data <- read.csv("C:/Users/Justin/Box/Masks_Vaccines/Data/0330_UnitedStates_Vaccine_Data.csv")


###### Initial Population Size (US Population) ######
#https://www.census.gov/quickfacts/fact/table/US/PST045219
            n <- 328239523        

###### Baseline Disease Characteristics ######
        ###RO
        #https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
            r0 <- 2.5
        
        ###R0 - B.1.1.7        
        #https://www.medrxiv.org/content/10.1101/2020.12.24.20248822v3    
            increase <- 0.5
            #b117_prop <- 0.5
            r0_b117 <- (r0*increase)+r0 
            #pop_r0 <- (r0_b117*b117_prop)+(r0*(1-b117_prop))

        ###Incubation Period
        #https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
            e.dur <- 6
            l <- 1/e.dur
        
        ###Symptomatic period
        #https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
            i.dur <- 10
            #Recovery Rate
            p <- 1/i.dur
        
        ###Duration of infectiousness (E+[A or I])
            infectious_duration <- i.dur + e.dur
        
        ###Proportion asymptomatic
        #https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
            p_asymp <- 0.3
            p_symp <- 1 - p_asymp
        
        ###Beta
            beta <- r0 *(1/infectious_duration)
            beta_117 <- r0_b117 *(1/infectious_duration)
            
        
        ###Infectiousness - Relative to Symptomatic Disease
        #https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
            #Presymptomatic Infectivity
            presymp_infect <- 0.5 
            
            #Asymptomatic Infectivity
            asymp_infect <- 0.75
        
        ### Hospitalization
        #https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30566-3/fulltext
        #https://www.cdc.gov/coronavirus/2019-ncov/hcp/clinical-guidance-management-patients.html
            #Hospitalization Rate
            h <- 1/11
            
            #Hospitalization Recovery Rate
            g <- 1/12
            
            #CFR
            #https://ourworldindata.org/mortality-risk-covid#the-case-fatality-rate
            cfr <- 0.00195
            

###### Interventions ######
        ### Mask Efficacy
        #https://www.healthdata.org/sites/default/files/files/Projects/COVID/Estimation_update_062520.pdf
            mask_eff <- 0.3
            
        ### Mask Proportion Multiplier
            m_v <- 1
            m_u <- 1
        
        ### Vaccine Effectiveness - day 14+ after second dose
            # Disease - Symptomatic and Asymptomatic 
            # https://www.cdc.gov/mmwr/volumes/70/wr/mm7013e3.htm
            vax_eff <- 0.9 
            
            # Hospitalization
            # https://www.nejm.org/doi/full/10.1056/NEJMoa2101765

            hosp_eff <- 0.87
            hv <- h * (1 - hosp_eff) # IV -> HV
            pv <- p *(1 + hosp_eff)  # IV -> R
        
            # Effectiveness of vaccine against death
            # Assumed from "Severe Disease" - no data on deaths 
            # https://www.nejm.org/doi/full/10.1056/NEJMoa2101765
            death_v <- 0.92 
            cfr_v <- cfr * (1 - death_v) # HV -> M
            gv <- g * (1 + death_v) # HV -> R
            
        ### Vaccination Rate
        # https://ourworldindata.org/covid-vaccinations
            lm_vax_rate <- lm(Difference ~ n, data=vax_data)   #Finding linear equation based on real-world data
            vax_first74 <- round(((coef(lm_vax_rate)["n"])*c(1:74)) + (coef(lm_vax_rate)["(Intercept)"]))
            
        
###### Initial Population Characteristics ######
#https://www.cdc.gov/nchs/covid19/covid-19-mortality-data-files.htm

            E <- 1318223            # Incubating;
            I <- 2601595 * p_symp   # Symptomatically Infectious; 
            A <- I * p_asymp        # Asymptomatically Infectious;
            H <- 128947             # Hospitalized;
            V <- 0                  # Vaccinated;
            EV <- 0                 # Incubating - Vaccinated;
            IV <- 0                 # Symptomatically Infectious - Vaccinated;
            AV <- 0                 # Asymptomatically Infectious - Vaccinated;
            HV <- 0                 # Hospitalized;
            R <- 19558126           # Recovered;
            M <- 391607             # Deceased;
            S <- n - (E+I+A+H+V+EV+IV+AV+HV+R+M)

    