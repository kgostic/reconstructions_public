rm(list = ls())

##  ____________________________________________
##  A. FIT MODELS ASSUMING NO EXPOSURE
##  ____________________________________________
# Tell the input script to ignore data on rates of age-specific poultry exposure ("E" in models)
exposure.switch = FALSE
source('Model_inputs.R')
# Imports the following data as global variables:
#    H1N1 weights - master.H1p_11, master.H1p_10, master.H1p_01, master.H1p_00
#    H2N2 weights - master.H2p_11, master.H2p_10, master.H2p_01, master.H2p_00
#    H3N2 weights - master.H3p_11, master.H3p_10, master.H3p_01, master.H3p_00
#    naive weights - master.naive
#    age groups - Uc.master (0-4), Um.master (5-64), Ue.master (65+)
#    case incidence - data.master
#    case fatalities - deaths.master
#    demography - p0.y
# See the top of Model_Comp_LikelihoodFunctions for a better description of these variables.
# In general, they are all matrices where each row represents a unique country and year of data observation, and each column represents a birth year.


# source likelihood functions, which are stored in another script.
# you will have to edit these to include your new models
#source('Model_Comp_LikelihoodFunctions_Final.R')
source('Model_Comp_LikelihoodFunctions.R')

#use only g1 to avoid double counting
weights.master.g1 = master.g1p_10 + master.g1p_11 + master.g1p_1n
weights.master.g2 = master.g1p_01 + master.g1p_00 + master.g1p_0n
#weights.master.gx is used for older models
nested_wmm = master.g1p_11
nested_wmo = master.g1p_10 + master.g1p_1n
nested_wom = master.g1p_01
nested_woo = master.g1p_00 + master.g1p_0n
nested_wnn = master.p_nn

## ___________________________
## 1. D - demography only
## ---------------------------
# Since this model has no parameters we don't need the optim function.
# We set Hm = 1, which means Hm has no effect (protected susceptibility = unprotected susceptibility. The risk ratio is 1.)
lk.null = lk.D.H(pars = c(Hm = 1), wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, p0 = p0.master, xx = data.master)
# Call function lk.D.H
# Give inputs:
#       pars = c(Hm = 1)     tells the function that parameter Hm takes value 1, which means no protective effect. (We set no effect because this model assumes no effect of imprinting!)
#       wm = fraction of each bith year with Matched imprinting to H5N1 (i.e. fraction with imprinting to H1 or H2, which are in the same group as H5)
#       wo = fraction of each birth year with Other imprinting (not matched to H5)
#       wn = naive fraction
#       p0 = matrix whose rows give the demographic vector from the country and year of interest
#       xx = matrix of observed cases
#
# Outputs: The negative log likelihood of the parameters given the data.

## ___________________________
## 2. DH - hemagglutinin history
## ---------------------------
# Here, we use pretty much the same imputs, except now we have to optimize the Hm parameter to find the best fit value. Hm gives the risk of infection for protected individuals relative to unprotected indiiduals. It takes values from 0 (no risk) to 1 (same risk as unprotected individuals.)
lk.DH = optim(par = c(Hm = .2), fn = lk.D.H, wm = (weights.master.g1), wo = weights.master.g2, wn = master.p_nn, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001), upper = c(1)); lk.DH
# See ? optim for a description of inputs
# Basically :
# optim([Parameters to be fit], [name of function we are optimizing], [inputs for the function we are optimizing], [more inputs....], ... [optimization arguments like method, and the lower and upper parameter bounds, if necessary])


## ___________________________
## 3. DN - neuraminidase history UNCHANGED
## ---------------------------
# Repeat the above strategy for different models, which may require different inputs or weights.
lk.DN = optim(par = c(Nm = .2), fn = lk.D.N, wm = (weights.master.1), wo = (weights.master.2+weights.master.3), wn = weights.master.naiive, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001), upper = c(1)); lk.DN

## ___________________________
## 4. DA - age risk
## ---------------------------
lk.DA = optim(par = c(Ac = 2, Ae = 2), fn = lk.D.A, uc = Uc.master, um = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(1,1), upper = c(100, 100)); lk.DA

## _________________________________________
## 5. DHN - Hemagglutinin and neuraminidase history UNCHANGED
## -----------------------------------------
lk.DHN = optim(par = c(Hm = .2, Nm = .2), fn = lk.D.H.N, wmm = weights.master.1, wmo = weights.master.2, woo = weights.master.3, wnn = weights.master.naiive, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001), upper = c(1, 1)); lk.DHN

## _________________________________________
## 6. DAH - age risk and hemagglutinin history
## -----------------------------------------
lk.DAH = optim(par = c(Hm = .2, Ac = 2, Ae = 2), fn = lk.D.A.H, wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 1, 1), upper = c(1, Inf, Inf)); lk.DAH

## _________________________________________
## 7. DAN - age risk and neuraminidase history SKIP
## -----------------------------------------
lk.DAN = optim(par = c(Nm = .2, Ac = 2, Ae = 2), fn = lk.D.A.N, wm = (weights.master.1), wo = (weights.master.2+weights.master.3), wn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 1, 1), upper = c(1, Inf, Inf)); lk.DAN

## _________________________________________
## 8. DAHN - age risk, hemagglutinin and neuraminidase history SKIP
## -----------------------------------------
lk.DAHN = optim(par = c(Hm = .2, Nm = .2,  Ac = 2, Ae = 2), fn = lk.D.A.H.N, wmm = weights.master.1, wmo = weights.master.2, woo = weights.master.3, wnn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 1, 1), upper = c(1, 1, 100, 100)); lk.DAHN

## _________________________________________
## 9. DHmep - Hemagglutinin e1 and e2 history
## -----------------------------------------
lk.DHmep = optim(par = c(Hm = .2, e1 = .2, e2 = .7), fn = lk.D.H.ep, wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 0.001), upper = c(1, 1, 1)); lk.DHmep

## _________________________________________
## 10. DAHmep - age risk, hemagglutinin e1 and e2 history
## -----------------------------------------
lk.DAHmep = optim(par = c(Hm = .2, e1 = .2, e2 = .7, Ac = 1.5, Ae = 1.5), fn = lk.D.A.H.ep, wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 0.001, 1, 1), upper = c(1, 1, 1, 100, 100)); lk.DAHmep



##  ____________________________________________
##  B. EXPOSURE
##  ____________________________________________
exposure.switch = TRUE
source('Model_inputs.R')
# Get:
#    H1N1 weights - weights.master.1
#    H2N2 weights - weights.master.2
#    H3N2 weights - weights.master.3
#    naive weights - weights.master.naive
#    age groups - Uc.master (0-4), Um.master (5-64), Ue.master (65+)
#    case incidence - data.master
#    case fatalities - deaths.master
#    demography - p0.y

## ___________________________
## 11. DE - run any model with all pars fixed at 0
## ---------------------------
lk.null.E = lk.D.H(pars = c(Hm = 1), wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, p0 = p0.master, xx = data.master)

## ___________________________
## 12. DEH - hemagglutinin history
## ---------------------------
lk.DEH = optim(par = c(Hm = .2), fn = lk.D.H, wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001), upper = c(1)); lk.DEH

## ___________________________
## 13. DEN - neuraminidase history SKIP
## ---------------------------
lk.DEN = optim(par = c(Nm = .2), fn = lk.D.N, wm = (weights.master.1), wo = (weights.master.2+weights.master.3), wn = weights.master.naiive, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001), upper = c(1)); lk.DEN

## ___________________________
## 14. DEA - age risk
## ---------------------------
lk.DEA = optim(par = c(Ac = 2, Ae = 2), fn = lk.D.A, uc = Uc.master, um = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(1,1), upper = c(Inf, Inf)); lk.DEA

## _________________________________________
## 15. DEHN - Hemagglutinin and neuraminidase history SKIP
## -----------------------------------------
lk.DEHN = optim(par = c(Hm = .2, Nm = .2), fn = lk.D.H.N, wmm = weights.master.1, wmo = weights.master.2, woo = weights.master.3, wnn = weights.master.naiive, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001), upper = c(1, 1)); lk.DEHN


## _________________________________________
## 16. DEAH - age risk and hemagglutinin history
## -----------------------------------------
lk.DEAH = optim(par = c(Hm = .2, Ac = 2, Ae = 2), fn = lk.D.A.H, wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 1, 1), upper = c(1, Inf, Inf)); lk.DEAH


## _________________________________________
## 17. DEAN - age risk and neuraminidase history SKIP
## -----------------------------------------
lk.DEAN = optim(par = c(Nm = .2, Ac = 2, Ae = 2), fn = lk.D.A.N, wm = (weights.master.1), wo = (weights.master.2+weights.master.3), wn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 1, 1), upper = c(1, Inf, Inf)); lk.DEAN

## _________________________________________
## 18. DEAHN - age risk, hemagglutinin and neuraminidase history SKIP
## -----------------------------------------
lk.DEAHN = optim(par = c(Hm = .2, Nm = .2, Ac = 2, Ae = 2), fn = lk.D.A.H.N, wmm = weights.master.1, wmo = weights.master.2, woo = weights.master.3, wnn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 1, 1), upper = c(1, 1, Inf, Inf)); lk.DEAHN

## _________________________________________
## 19. DEHmep - Hemagglutinin e1 and e2 history
## -----------------------------------------
lk.DEHmep = optim(par = c(He1m = .2, He2m = .2), fn = lk.D.He1.He2, wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001), upper = c(1, 1)); lk.DEHe1He2

## _________________________________________
## 20. DEAHmep - age risk, hemagglutinin e1 and e2 history
## -----------------------------------------
lk.DEAHmep = optim(par = c(He1m = .2, He2m = .2,  Ac = 2, Ae = 2), fn = lk.D.A.He1.He2, wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 1, 1), upper = c(1, 1, Inf, Inf)); lk.DEAHe1He2





## AIC compares model fits. The best AIC --> the best model.
##  ____________________________________________
##  C. AIC values
##  ____________________________________________
models = c('D', 'DH', 'DN', 'DA', 'DHN', 'DAH', 'DAN', 'DAHN', 'DHmep', 'DAHmep', 'DE', 'DEH', 'DEN', 'DEA', 'DEHN', 'DEAH', 'DEAN', 'DEAHN', 'DEHmep', 'DEAHmep')
# Calculate AIC values for each model and store in a named vector
AIC.vec = rep(0, length(models)); names(AIC.vec) = models
AIC.vec['D'] = AIC(neg.log.lk = lk.null, n.pars = 0) #1
AIC.vec['DH'] = AIC(neg.log.lk = lk.DH$value, n.pars = 1) #2
AIC.vec['DN'] = AIC(neg.log.lk = lk.DN$value, n.pars = 1) #3
AIC.vec['DA'] = AIC(neg.log.lk = lk.DA$value, n.pars = 2) #4
AIC.vec['DHN'] = AIC(neg.log.lk = lk.DHN$value, n.pars = 2) #5
AIC.vec['DAH'] = AIC(neg.log.lk = lk.DAH$value, n.pars = 3) #6
AIC.vec['DAN'] = AIC(neg.log.lk = lk.DAN$value, n.pars = 3) #7
AIC.vec['DAHN'] = AIC(neg.log.lk = lk.DAHN$value, n.pars = 4) #8
AIC.vec['DHmep'] = AIC(neg.log.lk = lk.DHmep$value, n.pars = 3) #9
AIC.vec['DAHmep'] = AIC(neg.log.lk = lk.DAHmep$value, n.pars = 4) #10
AIC.vec['DE'] = AIC(neg.log.lk = lk.null.E, n.pars = 0) #11
AIC.vec['DEH'] = AIC(neg.log.lk = lk.DEH$value, n.pars = 1) #12
AIC.vec['DEN'] = AIC(neg.log.lk = lk.DEN$value, n.pars = 1) #13
AIC.vec['DEA'] = AIC(neg.log.lk = lk.DEA$value, n.pars = 2) #14
AIC.vec['DEHN'] = AIC(neg.log.lk = lk.DEHN$value, n.pars = 2) #15
AIC.vec['DEAH'] = AIC(neg.log.lk = lk.DEAH$value, n.pars = 3) #16
AIC.vec['DEAN'] = AIC(neg.log.lk = lk.DEAN$value, n.pars = 3) #17
AIC.vec['DEAHN'] = AIC(neg.log.lk = lk.DEAHN$value, n.pars = 4) #18
AIC.vec['DEHmep'] = AIC(neg.log.lk = lk.DEHmep$value, n.pars = 3) #19
AIC.vec['DEAHmep'] = AIC(neg.log.lk = lk.DEAHmep$value, n.pars = 4) #20

#Sort
AIC.vec = AIC.vec[order(AIC.vec)]; AIC.vec

# Find delta AIC (best - AIC.vec)
del.AIC = AIC.vec-AIC.vec[1]; del.AIC
models = names(AIC.vec)





# ##  ____________________________________________
# ##  C. Profiles
# ##  ____________________________________________
exposure.switch = FALSE # re-load non-exposure inputs
source('Model_inputs.R')

## Profile for DAH:
# 1. Generate a sequence of Hm values to test:
# Choose a min and max value that span the range around the best estimate.
# Here, the best estimate is 0.268, so we are going to go from close to 0 up to .5
# Choose a reasonably small step size, but too small will make the code run forever.
Hms = seq(0.005, .5, by = .005)

#2. Initialize an empty vector to store neg log lk values at each Hm value
prof_DAH_Hm_nlls = numeric(length(Hms))

#3. Fix the Hm input value and optimize all other parameters.
#   Repeat for each Hm value in the above vector
#   Store the resulting nll in the appropriate vector
for(ii in 1:length(Hms)){
  # Note that this is the same as the likelihood optimization above with a few key changes:
  #   - Optimize lk.D.A.H.prof instead of lk.D.A.H
  #   - Input "Hm" as the prof.par
  #   - Input Hms[ii] (one of the Hm values above) as the "prof.par.val"
  #   - *Make sure you delete the boundaries corresponding to Hm from the lower and upper options. Otherwise, you'll get weird answers and your profile likelihoods will be way too high!
  #   - Output the value only (note that I added $value to the end of the whole big line of code.)
  prof_DAH_Hm_nlls[ii] = optim(par = c(Ac = 2, Ae = 2), fn = lk.D.A.H.prof, prof.par = 'Hm', prof.par.value = Hms[ii], wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(1, 1), upper = c(Inf, Inf))$value
}


## 4. Plot the nll values against the Hm vals
plot(Hms, prof_DAH_Hm_nlls, ylab = 'neg log lk', main = 'Profile on Hm\nModel DAH')
abline(h = lk.DAH$value+2, col = 'red') # add a red line 2 nll units above the best answer


### Ac for model DAH: best estimate at 2.0613034
Acs = seq(1.8, 2.3, by = .005) # values to test
prof_DAH_Ac_nlls = numeric(length(Acs)) # Empty vector
for(ii in 1:length(Acs)){
  prof_DAH_Ac_nlls[ii] = optim(par = c(Hm = .2, Ae = 2), fn = lk.D.A.H.prof, prof.par = 'Ac', prof.par.value = Acs[ii], wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(.001, 1), upper = c(1, Inf))$value
}
plot(Acs, prof_DAH_Ac_nlls, ylab = 'neg log lk', main = 'Profile on Ac\nModel DAH')
abline(h = lk.DAH$value+2) # add a red line 2 nll units above the best answer

# Ae for model DAH: best estimate at 1.000
Aes = seq(0.3, .8, by = .005) # values to test
prof_DAH_Ae_nlls = numeric(length(Aes)) # Empty vector
for(ii in 1:length(Aes)){
  prof_DAH_Ae_nlls[ii] = optim(par = c(Hm = .2, Ac = 2), fn = lk.D.A.H.prof, prof.par = 'Ae', prof.par.value = Aes[ii], wm = (weights.master.g1), wo = weights.master.g2, wn = weights.master.naiive, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(.001, 1), upper = c(1, Inf))$value
}
plot(Aes, prof_DAH_Ae_nlls, ylab = 'neg log lk', main = 'Profile on Ae\nModel DAH')
abline(h = lk.DAHmep$value+2) # add a red line 2 nll units above the best answer

####################################
## Profile for DAHmep

### Hm for DAHmep: best estimate at 1.0000
Hms = seq(0.01, 1, by = .005)
prof_DAHmep_Hm_nlls = numeric(length(Hms))
for(ii in 1:length(Hms)){
  prof_DAHmep_Hm_nlls[ii] = optim(par = c(e1 =.5, e2 = .5, Ac = 1.5, Ae = 1.5), fn = lk.D.A.H.ep.prof, prof.par = 'Hm', prof.par.value = Hms[ii], wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 1, 1), upper = c(1,1,Inf, Inf))$value
}
plot(Hms, prof_DAHmep_Hm_nlls, ylab = 'neg log lk', main = 'Profile on Hm\nModel DAHmep')
abline(h = lk.DAHmep$value+2, col = 'red') # add a red line 2 nll units above the best answer

### e1 for DAHmep: best estimate at 1.00000
e1s = seq(0.01, 1, by = .005)
prof_DAHmep_e1_nlls = numeric(length(e1s))
for(ii in 1:length(e1s)){
  prof_DAHmep_e1_nlls[ii] = optim(par = c(Hm = .2, e2 = .5, Ac = 1.5, Ae = 1.5), fn = lk.D.A.H.ep.prof, prof.par = 'e1', prof.par.value = e1s[ii], wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001,0.001,1, 1), upper = c(1, 1, Inf, Inf))$value
}
plot(e1s, prof_DAHmep_e1_nlls, ylab = 'neg log lk', main = 'Profile on e1\nModel DAHmep')
abline(h = lk.DAHmep$value+2, col = 'red') # add a red line 2 nll units above the best answer

### e2 for DAHmep: best estimate at 0.1987377
e2s = seq(0.005, .5, by = .005)
prof_DAHmep_e2_nlls = numeric(length(e2s))
for(ii in 1:length(e2s)){
  prof_DAHmep_e2_nlls[ii] = optim(par = c(Hm = .2, e1 = .5, Ac = 1.5, Ae = 1.5), fn = lk.D.A.H.ep.prof, prof.par = 'e2', prof.par.value = e2s[ii], wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(0.001, 0.001, 1, 1), upper = c(1, 1, Inf, Inf))$value
}
plot(e2s, prof_DAHmep_e2_nlls, ylab = 'neg log lk', main = 'Profile on e2\nModel DAHmep')
abline(h = lk.DAHmep$value+2, col = 'red') # add a red line 2 nll units above the best answer

### Ac for model DAHmep: best estimate at 1.9681725
Acs = seq(1.5,2.5, by = .005) # values to test
prof_DAHmep_Ac_nlls = numeric(length(Acs)) # Empty vector
for(ii in 1:length(Acs)){
  prof_DAHmep_Ac_nlls[ii] = optim(par = c(Hm = .2, e1 = .5, e2 = .5, Ae = 1.5), fn = lk.D.A.H.ep.prof, prof.par = 'Ac', prof.par.value = Acs[ii], wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(.001, 0.001, 0.001, 1), upper = c(1, 1, 1, Inf))$value
}
plot(Acs, prof_DAHmep_Ac_nlls, ylab = 'neg log lk', main = 'Profile on Ac\nModel DAHmep')
abline(h = lk.DAHmep$value+2) # add a red line 2 nll units above the best answer

# Ae for model DAHmep: best estimate at 1.000000
Aes = seq(0.5, 2, by = .005) # values to test
prof_DAHmep_Ae_nlls = numeric(length(Aes)) # Empty vector
for(ii in 1:length(Aes)){
  prof_DAHmep_Ae_nlls[ii] = optim(par = c(Hm = .2, e1 = .5, e2 = .5, Ac = 1.5), fn = lk.D.A.H.ep.prof, prof.par = 'Ae', prof.par.value = Aes[ii], wmm = nested_wmm, wmo = nested_wmo, wom = nested_wom, woo = nested_woo, wnn = nested_wnn, uc = Uc.master, ua = Um.master, ue = Ue.master, p0 = p0.master, xx = data.master, method = 'L-BFGS-B', lower = c(.001, 0.001,0.001, 1), upper = c(1, 1, 1, Inf))$value
}
plot(Aes, prof_DAHmep_Ae_nlls, ylab = 'neg log lk', main = 'Profile on Ae\nModel DAHmep')
abline(h = lk.DAHmep$value+2) # add a red line 2 nll units above the best answer


