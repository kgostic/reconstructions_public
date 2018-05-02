#############################################
# -----
#####     LIKELIHOOD FUNCTION - Broad immune effects
# -----
#############################################
#calculate the likelihood given demography, exposure and immune effects

##          ----- INPUTS ----
# wm --  matrix describing the fraction of matched (protective) HA imprinting in 
# -----  each birth cohort (columns) and country-year in which the challenge strain has been observed (rows)
# wn --  matrix describing the fraction of other (non-protective) HA imprinting 
# -----  first IAV exposures in each birth cohort (columns) and country-year in which the challenge strain has been observed (rows)
# wn --  matrix describing the fraction of naive children in each birth cohort (columns) and country-year 
# -----  in which the challenge strain has been observed (rows)
# p0 --  matrix describing demographic age distribution in each unique country-year (rows) of case 
# -----  observation, OR in models that consider poultry exposure, p0 describes the demographic 
# -----  age distribution scaled by age-specific poultry exposure risk
# xx --  matrix describing the number of observed cases or deaths in each unique country-year of 
# -----  case observation (rows), and in each birth year (columns)
# uc, ue, um --  matrices indicating membership in a specific age-risk group (uc = children ages 0-4, 
# -----  um = ages 5-64 and ue = elderly ages 65+). Rows describe unique country years of case 
# -----  observaiton and columns indicate birth year.
# pars -- a named vector of parameters, including Hm, Nm, Ac and Ae, as described in Extended Data Table 1.
# ------------ Hm indicates relative risk for those with matched HA imprinting
# ------------ He1m indicates relative risk for those with matched HA e1 imprinting
# ------------ He2m indicates relative risk for those with matched HA e1 imprinting
# ------------ Ac indicates relative risk for children ages 0-4
# ------------ Ae indicates relative risk for the elderly ages 65+


##          ----- OUTPUTS ----
# Negative log likelihood of the model with specified parameters. 
# optim() can use these functions to compute maximum likelihood estimates

## ___________________________
## 1. D - run any model with all pars fixed at 0
## ---------------------------


## ___________________________
## 2. DH - Hemagglutinin history only
## ---------------------------
lk.D.H = function(pars, wm, wo, wn, p0, xx){
  
  # 1. Assign parameters to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H first exposure
  #Ho = pars['Ho']  #Relative susceptibility of the naiive
  
  # 2. calculate p_i as a function of the parameters:
  # This step gives the model prediction
  pp = p0*(Hm*wm+wo+wn)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE) #This line returns the log multinomial density of the observed data, with expected probabilities governed by model predictions.
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}


## ___________________________
## 3. DN - Neuraminidase history only
## ---------------------------
lk.D.N = function(pars, wm, wo, wn, p0, xx){
  
  # 1. Assign parameters to be fit
  Nm = pars['Nm']    #Relative susceptibility atrributable to H first exposure

  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Nm*wm+wo+wn)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}



## ___________________________
## 4. DA - Age risk group only
## ---------------------------
lk.D.A = function(pars, uc, um, ue, p0, xx){
  
  # 1. Assign parameters to be fit
  Ac = pars['Ac']   #Relative susceptibility of young children
  Ae = pars['Ae']   #Relative susceptibility of the elderly
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Ac*uc+um+Ae*ue) 
  
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  
  lk
  # end function 
}



## _________________________________________
## 5. DHN - Hemagglutinin and neuraminidase history
## -----------------------------------------
lk.D.H.N = function(pars, wmm, wmo = 0, wom = 0, woo, wnn, p0, xx){
  ##            ---- INPUTS ----
  # wmm -- fraction of each birth cohort with matched imprinting to H and N
  # wmo -- fraction of each birth cohort with matched imprinting to H and mismatched imprinting to N
  # wom --fraction of each birth cohort with mismatched imprinting to H and matched imprinting to N
  # woo -- fraction of each birth cohort with mismatched imprinting to H and N
  # wnn -- fraction of each birth cohort naiive to influenza A
  # p0  -- demographic age distributions, or demographic age distributions scaled by poultry exposure in different age groups
  # xx -- the number of cases or deaths in each birth year
  
  # 1. Assign parameters to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H first exposure
  Nm = pars['Nm']    #Ditto for first N exposure
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Hm*Nm*wmm + Hm*wmo + Nm*wom + woo + wnn)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}


## _________________________________________
## 6. DAH - Age risk and hemagglutinin history
## -----------------------------------------
lk.D.A.H = function(pars, wm, wo, wn, uc, ua, ue, p0, xx){
  
  # 1. Assign parameters to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H first exposure
  Ac = pars['Ac']    #Relative susceptibility of children 0-4
  Ae = pars['Ae']    #Relative susceptibility of the elderly 65+
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Hm*wm+wo+wn)*(Ac*uc+ua+Ae*ue)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}



## _________________________________________
## 7. DAN - Age risk and neuraminidase history
## -----------------------------------------
lk.D.A.N = function(pars, wm, wo, wn, uc, ua, ue, p0, xx){
  
  # 1. Assign parameters to be fit
  Nm = pars['Nm']    #Relative susceptibility atrributable to H first exposure
  Ac = pars['Ac']    #Relative susceptibility of children 0-4
  Ae = pars['Ae']    #Relative susceptibility of the elderly 65+
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Nm*wm+wo+wn)*(Ac*uc+ua+Ae*ue)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}


## _________________________________________
## 8. DAHN - Age risk, hemagglutinin and neuraminidase history
## -----------------------------------------
lk.D.A.H.N = function(pars, wmm, wmo = 0, wom = 0, woo, wnn, uc, ua, ue, p0, xx){
  ##            ---- INPUTS ----
  # wmm -- fraction of each birth cohort with matched imprinting to H and N
  # wmo -- fraction of each birth cohort with matched imprinting to H and mismatched imprinting to N
  # wom --fraction of each birth cohort with mismatched imprinting to H and matched imprinting to N
  # woo -- fraction of each birth cohort with mismatched imprinting to H and N
  # wnn -- fraction of each birth cohort naiive to influenza A
  # p0  -- demographic age distributions, or demographic age distributions scaled by poultry exposure in different age groups
  # xx -- the number of cases or deaths in each birth year
  
  # 1. Assign parameters to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H first exposure
  Nm = pars['Nm']    #Relative susceptibility of those with matched first N exposure
  Ac = pars['Ac']    #Relative susceptibility of children 0-4
  Ae = pars['Ae']    #Relative susceptibility of the elderly, 65+
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Hm*Nm*wmm + Hm*wmo + Nm*wom + woo + wnn)*(Ac*uc+ua+Ae*ue)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}

## _________________________________________
## 9. DHmep - Hemagglutinin exposure 1 and exposure 2 history
## -----------------------------------------
lk.D.H.ep = function(pars, wmm, wmo, wom, woo, wnn, p0, xx){
  ##            ---- INPUTS ----
  # wmm -- fraction of each birth cohort with matched imprinting to H e1 and e2
  # wmo -- fraction of each birth cohort with matched imprinting to H e1 and mismatched imprinting to e2
  # wom --fraction of each birth cohort with mismatched imprinting to H e1 and matched imprinting to e2
  # woo -- fraction of each birth cohort with mismatched imprinting to H e1 and e2
  # wnn -- fraction of each birth cohort naiive to influenza A
  # p0  -- demographic age distributions, or demographic age distributions scaled by poultry exposure in different age groups
  # xx -- the number of cases or deaths in each birth year
  
  # 1. Assign parameters to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to any H exp
  e1 = pars['e1']    #Additonal bonus for 1st exp matched
  e2 = pars['e2']    #Additional bonus for wmm (consecutive match)
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Hm*e1*e2*wmm + Hm*e1*wmo + Hm*wom + (woo + wnn))
  #pp = pp/rowSums(pp)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}

## _________________________________________
## 10. DAHmep - Age risk, hemagglutinin exposure 1 and exposure 2 history
## -----------------------------------------
lk.D.A.H.ep = function(pars, wmm, wmo, wom, woo, wnn, uc, ua, ue, p0, xx){
  ##            ---- INPUTS ----
  # wmm -- fraction of each birth cohort with matched imprinting to H e1 and e2
  # wmo -- fraction of each birth cohort with matched imprinting to H e1 and mismatched imprinting to e2
  # wom --fraction of each birth cohort with mismatched imprinting to H e1 and matched imprinting to e2
  # woo -- fraction of each birth cohort with mismatched imprinting to H e1 and e2
  # wnn -- fraction of each birth cohort naiive to influenza A
  # p0  -- demographic age distributions, or demographic age distributions scaled by poultry exposure in different age groups
  # xx -- the number of cases or deaths in each birth year
  
  # 1. Assign parameters to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H
  e1 = pars['e1']    #Additonal bonus for 1st exp matched
  e2 = pars['e2']    #Additional bonus for wmm (consecutive match)
  Ac = pars['Ac']    #Relative susceptibility of children 0-4
  Ae = pars['Ae']    #Relative susceptibility of the elderly, 65+
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Hm*e1*e2*wmm + Hm*e1*wmo + Hm*wom + (woo + wnn))*(Ac*uc+ua+Ae*ue)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}


#### We will talk about the stuff below later!!
# 
# 
# 
# ##############
# ## Modified likelihood functions for calculating likelihood profiles
# ##############
# 
# ## ___________________________
# ## 2. DH - Hemagglutinin history only
# ## ---------------------------
# lk.D.H.prof = function(pars, prof.par, prof.par.value, wm, wo, wn, p0, xx){
# 
#   # 1. Assign parameters to be fit
#   Hm = ifelse('Hm' == prof.par, prof.par.value, pars['Hm'])
# 
#   # 2. calculate p_i as a function of the parameters:
#   pp = p0*(Hm*wm+wo+wn)
# 
#   #  3. Likelihood is based on the multinomial density
#   if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
#     lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
#   }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
#     storage = vector('numeric', dim(xx)[1])
#     for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
#       storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
#     }
#     lk = sum(storage)
#   }
#   lk # end function
# }
# 
# 
# ## ___________________________
# ## 3. DN - Neuraminidase history only
# ## ---------------------------
# lk.D.N.prof = function(pars, prof.par, prof.par.value, wm, wo, wn, p0, xx){
#   
#   # 1. Assign parameters to be fit
#   Nm = ifelse(prof.par == 'Nm', prof.par.value, pars['Nm'])
#   
#   # 2. calculate p_i as a function of the parameters:
#   pp = p0*(Nm*wm+wo+wn)
#   
#   #  3. Likelihood is based on the multinomial density
#   if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
#     lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
#   }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
#     storage = vector('numeric', dim(xx)[1])
#     for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
#       storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
#     }
#     lk = sum(storage) 
#   }
#   lk # end function 
# }
# 
# 
# 
# ## ___________________________
# ## 4. DA - Age risk group only
# ## ---------------------------
# lk.D.A.prof = function(pars, prof.par, prof.par.value, uc, um, ue, p0, xx){
#   
#   # 1. Assign parameters to be fit
#   Ac = ifelse(prof.par == 'Ac', prof.par.value, pars['Ac'])   
#   Ae = ifelse(prof.par == 'Ae', prof.par.value, pars['Ae'])   
#   
#   # 2. calculate p_i as a function of the parameters:
#   pp = p0*(Ac*uc+um+Ae*ue) 
# 
#   
#   
#   #  3. Likelihood is based on the multinomial density
#   if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
#     lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
#   }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
#     storage = vector('numeric', dim(xx)[1])
#     for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
#       storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
#     }
#     lk = sum(storage) 
#   }
#   
#   lk
#   # end function 
# }
# 
# 
# 
# ## _________________________________________
# ## 5. DHN - Hemagglutinin and neuraminidase history
# ## -----------------------------------------
# lk.D.H.N.prof = function(pars, prof.par, prof.par.value, wmm, wmo = 0, wom = 0, woo, wnn, p0, xx){
#   
#   # 1. Assign parameters to be fit
#   Hm = ifelse(prof.par == 'Hm', prof.par.value, pars['Hm'])
#   Nm = ifelse(prof.par == 'Nm', prof.par.value, pars['Nm'])
# 
#   
#   # 2. calculate p_i as a function of the parameters:
#   pp = p0*(Hm*Nm*wmm + Hm*wmo + Nm*wom + woo + wnn)
#   
#   #  3. Likelihood is based on the multinomial density
#   if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
#     lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
#   }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
#     storage = vector('numeric', dim(xx)[1])
#     for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
#       storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
#     }
#     lk = sum(storage) 
#   }
#   lk # end function 
# }
# 
# 
## _________________________________________
## 6. DAH - Age risk and hemagglutinin history
## -----------------------------------------
lk.D.A.H.prof = function(pars, prof.par, prof.par.value, wm, wo, wn, uc, ua, ue, p0, xx){
  
  # 1. Assign parameters to be fit
  Hm = ifelse(prof.par == 'Hm', prof.par.value, pars['Hm'])   #Relative susceptibility atrributable to H first exposure
  Ac = ifelse(prof.par == 'Ac', prof.par.value, pars['Ac'])
  Ae = ifelse(prof.par == 'Ae', prof.par.value, pars['Ae'])
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Hm*wm+wo+wn)*(Ac*uc+ua+Ae*ue)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}

# 
# ## _________________________________________
# ## 7. DAN - Age risk and neuraminidase history
# ## -----------------------------------------
# lk.D.A.N.prof = function(pars, prof.par, prof.par.value, wm, wo, wn, uc, ua, ue, p0, xx){
#   
#   # 1. Assign parameters to be fit
#   Nm = ifelse(prof.par == 'Nm', prof.par.value, pars['Nm'])
#   Ac = ifelse(prof.par == 'Ac', prof.par.value, pars['Ac'])
#   Ae = ifelse(prof.par == 'Ae', prof.par.value, pars['Ae'])
#   
#   # 2. calculate p_i as a function of the parameters:
#   pp = p0*(Nm*wm+wo+wn)*(Ac*uc+ua+Ae*ue)
#   
#   #  3. Likelihood is based on the multinomial density
#   if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
#     lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
#   }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
#     storage = vector('numeric', dim(xx)[1])
#     for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
#       storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
#     }
#     lk = sum(storage) 
#   }
#   lk # end function 
# }
# 
# 
# 
# ## _________________________________________
# ## 8. DAHN - Age risk, hemagglutinin and neuraminidase history
# ## -----------------------------------------
# lk.D.A.H.N.prof = function(pars, prof.par, prof.par.value, wmm, wmo = 0, wom = 0, woo, wnn, uc, ua, ue, p0, xx){
#   
#   # 1. Assign parameters to be fit
#   Hm = ifelse(prof.par == 'Hm', prof.par.value, pars['Hm'])    #Relative susceptibility atrributable to H first exposure
#   Nm = ifelse(prof.par == 'Nm', prof.par.value, pars['Nm'])   #Relative susceptibility of those with matched first N exposure
#   Ac = ifelse(prof.par == 'Ac', prof.par.value, pars['Ac'])    #Relative susceptibility of children 0-4
#   Ae = ifelse(prof.par == 'Ae', prof.par.value, pars['Ae'])    #Relative susceptibility of the elderly, 65+
#   
#   # 2. calculate p_i as a function of the parameters:
#   pp = p0*(Hm*Nm*wmm + Hm*wmo + Nm*wom + woo + wnn)*(Ac*uc+ua+Ae*ue)
#   
#   #  3. Likelihood is based on the multinomial density
#   if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
#     lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
#   }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
#     storage = vector('numeric', dim(xx)[1])
#     for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
#       storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
#     }
#     lk = sum(storage) 
#   }
#   lk # end function 
# }
# 
# 
## _________________________________________
## 10. DAHmep - Age risk, hemagglutinin exposure 1 and exposure 2 history
## -----------------------------------------
lk.D.A.H.ep.prof = function(pars, prof.par, prof.par.value, wmm, wmo, wom, woo, wnn, uc, ua, ue, p0, xx){
  
  # 1. Assign parameters to be fit
  Hm = ifelse(prof.par == 'Hm', prof.par.value, pars['Hm'])
  e1 = ifelse(prof.par == 'e1', prof.par.value, pars['e1'])
  e2 = ifelse(prof.par == 'e2', prof.par.value, pars['e2'])
  Ac = ifelse(prof.par == 'Ac', prof.par.value, pars['Ac'])
  Ae = ifelse(prof.par == 'Ae', prof.par.value, pars['Ae'])
  
  # 2. calculate p_i as a function of the parameters:
  pp = p0*(Hm*e1*e2*wmm + Hm*e1*wmo + Hm*wom + (woo + wnn))*(Ac*uc+ua+Ae*ue)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}

# 
# 
# 
# 
# 
# 
# 
# 
# 
# #############################################
# # -----
# #####     FUNCTIONS FOR CALCULATING CIs FROM LIKELIHOOD SLICES
# # -----
# #############################################
# 
# 
# 
# 
# ## INPUTS
# #      NLL_best is a numeric input that gives the negative log likelihood of the unconstrained model using maximum likelihood parameter estimates
# #  df is the number of degrees of freedom, or the number of paramters constrained in profile calculation. For likelihood slices, df = 1
# ## OUTPUTS
# #     threshold gives the threshold likelihood value. Paramter values falling below this threshold in the likelihood slice or profile will be contained in the 95% CI
# LR.Threshold = function(NLL_best, df){
#   threshold = qchisq(.95, df)/2+NLL_best
#   threshold
# }
# 
# 
# 
# ## INPUTS
# #      threshold is a numeric input, described above
# #      nll.vec represents the negative log likelihood values in the likelihood slice or profile 
# #      pars.vec is a vector of parameter values corresponding to the likelihood slice values in nll.vec
# ## OUTPUTS
# #     CI gives the 95% CI on the constrained parameter
# LR.CI = function(threshold, nll.vec, pars.vec){
#   if(any(which(nll.vec > threshold) < which.min(nll.vec))  &  any(which(nll.vec > threshold) > which.min(nll.vec)) ){ #If the minimum is not an end point
#     #Find string before and after the min value
#     lower = nll.vec[1:which.min(nll.vec)]
#     upper = nll.vec[(which.min(nll.vec)-1):length(nll.vec)]
#     #Extract low CI from first string, and upper CI from upper string
#     CI = c(pars.vec[which.min(abs(lower-threshold))], pars.vec[length(lower)+which.min(abs(upper-threshold))])
#   }else{
#     #If the first value is the minimum
#     if(any(which(nll.vec > threshold) < which.min(nll.vec)) == FALSE){ CI = c(pars.vec[1], pars.vec[which.min(abs(nll.vec-threshold))]) }
#     if(any(which(nll.vec > threshold) > which.min(nll.vec)) == FALSE){ CI = c(pars.vec[which.min(abs(nll.vec-threshold))], pars.vec[length(nll.vec)])}
#   }
#   CI
# }
# 
