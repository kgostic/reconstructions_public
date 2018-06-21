source('gete_ij_twohit.R')
## Inputs:
##    -> b_yr: Scalar. Birth year of interest
##    -> i_yr: Scalar. Year in which case was observed
##    -> flu_circ_ch: ???? Data on fraction of influenza A circulation caused by H1N1, H2N2 or H3N2 in the year of interest.


## Outputs:
## One list with the following elements
##  -> H1p_11: probs of exposure to any H1 virus @ 1st and second exposures
##  -> H1p_10: probs of exposure to any H1 virus @ 1st exposure, some other virus @ 2nd
##  -> H1p_01: probs of exposure to any H1 virus @ 2nd exposure, some other virus @ 1st
##  -> H1p_00: probs of exposure to a virus other than H1 @ 1st and 2nd exposures
##  -> H1p_1n: probs of exposure to any H1 virus @ 1st exposure and 2nd exposure hasn't yet occurred
##  -> H1p_on: probs of exposure to any virus other than H1 @ 1st exposure and 2nd exposure hasn't yet occurred
##  -> H2p_11: ... same as above for any H2 virus
##  -> H2p_10
##  -> H2p_01
##  -> H2p_00
##  -> H2p_1n
##  -> H2p_on
##  -> H3p_11 ... same as above for any H3 virus
##  -> H3p_10
##  -> H3p_01
##  -> H3p_00
##  -> H3p_1n
##  -> H3p_on
##  -> g1p_11 ... same as above for any HA group 1 virus (not repeated for g2, because H3 == g2)
##  -> g1p_10
##  -> g1p_01
##  -> g1p_00
##  -> g1p_1n
##  -> g1p_on
##  -> p_nn:  probs of no first or second exposure yet (naive to influenza) 


get_proportion = function(b_yr, i_yr, flu_circ_ch){
  func_result= get.e_ij(birth.year = b_yr, incidence.year = i_yr) # Outputs the probabilities of a first exposure in year i and second exposure in year j, as well as the probabilites of no first exposure, and of a first exposure in year i, but no second exposure.
  prob.matrix = func_result$total_mat # matrix: row = year of 1st exposure, column = year of 2nd exposure
  p_nn = func_result$totalPof_noexp1 # scalar
  no_exp2 = func_result$totalPof_yesexp1_and_noexp2 # vector: element = year of 1st exposure.
  
  nn = length(b_yr:min(b_yr+17, i_yr)) ## total number of years of possible first exposure years
  
  # Store probs of naive in 2nd exposure
  H1p_1n = H1p_0n = H2p_1n = H2p_0n = H3p_0n = H3p_1n = g1p_1n = g1p_0n = vector('numeric', nn)
  over_max_age = 0
  
  if(i_yr - b_yr > 17){ # If the birth cohorts is over age 17 at the time of observation
    # Assume everyone has had both a first and second exposure
    # Normalize the matrix so that probabilities sum to 1, this implicitly sets naive probs to 0
    prob.matrix = prob.matrix/sum(prob.matrix, na.rm = TRUE)
    over_max_age = 18
  } # Else, assume some individuals are still missing an first or second exposure, and calculate the naive fraction below
  
  ## Calculate probs of first and second exposure for H1N1
    H1p_11 = H1p_10 = H1p_01 = H1p_00 = matrix(nrow=nn, ncol=nn) # empty matrices
    for(exp1 in 1:nn){ # point at a specific year index in vector exp1
      mm = nn - exp1 # length of all poss exp2 years
      e1_H1 = flu_circ_ch[1,toString(exp1+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of first exposure
      if(mm>0){
        for(exp2 in (1+exp1):(mm+exp1)) { # Point at specific year index in vector exp2
          e2_H1 = flu_circ_ch[1,toString(exp2+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of second exposure
          #prob.matrix gives prob of any [exp1, exp2] combination of years
          #multiply by the probabilities of group 1 exposure in those years to calculate the probability of first and/or second exposure to a specific HA group
          H1p_11[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_H1) * (e2_H1)
          H1p_10[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_H1) * (1-e2_H1)
          H1p_01[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_H1) * (e2_H1)
          H1p_00[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_H1) * (1-e2_H1)
        }
      }
      if(over_max_age == 0){
        H1p_1n[exp1] = (e1_H1 * (no_exp2[exp1])) # wmn overall prob of H1 imprinting on 1st exp, but no second exposure
        H1p_0n[exp1] = ((1-e1_H1) * (no_exp2[exp1])) # won overall prob of H1 imprinting on 1st exp, but no second exposure
      }
      else{
        H1p_1n[exp1] = 0
        H1p_0n[exp1] = 0
        p_nn = 0
      }
    } # End loop for H1N1
   
    
    
    ########## Repeat for H2N2 childhood exposures
    # Store probs of naive in 2nd exposure
    over_max_age = 0
    
    if(i_yr - b_yr > 17){ # If the birth cohorts is over age 17 at the time of observation
      # Assume everyone has had both a first and second exposure
      # Normalize the matrix so that probabilities sum to 1, this implicitly sets naive probs to 0
      prob.matrix = prob.matrix/sum(prob.matrix, na.rm = TRUE)
      over_max_age = 18
    } # Else, assume some individuals are still missing an first or second exposure, and calculate the naive fraction below
    
    ## Calculate probs of first and second exposure for H2N1
    H2p_11 = H2p_10 = H2p_01 = H2p_00 = matrix(nrow=nn, ncol=nn) # empty matrices
    for(exp1 in 1:nn){ # point at a specific year index in vector exp1
      mm = nn - exp1 # length of all poss exp2 years
      e1_H2 = flu_circ_ch[2,toString(exp1+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of first exposure
      if(mm>0){
        for(exp2 in (1+exp1):(mm+exp1)) { # Point at specific year index in vector exp2
          e2_H2 = flu_circ_ch[2,toString(exp2+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of second exposure
          #prob.matrix gives prob of any [exp1, exp2] combination of years
          #multiply by the probabilities of group 1 exposure in those years to calculate the probability of first and/or second exposure to a specific HA group
          H2p_11[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_H2) * (e2_H2)
          H2p_10[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_H2) * (1-e2_H2)
          H2p_01[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_H2) * (e2_H2)
          H2p_00[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_H2) * (1-e2_H2)
        }
      }
      if(over_max_age == 0){
        H2p_1n[exp1] = (e1_H2 * (no_exp2[exp1])) # wmn overall prob of H2 imprinting on 1st exp, but no second exposure
        H2p_0n[exp1] = ((1-e1_H2) * (no_exp2[exp1])) # won overall prob of H2 imprinting on 1st exp, but no second exposure
      }
      else{
        H2p_1n[exp1] = 0
        H2p_0n[exp1] = 0
        p_nn = 0
      }
    } # End loop for H2N2
    
    
    
    
    
    # Store probs of naive in 2nd exposure
    over_max_age = 0
    
    if(i_yr - b_yr > 17){ # If the birth cohorts is over age 17 at the time of observation
      # Assume everyone has had both a first and second exposure
      # Normalize the matrix so that probabilities sum to 1, this implicitly sets naive probs to 0
      prob.matrix = prob.matrix/sum(prob.matrix, na.rm = TRUE)
      over_max_age = 18
    } # Else, assume some individuals are still missing an first or second exposure, and calculate the naive fraction below
    
    ## Calculate probs of first and second exposure for H3N2
    H3p_11 = H3p_10 = H3p_01 = H3p_00 = matrix(nrow=nn, ncol=nn) # empty matrices
    for(exp1 in 1:nn){ # point at a specific year index in vector exp1
      mm = nn - exp1 # length of all poss exp2 years
      e1_H3 = flu_circ_ch[3,toString(exp1+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of first exposure
      if(mm>0){
        for(exp2 in (1+exp1):(mm+exp1)) { # Point at specific year index in vector exp2
          e2_H3 = flu_circ_ch[3,toString(exp2+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of second exposure
          #prob.matrix gives prob of any [exp1, exp2] combination of years
          #multiply by the probabilities of group 1 exposure in those years to calculate the probability of first and/or second exposure to a specific HA group
          H3p_11[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_H3) * (e2_H3)
          H3p_10[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_H3) * (1-e2_H3)
          H3p_01[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_H3) * (e2_H3)
          H3p_00[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_H3) * (1-e2_H3)
        }
      }
      if(over_max_age == 0){
        H3p_1n[exp1] = (e1_H3 * (no_exp2[exp1])) # wmn overall prob of H3 imprinting on 1st exp, but no second exposure
        H3p_0n[exp1] = ((1-e1_H3) * (no_exp2[exp1])) # won overall prob of H3 imprinting on 1st exp, but no second exposure
      }
      else{
        H3p_1n[exp1] = 0
        H3p_0n[exp1] = 0
        p_nn = 0
      }
    } # End loop for H3N2
    
    
    # Store probs of naive in 2nd exposure
    over_max_age = 0
    
    if(i_yr - b_yr > 17){ # If the birth cohorts is over age 17 at the time of observation
      # Assume everyone has had both a first and second exposure
      # Normalize the matrix so that probabilities sum to 1, this implicitly sets naive probs to 0
      prob.matrix = prob.matrix/sum(prob.matrix, na.rm = TRUE)
      over_max_age = 18
    } # Else, assume some individuals are still missing an first or second exposure, and calculate the naive fraction below
    
    ## Calculate probs of first and second exposure for g1N2
    g1p_11 = g1p_10 = g1p_01 = g1p_00 = matrix(nrow=nn, ncol=nn) # empty matrices
    for(exp1 in 1:nn){ # point at a specific year index in vector exp1
      mm = nn - exp1 # length of all poss exp2 years
      e1_g1 = flu_circ_ch[4,toString(exp1+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of first exposure
      if(mm>0){
        for(exp2 in (1+exp1):(mm+exp1)) { # Point at specific year index in vector exp2
          e2_g1 = flu_circ_ch[4,toString(exp2+b_yr-1)] # Extract fraction of all flu A circulation caused by group 1 in the year of second exposure
          #prob.matrix gives prob of any [exp1, exp2] combination of years
          #multiply by the probabilities of group 1 exposure in those years to calculate the probability of first and/or second exposure to a specific HA group
          g1p_11[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_g1) * (e2_g1)
          g1p_10[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_g1) * (1-e2_g1)
          g1p_01[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_g1) * (e2_g1)
          g1p_00[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_g1) * (1-e2_g1)
        }
      }
      if(over_max_age == 0){
        g1p_1n[exp1] = (e1_g1 * (no_exp2[exp1])) # wmn overall prob of g1 imprinting on 1st exp, but no second exposure
        g1p_0n[exp1] = ((1-e1_g1) * (no_exp2[exp1])) # won overall prob of g1 imprinting on 1st exp, but no second exposure
      }
      else{
        g1p_1n[exp1] = 0
        g1p_0n[exp1] = 0
        p_nn = 0
      }
    } # End loop for H3N2
    
    
# return
    result <- list("H1p_00" = H1p_00, "H1p_01" = H1p_01, "H1p_0n" = H1p_0n, "H1p_10" = H1p_10, "H1p_11" = H1p_11, "H1p_1n" = H1p_1n, "H2p_00" = H2p_00, "H2p_01" = H2p_01, "H2p_0n" = H2p_0n, "H2p_10" = H2p_10, "H2p_11" = H2p_11, "H2p_1n" = H2p_1n, "H3p_00" = H3p_00, "H3p_01" = H3p_01, "H3p_0n" = H3p_0n,  "H3p_10" = H3p_10, "H3p_11" = H3p_11, "H3p_1n" = H3p_1n, "g1p_00" = g1p_00, "g1p_01" = g1p_01, "g1p_0n" = g1p_0n, "g1p_10" = g1p_10, "g1p_11" = g1p_11, "g1p_1n" = g1p_1n, "p_nn" = p_nn)
    result
}




# Test
setwd('~/Dropbox/R/Reconstructions/')
source('Regina two hit/country_data_import.R')
flu_circ_input = import.country.dat('Egypt')
## Test for a birth year and observation year less than 17 years apart
# Tests below should sum to 1
test = get_proportion(b_yr = 2000, i_yr = 2009, flu_circ_ch = flu_circ_input)
attach(test)
sum(H1p_11, na.rm = T)+sum(H1p_10, na.rm = T)+sum(H1p_01, na.rm = T)+sum(H1p_00, na.rm = T)+sum(H1p_1n, na.rm = T)+sum(H1p_0n, na.rm = T)+sum(p_nn)
sum(H2p_11, na.rm = T)+sum(H2p_10, na.rm = T)+sum(H2p_01, na.rm = T)+sum(H2p_00, na.rm = T)+sum(H2p_1n, na.rm = T)+sum(H2p_0n, na.rm = T)+sum(p_nn)
sum(H3p_11, na.rm = T)+sum(H3p_10, na.rm = T)+sum(H3p_01, na.rm = T)+sum(H3p_00, na.rm = T)+sum(H3p_1n, na.rm = T)+sum(H3p_0n, na.rm = T)+sum(p_nn)
sum(g1p_11, na.rm = T)+sum(g1p_10, na.rm = T)+sum(g1p_01, na.rm = T)+sum(g1p_00, na.rm = T)+sum(g1p_1n, na.rm = T)+sum(g1p_0n, na.rm = T)+sum(p_nn)
detach(test)

## Test for a birth year and observation year more than 17 years apart
# Tests below should sum to 1
test = get_proportion(b_yr = 1981, i_yr = 2009, flu_circ_ch = flu_circ_input)
attach(test)
sum(H1p_11, na.rm = T)+sum(H1p_10, na.rm = T)+sum(H1p_01, na.rm = T)+sum(H1p_00, na.rm = T)+sum(H1p_1n, na.rm = T)+sum(H1p_0n, na.rm = T)+sum(p_nn)
sum(H2p_11, na.rm = T)+sum(H2p_10, na.rm = T)+sum(H2p_01, na.rm = T)+sum(H2p_00, na.rm = T)+sum(H2p_1n, na.rm = T)+sum(H2p_0n, na.rm = T)+sum(p_nn)
sum(H3p_11, na.rm = T)+sum(H3p_10, na.rm = T)+sum(H3p_01, na.rm = T)+sum(H3p_00, na.rm = T)+sum(H3p_1n, na.rm = T)+sum(H3p_0n, na.rm = T)+sum(p_nn)
sum(g1p_11, na.rm = T)+sum(g1p_10, na.rm = T)+sum(g1p_01, na.rm = T)+sum(g1p_00, na.rm = T)+sum(g1p_1n, na.rm = T)+sum(g1p_0n, na.rm = T)+sum(p_nn)
detach(test)
