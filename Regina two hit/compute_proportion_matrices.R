source('regina_pt3code.R')

get_proportion = function(b_yr, i_yr, flu_circ_ch){
   # flu_circ_ch = read.csv('China_flu_circulation.csv'); colnames(flu_circ_ch) = c('',2017:1901); rownames(flu_circ_ch) = c('H1','H2','H3','Group1','Group2');
   # rownames(flu_circ_ch) <- flu_circ_ch[,1]
   # flu_circ_ch <- flu_circ_ch[,-1]
  # 
  ## Use existing function to output probabilities that people born in 1995 had their first and second exposures in 1995, 1996, 1997...etc.
  func_result= get.e_ij(birth.year = b_yr, incidence.year = i_yr) # Tuple
  prob.matrix = func_result$total_mat
  p_nn = func_result$totalPof_noexp1 # same as no_exp1
  no_exp2 = func_result$totalPof_yesexp1_and_noexp2
  
  nn = length(b_yr:min(b_yr+17, i_yr))
  naive = matrix(0L,nrow=nn, ncol=nn) #incident year successfully starts with naive = 1
  
  # Store probs of naive in 2nd exposure
  g1p_1n = g1p_0n = g2p_1n = g2p_0n = vector('numeric', nn)
  over_max_age = FALSE
  
  ## Normalize so prob.matrix sums to 1
  if(i_yr - b_yr > 17){ # This birth cohorts is currently over the max age
    prob.matrix = prob.matrix/sum(prob.matrix, na.rm = TRUE)
    over_max_age = 18
  #}else{ # This birth cohort is currently under the minumum age
    #naive = 1 - sum(prob.matrix, na.rm = TRUE)
  }
  ## Calculate probs of first and second exposure for H1N1 and H2N2 (Group 1)
    g1p_11 = g1p_10 = g1p_01 = g1p_00 = matrix(nrow=nn, ncol=nn)
    
    for(exp1 in 1:nn){ # point at a specific year index in vector exp1
      mm = nn - exp1 # length of all poss exp2 years
      
      e1_g1 = flu_circ_ch[4,toString(exp1+b_yr-1)] # Probability of that year for H1 AND H2
      if(mm>0){
        for(exp2 in (1+exp1):(mm+exp1)) { # Point at specific year index in vector exp2
          e2_g1 = flu_circ_ch[4,toString(exp2+b_yr-1)]
          
          #prob.matrix supposed to give exp to ANY strain
          g1p_11[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_g1) * (e2_g1)
          g1p_10[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_g1) * (1-e2_g1)
          g1p_01[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_g1) * (e2_g1)
          g1p_00[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_g1) * (1-e2_g1)
        }
      }
      if(over_max_age == 0){
        g1p_1n[exp1] = (e1_g1 * (no_exp2[exp1])) # wmn overall prob of g1 imprinting on 1st exp
        g1p_0n[exp1] = ((1-e1_g1) * (no_exp2[exp1])) # won overall prob of g1 imprinting on 1st exp
      }
      else{
        over_max_age= over_max_age - 1; # decrease over_max_age
        g1p_1n[exp1] = 0
        g1p_0n[exp1] = 0
        p_nn = 0
      }
    } # End loop over years of 1st exposure
   
  #########################################################################
    
  ## Calculate probs of first and second exposure for H2N2 (updated version already implies this within Group1)
    # 
    # H2p_11 = H2p_10 = H2p_01 = H2p_00 = matrix(nrow=nn, ncol=nn)
    # for(exp1 in 1:nn){ # point at a specific year index in vector exp1
    #   mm = nn - exp1 #length of all poss exp2 years
    # 
    #   e1_h2n2 = flu_circ_ch[2,toString(exp1+b_yr-1)] #constant, obtain prob of that year
    #   if(mm>0){
    #     for(exp2 in (1+exp1):(mm+exp1)) { # point at a specific year index in vector exp2
    #       e2_h2n2 = flu_circ_ch[2,toString(exp2+b_yr-1)]
    # 
    #       #prob.matrix supposed to give exp to ANY strain
    #       H2p_11[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_h2n2) * (e2_h2n2)
    #       H2p_10[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_h2n2) * (1-e2_h2n2)
    #       H2p_01[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_h2n2) * (e2_h2n2)
    #       H2p_00[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_h2n2) * (1-e2_h2n2)
    #     }
    #   }
    # }

  ##########################################################################
  ## Calculate probs of first and second exposure for H3N2 (Group 2)
    # reset no_exp1 and counter
    p_nn = func_result$totalPof_noexp1
    if(i_yr - b_yr > 17){ # This birth cohorts is currently over the max age
      prob.matrix = prob.matrix/sum(prob.matrix, na.rm = TRUE)
      over_max_age = 18
    }
    
    g2p_11 = g2p_10 = g2p_01 = g2p_00 = matrix(nrow=nn, ncol=nn)
    
    for(exp1 in 1:nn){ # point at a specific year index in vector exp1
      mm = nn - exp1 # length of all poss exp2 years
      e1_g2 = flu_circ_ch[5,toString(exp1+b_yr-1)] # constant, obtain prob of that year
      if(mm>0){
        for(exp2 in (1+exp1):(mm+exp1)) { # point at a specific year index in vector exp2
          e2_g2 = flu_circ_ch[5,toString(exp2+b_yr-1)]
          
          g2p_11[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_g2) * (e2_g2)
          g2p_10[exp1,exp2] = prob.matrix[exp1, exp2] * (e1_g2) * (1-e2_g2)
          g2p_01[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_g2) * (e2_g2)
          g2p_00[exp1,exp2] = prob.matrix[exp1, exp2] * (1-e1_g2) * (1-e2_g2)
        } # End loop over years of second exposure
      }
      if(over_max_age == 0){
        g2p_1n[exp1] = (e1_g2 * (no_exp2[exp1])) # wmn overall prob of g2 imprinting on 1st exp
        g2p_0n[exp1] = ((1-e1_g2) * (no_exp2[exp1])) # won overall prob of g2 imprinting on 1st exp
      }
      else{
        over_max_age= over_max_age - 1; # decrease over_max_age
        g2p_1n[exp1] = 0
        g2p_0n[exp1] = 0
        p_nn = 0
      }
    } # End loop over years of 1st exposure
    
    
    
    #result <- list("naive" = naive, "H1p_11" = H1p_11, "H1p_10" = H1p_10, "H1p_01" = H1p_01, "H1p_00" = H1p_00, "H2p_11" = H2p_11, "H2p_10" = H2p_10, "H2p_01" = H2p_01, "H2p_00" = H2p_00, "H3p_11" = H3p_11, "H3p_10" = H3p_10, "H3p_01" = H3p_01, "H3p_00" = H3p_00)
    #result <- list("naive" = naive, "g1p_11" = g1p_11, "g1p_10" = g1p_10, "g1p_01" = g1p_01, "g1p_00" = g1p_00, "g2p_11" = g2p_11, "g2p_10" = g2p_10, "g2p_01" = g2p_01, "g2p_00" = g2p_00)
    result <- list("g1p_11" = g1p_11, "g1p_10" = g1p_10, "g1p_01" = g1p_01, "g1p_00" = g1p_00, "g1p_1n" = g1p_1n, "g1p_0n" = g1p_0n, "g2p_11" = g2p_11, "g2p_10" = g2p_10, "g2p_01" = g2p_01, "g2p_00" = g2p_00, "g2p_1n" = g2p_1n, "g2p_0n" = g2p_0n, "p_nn" = p_nn)
    result
}


#sum(p_11, na.rm = TRUE)+sum(p_10, na.rm = TRUE) +sum(p_01, na.rm = TRUE) +sum(p_00, na.rm = TRUE)
#sum(prob.matrix, na.rm = TRUE)=.8650201