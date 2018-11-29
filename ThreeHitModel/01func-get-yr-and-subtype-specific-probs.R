## Load function to calculate the probabilities of first, second and third exposures in years
## [i,j,k], given birth year and observation year:
source('ThreeHitModel/00func-get-yr-specific-probs.R')



### Write a new funtion: -----------------------------
## THIS FUNCTION INPUTS A BIRTH YEAR AND OBSERVATION YEAR OF INTEREST
## THIS FUNCTION OUTPUTS A LIST OF PROBABILITIES OF 1ST, 2ND and 3RD EXPOSURE
##  TO A PARTICULAR COMBINATION OF INFLUENZA SUBTYPES, GIVEN THE
##  INPUT BIRTH YEAR AND OBSERVATION YEAR.

## Inputs:
##    -> b_yr: Scalar. Birth year of interest
##    -> i_yr: Scalar. Year in which case was observed
##    -> circ_history: Data on fraction of influenza A circulation caused by H1N1, H2N2 or H3N2 in the specific years and countries of interest. These are imported using function import.country.dat(), defined in the script ../0func-country-data-import.R (See bottom of script for example import of these inputs)

## Outputs:
## One list with the following elements
## pH1H1 --> Probability of 1st and second exposure to H1
## pH1H2 --> Probability of 1st exposure to H1, and second exposure to H2
## pH1H3 --> Probability of 1st exposure to H1, and second exposure to H3
## pH1n --> Probability of 1st exposure to H1, and no second exposure (may be nonzero for ages <= 17)
## pH2H1 --> ...etc:
## pH2H2
## pH2H3
## pH2n
## pH3H1
## pH3H2
## pH3H3
## pH3n
## pnn
partition_to_subtype = function(b_yr, i_yr, circ_history){
  func_result= get.e_ij(birth.year = b_yr, incidence.year = i_yr) # Outputs the probabilities of a first exposure in year i, second exposure in year j, and third exposure in year k, as well as the probabilites of no first exposure, and of a first exposure in year i, but no second exposure.
  
  #3D array of probs of first, second and third exposure in years [i,j,k]
  three.hit = func_result$three.hit 
  naivex1 = func_result$naivex1 #2D array of probs that the third exposure has not yet occurred (one missing), given first and second exposure in years [i,j]. 
  naivex2 = func_result$naivex2 #vector of probs that the second AND third exposures are missing, given first exposure in year [i]
  naivex3 = func_result$naivex3 #scalar, prob all three exposures have not yet occurred

  nn = dim(three.hit)[1] ## n years exp1 possible
  mm = dim(three.hit)[2] ## n years exp2 possible
  ll = dim(three.hit)[3] ## n years exp3 possible
  naive.possible = TRUE  ## Default setting for logical switch

  if(i_yr - b_yr > 25){ # If the birth cohorts is over age 25 at the time of observation
    # Assume everyone has had at least three exposures
    # Normalize the three.hit matrix so that probabilities sum to 1, this implicitly sets naive probs to 0
    # Normalization is necessary because calculated sum of probs will approach 1, but is not 1 exactly unless we normalize
    three.hit = three.hit/sum(three.hit, na.rm = TRUE)
    # Assume it is not possible to be naive at the first, second or third position, so set these to 0
    naivex3 = naivex2 = naivex1 = 0
    naive.possible = FALSE
  } # Else, assume some individuals are young enoung to still missing an first, second or third exposure, and calculate the naive fraction below
  
  # List of subtypes that circulated in each year since 1918. 1 = H1N1, 2 = H2N2, 3 = H3N3, 4 = naive
  # if(naive.possible == TRUE){
  #   # Include naive
  #   valid.exp.types = apply(rbind(circ_history[1:3, ], rep(1, ncol(circ_history))), MARGIN = 2, FUN = function(xx) which(xx>0)) 
  # }else{
  #   # Don't include naive as a possibility
    valid.exp.types = apply(circ_history[1:3, ], MARGIN = 2, FUN = function(xx) which(xx>0))
  # }
  
  
  ## Creat a matrix representing all possible combinations of subtypes encountered at exp1 (row1), exp2 (row2), exp3 (row3)
  #1st row = first exp. subtype, 2nd row = 2nd exp subtype
  # 1 = H1N1, 2 = H2N2, 3 = H3N2, 4 = naive
  lindex = cbind(rbind(rep(1:3, each = 12), rep(rep(1:3, each = 4), 3), rep(1:4, times = 9)), c(1,4,4), c(2,4,4), c(3,4,4), c(4,4,4))
  
  ## Generate a list of empty matrices. One matrix entry per subtype-specific exposure permuation. This will be used to store output probabilities
  outlist = lapply(1:ncol(lindex), FUN = function(xx){array(0, dim = c(nn, mm, ll))})
  ## Add list names
  ch.names = cbind(rbind(rep(c('H1', 'H2', 'H3'), each = 12), rep(rep(c('H1', 'H2', 'H3'), each = 4), 3), rep(c('H1', 'H2', 'H3', 'n'), times = 9)), c('H1', 'n', 'n'), c('H2', 'n', 'n'), c('H3', 'n', 'n'), c('n', 'n', 'n'))
  for(ii in 1:length(outlist)){
    names(outlist)[ii] = paste(ch.names[,ii], collapse = '_')
  }

  
  ### First, track everyone who has had three exposures.
  ### This is only possible if nn >= 3 (3 possible years of exposure, i.e. birth cohort is at least 2 years old)
  if(nn >= 3){
  for(exp1 in 1:(nn-2)){ # FOR EACH YEAR OF 1st EXPOSURE
      for(exp2 in (exp1+1):(nn-1)){ # AND FOR EACH YEAR OF 2ND EXPOSURE
        for(exp3 in (exp2+1):nn){ # AND FOR EACH YEAR OF 3rd EXPOSURE
        # Extract possible subtypes causing the first and 2nd exposure, given the year of 1st and 2nd exposure indicated by our position in the loop
        pos.e1 = valid.exp.types[[as.character(b_yr+exp1-1)]]
        pos.e2 = valid.exp.types[[as.character(b_yr+exp2-1)]]
        pos.e3 = valid.exp.types[[as.character(b_yr+exp3-1)]]
        # Loop over possible exposure types
        for(i1 in pos.e1){
          for(i2 in pos.e2){
            for(i3 in pos.e3){
            ## Partition by subtype
            valid = which(lindex[1,] == i1 & lindex[2,] == i2 & lindex[3,]  == i3) #Which list element to fill in?
            f1 = circ_history[i1, as.character(b_yr+exp1-1)] # Fraction of circulation in year of 1st exposure caused by subtype i1
            f2 = circ_history[i2, as.character(b_yr+exp2-1)] # Fraction of circulation in year of 2nd exposure caused by subtype i2
            f3 = circ_history[i3, as.character(b_yr+exp3-1)] # Fraction of circulation in year of 2nd exposure caused by subtype i3
            ## prob. 1st, 2nd, 3rd exposure in years exp1, exp2, exp3, and to subtypes f1, f2, f3 is given by p(exp1, exp2, exp3)*f1*f2*f3
              outlist[[valid]][exp1, exp2, exp3] = three.hit[exp1, exp2, exp3]*f1*f2*f3
           }  # Close loop over possible exp3 subtypes
          } # Close loop over possible exp2 subtypes
        } # Close loop over possible exp1 subtypes
      } # Close loop over years of 3rd exposure
     } # Close loop over years of 2nd exposure
    } # Close loop over years of 1st exposure
  } # Close loop over if
    
    
  ## At this point, we have our outputs for everyone who is old enough to have had three exposures.
  ## If the cohort was under age 25 in the year of case observation, we need to calculate the fraction
  ## that is still naive at the first, second or third position.
  
  ## First, let's deal with inidividuals who are only naive at the third position.
  ##  We still care about what years their first and second exposures occurred in, and to what
  ##  subtypes they were exposed
  if(naive.possible == TRUE){
    # Keep track of everyone that has exp1 and exp2 but lacks exp3
    # this is only possible if nn >= 2
    if(nn >= 2){
      for(exp1 in 1:(nn-1)){ # FOR EACH YEAR OF 1st EXPOSURE
        for(exp2 in (exp1+1):(nn)){ # AND FOR EACH YEAR OF 2ND EXPOSURE
          exp3 = nn # dummy setting
            # Extract possible subtypes causing the first and 2nd exposure, given the year of 1st and 2nd exposure indicated by our position in the loop
            pos.e1 = valid.exp.types[[as.character(b_yr+exp1-1)]]
            pos.e2 = valid.exp.types[[as.character(b_yr+exp2-1)]]
            # Loop over possible exposure types
            for(i1 in pos.e1){
              for(i2 in pos.e2){
                  ## Partition by subtype
                  valid = which(lindex[1,] == i1 & lindex[2,] == i2 & lindex[3,]  == 4) #Which list element to fill in?
                  f1 = circ_history[i1, as.character(b_yr+exp1-1)] # Fraction of circulation in year of 1st exposure caused by subtype i1
                  f2 = circ_history[i2, as.character(b_yr+exp2-1)] # Fraction of circulation in year of 2nd exposure caused by subtype i2
                  outlist[[valid]][exp1, exp2, exp3] = naivex1[exp1, exp2]*f1*f2
              } # Close loop over possible exp2 subtypes
            } # Close loop over possible exp1 subtypes
        } # Close loop over years of 2nd exposure
      } # Close loop over years of 1st exposure
    } # Close if
    
    
    # Now, keep track of everyone that has exp1 but lacks exp2 and exp3
    # We still care about the year in which they had exp1, and the subtype responsible for exp1
      for(exp1 in 1:(nn)){ # FOR EACH YEAR OF 1st EXPOSURE
          exp2 = nn # dummy setting. Exp2 and 3 haven't actually happened, but we need to store them somewhere
          exp3 = nn # dummy setting
          # Extract possible subtypes causing the first and 2nd exposure, given the year of 1st and 2nd exposure indicated by our position in the loop
          pos.e1 = valid.exp.types[[as.character(b_yr+exp1-1)]]
          # Loop over possible exposure types
          for(i1 in pos.e1){
              ## Partition by subtype
              valid = which(lindex[1,] == i1 & lindex[2,] == 4 & lindex[3,]  == 4) #Which list element to fill in?
              f1 = circ_history[i1, as.character(b_yr+exp1-1)] # Fraction of circulation in year of 1st exposure caused by subtype i1
              outlist[[valid]][exp1, exp2, exp3] = naivex2[exp1]*f1
          } # Close loop over possible exp1 subtypes
      } # Close loop over years of 1st exposure

    # Keep track of everyone that is totally naive
    valid = which(lindex[1,] == 4 & lindex[2,] == 4 & lindex[3,]  == 4) #Which list element to fill in?
    outlist[[valid]][nn, nn, nn] = naivex3
    
  }#Close if naive.possible == TRUE
    
 return(outlist)
}
  



# # # Test
# setwd('~/Dropbox/R/Reconstructions/')
# source('0func-country_data_import.R')
# circ_history = import.country.dat('USA') ## See within 0func-fountry_data_import for a list of available countries
# ## ../readme.txt describes workflow for adding new countries
# out = partition_to_subtype(b_yr = 1957, i_yr = 2000, circ_history = circ_history)
# sapply(out, sum) ## For a birth year of 1957, only H2H2H2, H2H2H3, H2H2H1, H2H3H3, H2H3H1, H2H1H1, H1H1H1, H3H3H3 should be nonzero
# sum(sapply(out, sum)) ## All categories should sum to 1
# 
# out = partition_to_subtype(b_yr = 1939, i_yr = 2000, circ_history = circ_history)
# sapply(out, sum) ## For a birth year of 1957, only H2H2H2, H2H2H3, H2H2H1, H2H3H3, H2H3H1, H2H1H1, H1H1H1, H3H3H3 should be nonzero
# sum(sapply(out, sum)) ## All categories should sum to 1
# 
# 
# 
# 
# 
# ## Test when naive individuals should be included
# out = partition_to_subtype(b_yr = 2011, i_yr = 2012, circ_history = circ_history)
# sapply(out, sum) ## For a birth year of 1957, only H2H2H2, H2H2H3, H2H2H1, H2H3H3, H2H3H1, H2H1H1, H1H1H1, H3H3H3 should be nonzero
# raw = get.e_ij(birth.year = 2011, incidence.year = 2012)
# # No exposures, should be equal
# raw$naivex3
# sum(out$n_n_n)
# 
# # One exposure, should be equal in sum
# raw$naivex2 #p(missing two | exp1)
# n2 = out$H1_n_n + out$H2_n_n + out$H3_n_n
# rowSums(n2) #p(missing two | exp1)
# 
# # One exposure, should be equal in sum
# raw$naivex1 #p(missing one | exp1, exp2)
# n3 = out$H1_H1_n + out$H1_H2_n + out$H1_H3_n + out$H2_H1_n + out$H2_H2_n + out$H2_H3_n + out$H3_H1_n + out$H3_H2_n + out$H3_H3_n
# rowSums(n3, dims = 2) #p(missing one | exp1)
# sum(round(rowSums(n3, dims = 2)  - raw$naivex1, digits = 10)) # Total difference should be 0
# 
# sum(sapply(out, sum)) ## All categories should sum to 1




