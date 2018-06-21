source('~/Dropbox/R/Reconstructions/TwoHitModel/gete_ij_twohit.R')

## THIS FUNCTION INPUTS A BIRTH YEAR AND OBSERVATION YEAR OF INTEREST
## THIS FUNCTION OUTPUTS A LIST OF PROBABILITIES OF 1ST AND 2ND EXPOSURE
##  TO A PARTICULAR COMBINATION OF INFLUENZA SUBTYPES, GIVEN THE
##  INPUT BIRTH YEAR AND OBSERVATION YEAR.

## Inputs:
##    -> b_yr: Scalar. Birth year of interest
##    -> i_yr: Scalar. Year in which case was observed
##    -> circ_history: Data on fraction of influenza A circulation caused by H1N1, H2N2 or H3N2 in the specific years and countries of interest. (See function test at bottom of script for import of this data.)

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
  func_result= get.e_ij(birth.year = b_yr, incidence.year = i_yr) # Outputs the probabilities of a first exposure in year i and second exposure in year j, as well as the probabilites of no first exposure, and of a first exposure in year i, but no second exposure.
  twoexp = func_result$total_mat # matrix: row = year of 1st exposure, column = year of 2nd exposure
  oneexp = func_result$totalPof_yesexp1_and_noexp2 # vector: element = year of 1st exposure, 2nd exp has not yet happened
  noexp = func_result$totalPof_noexp1 # scalar. P(completely naive to flu | birth year)

  nn = length(b_yr:min(b_yr+17, i_yr)) ## total number of years of possible first exposure years, given birth year and observation year
  naive.possible = TRUE
  valid.exp.types = apply(rbind(circ_history[1:3, ], rep(1, ncol(circ_history))), MARGIN = 2, FUN = function(xx) which(xx>0)) # List of subtypes that circulated in each year since 1918. 1 = H1N1, 2 = H2N2, 3 = H3N3, 4 = naive

  if(i_yr - b_yr > 17){ # If the birth cohorts is over age 17 at the time of observation
    # Assume everyone has had both a first and second exposure
    # Normalize the twoexp matrix so that probabilities sum to 1, this implicitly sets naive probs to 0
    twoexp = twoexp/sum(twoexp, na.rm = TRUE)
    oneexp = noexp = 0
    #redo valid.exp.types so that "naive" isn't an option
    naive.possible = FALSE
  } # Else, assume some individuals are still missing an first or second exposure, and calculate the naive fraction below
  
  # List of subtypes that circulated in each year since 1918. 1 = H1N1, 2 = H2N2, 3 = H3N3, 4 = naive
  if(naive.possible == TRUE){
    # Include naive
    valid.exp.types = apply(rbind(circ_history[1:3, ], rep(1, ncol(circ_history))), MARGIN = 2, FUN = function(xx) which(xx>0)) 
  }else{
    # Don't include naive as a possibility
    valid.exp.types = apply(circ_history[1:3, ], MARGIN = 2, FUN = function(xx) which(xx>0))
  }
  
  
  ## Initialize master list
  outlist = lapply(1:13, FUN = function(xx){matrix(0, nrow=nn, ncol=nn)})
  names(outlist) = c('H1H1', 'H1H2', 'H1H3', 'H1n', 
                     'H2H1', 'H2H2', 'H2H3', 'H2n',
                     'H3H1', 'H3H2', 'H3H3', 'H3n', 'nn')
  lindex = rbind(c(rep(c(1,2,3), each = 4), 4), c(rep(c(1,2,3,4), 3), 4)) #1st row = first exp. subtype, 2nd row = 2nd exp subtype
  
  for(exp1 in 1:nn){ # FOR EACH YEAR OF 1st EXPOSURE
    mm = nn-exp1 # n possible exp2 years, assume exp1 and exp2 never happen in the same year.
    if(mm>0){
      for(exp2 in exp1+(1:mm)){ # AND FOR EACH YEAR OF 2ND EXPOSURE
        # Extract possible subtypes causing the first and 2nd exposure, given the year of 1st and 2nd exposure indicated by our position in the loop
        pos.e1 = valid.exp.types[[as.character(b_yr+exp1-1)]]
        pos.e2 = valid.exp.types[[as.character(b_yr+exp2-1)]]
        for(i1 in pos.e1){
          for(i2 in pos.e2){
            ## Partition by subtype
            if(i1 == 4 & i2 != 4){
            }else{
            valid = which(lindex[1,] == i1 & lindex[2,] == i2) #Which list element to fill in?
            f1 = circ_history[i1, as.character(b_yr+exp1-1)] # Fraction of circulation in year of 1st exposure caused by subtype i1
            f2 = circ_history[i2, as.character(b_yr+exp2-1)] # Fraction of circulation in year of 2nd exposure caused by subtype i2
            if(i1 == 4){ # Look in estimates of naive probs if i1 or i2 = 4, which corresponds to no exposure
              outlist[[valid]][exp1, exp2] = noexp/(nn*(nn-1)/2) # Divide the total prob of no exp evenly across the lower triangle of the matrix.
            }else if(i2 == 4){ # If one exposure, but no second
              outlist[[valid]][exp1, exp2] = oneexp[exp1]*f1/mm
            }else{
              outlist[[valid]][exp1, exp2] = twoexp[exp1,exp2]*f1*f2
            }
            }
          } # Close loop over possible exp2 subtypes
        } # Close loop over possible exp1 subtypes
      } # Close loop over years of 2nd exposure
    } # Close if
  } # Close loop over years of 1st exposure
  if(naive.possible == TRUE){
    exp1 = nn # account for the naive fraction from the last possible year of exposure
    pos.e1 = valid.exp.types[[as.character(b_yr+exp1-1)]]
    pos.e1 = pos.e1[-length(pos.e1)]
    for(i1 in pos.e1){
      valid = which(lindex[1,] == i1 & lindex[2,] == 4) #Which list element to fill in?
      f1 = circ_history[i1, as.character(b_yr+exp1-1)]
      outlist[[valid]][exp1, exp1] = f1*oneexp[exp1]
    }
  }
  outarray = simplify2array(outlist)
 return(outarray)
}
  




# Test
setwd('~/Dropbox/R/Reconstructions/')
source('TwoHitModel/country_data_import.R')
circ_history = import.country.dat('Egypt')
out = partition_to_subtype(b_yr = 1957, i_yr = 2000, circ_history = circ_history)
apply(out, MARGIN = 3, FUN = sum) ## For a birth year of 1957, only H2H2, H2H3 and H3H3 should be nonzero
sum(out) ## All categories should sum to 1

## Test when naive individuals should be included
out = partition_to_subtype(b_yr = 2000, i_yr = 2012, circ_history = circ_history)
apply(out, MARGIN = 3, FUN = sum) ## For a birth year of 2000, H1H1, H1H3, H3H1, H3H3, H1n, H3n and nn are all possible
sum(out) ## All categories should sum to 1



