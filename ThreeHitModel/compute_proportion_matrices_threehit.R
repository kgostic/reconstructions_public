source('~/Dropbox/R/Reconstructions/ThreeHitModel/gete_ij_threehit.R')

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
  three.hit = func_result$three.hit #3d array of years of exp1, exp2, exp3
  naivex1 = func_result$naivex1 #2d array of years of exp1, exp2
  naivex2 = func_result$naivex2 #vector of years of exp1
  naivex3 = func_result$naivex3 #scalar, prob no exposure

  nn = dim(three.hit)[1] ## n years exp1 possible
  mm = dim(three.hit)[2] ## n years exp2 possible
  ll = dim(three.hit)[3] ## n years exp3 possible
  naive.possible = TRUE  ## Default setting for logical switch

  if(i_yr - b_yr > 25){ # If the birth cohorts is over age 20 at the time of observation
    # Assume everyone has had three exposures
    # Normalize the twoexp matrix so that probabilities sum to 1, this implicitly sets naive probs to 0
    three.hit = three.hit/sum(three.hit, na.rm = TRUE)
    naivex3 = naivex2 = naivex1 = 0
    naive.possible = FALSE
  } # Else, assume some individuals are still missing an first, second or third exposure, and calculate the naive fraction below
  
  # List of subtypes that circulated in each year since 1918. 1 = H1N1, 2 = H2N2, 3 = H3N3, 4 = naive
  if(naive.possible == TRUE){
    # Include naive
    valid.exp.types = apply(rbind(circ_history[1:3, ], rep(1, ncol(circ_history))), MARGIN = 2, FUN = function(xx) which(xx>0)) 
  }else{
    # Don't include naive as a possibility
    valid.exp.types = apply(circ_history[1:3, ], MARGIN = 2, FUN = function(xx) which(xx>0))
  }
  
  
  ## Initialize master list of all possible combinations of exp1, exp2, exp3
  #1st row = first exp. subtype, 2nd row = 2nd exp subtype
  lindex = cbind(rbind(rep(1:3, each = 12), rep(rep(1:3, each = 4), 3), rep(1:4, times = 9)), c(1,4,4), c(2,4,4), c(3,4,4), c(4,4,4))
  
  outlist = lapply(1:ncol(lindex), FUN = function(xx){array(0, dim = c(nn, mm, ll))})
  ch.names = cbind(rbind(rep(c('H1', 'H2', 'H3'), each = 12), rep(rep(c('H1', 'H2', 'H3'), each = 4), 3), rep(c('H1', 'H2', 'H3', 'n'), times = 9)), c('H1', 'n', 'n'), c('H2', 'n', 'n'), c('H3', 'n', 'n'), c('n', 'n', 'n'))
  for(ii in 1:40){
    names(outlist)[ii] = paste(ch.names[,ii], collapse = '_')
  }

  
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
            if( (i1 == 4 & (i2 != 4 | i3 !=4)) | (i2 == 4 & i3 !=4)){
              # Combinations in which a valid exposure follows a naive exposure
              # are not biologically possible 
              # (e.g. can't have a 3rd exposure if you don't yet have a 2nd)
              # Skip these cases
            }else{
            valid = which(lindex[1,] == i1 & lindex[2,] == i2 & lindex[3,]  == i3) #Which list element to fill in?
            f1 = circ_history[i1, as.character(b_yr+exp1-1)] # Fraction of circulation in year of 1st exposure caused by subtype i1
            f2 = circ_history[i2, as.character(b_yr+exp2-1)] # Fraction of circulation in year of 2nd exposure caused by subtype i2
            f3 = circ_history[i3, as.character(b_yr+exp3-1)] # Fraction of circulation in year of 2nd exposure caused by subtype i2
            if(i1 == 4){ # Look in estimates of naive probs if i1 or i2 = 4, which corresponds to no exposure
              outlist[[valid]][exp1, exp2, exp3] = naivex3/(nn*(nn-1)*(nn-2)/6) # Divide the total prob of no exp evenly across the lower triangle of the matrix.
            }else if(i2 == 4){ # If one exposure, but no second or third,
              mm = nn-exp1
              outlist[[valid]][exp1, exp2, exp3] = naivex2[exp1]*f1/((mm-1)*mm/2)
            }else if(i3 == 4){
              mm = nn-exp2
              outlist[[valid]][exp1, exp2, exp3] = naivex1[exp1, exp2]*f1*f2/mm
              }else{
              outlist[[valid]][exp1, exp2, exp3] = three.hit[exp1, exp2, exp3]*f1*f2*f3
              }
            }
           }  # Close loop over possible exp3 subtypes
          } # Close loop over possible exp2 subtypes
        } # Close loop over possible exp1 subtypes
      } # Close loop over years of 3rd exposure
     } # Close loop over years of 2nd exposure
    } # Close loop over years of 1st exposure
  
  if(naive.possible == TRUE){
    # account for second exposures that happened in the last possible year
    for(exp1 in 1:(nn-1)){ # FOR EACH YEAR OF 1st EXPOSURE
      exp2 = nn
          pos.e1 = valid.exp.types[[as.character(b_yr+exp1-1)]]
          pos.e2 = valid.exp.types[[as.character(b_yr+exp2-1)]]
          # Loop over possible exposure types
          for(i1 in pos.e1){
            if(i1 == 4){ # skip this
            }else{
            for(i2 in pos.e2){
              if(i2 == 4){
                # skip
              }else{
              i3 == 4
              valid = which(lindex[1,] == i1 & lindex[2,] == i2 & lindex[3,]  == 4) #Which list element to fill in?
              f1 = circ_history[i1, as.character(b_yr+exp1-1)] # Fraction of circulation in year of 1st exposure caused by subtype i1
              f2 = circ_history[i2, as.character(b_yr+exp2-1)] # Fraction of circulation in year of 2nd exposure caused by subtype i2
              outlist[[valid]][exp1, exp2, nn] = naivex1[exp1, exp2]*f1*f2
            }}}}}
    #account for first exposures that happened in the last or penultimate possible year, assuming only one exposure happenend
    for(exp1 in nn-1:0){ # FOR EACH YEAR OF 1st EXPOSURE
      exp2 = nn
      exp3 = nn
      pos.e1 = valid.exp.types[[as.character(b_yr+exp1-1)]]
      # Loop over possible exposure types
      for(i1 in pos.e1){
        if(i1 == 4){ # skip this
        }else{
          i2 = i3 = 4
          valid = which(lindex[1,] == i1 & lindex[2,] == i2 & lindex[3,]  == 4) #Which list element to fill in?
          f1 = circ_history[i1, as.character(b_yr+exp1-1)] # Fraction of circulation in year of 1st exposure caused by subtype i1
          outlist[[valid]][exp1, nn, nn] = naivex2[exp1]*f1
        }}}
    }
    
 return(outlist)
}
  




# Test
setwd('~/Dropbox/R/Reconstructions/')
source('ThreeHitModel/country_data_import.R')
circ_history = import.country.dat('Egypt')
out = partition_to_subtype(b_yr = 1957, i_yr = 2000, circ_history = circ_history)
sapply(out, sum) ## For a birth year of 1957, only H2H2H2, H2H2H3, H2H2H1, H2H3H3, H2H3H1, H2H1H1, H1H1H1, H3H3H3 should be nonzero
sum(sapply(out, sum)) ## All categories should sum to 1


## Test when naive individuals should be included
out = partition_to_subtype(b_yr = 2000, i_yr = 2012, circ_history = circ_history)
sapply(out, sum) ## For a birth year of 1957, only H2H2H2, H2H2H3, H2H2H1, H2H3H3, H2H3H1, H2H1H1, H1H1H1, H3H3H3 should be nonzero
raw = get.e_ij(birth.year = 2000, incidence.year = 2012)
# No exposures, should be equal
raw$naivex3    
sum(out$n_n_n)

# One exposure, should be equal in sum
raw$naivex2 #p(missing two | exp1)
n2 = out$H1_n_n + out$H2_n_n + out$H3_n_n
rowSums(n2) #p(missing two | exp1)

# One exposure, should be equal in sum
raw$naivex1 #p(missing one | exp1, exp2)
n3 = out$H1_H1_n + out$H1_H2_n + out$H1_H3_n + out$H2_H1_n + out$H2_H2_n + out$H2_H3_n + out$H3_H1_n + out$H3_H2_n + out$H3_H3_n 
rowSums(n3, dims = 2) #p(missing one | exp1)
sum(round(rowSums(n3, dims = 2)  - raw$naivex1, digits = 10)) # Total difference should be 0

sum(sapply(out, sum)) ## All categories should sum to 1




