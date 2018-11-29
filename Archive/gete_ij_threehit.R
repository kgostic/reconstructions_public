## Function inputs the year in which an individual was born (birth.year), and in which the individual became infected with bird flu (incidence year)

## Function outputs a matrix of probabilities, the first representing the probability of first flu infection in the first year of life (age 0), the second representing the probability of first flu infection in the second year of life (age 1), and so on up to the 18th year of life (age 17)
# Output a list, where the first entry gives the original total matrix, the second entry gives the total prob of having no first exposure and the third entry gives the total prob of a first exposure, but no second exposure.

## Function is called by another script which calculates normalized probabilities (see denominators in equations 3 and 4), and calculates probability of first exposure to a specific HA subtype.

get.e_ij = function(birth.year, incidence.year){
  ## Inputs
  #Load saved data on intensity of influenza circulation in specific years of first infection
  intensities = read.csv('Intensitymatser.csv', col.names = c('Year', 'Intensity')); rownames(intensities) = 1911:2017
 
  # Weighted attack rate = annual prob infection scaled by historical data on year-specific circulation intensity (from Intensitymaster.csv)
  weighted.attack.rate = 0.28*(intensities$Intensity); names(weighted.attack.rate) = 1911:2017
  
  ## Calculations (vector, constant, vector)
  jjs = birth.year:min(birth.year+25, incidence.year) #Possible years of first, second and third infection (ages 0-25)
  ajs = weighted.attack.rate[as.character(jjs)] #Get weighted attack rates "p" corresponding to possible years of first infection
  
  nn = length(jjs) # Number of possible years of 1st inf.
  
  
  ## Initilaize storage
  ## entry [i,j,k] gives the probability of having a first exposure in year i, a second exposure in year j, and a third exposure in year k
  one.hit = vector(mode = 'numeric', length = nn) # Stores the probabilities of first exposure [i]
  two.hit = array(data = 0, dim = c(nn, nn)) # Stores the probabilities of first x 2nd exposure [i,j]
  three.hit = array(data = 0, dim = c(nn, nn, nn)) # Stores the probabilities of first x 2nd x3rd exposure [i,j,k]
  
  ####--Calculate one-hit probs--####
  #### Overall formula is (1-a_0)*...(1-a_[i-1])*a_i, where a_0 through a_[i-1] represent the hazard of infection starting in the year of birth, up to the year before first infection
  for(exp1 in 1:nn){
    if(exp1 == 1){
      wait.to.first = 1 # If first infection happens in the year of birth, no (1-ax) factor. Multiply by 1.
    }else{
      wait.to.first = prod((1-ajs)[1:(exp1-1)]) # Else, multiply across probabilities of no infection for the years prior to first exposure
    }
    one.hit[exp1] = wait.to.first*ajs[exp1]
  }
  #sum(one.hit, na.rm = TRUE)
  
  
####--Calculate two-hit probs--####
  if(nn >= 2){ # nn gives the number of possible years in which childhood exposures can occur, and may be low for recently born cohorts
    ## Specifically, n == 1 for newborns (only one possible year of first infection)
    ## Because we assume two infections are not possible in the same year, a second infection is not possible if nn == 1
  for(exp1 in 1:(nn-1)){ # Loop over all years of first exposure
    for(exp2 in (exp1+1):(nn)){ # And all corresponding years of second exposure
      
      ## Calculate the probability of second exposure in year exp2
      if(exp2 == exp1+1){
        wait.to.second = 1 # If no wait to first exposure, just multiply by 1
      }else{
        wait.to.second = prod((1-ajs)[(exp1+1):(exp2-1)]) # Else, multiply across probabilities of no infection for the years prior to first exposure
      }
      ## The overall probability of first and second exposure in years [exp1, exp2] is given by the product: 
      ##                    ....p(exp1)..*........p(exp2).........
      two.hit[exp1, exp2] = one.hit[exp1]*wait.to.second*ajs[exp2]
    }
  }
  }
  #sum(two.hit, na.rm = TRUE)
  
  
  

####--Calculate three-hit probs--####
  if(nn >= 3){
  for(exp1 in 1:(nn-2)){ # Loop over exp1 years
    for(exp2 in (exp1+1):(nn-1)){ # And corresponding years of 2nd exposure
      for(exp3 in (exp2+1):nn){ # And corresponding years of 3rd exposure
        # Calculate probability of third exopsure in year exp3
        if(exp3 == exp2+1){
          wait.to.third = 1 # If no wait to first exposure, just multiply by 1
        }else{
          wait.to.third = prod((1-ajs)[(exp2+1):(exp3-1)]) # Else, multiply across probabilities of no infection for the years prior to first exposure
        }
        ## Overall prob of first, second, third exposure in years [exp1, exp2, exp3] is given by:
        ##                           ...p(exp1 n exp2)...*....p(exp3)............
        three.hit[exp1, exp2, exp3] = two.hit[exp1, exp2]*wait.to.third*ajs[exp3]
        
      } # Close loop over exp3
      } # Close loop over exp2 
    } # Close loop over exp1 
  }
#sum(three.hit, na.rm = TRUE)

## Calculate naive fractions
p_noexp1 = 1-sum(one.hit) # Scalar, probability of no first exposure
p_exp1_noexp2 = one.hit-rowSums(two.hit, na.rm = T) # Vector, probabilitiy of no second exposure, given first exposure in year exp1. Vector entries correspond to exp1.
p_exp1_exp2_noexp3 = two.hit - rowSums(three.hit, na.rm = T, dims = 2) # Matrix, probability of no third exposure, given exp1 and exp2. Matrix entries indicate years of first and second exposure

return(list('one.hit' = one.hit, 'two.hit' = two.hit, 'three.hit' = three.hit, 'naivex3' = p_noexp1, 'naivex2' = p_exp1_noexp2, 'naivex1' = p_exp1_exp2_noexp3))
} # End function

#
# # # ## check:
# get.e_ij(1990, 2017)
# test = get.e_ij(1990, 2017)
# sum(test$one.hit)
# sum(test$two.hit)
# sum(test$three.hit)
# sum(test$naivex1)
# sum(test$naivex2)
# sum(test$naivex3)
# sum(test$three.hit, test$naivex1, test$naivex2, test$naivex3)
# 
# 
# test = get.e_ij(2009, 2009)
# # Sum should equal 1
# sum(test$one.hit)
# sum(test$two.hit)
# sum(test$three.hit)
# sum(test$naivex1)
# sum(test$naivex2)
# sum(test$naivex3)
# sum(test$three.hit, test$naivex1, test$naivex2, test$naivex3)
