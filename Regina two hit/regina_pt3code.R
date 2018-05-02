#setwd("/Users/reginalee/Dropbox/Classical Disease Papers/Research scripts/Regina_research/")
#TIP: press TAB for setwd()
#Check that you set your working directory correctly
getwd()

## Note - Intensitymaster.csv contains data loaded into this function

## Function inputs the year in which an individual was born (birth.year), and in which the individual became infected with bird flu (incidence year)

## Function outputs a vector of 13 probabilities, the first representing the probability of first flu infection in the first year of life (age 0), the second representing the probability of first flu infection in the second year of life (age 1), and so on up to the 18th year of life (age 17)

## Outputs correspond to the numerator of equation 3 in the supplementary methods

## Function is called by another script (Infection_age_structure.R), which calculates normalized probabilities (see denominators in equations 3 and 4), and calculates probability of first exposure to a specific HA subtype.

get.e_ij = function(birth.year, incidence.year){
  ## Inputs
  #Load saved data on intensity of influenza circulation in specific years of first infection
  intensities = read.csv('Intensitymatser.csv', col.names = c('Year', 'Intensity')); rownames(intensities) = 1911:2017
  load('pest.RData') # Load the annual probability of first infection, estimated from serological data (see two papers by Sauerbrei et al.)
  #TIP: save results as lists in Rdata file, check when workspace is empty
  
  # Weighted attack rate = annual prob infection weighted by circulation intensity
  weighted.attack.rate = p.est*(intensities$Intensity); names(weighted.attack.rate) = 1911:2017
  
  ## Calculations (vector, constant, vector)
  jjs = birth.year:min(birth.year+17, incidence.year) #Possible years of first infection (ages 0-17)
  nn = length(jjs) # Number of possible years of 1st inf.
  ajs = weighted.attack.rate[as.character(jjs)] #Get weighted attack rates "p" corresponding to possible years of first infection
  #names(e_ij) = jjs; names(naiive) = jjs
  
  
  #### BEGIN CALCULATION OF PROBABILITIES OF EXPOSURE 1 AND EXPOSURE 2 AT A TIME 0, 1, 2, ... 17 YEARS AFTER BIRTH.
  ## Let i represent row index and j represent column index.
  ## entry total_ij gives the total probability of having a first exposure in year i and a second exposure in year j
  total <- matrix(NA, nrow=nn, ncol=nn); #want to be a MATRIX
  ## KG addition:
  # all_exp1_probs is a vector that stores the total probability of 1st exposure in year exp1, regardless of exp2 status.
  all_exp1_probs = vector('numeric', nn)
  
  # OUTER LOOP--CALCULATE THE PROBABILITY OF FIRST EXPOSURE IN EACH POSSIBLE YEAR FROM 0-17 YEARS AFTER BIRTH
  for(exp1 in 1:nn){ #year is, increment "0-1" "0-2" "0-3 etc", CONSTANT length nn
    mm = nn - exp1 # Number of possible years of 2nd inf.
    if(exp1 > 1){ #IF exp1 NOT in first poss year
      years.not.infected1 = prod(1-ajs[1:(exp1-1)]) #probability p of NOT getting infected from year 0 - exp1: "(1-p)^(x-1)"
    }else{ #IF exp1 OCCURRED in very first poss year
      # Then we don't need any factors of (1-p). Store a 1, so we can multiply by 1*p below.
      years.not.infected1 = 1 #1-ajs[1] (multiplication purposes)
    }
    # 2. Calculate and store the total probabilitiy of first exposure in year exp1
    all_exp1_probs[exp1] = years.not.infected1*ajs[exp1]
    
    
    # INNER LOOP--CALCULATE THE PROBABILITY OF SECOND EXPOSURE IN EACH POSSIBLE YEAR, GIVEN FIRST EXPOSURE IN YEAR EXP1
    start_2 = exp1 + 1 # Set the first possible year of 2nd exposure as the year after exp1
    if(mm > 0){ # If there is more than one possible year of 2nd exposure, loop through every possible year of second exposure.
      for(exp2 in start_2:(start_2+mm-1)){ #exp2 interval begins AFTER exp1
        if(exp2 > start_2){ #If exp2 NOT in first poss year, calculate the probability of no 2nd exposure from year start_2:(exp2-1)
          years.not.infected2 = prod(1-ajs[start_2:(exp2-1)])        
          }else{ #exp2 in very first poss year
          # Then we don't need any factors of (1-p). Store 1, so we can multiply by 1*p below
          years.not.infected2 = 1 #1-ajs[1] (multiplication purposes)
        }
        
        ## For each combo of exp1 and exp2, fill in the probability of 1st exposure in year 1 AND 2nd exposure in year 2
        total[exp1,exp2] <- (years.not.infected1*ajs[exp1])*(years.not.infected2*ajs[exp2])
        
      }## Exit if statement for >0 possible years of 2nd exposure
    }## Exit inner loop: 2nd exposure calculation
  }## Exit outter loop: 1st exposure calculation
  ## Return total
  #return(list('total_mat' = total, 'totalPof_noexp1' = 1-sum(all_exp1_probs), 'totalPof_yesexp1_and_noexp2' = sum(all_exp1_probs)-sum(total, na.rm = TRUE)))
  return(list('total_mat' = total, 'totalPof_noexp1' = 1-sum(all_exp1_probs), 'totalPof_yesexp1_and_noexp2' = all_exp1_probs-rowSums(total, na.rm = TRUE)))
} # End function

## KG - 
# Now the function outputs a list, where the first entry gives the original total matrix, the second entry gives the total prob of having no first exposure and the third entry gives the total prob of a first exposure, but no second exposure.


#0.23235532 0.10578127 0.33370429 0.05259924 0.08451051 0.05644345 from unnested version
