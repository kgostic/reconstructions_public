setwd("/Users/reginalee/Dropbox/Classical Disease Papers/Research scripts/Regina_research/")
#TIP: press TAB for setwd()
#Check that you set your working directory correctly
getwd()

## Note - Intensitymaster.csv contains data loaded into this function

## Function inputs the year in which an individual was born (birth.year), and in which the individual became infected with bird flu (incidence year)

## Function outputs a vector of 18 probabilities, the first representing the probability of first flu infection in the first year of life (age 0), the second representing the probability of first flu infection in the second year of life (age 1), and so on up to the 18th year of life (age 17)

## Outputs correspond to the numerator of equation 3 in the supplementary methods

## Functions called by another script (Infection_age_structure.R), which calculates normalized probabilities (see denominators in equations 3 and 4), and calculates probability of first exposure to a specific HA subtype.

first_get.e_ij = function(birth.year, incidence.year){
  ## Inputs
  #Load saved data on intensity of influenza circulation in specific years of first infection
  intensities = read.csv('Intensitymatser.csv', col.names = c('Year', 'Intensity')); rownames(intensities) = 1911:2017
  load('pest.RData') # Load the annual probability of first infection, estimated from serological data (see two papers by Sauerbrei et al.)
  #TIP: save results as lists in Rdata file, check when workspace is empty
  
  # Weighted attack rate = annual prob infection weighted by circulation intensity
  weighted.attack.rate = p.est*(intensities$Intensity); names(weighted.attack.rate) = 1911:2017
  
  ## Calculations (vector, constant, vector)
  jjs = birth.year:min(birth.year+17, incidence.year) #Possible years of first infection (ages 0-17)
  nn = length(jjs) # How many possible years of infection?
  ajs = weighted.attack.rate[as.character(jjs)] #Get weighted attack rates corresponding to possible years of first infection
  #names(e_ij) = jjs; names(naiive) = jjs
  
  total = numeric(nn); #want to be a vector
  #DEALING WITH SINGLE EXPOSURE ONLY
  for(occ in 1:nn){ #year is, increment "0-1" "0-2" "0-3 etc", length nn
    if(occ > 1){ #for occurences past year zero, min 1:1
      years.not.infected = prod(1-ajs[1:(occ-1)]) #probability p of NOT getting infected: "(1-p)^(x-1)"
    }
    else{
      years.not.infected = 1 #1-ajs[1] (multiplication purposes)
    }
    total[occ] = years.not.infected * ajs[occ] #combine with infected year
    }
 total
}
