## Note - Intensitymaster.csv contains data loaded into this function

## Function inputs the year in which an individual was born (birth.year), and in which the individual became infected with bird flu (incidence year)

## Function outputs a vector of 13 probabilities, the first representing the probability of first flu infection in the first year of life (age 0), the second representing the probability of first flu infection in the second year of life (age 1), and so on up to the 13th year of life (age 12)

## Outputs correspond to the numerator of equation 3 in the supplementary methods

## Function is called by another script (Infection_age_structure.R), which calculates normalized probabilities (see denominators in equations 3 and 4), and calculates probability of first exposure to a specific HA subtype.

get.e_ij = function(birth.year, incidence.year){
  ## Inputs
  #Load saved data on intensity of influenza circulation in specific years of first infection
  intensities = read.csv('Intensitymatser.csv', col.names = c('Year', 'Intensity')); rownames(intensities) = 1911:2017
  load('pest.RData') # Load the annual probability of first infection, estimated from serological data (see two papers by Sauerbrei et al.)
  
  # Weighted attack rate = annual prob infection weighted by circulation intensity
  weighted.attack.rate = p.est*(intensities$Intensity); names(weighted.attack.rate) = 1911:2017
  
  ## Calculations
  jjs = birth.year:min(birth.year+12, incidence.year) #Possible years of first infection (ages 0-12)
  
  nn = length(jjs) # How many possible years of infection?
  ajs = weighted.attack.rate[as.character(jjs)] #Get weighted attack rates corresponding to possible years of first infection
  #names(e_ij) = jjs; names(naiive) = jjs
  
  ## Create matrices of 0s and 1s, which will be used below to vectorize and speed calculations
  ii = not_ii = matrix(0, nn, nn) #create two square matrices of dim nn
  #Fill in diagonal of one with 1s. Diagonal represents the year of first infection. 
  diag(ii) = 1
  #Fill in sub-diagonal for all the years since birth in which the individual escaped infection
  not_ii[lower.tri(not_ii)] = 1
  # Below, multiplying across rows of this matrix will calculate overall probability of first infection at different ages (see equation 3 in the supplement of the Science paper)

  #Create a matrix that takes (1-aj) for all non-infection years and aj for all infection years
  prod.mat = ajs*ii+matrix(rep(1-ajs, nn), nn, nn, byrow = T)*not_ii
  #Fill in upper triangle with 1s to make multiplication possible
  prod.mat[upper.tri(prod.mat)] = 1
  
  #Take product across rows
  e_ij = apply(prod.mat, 1, prod)
  e_ij
}
  