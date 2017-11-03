# #Load year-of-infection intensities
# intensities = read.csv('Intensitymatser.csv', col.names = c('Year', 'Intensity')); rownames(intensities) = 1911:2015
# #Weight the annual probability of infection by intensity
# load('pest.RData')
# 
# 
# weighted.attack.rate = p.est*(intensities$Intensity); names(weighted.attack.rate) = 1911:2015


get.type.weights.AB = function(years.out, Countries.out, type = 5){
  
  # This script imports a .csv that holds data on which subtypes (H1N1, H3N2 and influenza B) were observed in specific countries, and in specific years.
  # This script also sources two functions, get.cocirculation.ref and get.country.data, which we will use below to format a matrix that tells what fraction of circulation was driven by H1N1 vs. H3N2 in a given season
  source('Clean_CocirculationImport_2017.R')
  head(cocirculation) # Data frame that records which subtypes circulated in each country and year
  
  # Set years of birth (we will estimate probabilities of imprinting to H1N1, H2N2 and H3N2 for each of these birth years)
  birth.years = 1918:2017
  infection.years = birth.years
  
  

  
  #Initialize master weight matrix
  #Rows - which country and year are we doing the reconstruction from the persepctive of?
  #Cols - what birth year are we estimating imprinting probabilities for?
  weights.master.1 = matrix(NA, nrow = length(Countries.out)*length(years.out), ncol = length(birth.years), dimnames = list(paste(rep(years.out, length(Countries.out)), rep(Countries.out, each = length(years.out)), sep = ''), rev(birth.years))) # Initialize a matrix in which to store probs of H1N1 imprinting
  weights.master.2 = weights.master.naiive = weights.master.3 = weights.master.1
  # weights.master.1 stores probs of H1N1 imprinting
  #             ...2 stores probs of H2N2 imrpinting
  #             ...3 stores probs of H3N2 imprinting
  #             ...naiive stores probs of no influenza exposure (no imprintng). This will only take nonzero values for very young birth cohorts.
  
  # The code below fills in weights.master.1, etc.
  #   Repeat the loop for each country of interest
  for(cc in 1:length(Countries.out)){ 
    
    # Extract country-sepcific data on which subtypes circulated in past years
    if(Countries.out[cc] == 'China'){cocirculation.dat = get.cocirculation.ref('China')}
    if(Countries.out[cc] == 'Egypt'){cocirculation.dat = get.cocirculation.ref('Egypt')}
    if(Countries.out[cc] == 'Indonesia'){cocirculation.dat = get.cocirculation.ref('Indonesia')}
    if(Countries.out[cc] == 'Cambodia'){cocirculation.dat = get.cocirculation.ref('Cambodia')}
    if(Countries.out[cc] == 'Vietnam'){cocirculation.dat = get.cocirculation.ref('Vietnam')}
    if(Countries.out[cc] == 'Thailand'){cocirculation.dat = get.cocirculation.ref('Thailand')}
    if(Countries.out[cc] == 'USA'){cocirculation.dat = get.cocirculation.ref('USA')}
    if(Countries.out[cc] == 'UK'){cocirculation.dat = get.cocirculation.ref('UK')}
    if(Countries.out[cc] == 'Turkey'){cocirculation.dat = get.cocirculation.ref('Turkey')}
    if(Countries.out[cc] == 'Iraq'){cocirculation.dat = get.cocirculation.ref('Iraq')}
    if(Countries.out[cc] == 'Azerbaijan'){cocirculation.dat = get.cocirculation.ref('Azerbaijan')}
    if(Countries.out[cc] == 'Pakistan'){cocirculation.dat = get.cocirculation.ref('Pakistan')}
    if(Countries.out[cc] == 'Nigeria'){cocirculation.dat = get.cocirculation.ref('Nigeria')}

    
     
    #Extract and re-order relevant data
    #These describe the fraction of circulating influenza viruses isolated in a given year that were of subtype H1N1 (type1), H2N2 (type2), or H3N2 (type3)
    type1.dat = (cocirculation.dat['H1', as.character(birth.years)]) 
    type2.dat = (cocirculation.dat['H2', as.character(birth.years)]) 
    type3.dat = (cocirculation.dat['H3', as.character(birth.years)]) 

## THIS FUNCTION RETURNS A LIST OF:
# 1. The annual probabilities of infection for a given birth year, from the perspective of a given incidence year
# 2. The naiive fraction if the birth year is within 12 years of the incidence year

## Calculate e_ij, prob of any infection in a given year given birth in year i
#Output a vector of probabilities of being first infected in any year, j given birth in year i
source('gete_ij_2017.R')

### This function fills in a matrix of exposure weights (equations 2 and 3 in manuscript), for each year of reference
# This function is a wrapper for  function get.e_ij
# OUTPUT - A matrix of all probabilites, e_ij with birth years on rows, and years of first infection on columns
# INPUT - From the perspective of what year are we estimating imprinting patterns.
    
get.infection.weights = function(years.out){
  H1.mat = matrix(0, nrow = length(years.out), ncol = length(birth.years), dimnames = list((years.out), (birth.years)))
  naiive.mat = H2.mat = H3.mat = H1.mat
  for(jj in 1:length(years.out)){
    for(ii in 1:(years.out[jj]-1918+1)){ #for all relevant birth years
      n.inf.years = min(12, years.out[jj]-birth.years[ii])
      inf.probs = get.e_ij(birth.years[ii], years.out[jj])
      
      #If all 13 possible years of infection have passed, normalize
      if(length(inf.probs) == 13) inf.probs = inf.probs/sum(inf.probs)
      
      #Fill in the appropriate row (observation year) and column (birth year) of the output matrix
      H1.mat[jj, ii] = sum(inf.probs*type1.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
      H2.mat[jj, ii] = sum(inf.probs*type2.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
      H3.mat[jj, ii] = sum(inf.probs*type3.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
      naiive.mat[jj, ii] = round(1-sum(inf.probs), digits = 8) #Rounds to the nearest 8 to avoid machine 0 errors
    }
  }
  #return the output in order of 2015:2918
  H1.mat = H1.mat[,as.character(max(birth.years):min(birth.years))]
  H2.mat = H2.mat[,as.character(max(birth.years):min(birth.years))]
  H3.mat = H3.mat[,as.character(max(birth.years):min(birth.years))]
  naiive.mat = naiive.mat[,as.character(max(birth.years):min(birth.years))]
  
  
#   if(type == 1){return(H1.mat)}
#   if(type == 2){return(H2.mat)}
#   if(type == 3){return(H3.mat)}
#   if(type == 4){return(naiive.mat)}
  return(list(H1.mat = H1.mat, H2.mat = H2.mat, H3.mat = H3.mat, naiive.mat = naiive.mat))
}
    
  
    
    #Fill in the appropriate country-specific rows of the master matrix
wts = get.infection.weights(years.out = years.out) #Get all weight types for the year of interest

weights.master.1[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H1.mat #Fill in the appropriate country-year row of the weights.master.## matrix
weights.master.2[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H2.mat
weights.master.3[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H3.mat
weights.master.naiive[((cc-1)*length(years.out))+1:length(years.out), ] = wts$naiive.mat
  }
  
  if(type == 1){ return(weights.master.1)}
  if(type == 2){ return(weights.master.2)}
  if(type == 3){return(weights.master.3)}
  if(type == 4){ return(weights.master.naiive)}
if(type == 5){return(list(weights.master.1, weights.master.2, weights.master.3, weights.master.naiive))}
  
}


