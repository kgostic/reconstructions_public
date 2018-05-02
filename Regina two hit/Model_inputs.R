##################################################
### --
#####        Source this script to import inputs
#####        for H5N1 and H7N9 model comparison
### --
##################################################
#setwd("/Users/reginalee/Dropbox/Classical Disease Papers/Research scripts/Regina_research/")
max.yr = 2017
print(paste('max year =', max.yr))

## Check that "exposure' input has been defined
if(is.logical(exposure.switch) == FALSE) stop('Specify a logical variable, exposure')
##  KG note - I set up this script so that before your source this script you should define a boolean variable called exposure.switch
##            exposure.switch = FALSE ignores data on the rate of poultry exposure in different age groups
##            exposure.switch = TRUE uses age-specific poultry exposure data to modulate the demography vector
##  I did this because it makes it easier to code the models with/without "E" in the ModelComparison scripts.
##  If you forget to define "exposure.switch" you'll get an error if you try to source the inputs.


##################################################
#####        Import Exposure
##################################################
#define years for which to import data
years = use.years =  c(1997, 2003:2017)

source('Exposure_Inputs_H7N9.R')
exposure.by.year = get.exposure.data()
# Returns a 98-vector.
# Entires give the fraction of total poultry exposures occurring in people of age 0, 1, 2, 3, ... 98

# Use the above vector to fill in a matrix whose columns each represent a birth year and whose rows represent each 
# unique country-year in which H5N1 or H7N9 cases were observed.
exposure.master = matrix(nrow = length(years), ncol = length(max.yr:1918), dimnames = list(years, max.yr:1918))
# Declare an empty matrix with named rows and columns
for(ii in 1:length(years)){
  skips = max.yr - years[ii]
  exposure.master[ii, ] = c(rep(0, skips), exposure.by.year[1:(length(max.yr:1918)-skips)])
}

# If you view exposure master, you'll see that we've filled in a 0 in all birth years after the year of case observation (e.g. if we collected data in 1997, then birth years 1998, 1999, 2000... had not yet occurred. So we will in a 0). We start inputting the exposure.by.year vector in the birth year that matches the year of data collection, where people are age 0




##################################################
#####        Import Demography
##################################################
## KG note - I'm sorry this is such a mess. Don't worry too much about this code. It's designed to take demographic data and to put it into the master matrix row that corresponds to the country and year of case observation.
source('census_sort.R')

## China - import
china.raw = read.csv('census_data_China_1997_2018.csv', skip = 1, header = T)
# ----- Sort by year
china.demography = matrix(NA, nrow = length(years), ncol = length(1918:max.yr), dimnames = list(years, rev(1918:max.yr)))
for(ii in 1:length(years)){
  year.pop = single.year.sort(china.raw[china.raw$Year == years[ii],])
  #For years < max.yr, fill max.yr:actual year with 0s and remove the same number of entries from end
  skips = max.yr-years[ii]
  china.demography[ii, ] = c(rep(0, skips), year.pop[1:(length(1918:max.yr)-skips)])
}

## Egypt - import
egypt.raw = read.csv('census_data_Egypt_1997_2018.csv', skip = 1, header = T)
# ----- Sort by year
egypt.demography = matrix(NA, nrow = length(years), ncol = length(1918:max.yr), dimnames = list(years, rev(1918:max.yr)))
for(ii in 1:length(years)){
  year.pop = single.year.sort(egypt.raw[egypt.raw$Year == years[ii],])
  #For years < max.yr, fill max.yr:actual year with 0s and remove the same number of entries from end
  skips = max.yr-years[ii]
  egypt.demography[ii, ] = c(rep(0, skips), year.pop[1:(length(1918:max.yr)-skips)])
}

## Indonesia - import
indonesia.raw = read.csv('census_data_Indonesia_1997_2018.csv', skip = 1, header = T)
# ----- Sort by year
indonesia.demography = matrix(NA, nrow = length(years), ncol = length(1918:max.yr), dimnames = list(years, rev(1918:max.yr)))
for(ii in 1:length(years)){
  year.pop = single.year.sort(indonesia.raw[indonesia.raw$Year == years[ii],])
  #For years < max.yr, fill max.yr:actual year with 0s and remove the same number of entries from end
  skips = max.yr-years[ii]
  indonesia.demography[ii, ] = c(rep(0, skips), year.pop[1:(length(1918:max.yr)-skips)])
}

## vietnam - import
vietnam.raw = read.csv('census_data_Vietnam_1997_2018.csv', skip = 1, header = T)
# ----- Sort by year
vietnam.demography = matrix(NA, nrow = length(years), ncol = length(1918:max.yr), dimnames = list(years, rev(1918:max.yr)))
for(ii in 1:length(years)){
  year.pop = single.year.sort(vietnam.raw[vietnam.raw$Year == years[ii],])
  #For years < max.yr, fill max.yr:actual year with 0s and remove the same number of entries from end
  skips = max.yr-years[ii]
  vietnam.demography[ii, ] = c(rep(0, skips), year.pop[1:(length(1918:max.yr)-skips)])
}

## cambodia - import
cambodia.raw = read.csv('census_data_Cambodia_1997_2018.csv', skip = 1, header = T)
# ----- Sort by year
cambodia.demography = matrix(NA, nrow = length(years), ncol = length(1918:max.yr), dimnames = list(years, rev(1918:max.yr)))
for(ii in 1:length(years)){
  year.pop = single.year.sort(cambodia.raw[cambodia.raw$Year == years[ii],])
  #For years < max.yr, fill max.yr:actual year with 0s and remove the same number of entries from end
  skips = max.yr-years[ii]
  cambodia.demography[ii, ] = c(rep(0, skips), year.pop[1:(length(1918:max.yr)-skips)])
}

## Thailand - import
thailand.raw = read.csv('census_data_Thailand_1997_2018.csv', skip = 1, header = T)
# ----- Sort by year
thailand.demography = matrix(NA, nrow = length(years), ncol = length(1918:max.yr), dimnames = list(years, rev(1918:max.yr)))
for(ii in 1:length(years)){
  year.pop = single.year.sort(thailand.raw[thailand.raw$Year == years[ii],])
  #For years < max.yr, fill max.yr:actual year with 0s and remove the same number of entries from end
  skips = max.yr-years[ii]
  thailand.demography[ii, ] = c(rep(0, skips), year.pop[1:(length(1918:max.yr)-skips)])
}


##################################################
#####        Import Data
##################################################
# ## 1997
# H5.1997 = read.csv('H5N1_Linelist_1997cases.csv', header = T, stringsAsFactors = FALSE)
# 
# ## 2006-2010
# H5.RKI = read.csv('H5N1_Final_Linelist_RKI_2006_2010cases.csv', header = T, stringsAsFactors = FALSE)
# 
# ## all others
# H5.new = read.csv('H5N1_linelist_KG_MRA.csv', header = T, stringsAsFactors = FALSE)
# 
# full.data = (data.frame(country = c(H5.1997$country, H5.new$country, H5.RKI$country), age = c(H5.1997$age, H5.new$age, H5.RKI$age), year = c(H5.1997$year, H5.new$year, H5.RKI$year), outcome = c(H5.1997$outcome, H5.new$outcome, H5.RKI$outcome), confirmed = c(H5.1997$confirmed, H5.new$confirmed, H5.RKI$confirmed), province = c(H5.1997$province, H5.new$province, H5.RKI$province)))



#Define birth years of interest
birth.yrs = rev(1918:max.yr)

#Format imported data into useful data frame
source('Get_data_by_country.R')
Egypt.cases.by = get.country.data.H5N1('Egypt')
China.cases.by = get.country.data.H5N1('China')
Vietnam.cases.by = get.country.data.H5N1('Vietnam')
Cambodia.cases.by = get.country.data.H5N1('Cambodia')
Indonesia.cases.by = get.country.data.H5N1('Indonesia')
Thailand.cases.by = get.country.data.H5N1('Thailand')


source('Get_country_deaths.R')
#Mortality data
Egypt.deaths.by = get.country.deaths.H5N1('Egypt')
China.deaths.by = get.country.deaths.H5N1('China')
Vietnam.deaths.by = get.country.deaths.H5N1('Vietnam')
Cambodia.deaths.by = get.country.deaths.H5N1('Cambodia')
Indonesia.deaths.by = get.country.deaths.H5N1('Indonesia')
Thailand.deaths.by = get.country.deaths.H5N1('Thailand')




#######################################
# -----
#####     GENERATE YOUNG/OLD WEIGHTS BY YEAR
# -----
######################################
Uy= c(rep(1, 5), rep(0, 100)) #age 0-4 high risk
## MID
Um= c(rep(0,5), rep(1, 60), rep(0, 40)) #age 5-64 low risk
## Old
Uo = c(rep(0, 65), rep(1, 40))

nr = length(1918:max.yr)
if(nr > 105){error("Need longer Uy, Um and Uo vectors!")}

Uo.by.year = Um.by.year = Uy.by.year = matrix(NA, nrow = length(use.years), ncol = nr, dimnames = list(use.years, max.yr:1918))
for(ii in 1:length(use.years)){
  skips = max.yr-use.years[ii]
  Uy.by.year[ii, ] = c(rep(0, skips), Uy[1:(nr-skips)])
  Um.by.year[ii, ] = c(rep(0, skips), Um[1:(nr-skips)])
  Uo.by.year[ii, ] = c(rep(0, skips), Uo[1:(nr-skips)])
}

Uy.master = rbind(Uy.by.year, Uy.by.year, Uy.by.year, Uy.by.year, Uy.by.year, Uy.by.year); rownames(Uy.master) = paste(rep(use.years, 6), rep(c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand'), each = length(use.years)), sep = '') #for 5 countries
Uo.master = rbind(Uo.by.year, Uo.by.year, Uo.by.year, Uo.by.year, Uo.by.year, Uo.by.year); rownames(Uo.master) = paste(rep(use.years, 6), rep(c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand'), each = length(use.years)), sep = '') #for 5 countries
Um.master = rbind(Um.by.year, Um.by.year, Um.by.year, Um.by.year, Um.by.year, Um.by.year); rownames(Um.master)  = paste(rep(use.years, 6), rep(c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand'), each = length(use.years)), sep = '') #for 5 countries

###################################
# -----
#####     Likelihood and AIC functions
# -----
##################################
#source('Model_Comp_LikelihoodFunctions_H2H1_H5N1_revised.R')

AIC = function(neg.log.lk, n.pars){
  2*n.pars + 2*neg.log.lk
}
source('LR.test.R')


###################################
# -----
#####     ANALYSIS
# -----
##################################

#calculate the likelihood given demography, exposure and immune effects
# estimate g1, the protective effect of g1 exposure to a g1 virus
# OR fix g1 = 1 to return the null model: demography and exposure only

## immune.group.lk = function(pars, g1 = 1, g1.weight, g2.weight, p0, xx)
#  estimate g2 only

#Set up a data matrix
data.master = rbind(China.cases.by, Cambodia.cases.by, Egypt.cases.by, Indonesia.cases.by, Vietnam.cases.by, Thailand.cases.by)
rownames(data.master) = paste(rep(use.years, 6), rep(c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand'), each = length(use.years)), sep = '')

deaths.master = rbind(China.deaths.by, Cambodia.deaths.by, Egypt.deaths.by, Indonesia.deaths.by, Vietnam.deaths.by, Thailand.deaths.by)
rownames(deaths.master) = paste(rep(use.years, 6), rep(c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand'), each = length(use.years)), sep = '')

#Set up demographic matrix, each row corresponds to a given country-year
demography.master = rbind(china.demography, cambodia.demography, egypt.demography, indonesia.demography, vietnam.demography, thailand.demography)
rownames(demography.master) = paste(rep(use.years, 6), rep(c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand'), each = length(use.years)), sep = '')


#Set up p0 - use demography only since it has better AIC than exposure
if(exposure.switch == TRUE){
  p0.master =  demography.master*rbind(exposure.master, exposure.master, exposure.master, exposure.master, exposure.master, exposure.master)
}else{
  if(exposure.switch == FALSE) p0.master = demography.master
}

exposure.master = rbind(exposure.master, exposure.master, exposure.master, exposure.master, exposure.master, exposure.master)


## This is where I generate weights, but because it takes a few minutes, I generally save the outputs and load a saved file instead of re-generating the weights each time.
# source('Infection.age.structure.AB.R')
# wts = get.type.weights.AB(years.out = use.years, Countries.out = c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand'), type = 5)
# weights.master.1 = wts[[1]]
# weights.master.2 = wts[[2]]
# weights.master.3 = wts[[3]]
# weights.master.naiive = wts[[4]]


#Load saved data
load('H5N1_weights_Jan31.RData')
load('H5N1_cleaned_9weights.RData')

#drop country-years with less than 1 case
drop.me = which(rowSums(data.master) < 1)
data.master = data.master[-drop.me,]
deaths.master = deaths.master[-drop.me,]
p0.master = p0.master[-drop.me,]
demography.master = demography.master[-drop.me,]
Uc.master = Uy.master[-drop.me, ]
Um.master = Um.master[-drop.me, ]
Ue.master = Uo.master[-drop.me, ]
exposure.master = exposure.master[-drop.me, ]
Ua.master = Um.master
rm(Uy.master)
rm(Uo.master)

# Only run this immediately after generating new rates
# weights.master.1 = weights.master.1[-drop.me, ]
# weights.master.2 = weights.master.2[-drop.me, ]
# weights.master.3 = weights.master.3[-drop.me, ]
# weights.master.naiive = weights.master.naiive[-drop.me, ]
# save(weights.master.1, weights.master.2, weights.master.3, weights.master.naiive, file = 'H5N1_weights_77-78-H1dominant.RData')
#

#
# master.g1p_11 = master.g1p_11[-drop.me, ]
# master.g1p_10 = master.g1p_10[-drop.me, ]
# master.g1p_01 = master.g1p_01[-drop.me, ]
# master.g1p_00 = master.g1p_00[-drop.me, ]
# master.g1p_1n = master.g1p_1n[-drop.me, ]
# master.g1p_0n = master.g1p_0n[-drop.me, ]
# master.g2p_11 = master.g2p_11[-drop.me, ]
# master.g2p_10 = master.g2p_10[-drop.me, ]
# master.g2p_01 = master.g2p_01[-drop.me, ]
# master.g2p_00 = master.g2p_00[-drop.me, ]
# master.g2p_1n = master.g2p_1n[-drop.me, ]
# master.g2p_0n = master.g2p_0n[-drop.me, ]
# master.p_nn = master.p_nn[-drop.me, ]
# master.g1p_11[which(is.na(master.g1p_11))] = 0
# master.g1p_01[which(is.na(master.g1p_01))] = 0
# master.g1p_10[which(is.na(master.g1p_10))] = 0
# master.g1p_00[which(is.na(master.g1p_00))] = 0
# master.g1p_1n[which(is.na(master.g1p_1n))] = 0
# master.g1p_0n[which(is.na(master.g1p_0n))] = 0
# master.g2p_11[which(is.na(master.g2p_11))] = 0
# master.g2p_01[which(is.na(master.g2p_01))] = 0
# master.g2p_10[which(is.na(master.g2p_10))] = 0
# master.g2p_00[which(is.na(master.g2p_00))] = 0
# master.g2p_1n[which(is.na(master.g2p_1n))] = 0
# master.g2p_0n[which(is.na(master.g2p_0n))] = 0
# master.p_nn[which(is.na(master.p_nn))] = 0
# save(master.g1p_11, master.g1p_10, master.g1p_01, master.g1p_00, master.g1p_1n, master.g1p_0n, master.g2p_11, master.g2p_10, master.g2p_01, master.g2p_00, master.g2p_1n, master.g2p_0n, master.p_nn, file = 'H5N1_cleaned_9weights.RData')

## Save the data
# save(demography.master, file = 'Census_H5N1_by.RData')
# save(data.master, file = 'H5N1Incidence_by.RData')
# save(deaths.master, file = 'H5N1Mortality_by.RData')

