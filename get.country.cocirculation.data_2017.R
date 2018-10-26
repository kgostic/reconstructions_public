## Virological surveillance data from flu net needs to be reformatted, and we need to calculate the fraction of IAV circulation in a given country and year caused by a particular seasonal subtype.
## This function performs these tasks
##  
## INPUTS: country, string naming country of interest. Signals function to pull data from the country of interest from the master .csv. 
##         currently, choose form 

## Output - Matrix with possible subtypes on rows (A, H1 or H3), and with year of observation on columns

get.country.data = function(country){
  years = 1997:2017
  types = c('All', 'A', 'H1', 'H3')
  country.dat = matrix(NA, nrow = length(types), ncol = length(years), dimnames = list(types, years))
  for(ii in 1:length(years)){
    #Does the year exist for the given country
    #If yes, fill it in with the total/percents
    #If no, fill in NA
    # -- Total       ## Inserted NA because I no longer use "All" but don't want to change the output format
    country.dat['All', ii] = NA
    # -- A
    country.dat['A', ii] = ifelse(years[ii] %in% subset(cocirculation, Country == country & Subtype == 'A')$Year,
                                  subset(cocirculation, Country == country & Subtype == 'A' & Year == years[ii])$Observed, NA)
    # -- H1
    country.dat['H1', ii] = ifelse(years[ii] %in% subset(cocirculation, Country == country & Subtype == 'H1')$Year | 
                                     years[ii] %in% subset(cocirculation, Country == country & Subtype == 'H1p09')$Year,
                                   sum(subset(cocirculation, Country == country & Subtype %in% c('H1', 'H1p09') & Year == years[ii])$Observed)/
                                     sum(subset(cocirculation, Country == country & Subtype %in% c('H1', "H3", "H1p09") & Year == years[ii])$Observed), NA)
    # -- H3
    country.dat['H3', ii] = ifelse(years[ii] %in% subset(cocirculation, Country == country & Subtype == 'H3')$Year,
                                   sum(subset(cocirculation, Country == country & Subtype == 'H3' & Year == years[ii])$Observed)/
                                     sum(subset(cocirculation, Country == country & Subtype %in% c('H1', "H3", "H1p09") & Year == years[ii])$Observed), NA)
  }
  ## Return 
  country.dat
}