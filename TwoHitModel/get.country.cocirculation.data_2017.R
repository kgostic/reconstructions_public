get.country.data = function(country){
  years = 1997:2017
  types = c('All', 'A', 'H1', 'H3')
  country.dat = matrix(NA, nrow = length(types), ncol = length(years), dimnames = list(types, years))
  for(ii in 1:length(years)){
    #Does the year exist for the given country
    #If yes, fill it in with the total/percents
    #If no, fill in NA
    # -- Total       ## Commented these out because I don't think I use them
    country.dat['All', ii] = NA #ifelse(years[ii] %in% subset(cocirculation, Country == country & Subtype == 'All')$Year,
                                    #subset(cocirculation, Country == country & Subtype == 'All' & Year == years[ii])$Observed, NA)
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
  #country.dat[2:5,] = country.dat[2:5,]/country.dat[1,] #Return proportions
  country.dat
}