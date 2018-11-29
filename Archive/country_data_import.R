
################################################ FUNCTION 1 #######################################################
############################################# get.country.data ####################################################
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
############################################# END FUNCTION 1 ######################################################
############################################# get.country.data ####################################################




################################################ FUNCTION 2 #######################################################
############################################# get.cocirculation.ref ###############################################
## This function reads the output of get.country.data, and then checks that imported data contains a sufficient number
##    of reported viruses to calculate a reasonably accurate H1 and H3 proportion. If not enough cases observed, this
##    function substutites aggregate data form across all countries in the region, or from the US.
## This function also tabluates circulation fractions from 1918-2017, whereas get.country.data only tabulates
##    fractions from 1997:2017 (the years in which virological surveillance are reliably available from Flu Net)

## Input: country of interest
##     Reads and re-formats relevant data from master .csv files saved in directory
##     Region: "Asia" takes the average from all Asian countries as a fallback
##             "Euro" takes average from all european countries as a fallback
##              NA or any other input takes Thompson data as a fallback

## Output: matrix with years across columns, H1, H2, H3, G1, G2 down rows for 2017:1918

## Import master spreadsheet
cocirculation = read.csv('../CocirculationData.csv', header = TRUE)
years = as.character(seq(1997, 2017))
#source('get.country.cocirculation.data_2017.R')

get.cocirculation.ref = function(Country, region = 'default'){

  
  # 1. Load template from .csv
  #      Template includes all years from 1976:1918, which are fixed across countries
  template = as.data.frame(read.csv('../Cocirculation_template_2017.csv', header = T))
  rownames(template) = template[,1]; template = template[,-1]
  colnames(template) = 2017:1901
  
  # 2. Fill in 1996:1977 data from Thompson paper
  Thompson.data = read.csv('../Thompson_data.csv', header = FALSE, skip = 1, col.names = c('Year', 'Group1', 'Group2', 'Source', 'X'), row.names = as.character(1977:2017))[,1:3]
  
  template['Group1', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group1']
  template['H1', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group1']
  template['H2', as.character(1996:1977)] = 0
  template['Group2', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group2']
  template['H3', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group2']
  
  
  # 3. Fill in country-specific data from flu net for 2017:1997 where available
  #get.country.data is a Function that formats data to be input into template
  ## Format country -specific data
  country.data.1997_2017 = get.country.data(country = Country)
  ## Quality checks
  #    - remove data from years with < 50 specimens
  if(any(country.data.1997_2017['A', is.na(country.data.1997_2017['A', ]) == FALSE] < 50)) country.data.1997_2017[ , which(country.data.1997_2017['A', ] < 50)] = rep(NA, 4)
  
  ## Fill in template
  template['H1', as.character(2017:1997)] = country.data.1997_2017['H1', as.character(2017:1997)]
  template['H3', as.character(2017:1997)] = country.data.1997_2017['H3', as.character(2017:1997)]
  template['H2', as.character(2017:1997)] = 0
  template['Group1', as.character(2017:1997)] = colSums(template[c('H1', 'H2'), as.character(2017:1997)])
  template['Group2', as.character(2017:1997)] = colSums(template[c('H3'), as.character(2017:1997)])
  
  
  # 4. Fill in general data where not available
  #    Load fallback data
  #    fallback.data = Thompson.data[1997:2017]
  Asia.fallback = read.csv('../Asia_fallback.csv', header = FALSE, skip = 1, col.names = c('Year', 'Group1', 'Group2', 'Source'), row.names = as.character(1997:2017))[,1:3]
  Euro.fallback = read.csv('../Euro_fallback.csv', header = FALSE, skip = 1, col.names = c('Year', 'Group1', 'Group2', 'Source'), row.names = as.character(1997:2017))[,1:3]
  #    Figure out which years are still NA, and replace with fallback data.
  
  if(any(is.na(template['Group1', ]))){
    proxy.years = as.character(2017:1997)[which(is.na(template['Group1', as.character(2017:1997)]))]
    if(region == 'Asia'){ # If Asian country, substitute pooled data from asian countries in years where not enough samples exist for a single country
      template['Group1', proxy.years] = Asia.fallback[proxy.years, 'Group1']
      template['H1', proxy.years] = Asia.fallback[proxy.years, 'Group1'] #H1 was the only G1 virus circulating
      template['H2', proxy.years] = 0                                    #H2 was not circulating
      template['Group2', proxy.years] = Asia.fallback[proxy.years, 'Group2']
    }else if(region == 'Euro'){ # Same as above for European countries
      template['Group1', proxy.years] = Euro.fallback[proxy.years, 'Group1']
      template['H1', proxy.years] = Euro.fallback[proxy.years, 'Group1'] #H1 was the only G1 virus circulating
      template['H2', proxy.years] = 0                                    #H2 was not circulating
      template['Group2', proxy.years] = Euro.fallback[proxy.years, 'Group2']
      template['H3', proxy.years] = Euro.fallback[proxy.years, 'Group2'] #H3 was the only G2 virus circulating
    }else{ # Otherwise, fall back on data from the US (Thompson et al.)
      template['Group1', proxy.years] = Thompson.data[proxy.years, 'Group1']
      template['H1', proxy.years] = Thompson.data[proxy.years, 'Group1'] #H1 was the only G1 virus circulating
      template['H2', proxy.years] = 0                                    #H2 was not circulating
      template['Group2', proxy.years] = Thompson.data[proxy.years, 'Group2']
      template['H3', proxy.years] = Thompson.data[proxy.years, 'Group2'] #H3 was the only G2 virus circulating
    }
  }
  
  # Return output
  template
}
############################################## END FUNCTION 2 #####################################################
############################################# get.cocirculation.ref ###############################################









################################################ FUNCTION 3 #######################################################
############################################# import country data ###############################################
##  This wrapper combines outputs of the above two functions, and inputs only the name of the country of interest
##    
## Input - 
#         Country.out - string. Valid inputs are listed inside the function below.
#         region.in - string vector of the same length as 'countries in'. In cases where there is not enough data on which subtypes circulated in the country of interest in particular year, we are forced to substitute data from the surrounding region. This option tells the function which "fallback" data to draw from. Options are: 'default' - pulls data from the USA. 'Asia' (aggregate data from Asian countries) or 'Euro' (aggregate data from European countries). If no input is given, the function will automatically use data from the USA for all reconstructions.
## Output - Matrix of values showing the fraction of influenza viruses of different types that circulated in years from 1918:2017

#THIS REFORMATS COUNTRY DATA
import.country.dat = function(Country.out, region.in = 'default'){
  ######## 1 - IMPORT DATA ON WHICH SUBTYPES CIRCULATED IN PAST YEARS IN SPECIFIC COUNTRIES.  
  ############---------------------------------------------------------------------------  
  
  # Sourced script imports a .csv that holds data on which subtypes (H1N1, H3N2 and influenza B) were observed in specific countries, and in specific years.
  # This script also sources two functions, get.cocirculation.ref and get.country.data, which we will use below to format a matrix that tells what fraction of circulation was driven by H1N1 vs. H3N2 in a given season
  #source('~/Dropbox/R/Reconstructions/ThreeHitModel/Clean_CocirculationImport_2017.R')
    
    # Extract country-sepcific data on which subtypes circulated in past years
    if(Country.out == 'China'){cocirculation.dat = get.cocirculation.ref('China', region.in)}
    if(Country.out == 'Egypt'){cocirculation.dat = get.cocirculation.ref('Egypt', region.in)}
    if(Country.out == 'Indonesia'){cocirculation.dat = get.cocirculation.ref('Indonesia', region.in)}
    if(Country.out == 'Cambodia'){cocirculation.dat = get.cocirculation.ref('Cambodia', region.in)}
    if(Country.out == 'Vietnam'){cocirculation.dat = get.cocirculation.ref('Vietnam', region.in)}
    if(Country.out == 'Thailand'){cocirculation.dat = get.cocirculation.ref('Thailand', region.in)}
    if(Country.out == 'USA'){cocirculation.dat = get.cocirculation.ref('USA', region.in)}
    if(Country.out == 'UK'){cocirculation.dat = get.cocirculation.ref('UK', region.in)}
    if(Country.out == 'Turkey'){cocirculation.dat = get.cocirculation.ref('Turkey', region.in)}
    if(Country.out == 'Iraq'){cocirculation.dat = get.cocirculation.ref('Iraq', region.in)}
    if(Country.out == 'Azerbaijan'){cocirculation.dat = get.cocirculation.ref('Azerbaijan', region.in)}
    if(Country.out == 'Pakistan'){cocirculation.dat = get.cocirculation.ref('Pakistan', region.in)}
    if(Country.out == 'Nigeria'){cocirculation.dat = get.cocirculation.ref('Nigeria', region.in)}
    if(Country.out == 'Austria'){cocirculation.dat = get.cocirculation.ref('Austria', region.in)}
    if(Country.out == 'Belguim'){cocirculation.dat = get.cocirculation.ref('Belgium', region.in)}
    if(Country.out == 'Denmark'){cocirculation.dat = get.cocirculation.ref('Denmark', region.in)}
    if(Country.out == 'Estonia'){cocirculation.dat = get.cocirculation.ref('Estonia', region.in)}
    if(Country.out == 'Greece'){cocirculation.dat = get.cocirculation.ref('Greece', region.in)}
    if(Country.out == 'Norway'){cocirculation.dat = get.cocirculation.ref('Norway', region.in)}
    if(Country.out == 'Australia'){cocirculation.dat = get.cocirculation.ref('Australia', region.in)}
    if(Country.out == 'Poland'){cocirculation.dat = get.cocirculation.ref('Poland', region.in)}
    if(Country.out == 'Portugal'){cocirculation.dat = get.cocirculation.ref('Portugal', region.in)}
    if(Country.out == 'Spain'){cocirculation.dat = get.cocirculation.ref('Spain', region.in)}
    if(Country.out == 'Argentina'){cocirculation.dat = get.cocirculation.ref('Argentina', region.in)}
    if(Country.out == 'Peru'){cocirculation.dat = get.cocirculation.ref('Peru', region.in)}
    if(Country.out == 'Chile'){cocirculation.dat = get.cocirculation.ref('Chile', region.in)}
    if(Country.out == 'Germany'){cocirculation.dat = get.cocirculation.ref('Germany', region.in)}
    if(Country.out == 'Japan'){cocirculation.dat = get.cocirculation.ref('Japan', region.in)}
    if(Country.out == 'Singapore'){cocirculation.dat = get.cocirculation.ref('Singapore', region.in)}
  
  cocirculation.dat
}
### End function ###




# ## Examples
# import.country.dat('Vietnam', region.in = 'Asia')
# import.country.dat('USA')
    
    