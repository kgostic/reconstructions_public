## Cleaned up cocirculation import

## Import master spreadsheet
cocirculation = read.csv('../CocirculationData.csv', header = TRUE)
years = as.character(seq(1997, 2017))
source('get.country.cocirculation.data_2017.R')

get.cocirculation.ref = function(Country, region = 'default'){
## Output: matrix with years across columns, H1, H2, H3, G1, G2 down rows for 2017:1918
## Input: country of interest
##     Reads and re-formats relevant data from master .csv files saved in directory
##     Region: "Asia" takes the average from all Asian countries as a fallback
##             "Euro" takes average from all european countries as a fallback
##              NA or any other input takes Thompson data as a fallback

# 1. Load template from .csv
#      Template includes all years from 1976:1918, which are fixed across countries
template = as.data.frame(read.csv('Cocirculation_template_2017.csv', header = T))
rownames(template) = template[,1]; template = template[,-1]
colnames(template) = 2017:1901

# 2. Fill in 1996:1977 data from Thompson paper
Thompson.data = read.csv('Thompson_data.csv', header = FALSE, skip = 1, col.names = c('Year', 'Group1', 'Group2', 'Source', 'X'), row.names = as.character(1977:2017))[,1:3]

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
Asia.fallback = read.csv('Asia_fallback.csv', header = FALSE, skip = 1, col.names = c('Year', 'Group1', 'Group2', 'Source'), row.names = as.character(1997:2017))[,1:3]

Euro.fallback = read.csv('Euro_fallback.csv', header = FALSE, skip = 1, col.names = c('Year', 'Group1', 'Group2', 'Source'), row.names = as.character(1997:2017))[,1:3]

#    Figure out which years are still NA, and replace

if(any(is.na(template['Group1', ]))){
  proxy.years = as.character(2017:1997)[which(is.na(template['Group1', as.character(2017:1997)]))]
  if(region == 'Asia'){ # If Asian country, substitute pooled data from asian countries in years where not enough samples exist for a single country
    template['Group1', proxy.years] = Asia.fallback[proxy.years, 'Group1']
    template['H1', proxy.years] = Asia.fallback[proxy.years, 'Group1'] #H1 was the only G1 virus circulating
    template['H2', proxy.years] = 0                                    #H2 was not circulating
    template['Group2', proxy.years] = Asia.fallback[proxy.years, 'Group2']
    template['H3', proxy.years] = Asia.fallback[proxy.years, 'Group2'] #H3 was the only G2 virus circulating
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

