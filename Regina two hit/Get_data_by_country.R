get.country.data.H5N1 = function(Country, use.years = c(1997, 2003:max.yr), birth.yrs = max.yr:1918){
  
  ## 1997
  H5.1997 = read.csv('H5N1_Linelist_1997cases.csv', header = T, stringsAsFactors = FALSE)

  ## 2006-2010
  H5.RKI = read.csv('H5N1_Final_Linelist_RKI_2006_2010cases.csv', header = T, stringsAsFactors = FALSE)

  ## all others
  H5.new = read.csv('H5N1_linelist_KG_MRA.csv', header = T, stringsAsFactors = FALSE)
  

### EXTRACT AGE DATA, COUNTRY AND OUTCOME
raw.cases.out = data.frame(year = c(subset(H5.1997, country == Country, select = year)$year,
                                      subset(H5.RKI, country == Country, select = year)$year, 
                                      subset(H5.new, country == Country, select = year)$year),
                             age = c(subset(H5.1997, country == Country, select = age)$age,
                                     subset(H5.RKI, country == Country, select = age)$age,
                                     subset(H5.new, country == Country, select = age)$age),
                             outcome = c(subset(H5.1997, country == Country, select = outcome)$outcome,
                                         subset(H5.RKI, country == Country, select = outcome)$outcome,
                                         subset(H5.new, country == Country, select = outcome)$outcome))

raw.cases.out$birth.year = raw.cases.out$year-raw.cases.out$age


#get rid of NA ages
if(any(is.na(raw.cases.out$age) == TRUE)){
drop = which(is.na(raw.cases.out$age))
raw.cases.out = raw.cases.out[-drop, ]
}


cases.by = matrix(NA, nrow = length(use.years), ncol = (max.yr-1918+1), dimnames = list(use.years, rev(1918:max.yr)))
for(jj in 1:length(use.years)){ #For each year of incidence...tabluate birth years
  birth.years = subset(raw.cases.out, year == use.years[jj], select = birth.year)$birth.year
  for(ii in 1:(max.yr-1918+1)){cases.by[jj, ii] = sum(birth.years == birth.yrs[ii])}
}


cases.by
}