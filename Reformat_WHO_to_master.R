## Write a function to reformat WHO data
rm(list = ls())
setwd('~/Dropbox/R/Reconstructions/')

## Load data
raw = read.csv('~/Dropbox/R/Reconstructions/Belgium_cocirculation.csv', skip = 3, stringsAsFactors = FALSE)




# Read available countries
countries = sort(unique(raw$Country))


for(cc in 1:length(countries)){
# Extract data from country of interest
valid = subset(raw, Country == countries[cc])

# Read available years
yrs = sort(unique(valid$Year))  
  
# Initialize output data frame
c1 = rep(countries[cc], each = length(yrs)*4)
c2 = rep(yrs, each = 4)
c3 = rep(c('H1', 'H1p09', 'H3', 'A'), times = length(yrs))
c4 = integer(length(yrs*4))

# Fill in the number of samples for each year and subtype
for(yy in 1:length(yrs)){
  valid2 = subset(valid, Year == yrs[yy], select = c('AH1', 'AH1N12009', 'AH3'))
  yr.dat = c(colSums(valid2, na.rm = TRUE), sum(valid2, na.rm = TRUE))
  c4[(yy-1)*4 + 1:4] = as.integer(yr.dat)
}
# Format this country's results
output = data.frame(c1, c2, c3, c4)
# If only one country output the results
# Else, rbind results form each country and output
if(cc == 1){
  full.list = output
}else{
  full.list = rbind(full.list, output)
}

}

full.list

write.csv(full.list, file = 'Reformatted_output.csv')
