dat = read.csv('CocirculationData.csv', stringsAsFactors = FALSE)

# Extract data from European countries
valid = subset(dat, Country %in% c('Austria', 'Belgium', 'Denmark', 'Estonia', 'Greece', 'Norway', 'Poland', 'Portugal', 'Spain', 'Germany'))
#valid = subset(dat, Country %in% c('Bangladesh', 'Cambodia', 'China', 'Indonesia', 'Thailand', 'Vietnam', 'Japan'))

# Initalize
yrs = 1997:2017
output = matrix(NA, nrow = length(yrs), ncol = 2, dimnames = list(yrs, c('Group1', 'Group2')))


for(ii in 1:length(yrs)){
  total.g1 = sum(valid$Observed[valid$Year == yrs[ii] & valid$Subtype %in% c('H1', 'H1p09')])
  total.g2 = sum(valid$Observed[valid$Year == yrs[ii] & valid$Subtype == 'H3'])
  
  output[ii, ] = c(total.g1/(total.g2+total.g1), total.g2/(total.g2+total.g1))
  print(total.g1+total.g2)
}

write.csv(output, file = 'formatted_fallback.csv')

