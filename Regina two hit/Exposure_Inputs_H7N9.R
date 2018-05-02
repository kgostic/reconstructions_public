###########################################################
#####     Import exposure data - Cowling et al, 2013  #####
#####     LATER - Look at Rivers Data again and 
#####             see if you can incorporate it.
###########################################################
#This is the cowling data, stratified by location
exposures = read.csv('Exposures.csv', skip = 1)
attach(exposures)
n.old.reps = length(max.yr:1918)-sum(15, 10, 10, 10, 10, 10)
expand = function(x){rep(x, times = c(15, 10, 10, 10, 10, 10, n.old.reps))}

#Expand exposures data frame to repeat the rate once for each age group
exposures.by = as.data.frame(apply(exposures[,-c(1:3)], 2, expand)); rownames(exposures.by) = 0:(length(max.yr:1918)-1)

#Normalize each city
normalize = function(x){x/sum(x)}
exposures.by = apply(exposures.by, 2, normalize)

#Take average across rows
exposures.by = rowMeans(exposures.by)

detach(exposures)


#barplot(exposure, main = 'Final Distribution', names.arg = age.classes)

get.exposure.data = function(){exposures.by}