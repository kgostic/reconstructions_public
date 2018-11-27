## Three-hit reconstructions for Haley

### Code at the bottom plots the reconstruction outputs.
setwd('~/Dropbox/R/Reconstructions/ThreeHitModel/')
source('country_data_import.R')  ## Source functions to import data on which subtypes cicrculated in which seasons in different countries
source('compute_proportion_matrices_threehit.R') ## Source functions that compute the probability of first, second and third exposure at 
## age x, y, and z (i.e. first exposure occurring x, y and z years after birth.) Then cross-reference these probabilities with data on
## the subtypes circulating in years x, y and z, to calculate overall probabilities of first, second and third exposure to a particular subtype.

Countries.out = c('USA') ## All subjects from USA
Years.out = 2017 ## Reconstruct from the perspective of 2017. Subjects are not young enough to be influenza naive, so a choice of 2016 or 2015 would yield identical reconstructions.

# Data contain birth years from 1968:2000
birth.years = 1968:2017
######################################################################################
# 2.a Initialize master matrices
#Rows - which country and year are we doing the reconstruction from the persepctive of?
#Cols - what birth year are we estimating imprinting probabilities for?

# Initialize a matrix in which to store probs of HxNx imprinting
# master.Hxp_110 stores probs of 1st and 2nd exposure to Hx, and 3rd exposure to some other subtype
# master.Hxp_101 stores probs of 1st exposure to Hx, but 2nd exposure to some other subtype, and third exposure to Hx
# ...etc.
master.H1p_111 = master.H1p_110 = master.H1p_101 = master.H1p_100 = master.H1p_011 = master.H1p_010 = master.H1p_001 = master.H1p_000 =
  master.H3p_111 = master.H3p_110 = master.H3p_101 = master.H3p_100 = master.H3p_011 = master.H3p_010 = master.H3p_001 = master.H3p_000 =  
  master.H1_H1_H1 = master.H1_H1_H3 = master.H1_H3_H1 = master.H1_H3_H3 = master.H3_H1_H1 = master.H3_H1_H3 = master.H3_H3_H1 = master.H3_H3_H3 = matrix(0, nrow = length(Countries.out)*length(Years.out), ncol = length(birth.years), dimnames = list(paste(rep(Years.out, length(Countries.out)), rep(Countries.out, each = length(Years.out)), sep = ''), rev(birth.years)))


######################################################################################
# 2.b Fill in the master matrices:
for(cc in Countries.out){ # COUNTRY LOOP
  ## Note, this loop is unnecessary because we're only importing for USA, but I'll leave it
  ## in case you want to do other countries later.
  
  # Call the cocirculation data from the country of interest
  flu_circ_input = import.country.dat(cc) ## Function defined in "country_data_import.R"
  
  for(yy in Years.out){ #INCIDENCE YEAR LOOP
    ## Again, unnecessary because we are only using 2017, but I"ll leave this in case you
    ## want to try multiple years
    
    for(bb in birth.years){ #BIRTH YEAR LOOP
      print(bb)
      if(bb > yy){ # break if the birth year hasn't yet happened
        stop(paste('birth year', bb, '> incidence year', yy, '. Birth year has not yet occurred.'))
      }
      
      # Output probs of every possible 3-hit combo for H1N1, H2N2 and H3N2 exposures
      LISTofall = c(partition_to_subtype(b_yr = bb, i_yr = yy, flu_circ_input))
      
      
      ## Fill in the master matrices initialized above
      row.string = paste(yy, cc, sep = '') # This will yield a string, e.g. "2013Vietnam". Use that string to call the named row of the master matrix below.
      # Fill in correct row and column of each output matrix by taking the sum of the relevant probabilities
      master.H1p_111[row.string, as.character(bb)] = sum(LISTofall$H1_H1_H1, na.rm = TRUE) 
      master.H1p_110[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H1_H1_H?[23n]', x = names(LISTofall))])
      master.H1p_101[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H1_H?[23n]_H1', x = names(LISTofall))])
      master.H1p_100[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H1_H?[23n]_H?[23n]', x = names(LISTofall))])
      master.H1p_011[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[23n]_H1_H1', x = names(LISTofall))])
      master.H1p_010[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[23n]_H1_H?[23n]', x = names(LISTofall))])
      master.H1p_001[row.string, as.character(bb)] =  Reduce(sum, LISTofall[grep(pattern = 'H[23]_H[23]_H1', x = names(LISTofall))])
      master.H1p_000[row.string, as.character(bb)] =  Reduce(sum, LISTofall[grep(pattern = 'H?[23n]_H?[23n]_H?[23n]', x = names(LISTofall))])
      # ## Test, should all sum to 1 for the birth year
      # (master.H1p_111 + master.H1p_110 + master.H1p_101 + master.H1p_100 + master.H1p_011 + master.H1p_010 + master.H1p_001 + master.H1p_000)[row.string, as.character(bb)]
      
      
      
      
      # Fill in correct row and column of each output matrix by taking the sum of the relevant probabilities
      master.H3p_111[row.string, as.character(bb)] = sum(LISTofall$H3_H3_H3, na.rm = TRUE) 
      master.H3p_110[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H3_H3_H?[12n]', x = names(LISTofall))])
      master.H3p_101[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H3_H?[12n]_H3', x = names(LISTofall))])
      master.H3p_100[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H3_H?[1n]_H?[12n]', x = names(LISTofall))])
      master.H3p_011[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[12n]_H3_H3', x = names(LISTofall))])
      master.H3p_010[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[12n]_H3_H?[12n]', x = names(LISTofall))])
      master.H3p_001[row.string, as.character(bb)] =  Reduce(sum, LISTofall[grep(pattern = 'H[12]_H[12]_H3', x = names(LISTofall))])
      master.H3p_000[row.string, as.character(bb)] =  Reduce(sum, LISTofall[grep(pattern = 'H?[12n]_H?[12n]_H?[12n]', x = names(LISTofall))])
      # (master.H3p_111 + master.H3p_110 + master.H3p_101 + master.H3p_100 + master.H3p_011 + master.H3p_010 + master.H3p_001 + master.H3p_000)[row.string, as.character(bb)]
      
      master.H1_H1_H1[row.string, as.character(bb)] =  sum(LISTofall$H1_H1_H1)
      master.H1_H1_H3[row.string, as.character(bb)] =  sum(LISTofall$H1_H1_H3)
      master.H1_H3_H1[row.string, as.character(bb)] =  sum(LISTofall$H1_H3_H1)
      master.H1_H3_H3[row.string, as.character(bb)] =  sum(LISTofall$H1_H3_H3)
      master.H3_H1_H1[row.string, as.character(bb)] =  sum(LISTofall$H3_H1_H1)
      master.H3_H1_H3[row.string, as.character(bb)] =  sum(LISTofall$H3_H1_H3)
      master.H3_H3_H1[row.string, as.character(bb)] =  sum(LISTofall$H3_H3_H1)
      master.H3_H3_H3[row.string, as.character(bb)] =  sum(LISTofall$H3_H3_H3)
      
    } # CLOSE BIRTH YEAR LOOP
  }# CLOSE INCIDENCE YEAR LOOP
} # CLOSE COUNTRY LOOP


## Reformat into a table and save as a .csv
nms = ls(pattern = "master.\\w\\w+")
nms = nms[c(1:4, 13:16, 5:12, 17:24)]
out_table = sapply(nms, get)
out_table = as.data.frame(out_table)
colnames(out_table) = gsub(colnames(out_table), pattern = 'master.(H\\d)p_(\\d\\d\\d)', replacement = '\\1:\\2')
colnames(out_table) = gsub(colnames(out_table), pattern = 'master.(H\\d)_(H\\d)_(H\\d)', replacement = '(\\1, \\2, \\3)')
out_table$BirthYear = birth.years
out_table$SubjObservedIn = 'USA_2017'
write.csv(x = out_table, file = 'Haley_reconstructions.csv')
