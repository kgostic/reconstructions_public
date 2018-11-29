### This script populates and then saves master matrices for all imprinting cases in the two hit model
### TAKES A WHILE TO RUN IF MULTIPLE COUNTRIES, YEARS DESIRED. LOAD SAVED RESULTS IF POSSIBLE. SAVE RESULTS AFTER RUNNING.
### Code at the bottom plots the reconstruction outputs.
### This script is essentially a wrapper that repeatedly calls the function defined in ThreeHitModel/01func-get-yr-and....
### Outputs must be repeated for each unique country and year in the data
setwd('~/Dropbox/R/Reconstructions/')
source('0func-country_data_import.R')
source('ThreeHitModel/01func-get-yr-and-subtype-specific-probs.R')

## Define a vector of countries for which outputs are desired.
## See 0func-country-data-import in the main project folder for a list of available countries
## See readme in main project folder for workflow to add new countries
Countries.out = c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand')

## Define years of case observation. Include all years in which data was observed.
## Imprinting probabilities in children will change from year to year as children get older and have more exposure opportunities.
Years.out = as.integer(c(1997, 2003:2017))

# Define birth years of interest. We generally care about birth years from 1918:present
# Scripts are currently set up to do reconstructions up to 2017
birth.years = 1918:2017




######################################################################################
# 2.a Initialize master matrices
#Rows - which country and year are we doing the reconstruction from the persepctive of?
#Cols - what birth year are we estimating imprinting probabilities for?

# Initialize a matrix in which to store probs of HxNx imprinting
# master.Hxp_110 stores probs of 1st and 2nd exposure to group x, 3rd exposure to other group
# master.Hxp_100 stores probs of 1st exposure to group x, 2nd and 3rd exposure to other group
# ...etc.
master.g1p_111 = master.g1p_110 = master.g1p_101 = master.g1p_100 = master.g1p_011 = master.g1p_010 = master.g1p_001 = master.g1p_000 =
master.g2p_111 = master.g2p_110 = master.g2p_101 = master.g2p_100 = master.g2p_011 = master.g2p_010 = master.g2p_001 = master.g2p_000 =  matrix(0, nrow = length(Countries.out)*length(Years.out), ncol = length(birth.years), dimnames = list(paste(rep(Years.out, length(Countries.out)), rep(Countries.out, each = length(Years.out)), sep = ''), rev(birth.years)))

######################################################################################
# 2.b Fill in the master matrices:
for(cc in Countries.out){ # COUNTRY LOOP
  # Call the cocirculation data from the country of interest
  flu_circ_input = import.country.dat(cc)
  
  for(yy in Years.out){ #INCIDENCE YEAR LOOP
    for(bb in birth.years){ #BIRTH YEAR LOOP
      
      if(bb > yy){ # break if the birth year hasn't yet happened
        break
      }
      
      # Output probs of every possible 2-hit combo for H1N1, H2N2 and H3N2 exposures
      LISTofall = c(partition_to_subtype(b_yr = bb, i_yr = yy, flu_circ_input))
      
      
      ## Fill in the master matrices initialized above
      row.string = paste(yy, cc, sep = '') # This will yield a string, e.g. "2013Vietnam". Use that string to call the named row of the master matrix below.
      
      # Fill in correct row and column of each output matrix by taking the sum of the relevant probabilities
      master.g1p_111[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[12]_H[12]_H[12]', x = names(LISTofall))])
      master.g1p_110[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[12]_H[12]_H?[3n]', x = names(LISTofall))])
      master.g1p_101[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[12]_H?[3n]_H[12]', x = names(LISTofall))])
      master.g1p_100[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H[12]_H?[3n]_H?[3n]', x = names(LISTofall))])
      master.g1p_011[row.string, as.character(bb)] =  Reduce(sum, LISTofall[grep(pattern = 'H?[3n]_H[12]_H[12]', x = names(LISTofall))])
      master.g1p_010[row.string, as.character(bb)] =  Reduce(sum, LISTofall[grep(pattern = 'H?[3n]_H[12]_H?[3n]', x = names(LISTofall))])
      master.g1p_001[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H?[3n]_H?[3n]_H[12]', x = names(LISTofall))])
      master.g1p_000[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H?[3n]_H?[3n]_H?[3n]', x = names(LISTofall))])
      #(master.g1p_111 + master.g1p_110 + master.g1p_101 + master.g1p_100 + master.g1p_011 + master.g1p_010 + master.g1p_001 + master.g1p_000)[row.string, as.character(bb)]

      
      master.g2p_000[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H?[12n]_H?[12n]_H?[12n]', x = names(LISTofall))])
      master.g2p_001[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H?[12n]_H?[12n]_H3', x = names(LISTofall))])
      master.g2p_010[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H?[12n]_H3_H?[12n]', x = names(LISTofall))])
      master.g2p_011[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H?[12n]_H3_H3', x = names(LISTofall))])
      master.g2p_100[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H3_H?[12n]_H?[12n]', x = names(LISTofall))])
      master.g2p_101[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H3_H?[12n]_H3', x = names(LISTofall))])
      master.g2p_110[row.string, as.character(bb)] = Reduce(sum, LISTofall[grep(pattern = 'H3_H3_H?[12n]', x = names(LISTofall))])
      master.g2p_111[row.string, as.character(bb)] = Reduce(sum, LISTofall$H3_H3_H3)
      #(master.g2p_111 + master.g2p_110 + master.g2p_101 + master.g2p_100 + master.g2p_011 + master.g2p_010 + master.g2p_001 + master.g2p_000)[row.string, as.character(bb)]

      
    } # CLOSE BIRTH YEAR LOOP
  }# CLOSE INCIDENCE YEAR LOOP
} # CLOSE COUNTRY LOOP


## Last updated 7 Nov 2018
## Edited to verify accounting of naive individiuals
save(list = grep(pattern = "master.\\w+_\\w+", ls(), value = TRUE), file = "Three_hit_weights.RData")


## check
pdf('Threehit_reconstructions.pdf')
par(mfrow = c(2,1))
xx = barplot(rbind(colSums(master.g1p_111, na.rm = TRUE),
              colSums(master.g1p_110, na.rm = TRUE),
              colSums(master.g1p_101, na.rm = TRUE),
              colSums(master.g1p_011, na.rm = TRUE),
              colSums(master.g1p_100, na.rm = TRUE),
              colSums(master.g1p_010, na.rm = TRUE),
              colSums(master.g1p_001, na.rm = TRUE),
              colSums(master.g1p_000, na.rm = TRUE)),
              main = 'Group 1', 
              col = c('navy', 'royalblue', 'blue2', 'lightblue', 'green', 'limegreen', 'lightgreen', 'orange'), border = NA, space = 0)
#lines(xx, apply(X = master.g1p_00, MARGIN = 2, FUN = function(ff) sum(!is.na(ff))), col = 'hotpink')
legend('topright', fill = c('navy', 'royalblue', 'blue2', 'lightblue', 'green', 'limegreen', 'lightgreen', 'orange'), legend = c('111', '110', '101', '011', '100', '101', '001', '000'), bg = 'white', cex = .7)



xx = barplot(rbind(colSums(master.g2p_111, na.rm = TRUE),
                   colSums(master.g2p_110, na.rm = TRUE),
                   colSums(master.g2p_101, na.rm = TRUE),
                   colSums(master.g2p_011, na.rm = TRUE),
                   colSums(master.g2p_100, na.rm = TRUE),
                   colSums(master.g2p_010, na.rm = TRUE),
                   colSums(master.g2p_001, na.rm = TRUE),
                   colSums(master.g2p_000, na.rm = TRUE)),
             main = 'Group 2', 
             col = c('firebrick4', 'red3', 'red', 'tomato', 'khaki1', 'goldenrod1', 'yellow', 'yellow4'), border = NA, space = 0)
#lines(xx, apply(X = master.g1p_00, MARGIN = 2, FUN = function(ff) sum(!is.na(ff))), col = 'hotpink')
legend('topright', fill = c('firebrick4', 'red3', 'red', 'tomato', 'khaki1', 'goldenrod1', 'yellow', 'yellow4'), legend = c('111', '110', '101', '011', '100', '101', '001', '000'), bg = 'white', cex = .7)

dev.off()


## Test
mget(ls(pattern = "master.g1p_\\d\\d\\d")) -> allg1
array.sum = function(ll){ll[[1]]+ll[[2]]+ll[[3]]+ll[[4]]+ll[[5]]+ll[[6]]+ll[[7]]+ll[[8]]}
array.sum(allg1) -> totals
totals[which(!(round(totals, 5) == 0 | round(totals, 5) == 1))]
