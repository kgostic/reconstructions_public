### This script populates and then saves master matrices for all imprinting cases in the two hit model
### Load results from ""
### TAKES SEVERAL MINUTES TO RUN. LOAD SAVED RESULTS IF POSSIBLE.
### Code at the bottom plots the reconstruction outputs.
setwd('~/Dropbox/R/Reconstructions/ThreeHitModel/')
source('country_data_import.R')
source('compute_proportion_matrices_threehit.R')

# 1. To calculate outputs for each country and year of interest, use get_proportion(b_yr, i_yr, flu_circ)
#   Note that each time you call this function you will need to import data from each relevant country.
#   Use of flu_circ as a variable to input the proper matrix each time you want to calculate HXp_xx
#   i.e. flu_circ = import.country.dat('Vietnam', region.in = 'Asia')

# 2. Once you are able to use the pseudocode script to output H1p_11, H1p_10, ... etc., we can start filling in the master matrices.
#   Note that we have data from 1997, and 2003:2017, and from China, Cambodia, Egypt, Indonesia, Vietnam and Thailand. Set:
Countries.out = c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand')
#Countries.out = 'China'
Years.out = as.integer(c(1997, 2003:2017))
#Years.out = as.integer(2011)
#
# We also care about birth years from 1918:2017
birth.years = 1918:2017




######################################################################################
# 2.a Initialize master matrices
#Rows - which country and year are we doing the reconstruction from the persepctive of?
#Cols - what birth year are we estimating imprinting probabilities for?

# Initialize a matrix in which to store probs of HxNx imprinting
# master.Hxp_11 stores probs of 1st and 2nd exposure to Hx
# master.Hxp_10 stores probs of 1st exposure to Hx, but 2nd exposure to some other subtype
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
