#this is the most current code
source('country_data_import.R')
source('compute_proportion_matrices.R')

# 1. To calculate 9 outputs for each country and year of interest, use get_proportion(b_yr, i_yr, flu_circ)
#   Note that each time you call this function you will need to import data from each relevant country.
#   Use of flu_circ as a variable to input the proper matrix each time you want to calculate HXp_xx
#   i.e. flu_circ = import.country.dat('Vietnam', region.in = 'Asia')

# 2. Once you are able to use the pseudocode script to output H1p_11, H1p_10, ... etc., we can start filling in the master matrices.
#   Note that we have data from 1997, and 2003:2017, and from China, Cambodia, Egypt, Indonesia, Vietnam and Thailand. Set:
Countries.out = c('China', 'Cambodia', 'Egypt', 'Indonesia', 'Vietnam', 'Thailand')
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
# master.HXp_xx stores probs of HX imprinting
# master.naive stores probs of no 2nd influenza exposure (no imprintng). This will only take nonzero values for very young birth cohorts.
# Each country fills in diff row of matrix!!!
# master.H1p_11 = matrix(NA, nrow = length(Countries.out)*length(Years.out), ncol = length(birth.years), dimnames = list(paste(rep(Years.out, length(Countries.out)), rep(Countries.out, each = length(Years.out)), sep = ''), rev(birth.years))) 
# master.H1p_10 = master.H1p_01 = master.H1p_00 = master.H2p_11 = master.H2p_10 = master.H2p_01 = master.H2p_00 = master.H3p_11 = master.H3p_10 = master.H3p_01 = master.H3p_00 = master.naive = master.H1p_11
master.g1p_11 = matrix(NA, nrow = length(Countries.out)*length(Years.out), ncol = length(birth.years), dimnames = list(paste(rep(Years.out, length(Countries.out)), rep(Countries.out, each = length(Years.out)), sep = ''), rev(birth.years))) 
master.g1p_10 = master.g1p_01 = master.g1p_00 = master.g2p_11 = master.g2p_10 = master.g2p_01 = master.g2p_00 = master.g1p_1n = master.g1p_0n = master.g2p_1n = master.g2p_0n = master.p_nn = master.g1p_11 

######################################################################################
# 2.b Fill in the master matrices:

for(cc in Countries.out){ # COUNTRY LOOP
  # Call the cocirculation data from the country of interest
  flu_circ_input = import.country.dat(cc)
  
  for(yy in Years.out){ #INCIDENCE YEAR LOOP
    for(bb in birth.years){ #BIRTH YEAR LOOP
      
      if(bb > yy){ #check this, used to be bb >= yy but needs to be fixed
        break
      }
      #print(c(yy,bb))
      
      LISTof9 = c(get_proportion(b_yr = bb, i_yr = yy, flu_circ_input))
      
      ## Fill in the master matrices initialized above
      row.string = paste(yy, cc, sep = '') # This will yield a string, e.g. "2013Vietnam". you can then use that string to call the named row of the master matrix below.
      
      # Fill in correct row and column of each output matrix by taking the sum of the relevant probabilities
      master.g1p_11[row.string, as.character(bb)] = sum(LISTof9$g1p_11, na.rm = TRUE)
      master.g1p_10[row.string, as.character(bb)] = sum(LISTof9$g1p_10, na.rm = TRUE)
      master.g1p_01[row.string, as.character(bb)] = sum(LISTof9$g1p_01, na.rm = TRUE)
      master.g1p_00[row.string, as.character(bb)] = sum(LISTof9$g1p_00, na.rm = TRUE)
      master.g1p_1n[row.string, as.character(bb)] = sum(LISTof9$g1p_1n, na.rm = TRUE) 
      master.g1p_0n[row.string, as.character(bb)] = sum(LISTof9$g1p_0n, na.rm = TRUE) 
      #master.H2p_11[row.string, as.character(bb)] = sum(LISTof13$H2p_11, na.rm = TRUE)  
      #master.H2p_10[row.string, as.character(bb)] = sum(LISTof13$H2p_10, na.rm = TRUE)  
      #master.H2p_01[row.string, as.character(bb)] = sum(LISTof13$H2p_01, na.rm = TRUE)  
      #master.H2p_00[row.string, as.character(bb)] = sum(LISTof13$H2p_00, na.rm = TRUE)  
      master.g2p_11[row.string, as.character(bb)] = sum(LISTof9$g2p_11, na.rm = TRUE)  
      master.g2p_10[row.string, as.character(bb)] = sum(LISTof9$g2p_10, na.rm = TRUE)  
      master.g2p_01[row.string, as.character(bb)] = sum(LISTof9$g2p_01, na.rm = TRUE)  
      master.g2p_00[row.string, as.character(bb)] = sum(LISTof9$g2p_00, na.rm = TRUE)
      master.g2p_1n[row.string, as.character(bb)] = sum(LISTof9$g2p_1n, na.rm = TRUE) 
      master.g2p_0n[row.string, as.character(bb)] = sum(LISTof9$g2p_0n, na.rm = TRUE)
      master.p_nn[row.string, as.character(bb)] = sum(LISTof9$p_nn, na.rm = TRUE) 
      #master.naive[row.string, as.character(bb)] = sum(LISTof9$naive, na.rm = TRUE) 
    } # CLOSE BIRTH YEAR LOOP
  }# CLOSE INCIDENCE YEAR LOOP
} # CLOSE COUNTRY LOOP







