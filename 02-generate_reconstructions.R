rm(list = ls())
source('01-reconstruct-imprinting-functions.R')

## OUTPUT one-hit weights for analysis of H5N1 and H7N9 data
outfile1 = paste('exampleReconstructions_generated', Sys.Date(), '.RData', sep = '')
outfile1


# Outputs
yo = c(2011:2015)             ## Set vector of observation years (years in which imprinting status is observed, for all birth years from 1918:2017)
co = c('China', 'USA')        ## Set vector of countries. See 00-country_data_import_functions.R for countries currently available.
regions = c('Asia', NULL)     ## Set vector of regions. User is reponsible for ensuring this vector is the same length as the countries vector above. Options are "Asia", "Europe", or NULL (for all other regions). When the country of interest hasn't reported enough surveillance data to estimate the fraction of influenza A circulation caused by H1N1 or H3N2 in the post-1977 era, regions direct the code to fall back on data from Asia, from Europe, or from the USA (the default, corresponds to NULL input)
wts = get.weights.master(years.out = yo, Countries.out = co, region.in = regions)
attach(wts)
## wts is a list containing elements weights.master.H1N1, weights.master.H2N2, weights.master.H3N2, weights.master.naive
## entries in these matrices represent the fraction of individuals born in each year (columns) and observed in each year (rows) that imprinted to the subtype of interest, or remains naive.
save(weights.master.1, weights.master.2, weights.master.3, weights.master.naiive, file = outfile1)
