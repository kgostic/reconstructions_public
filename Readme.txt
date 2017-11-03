Readme

This repository contains files that can be used to reconstruct imprinting patterns.
Current data allows reoncstrucitons up to 2017.



Classes of .csv files:
-----------------------
.*cocirculation.csv contains raw data dowlonaded from WHO Flu Net for the countrie(s) and years of interest

.*summary.csv contains a reformatted version of the raw data in .*cocirculation.csv. At one point, I was creating these reformatted files and then copying the relevant parts into a master .csv file.

CocirculationData.csv is the master file that contains a summary of all virological data from all countries and all years up to 2017. This file is read by Clean_CocirculationImport_2017.R



Notes on .R files:
-------------------
Clean_CocirculationImport_2017.R    Imports virological surveillance data from CocirculationData.csv and adds this data to a standard cocirculation template. Sources functions that generate country-specific tables of H1N1, H2N2 and H3N2 dominance in a certain year.
- Pre-1977 we assume only one subtype circulated each year.
- From 1977-~1997 we use summary data from the United States
- From 1997 on, we download country-specific data when possible, or substitute summary data from the surrounding region.
