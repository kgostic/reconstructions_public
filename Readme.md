Readme

This repository contains files that can be used to reconstruct imprinting patterns.
Current data allows reconstructions up to 2017.

SUMMARY OF OVERALL WORKFLOW:
----------------------------
0
Import country-specific virological surveillance data from WHO Flu Net (details below)

0func-country_data_import.R
Define functions that read raw virological surveillance data into formatted tables that describe the fraction of seasonal circulation caused by H1N1, H3N2 or H2N2 in any year from 1918-present, in the country of interest

01-reconstruct-imprinting-histories.R
Define a function to calculate the probability of a first influenza exposure in year j, given birth in year i.
Then, define a wrapper function that calculates the probability of subtype-specific imprinting for all years and countries of interest. The wrapper combines virological surveillance data from above with repeated calculations of the probabilities of imprinting for each birth year.



0 continued...WORKFLOW for adding virological surveillance data for new countries or years:
----------------------------------
1. Go to WHO Flu Net, http://apps.who.int/flumart/Default?ReportNo=12
2. Display report for the country(ies) of interest starting at week 1 of the earliest possible year and continuing to the final week of the last year of interest
3. Click on the disk icon at the top middle of the page to save the output as a .csv. The file will automatically save to Downloads.
4. Move the file to ~/Dropbox/R/Reconstructions under the filename raw-flunet-outputs/COUNTRY_cocirculation.csv
5. Open Reformat_WHO_to_master.R
6. Import the new file from step 4 into the script and run the script to reformat. The script will save a reformatted file under "Temp-Reformatted_output.csv"
7. Copy the reformatted output into CocirculationData.csv and save
9 Add the countries of interest as options in the import.country.dat function in the 0func-country_data_import.R script.


.csv files:
-----------------------
COUNTRY_cocirculation.csv contains raw data on which subtypes circulated seasonally in which years, and was dowlonaded from WHO Flu Net for the countries and years of interest. See the workflow above to add new countries.

summary.csv contains a reformatted version of the raw data in .*cocirculation.csv. At one point, I was creating these reformatted files and then copying the relevant parts into a master .csv file, but now I have a script to do this, so this file is technically defunct.

CocirculationData.csv is the master file that contains a summary of all virological data from all countries and all years up to 2017. This file is read by Clean_CocirculationImport_2017.R




0func:
-------------------
0func-import-formatted-country-coriculation-data.R - input a country name, and this function will read in and format virological surveillance data from the countries and year of interest. Output is a table of the number of cases of each flutype reported in a given year, within the country of interest. NAs are added when country-specific data is not available. This function is called by 00-Master-Cocirculation-Import.



Notes on .R files:
-------------------
Clean_CocirculationImport_2017.R    Imports virological surveillance data from CocirculationData.csv and adds this data to a standard cocirculation template. Sources functions that generate country-specific tables of H1N1, H2N2 and H3N2 dominance in a certain year.
- Pre-1977 we assume only one subtype circulated each year.
- From 1977-~1997 we use summary data from the United States
- From 1997 on, we download country-specific data when possible, or substitute summary data from the surrounding region.

gete_ij_2017.R    Defines a function that estimates the probabiltiy of first exposure in year x, x+1, x+2, ... x+12 given birth in year x.

get.country.cocirculation.data_2017.R

Infection_age_structure_2017.R

Reformat_WHO_to_master.R    Loads a .csv using the raw WHO format (eg .*cocirculation.csv), and reformats the data to match the format of CocirculationData.csv


Notes on outputs (.RData):
--------------------------
old_Three_hit_weights.RData - an outdated version of the Three Hit Weights outputs from the ThreeHitModel/ directory. I had made a mistake in the way I was dealing with naive children originally. I'm keeping this on file only for debugging/work checking purposes.

