This repository contains files that can be used to reconstruct birth year-specific impriting to influenza A subtypes.
Current data allows reconstructions up to 2017.


The methods were originally described in: Gostic et al., Science, 2016. https://science.sciencemag.org/content/354/6313/722
Please cite as appropriate.

SUMMARY OF OVERALL WORKFLOW:
----------------------------
Import country-specific virological surveillance data from WHO Flu Net (details below) (00-Reformat_WHO_to_master.R)

Define functions that read raw virological surveillance data into formatted tables that describe the fraction of seasonal circulation caused by H1N1, H3N2 or H2N2 in any year from 1918-present, in the country of interest. (0func-country_data_import.R)


Define a function to calculate the probability of a first influenza exposure in year j, given birth in year i. Then, define a wrapper function that calculates the probability of subtype-specific imprinting for all years and countries of interest. The wrapper combines virological surveillance data from above with repeated calculations of the probabilities of imprinting for each birth year. (01-reconstruct-imprinting-histories.R)

Run code to output reconstructions for the countries and observation years of interest. (02-generate_reconstructions.R)



Workflow for adding virological surveillance data for new countries or years:
----------------------------------
This code can build reconstructions for many Asian and Euorpean countries for observation years up to 2017. This workflow is only necessary if you desire to reconstruct patterns for countries or observation years not yet supported. See import.country.dat function in the 0func-country_data_import.R script for a list of currently supported countries.

1. Go to WHO Flu Net, http://apps.who.int/flumart/Default?ReportNo=12
2. Display report for the country(ies) of interest starting at week 1 of the earliest possible year and continuing to the final week of the last year of interest
3. Click on the disk icon at the top middle of the page to save the output as a .csv. The file will automatically save to Downloads.
4. Move the file to this directory, using the filename raw-flunet-outputs/COUNTRY_cocirculation.csv
5. Open Reformat_WHO_to_master.R
6. Import the new file from step 4 into the script and run the script to reformat. The script will save a reformatted file under "Temp-Reformatted_output.csv"
7. Append the reformatted output into CocirculationData.csv and save. The new, reformatted data is now in the master .csv file.
9 Add the countries of interest as options in the import.country.dat function in the 0func-country_data_import.R script.

If extending to observation years beyond 2017, you will need to update the max_observation_year option in ##.



.csv files:
-----------------------
These files contain data on either the fraction of circulation caused by H1N1, H2N2 or H3N2 in a given year, or on the overall intensity of influenza A circulation.

Ideally, we use country-specific data to characterize the fraction of influenza A circulation caused by H1N1 or H3N2 in a given season. However, not all countries report a sufficient number of influenza A cases to WHO flu net in each year to support a reaonable estimate. When <50 samples are reported in a given country-year, we fall back on alternate data. As a first choice, we fall back on data from the region of interest. As a second choice (the only option for 1977-1996), we fall back on data from the US.

Asia_fallback.csv, Euro_fallback.csv -  These datasets are the first-choice fallback, and represent averages of reported data from specific regions.

USA_fallback.csv - In rows representing 1977-1996, no data is available on WHO flu net from any country, and so instead we input fractions derived from Table 1 in Thompson et al., JAMA, 2003 (https://jamanetwork.com/journals/jama/fullarticle/195750). Data in rows 1997-present is from WHO FLU NET (US-specific).

Cocirculation_template.csv - From 1918-1976, only one subtype dominated circulation in any given year. This spreadsheet provides a set of indicator variables that characterize circulation in each of those years.

CocirculationData.csv - This is a master data set containing country- and year-specific data on the fraction of influenza A circulation caused by H1N1 or H3N2. Data are obtained from WHO flu net. See above for details.

Intensitymatser.csv - This data set describes the relative annual intensity of influenza circulation from 1918-present. Methods are described in Gostic et al., 2016.




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

