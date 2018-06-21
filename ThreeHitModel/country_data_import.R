## Input - 
#         Country.out - string. Valid inputs are listed in the function below.
#         region.in - string vector of the same length as 'countries in'. In cases where there is not enough data on which subtypes circulated in the country of interest in particular year, we are forced to substitute data from the surrounding region. This option tells the function which "fallback" data to draw from. Options are: 'default' - pulls data from the USA. 'Asia' (aggregate data from Asian countries) or 'Euro' (aggregate data from European countries). If no input is given, the function will automatically use data from the USA for all reconstructions.


## Output - Matrix of values showing the fraction of influenza viruses of different types that circulated in different years. 
##        - Note, output is identical to the format of flu_circ_china

#THIS REFORMATS COUNTRY DATA
import.country.dat = function(Country.out, region.in = 'default'){
  ######## 1 - IMPORT DATA ON WHICH SUBTYPES CIRCULATED IN PAST YEARS IN SPECIFIC COUNTRIES.  
  ############---------------------------------------------------------------------------  
  
  # Sourced script imports a .csv that holds data on which subtypes (H1N1, H3N2 and influenza B) were observed in specific countries, and in specific years.
  # This script also sources two functions, get.cocirculation.ref and get.country.data, which we will use below to format a matrix that tells what fraction of circulation was driven by H1N1 vs. H3N2 in a given season
  source('~/Dropbox/R/Reconstructions/ThreeHitModel/Clean_CocirculationImport_2017.R')
    
    # Extract country-sepcific data on which subtypes circulated in past years
    if(Country.out == 'China'){cocirculation.dat = get.cocirculation.ref('China', region.in)}
    if(Country.out == 'Egypt'){cocirculation.dat = get.cocirculation.ref('Egypt', region.in)}
    if(Country.out == 'Indonesia'){cocirculation.dat = get.cocirculation.ref('Indonesia', region.in)}
    if(Country.out == 'Cambodia'){cocirculation.dat = get.cocirculation.ref('Cambodia', region.in)}
    if(Country.out == 'Vietnam'){cocirculation.dat = get.cocirculation.ref('Vietnam', region.in)}
    if(Country.out == 'Thailand'){cocirculation.dat = get.cocirculation.ref('Thailand', region.in)}
    if(Country.out == 'USA'){cocirculation.dat = get.cocirculation.ref('USA', region.in)}
    if(Country.out == 'UK'){cocirculation.dat = get.cocirculation.ref('UK', region.in)}
    if(Country.out == 'Turkey'){cocirculation.dat = get.cocirculation.ref('Turkey', region.in)}
    if(Country.out == 'Iraq'){cocirculation.dat = get.cocirculation.ref('Iraq', region.in)}
    if(Country.out == 'Azerbaijan'){cocirculation.dat = get.cocirculation.ref('Azerbaijan', region.in)}
    if(Country.out == 'Pakistan'){cocirculation.dat = get.cocirculation.ref('Pakistan', region.in)}
    if(Country.out == 'Nigeria'){cocirculation.dat = get.cocirculation.ref('Nigeria', region.in)}
    if(Country.out == 'Austria'){cocirculation.dat = get.cocirculation.ref('Austria', region.in)}
    if(Country.out == 'Belguim'){cocirculation.dat = get.cocirculation.ref('Belgium', region.in)}
    if(Country.out == 'Denmark'){cocirculation.dat = get.cocirculation.ref('Denmark', region.in)}
    if(Country.out == 'Estonia'){cocirculation.dat = get.cocirculation.ref('Estonia', region.in)}
    if(Country.out == 'Greece'){cocirculation.dat = get.cocirculation.ref('Greece', region.in)}
    if(Country.out == 'Norway'){cocirculation.dat = get.cocirculation.ref('Norway', region.in)}
    if(Country.out == 'Australia'){cocirculation.dat = get.cocirculation.ref('Australia', region.in)}
    if(Country.out == 'Poland'){cocirculation.dat = get.cocirculation.ref('Poland', region.in)}
    if(Country.out == 'Portugal'){cocirculation.dat = get.cocirculation.ref('Portugal', region.in)}
    if(Country.out == 'Spain'){cocirculation.dat = get.cocirculation.ref('Spain', region.in)}
    if(Country.out == 'Argentina'){cocirculation.dat = get.cocirculation.ref('Argentina', region.in)}
    if(Country.out == 'Peru'){cocirculation.dat = get.cocirculation.ref('Peru', region.in)}
    if(Country.out == 'Chile'){cocirculation.dat = get.cocirculation.ref('Chile', region.in)}
    if(Country.out == 'Germany'){cocirculation.dat = get.cocirculation.ref('Germany', region.in)}
    if(Country.out == 'Japan'){cocirculation.dat = get.cocirculation.ref('Japan', region.in)}
    if(Country.out == 'Singapore'){cocirculation.dat = get.cocirculation.ref('Singapore', region.in)}
  
  cocirculation.dat
}
### End function ###




## Examples
#import.country.dat('Vietnam', region.in = 'Asia') ??
#import.country.dat('USA')
    
    