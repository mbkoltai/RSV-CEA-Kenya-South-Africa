#############################################################################
# This file is part of the RSV modelling project.
# 
# => SCRIPT TO LOAD ALL PACKAGES AND FUNCTIONS OF MCMARCEL
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################


# load RSV-related functions
source('functions/RSV_config.R')
source('functions/RSV_get_burden.R')
source('functions/RSV_get_cost_data.R')
source('functions/RSV_get_incidence.R')
source('functions/RSV_get_life_table.R')
source('functions/RSV_plot_CEA_country.R')
source('functions/RSV_plot_CEAF_table.R')
source('functions/RSV_plot_maps.R')
source('functions/RSV_write_summary_tables.R')

# load help functions
source('functions/MISC_CEA.R')
source('functions/MISC_cli.R')
source('functions/MISC_file_management.R')
source('functions/MISC_parallel_admin.R')


# load all required packages for this project
#
# wpp2017       to obtain UN population data
# countrycode   to convert country names into iso3 codes
# scales        to use transparent colors 
# mgcv          to use the gam function
# rgdal         to plot results on Geospatial scales
# ggplot2       to make fancy plots
# dplyr         to aggregate simulation results (geospatial map)
# doParallel    to run foreach loops in parallel 
#
# ADDITIONAL:
# doParallel    to run foreach loops in parallel

all_packages <- c('wpp2017','countrycode','scales','mgcv','rgdal','ggplot2','dplyr')

# load the doParallel package seperatly so we can use 'all_packages' in a parallel foreach
all_packages_with_parallel <- c(all_packages,'doParallel')

# loop over the packages
for(package_i in all_packages_with_parallel){
  
  # if not present => install
  if(!package_i %in% rownames(installed.packages())){
    install.packages(package_i)
  }
  
  # load package (without messages)
  suppressPackageStartupMessages(
    library(package_i, 
            character.only=TRUE, 
            quietly = TRUE, 
            warn.conflicts = FALSE,
            verbose = FALSE)
  )
}

