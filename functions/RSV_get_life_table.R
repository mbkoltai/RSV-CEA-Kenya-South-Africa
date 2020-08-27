#############################################################################
# This file is part of the RSV modelling project.
# 
# => GET LIFE TABLE
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

# Get the life table with monthly life expectancy, based on mortality rates from the UN 
# f_country   country name
# f_year      number => is converted into string with 5 year interval, such as 2015-2020
# force_new   generate a new life table, irrespectively if there is a local copy (default = false)
get_life_table <- function(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
{
  
  life_table_filename   <- get_life_table_filename(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
  if(!file.exists(life_table_filename)) {
    generate_life_table(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
  } 
  
  load(life_table_filename)
  return(life_table) #end

}

get_life_table_filename <- function(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly){
 
  period <- get_year_category(f_year)[1]
  temp_output_folder   <- get_temp_output_folder(f_outputFileDir)
  file_name_base       <- file.path(temp_output_folder,paste0('life_table_',f_country_iso,'_',period,'_disc',gsub("\\.",'p',f_disc_rate_qaly)))
  life_table_filename <- paste0(file_name_base,'.RData')
  
  return(life_table_filename)
}

generate_life_table <- function(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
{

  # get period
  period <- get_year_category(f_year)[1]
  
  # filename
  life_table_filename   <- get_life_table_filename(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
  
  if(!file.exists(life_table_filename))
  {
    # convert the country iso3 code into the country name (corresponding the wpp2017 package)
    f_country <- get_UN_country_name(f_country_iso,f_outputFileDir)
    
    # load population data
    load_wpp2017_databases(f_outputFileDir)
  
    if(!f_country %in% sexRatio$name){
      cli_print("Unknown country: ", f_country)
      cli_print("Make sure that the country name starts with a capital letter")
      cli_print("All country names are provided by the function 'get_country_list()'")
      return(NA)
    } else {
      cli_print("Calculating the life table for:", f_country_iso,period,'disc',f_disc_rate_qaly)
    }
    
    # get column for year interval
    sel_col     <- names(mxM) == sub("_","-",period)
    
    # age: A character string representing an age interval (given by the starting age of the interval).
    # note: Data for ages 85-100+ are not the official UN data. While the published UN mortality datasets 
    # contain data only up to 85+, data for ages 85-100+ in this dataset were derived from UN published 
    # life table quantities.
    wpp_data <- data.frame(age_group_min = mxM$age[mxM$name%in%f_country],
                           mxM = mxM[mxM$name%in%f_country,sel_col],
                           mxF = mxF[mxF$name%in%f_country,sel_col],
                           sexRatio = sexRatio[sexRatio$name%in%f_country,sel_col])
    
    # Adjust ages
    wpp_data$age_group_mean <- wpp_data$age_group_min + diff(c(wpp_data$age_group_min,max(wpp_data$age_group_min+5)))/2
    
    # Account for gender balance
    # Estimates and projections of the sex ratio at birth derived as the number of female divided by the number of male.
    wpp_data$prop_female <- wpp_data$sexRatio / 2
    wpp_data$prop_male   <- 1 - wpp_data$prop_female
    
    # Get population nMx
    wpp_data$nMx <- wpp_data$mxM * wpp_data$prop_male + wpp_data$mxF * wpp_data$prop_female
    
    # Approximation under linear model
    #life_table_year<- approx(wpp_data$age_group_min,wpp_data$nMx/wpp_data$age_group_size,xout=seq(0,100,1),rule = 2:1)
    life_table_year<- approx(wpp_data$age_group_min,wpp_data$nMx,xout=seq(0,100,1),rule = 2:1)
    names(life_table_year) <- c('age','nMx')
    life_table_year<- data.frame(life_table_year) # store as data.frame
  
    # complete life table
    life_table_year<- create_full_life_table(life_table_year)
    
    # linear approximation for lx and life_expectancy
    num_months_year <- 12
    lx_month                   <- approx(life_table_year$age*num_months_year,life_table_year$lx,xout=seq(0,100,1))
    life_expectancy_month      <- approx(life_table_year$age*num_months_year,life_table_year$life_expectancy,xout=seq(0,100,1))
    nMx_month                  <- approx(life_table_year$age*num_months_year,life_table_year$nMx,xout=seq(0,100,1))
    life_table                 <- data.frame(lx_month)
    life_table$life_expectancy <- life_expectancy_month$y
    life_table$nMx             <- nMx_month$y
    names(life_table)          <- c('age','lx','life_expectancy','nMx')
    life_table$lx_rate         <- life_table$lx / life_table$lx[1]
  
    # discounted life expectancy
    life_table$life_expectancy_disc <-  NA
    for(i in 1:length(life_table$life_expectancy)){
      
      tmp_years <- 1:floor(life_table$life_expectancy[i])
      life_expectancy_disc_floor <- sum(1/((1+f_disc_rate_qaly)^tmp_years))
      
      tmp_years <- ceiling(life_table$life_expectancy[i])
      life_expectancy_disc_remaining <- sum(1/((1+f_disc_rate_qaly)^tmp_years)) * (life_table$life_expectancy[i] - floor(life_table$life_expectancy[i]))
      
      life_table$life_expectancy_disc[i] <- life_expectancy_disc_floor + life_expectancy_disc_remaining
    }
    
    if(exists('use_fixed_population') && use_fixed_population){
      life_table$lx       <- life_table$lx[1]
      life_table$lx_rate  <- life_table$lx_rate[1]
      life_table$nMx      <- 0
    }
    
    save(life_table,file=life_table_filename)
    
  }
  
}

get_country_list <- function()
{
  data(sexRatio, package = "wpp2017", envir = environment())
  return(levels(sexRatio$name))
}

get_UN_country_name <- function(country_iso3,f_outputFileDir){
  
  # get UN country name - ISO databse
  UNcountries <- get_UN_country_database(f_outputFileDir)
  
  # look for the given iso3 code in the country database
  flag <- country_iso3 == UNcountries$iso3 
  
  # return the country name
  return(UNcountries$name[flag])
}

get_country_iso3 <- function(country_name,f_outputFileDir){
  
  # get UN country name - ISO databse
  UNcountries <- get_UN_country_database(f_outputFileDir)
  
  # look for the given iso3 code in the country database
  flag <- country_name == UNcountries$name 
  
  # return the country name
  return(UNcountries$iso3[flag])
}

create_UN_country_database <- function(f_outputFileDir){
  
  # filename for the [country name - iso3 code] database
  temp_output_folder              <- get_temp_output_folder(f_outputFileDir)
  dictionary_UN_countrynames_file <- file.path(temp_output_folder,'dictionary_UN_countrynames.RData')
  
  # if this database does not exist => create this using the 'countrycode' library
  # computational expensive procedure! => store local copy
  if(!file.exists(dictionary_UN_countrynames_file)) {
    cli_print('Create UN country database')
    data("UNlocations",package='wpp2017')
    UNcountries      <- data.frame(name=UNlocations$name)
    UNcountries$iso3 <- countrycode(UNcountries$name,'country.name','iso3c',warn=F)
    UNcountries$iso3[is.na(UNcountries$iso3)] <- 0
    
    save(UNcountries,file=dictionary_UN_countrynames_file)
    cli_print('UN country database completed')
  }

}

get_UN_country_database <- function(f_outputFileDir){
  
  # filename for the [country name - iso3 code] database
  temp_output_folder              <- get_temp_output_folder(f_outputFileDir)
  dictionary_UN_countrynames_file <- file.path(temp_output_folder,'dictionary_UN_countrynames.RData')
  
  # if this database does not exist => create this using the 'countrycode' library
  # computational expensive procedure! => store local copy
  if(!file.exists(dictionary_UN_countrynames_file)) {
    create_UN_country_database(f_outputFileDir)
  }
  
  # load the [country name - iso3 code] database
  load(dictionary_UN_countrynames_file)
  
  return(UNcountries)
}

load_UN_country_database <- function(f_outputFileDir){
  
  # filename for the [country name - iso3 code] database
  temp_output_folder              <- get_temp_output_folder(f_outputFileDir)
  dictionary_UN_countrynames_file <- file.path(temp_output_folder,'dictionary_UN_countrynames.RData')
  
  # if this database does not exist => create this using the 'countrycode' library
  # computational expensive procedure! => store local copy
  if(!file.exists(dictionary_UN_countrynames_file)) {
    cli_print('create UN country database')
    data("UNlocations",package='wpp2017')
    UNcountries      <- data.frame(name=UNlocations$name)
    UNcountries$iso3 <- countrycode(UNcountries$name,'country.name','iso3c',warn=F)
    UNcountries$iso3[is.na(UNcountries$iso3)] <- 0
    
    save(UNcountries,file=dictionary_UN_countrynames_file)
    cli_print('UN country database completed')
  }
  
  # load the [country name - iso3 code] database
  load(dictionary_UN_countrynames_file)
  
  return(UNcountries)
}

load_wpp2017_databases <- function(f_outputFileDir)
{
  
  temp_output_folder <- get_temp_output_folder(f_outputFileDir)
  filename_wpp2017_values <- file.path(temp_output_folder,'wpp2017_values.RData')
  if(file.exists(filename_wpp2017_values)){
    load(filename_wpp2017_values)
  } else {
    data(mxM, package = "wpp2017")#, envir = environment())
    data(mxF, package = "wpp2017")#, envir = environment())
    data(sexRatio, package = "wpp2017")#, envir = environment())
    save(mxM,mxF,sexRatio,file=filename_wpp2017_values)
  }
  load(filename_wpp2017_values)
}

## FUNCTION TO COMPLETE LIFE TABLE BASED ON MORTALITY RATE
create_full_life_table <- function(f_life_table)
{

  # Copy life table
  # age - years
  # nMx - age-specific death rate between ages x and x+n
  life_table_out <- f_life_table
  
  # Convert rate to probability
  # nqx - probability of dying between ages x and x+n
  life_table_out$nqx=1-exp(-life_table_out$nMx)
  
  # Calculate:
  # - lx - number of people left alive at age x
  # - ndx - number of people dying between ages x and x+n
  # - nLx - person-years lived between ages x and x+n  
  life_table_out$lx  <- NA # init columns, otherwise we cannot use row indices
  life_table_out$ndx <- NA # init columns, otherwise we cannot use row indices
  life_table_out$nLx <- NA # init columns, otherwise we cannot use row indices
  for(i_age in 1:dim(life_table_out)[1]){
     if(i_age == 1) {
      # lx starts with 100 000 and every year it is reduced by the number of death in ndx
      life_table_out$lx[1]=100000 
    } else {
      life_table_out$lx[i_age]= life_table_out$lx[i_age-1] - life_table_out$ndx[i_age-1]  
    }
    life_table_out$ndx[i_age] = life_table_out$lx[i_age] * life_table_out$nqx[i_age]
    life_table_out$nLx[i_age] = life_table_out$ndx[i_age] / life_table_out$nMx[i_age] 
  }
  
  # calculate Tx - person-years lived above age x
  i_age = 1
  num_ages = dim(life_table_out)[1]
  for(i_age in 1:num_ages){
    life_table_out$Tx[i_age] = sum(life_table_out$nLx[i_age:num_ages])
  }
  
  # Calculate ex - expectation of life at age x  (by Ex=Tx/Lx)
  life_table_out$life_expectancy=life_table_out$Tx / life_table_out$lx
  
  return(life_table_out)
}


get_year_category <- function(f_year)
{
  year_cutoff <- (0:10)*5 + 2000
  year_index  <- max(which(f_year >= year_cutoff))
  period      <- paste(year_cutoff[c(year_index,year_index+1)],collapse='_')
  
  return(c(period,year_cutoff[c(year_index)]))
}

