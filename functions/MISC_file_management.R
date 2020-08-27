#############################################################################
# This file is part of the RSV modelling project.
# 
# => MISCELLANEOUS HELP FUNCTIONS FOR FILE MANAGEMENT
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

#############################################################################
# GET TEMPORARY OUTPUT FOLDER PATH (create folder if not present yet)
#############################################################################
get_temp_output_folder <- function(f_outputFileDir,f_subfolder=NA)
{
  
  # define "temp" folder, 
  temp_output_folder   <- file.path(f_outputFileDir,'temp')
  
  # option to create subfolder
  if(!is.na(f_subfolder)){
    temp_output_folder   <- file.path(temp_output_folder,f_subfolder)
  }
  
  # if it does not exist: create full path using the recursive option
  if(!file.exists(temp_output_folder)){
    cli_print('Create folder:',temp_output_folder)
    dir.create(temp_output_folder, recursive =  TRUE)
  }
  return(temp_output_folder)
}

#############################################################################
# GET PLOT OUTPUT FOLDER PATH ( create folder if not present yet)
#############################################################################
get_plot_output_folder <- function(sim_output,folder_name='4'){
  
  plot_output_folder <- paste0('./',unique(sim_output$outputFileDir),
                               '/plot_',file.path(folder_name),'/')

  # check if the output folder exist, and create the folder if not
  if(!dir.exists(plot_output_folder)){
    cli_print('Create folder:',plot_output_folder)
    dir.create(plot_output_folder,recursive = T)
  }
  
  return(plot_output_folder)
}

#############################################################################
# GET OUTPUT FOLDER PATH (create folder if not present yet)
#############################################################################
get_output_folder <- function(outputFileDir,folder_name=NA){
  
  if(is.na(folder_name)){
    output_folder <- file.path(outputFileDir)
  } else{
    output_folder <- file.path(outputFileDir,folder_name)  
  }
  
  # check if the output folder exist, and create the folder if not
  if(!dir.exists(output_folder)){
    cli_print('Create folder:',output_folder)
    dir.create(output_folder,recursive = T)
  }
  
  return(output_folder)
}

#############################################################################
# DERIVE THE RUN TAG FROM THE sim_output_filename
#############################################################################
get_run_tag <- function(sim_output_filename){
  
  return(basename(sim_output_filename))
}
