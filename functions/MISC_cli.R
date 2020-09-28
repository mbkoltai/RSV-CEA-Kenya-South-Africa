#############################################################################
# This file is part of the RSV modelling project.
# 
# => MISCELLANEOUS COMMAND LINE INTERFACE FUNCTIONS
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

#' @title Print message to command line interface
#'
#' @description Command line interface: print message
#'
#' @param ... (parts of the) message to print
#' @param WARNING boolean, to print the message in red
#'
#' @keywords internal
cli_print <- function(..., WARNING=F, FORCED=F) {
  
  # get function arguments
  function_arguments <- as.list(match.call(expand.dots=FALSE))$...

  # get function-call environment (to retrieve variable from that environment)
  pf <- parent.frame()
  
  #parse list => make character vector
  f_out <- ' '
for(i in 1:length(function_arguments)){ 
  f_out <- cbind(f_out,eval(unlist(function_arguments[[i]]),envir = pf)) 
}
  
  # add a space to each function arguments
  function_arguments <- paste(f_out,collapse = ' ')
  
  # set text color: black (default) or red (warning)
  web_color_black <- '\033[0;30m'
  web_color_red   <- '\033[0;31m'
  text_color      <- ifelse(WARNING,web_color_red,web_color_black)
  
  # print time + arguments (without spaces)
  cli_out <- paste0(c('echo "',text_color, '[',format(Sys.time(),'%H:%M:%S'),']',
                      function_arguments, web_color_black,'"'),collapse = '')
  
  # print if function is called by master-node or first slave
  if(!exists('par_nodes_info') || 
     Sys.getpid() == par_nodes_info$pid_master ||
     FORCED){
      system(cli_out)
  }  
 
  
  # add to R warnings
  if(WARNING){
    cli_warning <- paste0(c(text_color, '[',format(Sys.time(),'%H:%M:%S'),']',
                            function_arguments, web_color_black),collapse = '')
    warning(cli_warning,
            call. = FALSE, immediate.=FALSE)
  }
}

#########################################
## PROGRESS BAR
##########################################
cli_progress <- function(i_current,i_total,time_stamp_loop){
  
  # print if function is called by first node
  if(exists('par_nodes_info') && 
     Sys.getpid() == par_nodes_info$pid_slave1){
    
    # calculate progress
    progress_scen        <- floor(i_current/i_total*100)
    progress_time        <- round(difftime(Sys.time(),time_stamp_loop,units = "min"),digits=1)
    
    # estimate remaining time (after 15%)
    time_label <- ''
    if(progress_scen > 15 & progress_scen < 99) {
      estim_time <-  round(progress_time / progress_scen * (100-progress_scen),digits=1)
      if(estim_time<1) {estim_time <- '<1'}
      time_label <- paste0('[',estim_time,' min remaining]')
    }
    
    cli_print('RUNNING...',i_current,'/',i_total,time_label,FORCED=TRUE)
  }  
}
