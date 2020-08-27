#############################################################################
# This file is part of the RSV modelling project.
# 
# => MISCELLANEOUS HELP FUNCTIONS FOR THE PARALLEL WORKERS
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

################################################################
## START CLUSTER WITH PARALLEL WORKERS
################################################################
start_parallel_workers <- function()
{
  cli_print("START PARALLEL WORKERS")
  
  ## SETUP PARALLEL NODES
  # note: they will be removed after 600 seconds inactivity
  num_proc      <- detectCores()
  par_cluster   <- makeCluster(num_proc, cores=num_proc, timeout = 600) 
  registerDoParallel(par_cluster)
  
  # store the process id (pid) of the first slave
  pid_slave1 <- clusterEvalQ(par_cluster, { Sys.getpid() })[[1]]
  
  # CREATE GLOBAL VARIABLE
  par_nodes_info <<- list(par_cluster = par_cluster,
                          pid_master  = Sys.getpid(),
                          pid_slave1  = pid_slave1)
  
}

################################################################
## STOP CLUSTER WITH PARALLEL WORKERS
################################################################
stop_parallel_workers <- function()
{
  ## CLOSE NODES AND NODE INFO
  if(exists('par_nodes_info')){
    cli_print("STOP PARALLEL WORKERS")
    
    stopCluster(par_nodes_info$par_cluster); 
    rm(par_nodes_info,envir = .GlobalEnv) # REMOVE GLOBAL VARIABLE
  }
}

################################################################
## CHECK IF CLUSTER EXISTS AND START ONE IF NOT
################################################################
check_parallel_workers <- function(){
  if(!exists('par_nodes_info')){
    start_parallel_workers()
  } else if (!any(grepl(par_nodes_info$pid_slave1,system('ps -A',intern = T)))){
    start_parallel_workers()
  }
}

################################################################
## RESET CLUSTER
################################################################
# reset parallel workers
reset_parallel_workers <- function(){
  stop_parallel_workers()
  gc()
  start_parallel_workers()  
}
