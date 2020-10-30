run_tag  <- 'RSV_gavi72_basecase'  # 72 Gavi countries (basecase) 
# random number generater seed
rng_seed <- gsub('-','',Sys.Date()) # 20190118
# option to create geographical and country-specific plots # note: this might require substantial processing time
boolean_country_plots <- FALSE; boolean_map_plots <- FALSE # TRUE
# output directory postfix
output_dir_postfix <- paste0(run_tag,'_n',num_sim)
# add timestap to output directory name
output_dir <- paste0('output/',format(Sys.time(),'%m%d%H%M%S_'),output_dir_postfix)
# set seed
set.seed(rng_seed)
# config filename
config_filename <- paste0('./config/',run_tag,'.csv'); time_stamp_main <- Sys.time()
# always clear temporary results
cli_print('Clear all temporary output'); unlink(file.path(get_temp_output_folder(output_dir)),recursive = T)

sim_config_matrix <- read.table(config_filename,sep=',', dec='.',stringsAsFactors = F,header = T)
# set output file name prefix
sim_output_filename  <- file.path(output_dir,run_tag)
# add simulation details
sim_config_matrix$num_sim <- num_sim; sim_config_matrix$scenario_id <- 1:nrow(sim_config_matrix)
sim_config_matrix$rng_seed <- rng_seed; sim_config_matrix$outputFileDir <- get_output_folder(output_dir)
