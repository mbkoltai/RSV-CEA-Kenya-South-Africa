#############################################################################
# This file is part of the RSV modelling project.
# 
# => CREATE GEOSPATIAL FIGURES
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

################################################################
# Function header
################################################################

plot_maps <- function(sim_output_filename)
{

  ##################################################################################
  ## LOAD SIMULATION DATA
  ##################################################################################
  
  load(paste0(sim_output_filename,'.RData'))
  
  # get plot output directory 
  plot_output_folder <- get_plot_output_folder(sim_output,'maps')
  
  # fix for severity-specific effiacy results
  sim_output$efficacy_maternal[is.na(sim_output$efficacy_maternal)] <- 0
  
  # aggregate: get mean per country/scenario
  sim_data <- aggregate(. ~ config_tag + scenario + intervention + country_iso + outputFileDir, data=sim_output, mean, na.omit = NULL)
  
  # create additional output measure
  sim_data$total_DALY_lifeyear_disc_averted <- sim_data$total_DALY_disc_averted / sim_data$population * 1000
  sim_data$total_DALY_lifeyear_disc         <- sim_data$total_DALY_disc / sim_data$population * 1000
  
  ##################################################################################
  # Load map and ggplot settings
  ################################################################
  
  geo_map_dir <- './input/natural_earth_data'
  worldmap <- readOGR(dsn = file.path(geo_map_dir,'ne_50m_admin_0_countries'), layer = "ne_50m_admin_0_countries")
  #worldmap <- worldmap[worldmap$region_un != "Antarctica",]
  
  worldmap@data$country_iso <- as.character(worldmap@data$adm0_a3)
  worldmap@data$country_iso[worldmap@data$country_iso=="SDS"] = "SSD"
  
  # this is just a container to provide an outline for the world.
  bbox    <- readOGR(file.path(geo_map_dir,"ne_110m_graticules_all"), layer="ne_110m_wgs84_bounding_box") 
  bbox_df <- spTransform(bbox, CRS("+proj=robin"))  # reproject bounding box
  bbox_df <- fortify(bbox_df)
  
  themeopt <- theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                    axis.title.x = element_blank(), axis.title.y = element_blank(),
                    axis.ticks = element_blank(), # plot.title = element_text(size=12),
                    panel.background = element_rect(fill = NA)) 
  
  
  
  ##################################################################################
  ## PLOT SCENARIO DATA
  ##################################################################################
  
  # list scenario names to plot results on a map
  scenario_opt <- unique(sim_data$config_tag)

  # list output measures to plot with legend title: [column name, legend tag]
  plot_opt <- rbind(c('total_DALY_lifeyear_disc_averted','DALY averted\nper 1000PY (disc.)'),
                    c('total_DALY_lifeyear_disc','DALY per\n1000PY (disc.)'))
  
   # loop over the scenarios and create plots
  tag_scen <- scenario_opt[1]
  for(tag_scen in scenario_opt)
  {
    
    # select one scenario
    flag          <- sim_data$config_tag == tag_scen
    scenario_data <- sim_data[flag,]
    
    # add scenario data to worldmap variable
    scenario_worldmap      <- worldmap
    scenario_worldmap@data <- left_join(scenario_worldmap@data, scenario_data, by=c("country_iso"))
    
    # transform coordinates and reshape matrix
    scenario_worldmap    <- spTransform(scenario_worldmap, CRS("+proj=robin"))
    ggplotmap_f          <- fortify(scenario_worldmap)
    scenario_worldmap$id <- row.names(scenario_worldmap)
    ggplotmap_f          <- left_join(ggplotmap_f, scenario_worldmap@data,by=c("id"))
    
    # make plots
    for(i_p in 1:nrow(plot_opt))
    {
      # fix pdf name
      pdf_name    <- paste0(plot_output_folder,'map_',tag_scen,'_',plot_opt[i_p,1],'.pdf')
      
      # fix the color range per outcome over all interventions
      color_range <- range(sim_data[,plot_opt[i_p,1]])
      
      # fix the color scheme: yellow/red for burden and blue for averted burden
      plot_color_scale <- scale_fill_gradient(low="yellow", high="red", na.value="gray80",limits=color_range)
      if(grepl('averted',pdf_name)){
        plot_color_scale <- scale_fill_gradient(low="skyblue", high="navyblue", na.value="gray80",limits=color_range)
      }
      
      # create pdf and ggplot
      pdf(pdf_name,12,5)
      print(ggplot(data=bbox_df, aes(long,lat, group=group)) + 
            geom_path(colour="gray70") + 
            geom_polygon(fill="gray90", linetype="solid") + 
            geom_polygon(data=ggplotmap_f, aes_string('long', 'lat', group = 'group', fill = plot_opt[i_p,1])) + 
            coord_equal() + 
            themeopt +
            plot_color_scale +
            #scale_fill_gradient(low="skyblue", high="navyblue", na.value="gray80",limits=color_range) + 
            labs(fill = plot_opt[i_p,2]) ) 
        dev.off()
        print(pdf_name)
      
    } # end plotting for-loop
  } # end scenario for-loop
} # end "plot_maps" function
