#### Test myclim ####

folder_2<-file.path("C:","Program Files (x86)","Lolly","data","data_visit_2022_2")

data_file<-grep("data_9",list.files(folder_2,full.names = T),value=T)

library(myClim)
library(data.table)
library(sf)
library(raster)
meta_data
tms_data<-mc_read_files(data_file,dataformat_name = "TOMST")

local_table<- data.table(locality_id = coords$site_ID,elevation = coords$elevation,lon_wgs84= coords$x ,lat_wgs84=coords$y ,tz_offset= 120)

tms_data<-mc_read_data(files_table = data.table(path= data_file,locality_id = coords$site_ID[1:54], data_format ="TOMST" ),localities_table = local_table[1:54,])




tms.info <- mc_info_clean(tms_data) 





coords<-fread("Raw_data/loggers_microclim/Logger_gps_data_trimble.csv",dec=",")
coords[,sensor_ID:=NULL]
coords_sf<-st_as_sf(coords,coords =c("x","y") ,crs=st_crs(4326))
coords_sf_laea<-st_transform(coords_sf,st_crs(3035))
coords_sf<-st_transform(coords_sf,st_crs(2154))

coords[,xl93:=st_coordinates(coords_sf)[,1]]
coords[,yl93:=st_coordinates(coords_sf)[,2]]

coords$elevation<- extract(stack_all_variable_vallee$mnt_25_vosges,coords_sf)

infos<-mc_info(tms_data)

calibration_curves_T3<-fread("Raw_data/loggers_microclim/Calibration_T3.csv")
infos<-merge(infos,calibration_curves_T3,by.x="serial_number",by.y="sensor_ID",all.x=T)


prep_for_calib<-data.frame(serial_number = infos$serial_number,
           sensor_id = infos$sensor_id,
           datetime = as.POSIXct("2021-05-15",tz="UTC"),
           cor_factor = infos$offset,
           cor_slope = infos$slope-1)

tms.load <- mc_prep_calib_load(tms_data, prep_for_calib)

tms_data <- mc_prep_calib(tms.load, sensors = c("TM_T",
                                           "TMS_T1",
                                           "TMS_T2",
                                           "TMS_T3"))
mc_info_count(tms_data) #which returns the number of localities, loggers and sensors in myClim object
mc_info(tms_data)# returning data frame with summary per sensor
mc_info_meta(tms_data)# returning the data frame with locality metadata
mc_info_clean(tms_data) #returning the data frame with cleaning log

class(tms_data)

tms.plot <- mc_filter(tms_data, localities = c("L_S_1","H_N_2","cad_H_5","H_N_8"))
mc_plot_line(tms.plot,
             sensors = c("TMS_T3", "TMS_TMSmoisture"))

mc_plot_raster(tms.plot, sensors = c("TMS_T3"))

tms_day <- mc_agg(tms_data, fun = c("mean", "range", "coverage", "percentile"),
                  percentiles = 95, period = "day", min_coverage = 0.95)



tms.calc <- mc_calc_vwc(tms_data, soiltype = "loamy sand A")



