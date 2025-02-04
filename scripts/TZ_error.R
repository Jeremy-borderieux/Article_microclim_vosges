#### Error setting timezone offsets

library(myClim)
## preparing the locality table (aka metadata) for myclim, and load the path in the correct order
folder<-file.path("data","microclimate_loggers")

data_file<-grep("data_9",list.files(folder,full.names = T),value=T)
local_table<- data.table(locality_id = coords$site_ID,
                         elevation = coords$mnt_25_vosges,
                         lon_wgs84= coords$x ,
                         lat_wgs84=coords$y ,
                         tz_offset= -4*60)
local_table_total<-merge(local_table,meta_data[,c("site_ID","sensor_ID")],by.x="locality_id",by.y="site_ID")
local_table_total[,path_to_file:=file.path(folder,paste0("data_",sensor_ID,"_0.csv"))]
local_table_working_logs<-local_table_total[path_to_file%in%data_file,]

files_table = data.table(path= local_table_working_logs$path_to_file,
                         locality_id = local_table_working_logs$locality_id, 
                         data_format ="TOMST" )



##reading with myclim
tms_data<-mc_read_data(files_table = files_table)

tms_data <- mc_prep_meta_locality(tms_data, values = local_table_working_logs )

export = mc_reshape_long(tms_data)

head(export)## the first recording in utc is 2020-11-25 13:30:00
## however, since I set TZ = -4*60, it should be 2020-11-25 9:30:00
