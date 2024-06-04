#### packages ####

##installing the latest myclim version
install.packages("http://labgis.ibot.cas.cz/myclim/myClim_latest.tar.gz", repos=NULL, build_vignettes=TRUE)
library(myClim)
library(data.table)
library(sf)
library(raster)

source("Packages.R")

#### reding the meta data and covariables of the analysis #### 
## reading and preparing the metadata file
meta_data<-fread(file.path("data","metadata_loggers","meta_data_2022_points.csv"),dec = ",")
meta_data<-meta_data[site_ID!="",]
colnames(meta_data)[4]<-"date_setup"
meta_data[,date_setup:=dmy(date_setup)]
meta_data[,date_revisite_1:=ifelse(date_revisite_1=="",NA,date_revisite_1)]
meta_data[,date_revisite_1:=dmy(date_revisite_1)]
meta_data<-data.table::melt(meta_data,measure.vars = c("sensor_ID","sensor_ID_2"),value.name = "sensor_ID")
meta_data<-meta_data[order(sensor_ID),]
meta_data<-meta_data[!is.na(sensor_ID),]


## covariables
stack_all_variable_vallee<-stack(file.path("data","processed_raster",list.files(file.path("data","processed_raster"))))


# get copernicus data

tree_density<-raster(file.path("data","tree_density_raster","tree_cover_density_thur_2018_copernicus.tif"))
tree_density<-projectRaster(tree_density,crs=crs(stack_all_variable_vallee))

tree_density_reprojected<-projectRaster(tree_density,stack_all_variable_vallee)
stack_all_variable_vallee$tree_density_projected<-tree_density_reprojected


#stack_all_variable_vallee<-stack(stack_all_variable_vallee,tree_density_reprojected)

## precise coordinates with gnss trimble
coords<-fread(file.path("data","metadata_loggers","Logger_gps_data_trimble.csv"),dec=",")
coords[,sensor_ID:=NULL]
coords_sf<-st_as_sf(coords,coords =c("x","y") ,crs=st_crs(4326))
coords_sf_laea<-st_transform(coords_sf,st_crs(3035))
coords_sf<-st_transform(coords_sf,st_crs(2154))

coords[,xl93:=st_coordinates(coords_sf)[,1]]
coords[,yl93:=st_coordinates(coords_sf)[,2]]
## extraction of the covariables
env_var<- extract(stack_all_variable_vallee,coords_sf)
tree_density_2018<- extract(tree_density,coords_sf_laea)
env_var<-cbind(env_var,tree_density_2018);rm(tree_density_2018)

coords<-cbind(coords,env_var)
meta_data<-merge(meta_data,coords,by="site_ID")


## data from a meteological statino of grand ventron
Ventron_station<-fread(file.path("data","metadata_loggers","station_ventron_21_22.csv"))
Ventron_station[,time:=ymd_hms(time)]
Ventron_station[,time:=with_tz(time=time,"Europe/Paris")]
Ventron_station[,day:=floor_date(time,"day")]
Ventron_station[,month:=floor_date(time,"month")]
Ventron_station_day<-Ventron_station[,.(temp_ventr_station=mean(temperature),
                                        max_ventron=max(temperature),
                                        min_ventron=min(temperature),
                                        range_ventron=max(temperature)-min(temperature),
                                        n_range_ventron=(max(temperature)-min(temperature))/mean(temperature)),by=day]
Ventron_station_day[,date:=ymd(day)]
Ventron_station_day[,day:=NULL]


## canopy cover in a 25m radius data
cover_25m_sites<-fread(file.path("data","metadata_loggers","cover_25m_sites.csv"))
cover_25m_sites_aggr<-cover_25m_sites[,.(canopy_cover_25m=sum(free_cover),under_canopy_cover=sum(under_cano_cover),broadleaved_prop=(sum(ifelse(Broadleaved,free_cover,0))/sum(free_cover)  )),by=site_ID]
meta_data<-merge(meta_data,cover_25m_sites_aggr,by="site_ID",all.x=T)




#### read, clean, aggregate  microclimate data with myClim ####

folder_2<-file.path("C:","Program Files (x86)","Lolly","data","data_visit_2022_2")
## removing the intensive reads 
intensive_logs<-c("94195464","94195409")
path_intensive<-file.path(folder_2,paste0("data_",c("94195464","94195409"),"_0.csv"))

foreach(i = path_intensive)%do%{
  tmp<-fread(i)
  
  tmp[,V2:=ymd_hm(V2)]
  tmp[,time:=minute(V2)]
  tmp<-tmp[time%in%c(0,15,30,45),]
  tmp[,time:=NULL]
  tmp[,V2:=as.character(V2)]
  tmp[,V2:=substr(V2,1,nchar(V2)-3)]
  tmp[,V2:=str_replace_all(V2,"-",".")]
  write.table(tmp,i,sep=";",row.names = F,col.names = F,quote=F)
  tmp
}


## preparing the locality table (aka metadata) for myclim, and load the path in the correct order
data_file<-grep("data_9",list.files(folder_2,full.names = T),value=T)
local_table<- data.table(locality_id = coords$site_ID,elevation = coords$mnt_25_vosges,lon_wgs84= coords$x ,lat_wgs84=coords$y ,tz_offset= 120)
local_table_total<-merge(local_table,meta_data[,c("site_ID","sensor_ID")],by.x="locality_id",by.y="site_ID")
local_table_total[,path_to_file:=file.path(folder_2,paste0("data_",sensor_ID,"_0.csv"))]

local_table_working_logs<-local_table_total[path_to_file%in%data_file,]

##reading with myclim
tms_data<-mc_read_data(files_table = data.table(path= local_table_working_logs$path_to_file,locality_id = local_table_working_logs$locality_id, data_format ="TOMST" ),localities_table = local_table_working_logs)


infos<-mc_info(tms_data)

## calibration og the t3 sensors with curves obtained from lab calibration
calibration_curves_T3<-fread(file.path("data","metadata_loggers","Calibration_T3.csv"))
infos<-merge(infos,calibration_curves_T3,by.x="serial_number",by.y="sensor_ID",all.x=T)

prep_for_calib<-data.frame(serial_number = infos$serial_number,
                           sensor_id = infos$sensor_id,
                           datetime = as.POSIXct("2021-05-15",tz="UTC"),
                           cor_factor = infos$offset,
                           cor_slope = infos$slope-1)

tms.load <- mc_prep_calib_load(tms_data, prep_for_calib)
tms_data <- mc_prep_calib(tms.load, sensors = "TMS_T3")

rm(tms.load)
mc_info_count(tms_data) #which returns the number of localities, loggers and sensors in myClim object
mc_info(tms_data)# returning data frame with summary per sensor
mc_info_meta(tms_data)# returning the data frame with locality metadata
mc_info_clean(tms_data) #returning the data frame with cleaning log


tms.calc <- mc_calc_vwc(tms_data, soiltype = "universal")
## virtual sensor to estimate snow presence from 2 cm air temperature
tms.calc <- mc_calc_snow(tms.calc, sensor = "TMS_T2")
# vpd## virtual sensor with VPD
tms.calc <- mc_calc_vpd(tms.calc)


# min max mean temp ok, 
tms_error_flag <- mc_agg(tms.calc, fun = list(TMS_T1 = "error_flag"), period = "day",
                         custom_functions = list(error_flag = function(x){mean(na.rm=T,abs(x-c(NA,x[1:(length(x-1))])))}))

tms_error_flag <- data.table(mc_reshape_long(tms_error_flag))
# errflag ok


## summary data.frame of snow estimation
tms.snow <- mc_calc_snow_agg(tms.calc,period = 3)
tms.snow <-mc_agg(tms.calc,fun=list(snow="mean"),period = "day",min_coverage = 0.95)
tms.snow <- data.table(mc_reshape_long(tms.snow))

temp_env <- mc_agg(tms.calc, 
                   fun=list(TMS_T3=c("mean_t","percentile","range_t")),
                   period = "day", min_coverage = 0.95,percentiles = c(5,95),
                   custom_functions = list(mean_t=function(x) mean(x[x%between%  quantile(x,na.rm=T,probs=c(0.05,0.95))],na.rm=T) ,
                                           range_t=function(x) diff(range(x[x%between%  quantile(x,na.rm=T,probs=c(0.05,0.95))],na.rm=T) )))


temp_env_export<- data.table( mc_reshape_long(temp_env))

temp_env <- mc_env_temp(tms.calc, period = "day", min_coverage = 0.95)
temp_env<- data.table(temp_env)
temp_env<-temp_env[sensor_name%in%c("T.soil_8_cm.mean","T.air_15_cm.min5p","T.air_15_cm.mean","T.air_15_cm.max95p","T.air_15_cm.drange"),]

vpd_env <- mc_env_vpd(tms.calc, period = "day", min_coverage = 0.9)


## aggregated data= 
tms.calc<-mc_calc_gdd(tms.calc, sensor = "TMS_T3")
tms.calc<-mc_calc_fdd(tms.calc, sensor = "TMS_T3")

agg_growing_season<-mc_agg(tms.calc,
       fun=list(GDD5="sum"),
       period="custom",
       custom_start ="04-01" ,
       custom_end = "08-15")



agg_winter<-mc_agg(tms.calc,
                   fun=list(FDD0="sum",snow ="sum"),
                   period="custom",
                   custom_start ="10-01" ,
                   custom_end = "03-31")


all_tms_data_agg<-data.table(rbind(mc_reshape_long(agg_growing_season),mc_reshape_long(agg_winter)))
all_tms_data_agg[,height:=NULL]
all_tms_data_agg[,serial_number:=NULL]
colnames(all_tms_data_agg)[3]<-"date"
colnames(all_tms_data_agg)[4]<-"date_to"

all_tms_data_agg[,date:=ymd(date)]
all_tms_data_agg[,date_to:=ymd(date_to)]

all_tms_data_agg<-rbind(all_tms_data_agg[date==ymd("2022-04-01") & sensor_name=="GDD5_sum",],
                        all_tms_data_agg[date==ymd("2021-10-01") & sensor_name%in%c("FDD0_sum","snow_sum"),])


all_tms_data_agg[,date:=NULL]
all_tms_data_agg[,date_to:=NULL]
all_tms_data_agg[sensor_name=="snow_sum",value:=value/(4*24)]

all_tms_data_agg<-data.table::dcast(all_tms_data_agg, locality_id~sensor_name)

all_tms_data_agg<-merge(all_tms_data_agg,meta_data[!duplicated(site_ID),],by.x="locality_id",by.y="site_ID")


## temporal stability
tms_error_flag <- mc_agg(tms.calc, fun = list(TMS_T1 = "error_flag"), period = "day",
                         custom_functions = list(error_flag = function(x){mean(na.rm=T,abs(x-c(NA,x[1:(length(x-1))])))}))

tms_error_flag <- data.table(mc_reshape_long(tms_error_flag))



## export to a data.table
all_tms_data_day<-rbind(temp_env,tms.snow,tms_error_flag)
all_tms_data_day[,serial_number:=NULL]
all_tms_data_day<-data.table::dcast(all_tms_data_day, locality_id+datetime~sensor_name)
all_tms_data_day[,datetime:=ymd(datetime)]
colnames(all_tms_data_day)[2]<-"date"
colnames(all_tms_data_day)[3:7]<-c("T3_range","T3_max","T3_mean","T3_min","Tsoil_mean")




                 
all_tms_data_day<-rbind(temp_env_export,tms_error_flag)
all_tms_data_day[,serial_number:=NULL]
all_tms_data_day<-data.table::dcast(all_tms_data_day, locality_id+datetime~sensor_name)

colnames(all_tms_data_day)[2:7]<-c("date","TMS_T1_error_flag","T3_mean","T3_min","T3_max","T3_range")


all_tms_data_day[,n_T3_range:=T3_range/T3_mean]



unique(all_tms_data_day$sensor_name)


all_tms_data_day<-merge(all_tms_data_day,meta_data[!duplicated(site_ID),],by.x="locality_id",by.y="site_ID")
all_tms_data_day[,date:=ymd(date)]
all_tms_data_day<-merge(all_tms_data_day,Ventron_station_day,by="date",all.x=T)
all_tms_data_day<-all_tms_data_day[order(locality_id)]

### correcte the station temperature by elevation


#### out of soil sensor detection and management ####


all_tms_data_day[,rolling_mean_err:=frollapply(TMS_T1_error_flag,10,FUN = mean,align = "left"),by=locality_id]
treshold_out<-0.35

meta_data[,put_out_date_1:=NA]

all_tms_data_day[,treshold_crossed:=rolling_mean_err>treshold_out]
all_tms_data_day[,in_soil:=1]
for (sensor in unique(all_tms_data_day[,locality_id])){
  day_put_in<-unique(meta_data[site_ID==sensor,date_setup])
  print(day_put_in)
  day_crossing<-all_tms_data_day[treshold_crossed==T& locality_id==sensor& date>day_put_in,date][1]
  meta_data[,put_out_date_1:=ifelse(site_ID==sensor,as.character(day_crossing),put_out_date_1)]
  print(paste0("Sensor n° ",sensor," was put out of soil on ",day_crossing))
  if(is.na(day_crossing))next
  all_tms_data_day[,in_soil:=ifelse(locality_id==sensor,ifelse(date>=day_crossing,0,1),in_soil)]

  
  
}


meta_data[,put_out_date_1:=ymd(put_out_date_1)]
meta_data[, accessibility:=NULL]
meta_data[, Revisite:=NULL]
meta_data[, comments:=NULL]

meta_data[,back_in_soil_date:=date_revisite_1]


table(all_tms_data_day$in_soil,all_tms_data_day$locality_id)

ggplot(all_tms_data_day,aes(x=datetime,y=date,color=locality_id))+theme_bw()+geom_line()
ggplot(all_tms_data_day[locality_id=="L_N_7",],aes(x=date,y=TMS_T1_error_flag,color=locality_id))+theme_bw()+geom_line()
ggplot(all_tms_data_day[locality_id=="H_S_1",],aes(x=date,y=TMS_T1_error_flag,color=locality_id))+theme_bw()+geom_line()
ggplot(all_tms_data_day[locality_id=="H_N_2",],aes(x=date,y=TMS_T1_error_flag,color=locality_id))+theme_bw()+geom_line()


### two special case : 94207915 and  94207913
### 94207915 was in a building for several days, then it was planted again

plotly::ggplotly( ggplot(all_tms_data_day[locality_id=="H_N_2" ,],aes(x=date,y=TMS_T1_error_flag,color=locality_id))+theme_bw()+geom_line(lwd=1))#+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2))
date_out<-ymd(c("2021/08/27","2022/07/13"))
date_int<-ymd("2022/01/19")
all_tms_data_day[locality_id=="H_N_2"& date>date_out[1] ,in_soil:=0]
all_tms_data_day[locality_id=="H_N_2"& date>date_int[1] ,in_soil:=1]
all_tms_data_day[locality_id=="H_N_2"& date>date_out[2],in_soil:=0]
#plotly::ggplotly( ggplot(all_tms_data_day[locality_id=="H_S_1" ,],aes(x=date,y=TMS_T1_error_flag,color=locality_id))+theme_bw()+geom_line(lwd=1))#+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2))
date_out<-ymd(c("2021/06/03"))
date_int<-ymd("2022/06/04")
all_tms_data_day[locality_id=="H_S_1" ,in_soil:=1]
all_tms_data_day[locality_id=="H_S_1"& date>date_out[1] ,in_soil:=0]
all_tms_data_day[locality_id=="H_S_1"& date>date_int[1],in_soil:=1]


### end special case


#### modelling ####

#### daily scale ####

suspicious_loggers<-c("L_S_1","L_S_3","slope_H_1","L_S_7")
ggplot(all_tms_data_day[locality_id%in%suspicious_loggers,],aes(x=date,y=TMS_T1_error_flag,color=locality_id))+theme_bw()+geom_line()
ggplot(all_tms_data_day[locality_id%in%suspicious_loggers,],aes(x=date,y=T3_mean,color=locality_id))+theme_bw()+geom_line()


ggplot(all_tms_data_day_backup[locality_id%in%"H_N_6",],aes(x=date,y=snow_mean,color=locality_id))+theme_bw()+geom_line()

# 
# snow_data<-all_tms_data_day_backup[! locality_id %in% c(suspicious_loggers,new_loggers,broken_loggers) &
#                                      date%between% c(dmy("01/10/2021"),dmy("01/04/2022")),.(snow=mean(snow_mean,na.rm=T)),by=locality_id]
# 
# snow_data<-merge(snow_data,meta_data[variable=="sensor_ID",],by.x="locality_id",by.y="site_ID")
# snow_data_sf<-st_as_sf(snow_data,coords=c("xl93"  ,  "yl93"),crs=st_crs(2154))
# mapview::mapview(snow_data_sf,zcol="snow")
# 
# 
# ggplot(snow_data,aes(x=mnt_25_vosges,y=snow))+geom_point()+geom_smooth(method="lm")+theme_bw()
# 
# glm_snow<-glm(snow~ ipv_1000_25+Heat_load_index_25+mnt_25_vosges,data=snow_data,family = binomial())
# summary(glm_snow)
# tab_model(glm_snow)
# 
# 
mapview::mapview(coords_sf[coords_sf$site_ID%in% c("H_N_4","H_N_6","H_N_8_2","H_S_6","slope_H_5","slope_L_6"),])
coords_sf
all_tms_data_day_backup[,max(date),by=locality_id][V1==ymd("2022-08-17"),]


suspicious_loggers<-c("L_S_1","L_S_3","slope_H_1","L_S_7")
new_loggers<-c("new1")
broken_loggers<-c("cad_L_8")

all_tms_data_day_backup<-all_tms_data_day
all_tms_data_day<-all_tms_data_day_backup
all_tms_data_day<-all_tms_data_day[! locality_id %in% c(suspicious_loggers,new_loggers,broken_loggers)]

study_period<-c(dmy("01/04/2022"),dmy("21/08/2022"))
#all_tms_data_day<-all_tms_data_day[day%between%study_period,]
all_tms_data_day<-all_tms_data_day[date%between%study_period,]

all_tms_data_day[date %between% study_period,weight:=1/sum(!is.na(T3_mean)),by=locality_id]

## ##

ggplot(all_tms_data_day,aes(x=date,y=T3_mean,color=locality_id))+theme_bw()+geom_line()

plotly::ggplotly(ggplot(all_tms_data_day,aes(x=date,y=T3_mean,color=locality_id))+theme_bw()+geom_line())
ggplot(all_tms_data_day,aes(x=temp_ventr_station,y=T3_mean,color=mnt_25_vosges))+theme_bw()+geom_point(size=0.5,alpha=0.5)

##model
get_laps_rate<-lme4::lmer(T3_mean~mnt_25_vosges+temp_ventr_station+(1|locality_id)+(1|date),
                                      data=all_tms_data_day,
                                      weights = all_tms_data_day$weight)

summary(get_laps_rate)
tab_model(get_laps_rate)
laps_rate<-   -0.0064495  

all_tms_data_day[,temp_ventr_station_corrected:= temp_ventr_station + ( laps_rate*(mnt_25_vosges-1210))]
all_tms_data_day[,canopy_cover_glama:=str_replace_all(canopy_cover_glama,",",".")]
all_tms_data_day[,canopy_cover_glama:=as.numeric(canopy_cover_glama)]


#daily_model<-lme4::lmer(T3_mean~   temp_ventr_station*(mnt_25_vosges+Heat_load_index_25+tree_density_2018+ipv_1000_25)+
daily_model<-lme4::lmer(T3_mean~   temp_ventr_station*(mnt_25_vosges+Heat_load_index_25+tree_density_2018+ipv_1000_25+slope_25_vosges+distance_edge_25_vosges)+
                                                    
                          #daily_model<-lme4::lmer(T3_mean~  temp_ventr_station_corrected*(tree_density_2018+aspect_NS_25_vosges+ipv_1000_25)+
#daily_model<-lme4::lmer(T3_mean~  offset(temp_ventr_station_corrected)+ temp_ventr_station_corrected:(tree_density_2018+aspect_NS_25_vosges+ipv_1000_25)+tree_density_2018+aspect_NS_25_vosges+ipv_1000_25+
                        (1|locality_id)+(1|date),
                        data=all_tms_data_day,
                        weights = all_tms_data_day$weight,
                        )

daily_model_min<-lme4::lmer(T3_min~   min_ventron*(mnt_25_vosges+tree_density_2018+Heat_load_index_25+ipv_1000_25)+
                          #daily_model<-lme4::lmer(T3_mean~  offset(temp_ventr_station_corrected)+ temp_ventr_station_corrected:(tree_density_2018+aspect_NS_25_vosges+ipv_1000_25)+tree_density_2018+aspect_NS_25_vosges+ipv_1000_25+
                          (1|locality_id)+(1|date),
                        data=all_tms_data_day,
                        weights = all_tms_data_day$weight,
)

daily_model_max<-lme4::lmer(T3_max~   max_ventron*(mnt_25_vosges+tree_density_2018+Heat_load_index_25+ipv_1000_25)+
                          #daily_model<-lme4::lmer(T3_mean~  offset(temp_ventr_station_corrected)+ temp_ventr_station_corrected:(tree_density_2018+aspect_NS_25_vosges+ipv_1000_25)+tree_density_2018+aspect_NS_25_vosges+ipv_1000_25+
                          (1|locality_id)+(1|date),
                        data=all_tms_data_day,
                        weights = all_tms_data_day$weight,
)

daily_model_range<-lme4::lmer(T3_range~  range_ventron* (mnt_25_vosges+tree_density_2018+Heat_load_index_25+ipv_1000_25)+
                              #daily_model<-lme4::lmer(T3_mean~  offset(temp_ventr_station_corrected)+ temp_ventr_station_corrected:(tree_density_2018+aspect_NS_25_vosges+ipv_1000_25)+tree_density_2018+aspect_NS_25_vosges+ipv_1000_25+
                              (1|locality_id)+(1|date),
                            data=all_tms_data_day,
                            weights = all_tms_data_day$weight,
)




daily_meso_model<-lme4::lmer(T3_mean~  (temp_ventr_station_corrected) +
                        (1|locality_id)+(1|date),
                        data=all_tms_data_day,
                        weights = all_tms_data_day$weight,
                        )
summary(daily_model)

tab_model(daily_model)
tab_model(daily_model_min)
tab_model(daily_model_max)

tab_model(daily_model_range)

tab_model(daily_meso_model)


plot_model(daily_model,type="pred",terms=c("temp_ventr_station","aspect_NS_25_vosges"))
plot_model(daily_model,type="pred",terms=c("temp_ventr_station","Heat_load_index_25"))
plot_model(daily_model,type="pred",terms=c("Heat_load_index_25","temp_ventr_station"))

plot_model(daily_model,type="pred",terms=c("temp_ventr_station","slope_25_vosges"))
plot_model(daily_model,type="pred",terms=c("distance_edge_25_vosges","temp_ventr_station"))


plot_model(daily_model,type="pred",terms=c("temp_ventr_station","ipv_1000_25"))
plot_model(daily_model,type="pred",terms=c("temp_ventr_station","ipv_1000_25"))
plot_model(daily_model,type="pred",terms=c("temp_ventr_station","canopy_cover_glama"))

plot_model(daily_model,type="pred",terms=c("Heat_load_index_25","temp_ventr_station"))
plot_model(daily_model,type="pred",terms=c("ipv_1000_25","temp_ventr_station"))
plot_model(daily_model,type="pred",terms=c("tree_density_2018","temp_ventr_station"))


plot_model(daily_model_max,type="pred",terms=c("max_ventron","tree_density_2018"))
plot_model(daily_model_max,type="pred",terms=c("max_ventron","aspect_NS_25_vosges"))
plot_model(daily_model_max,type="pred",terms=c("max_ventron","Heat_load_index_25"))

plot_model(daily_model_min,type="pred",terms=c("min_ventron","ipv_1000_25"))
plot_model(daily_model_min,type="pred",terms=c("min_ventron","aspect_NS_25_vosges"))
plot_model(daily_model_min,type="pred",terms=c("min_ventron","Heat_load_index_25"))

plot_model(daily_model_range,type="pred",terms=c("range_ventron","ipv_1000_25"))
plot_model(daily_model_range,type="pred",terms=c("range_ventron","Heat_load_index_25"))
plot_model(daily_model_range,type="pred",terms=c("range_ventron","tree_density_2018"))


car::vif(daily_model)
car::vif(daily_model_max)
car::vif(daily_model_min)

#### aggregated scale ####

all_tms_data_day_agg_tmp<-all_tms_data_day[,.(T3_mean=mean(T3_mean)),by=locality_id]

all_tms_data_agg_backup<-all_tms_data_agg
all_tms_data_agg<-all_tms_data_agg[! locality_id %in% c(suspicious_loggers,new_loggers,broken_loggers),]
all_tms_data_agg<-merge(all_tms_data_agg,all_tms_data_day_agg_tmp,by="locality_id")

all_tms_data_agg[,]

lm_agg_gdd<-lm(GDD5_sum~mnt_25_vosges+Heat_load_index_25+tree_density_2018+ipv_1000_25 , data=all_tms_data_agg)

lm_agg_fdd<-lm(FDD0_sum~mnt_25_vosges+Heat_load_index_25+tree_density_2018+ipv_1000_25 , data=all_tms_data_agg)

lm_agg_snow<-lm(snow_sum~mnt_25_vosges+Heat_load_index_25+tree_density_2018+ipv_1000_25 , data=all_tms_data_agg)

lm_agg_mean<-lm(T3_mean~mnt_25_vosges+Heat_load_index_25+ipv_1000_25+tree_density_2018 , data=all_tms_data_agg)

lm_agg_mean<-lm(T3_mean~mnt_25_vosges_scale+Heat_load_index_25_scale*tree_density_2018_scale+ipv_1000_25_scale+ ipv_1000_25_scale:Heat_load_index_25_scale , data=all_tms_data_agg)


absolute_cover

lm_agg_mean<-lm(T3_mean~mnt_25_vosges , data=all_tms_data_agg)

all_tms_data_agg[,Heat_load_index_25_scale:=scale(Heat_load_index_25)]
all_tms_data_agg[,tree_density_2018_scale:=scale(tree_density_2018)]
all_tms_data_agg[,ipv_1000_25_scale:=scale(ipv_1000_25)]
all_tms_data_agg[,mnt_25_vosges_scale:=scale(mnt_25_vosges)]

lm_agg_mean<-lm(T3_mean~mnt_25_vosges+Heat_load_index_25+tree_density_2018+ipv_1000_25 , data=all_tms_data_agg)
summary(lm_agg_mean)


lm_agg_mean<-lm(T3_mean~mnt_25_vosges*(Heat_load_index_25+tree_density_2018+ipv_1000_25) , data=all_tms_data_agg)
lm_agg_mean<-lm(T3_mean~mnt_25_vosges_scale*(Heat_load_index_25_scale+tree_density_2018_scale+ipv_1000_25_scale) , data=all_tms_data_agg)

ggplot(all_tms_data_agg,aes(x=mnt_25_vosges,y=snow_sum))+geom_point()+geom_smooth(method="lm")

hist(all_tms_data_agg$snow_sum)

summary(lm_agg_mean)
summary(lm_agg_gdd)
summary(lm_agg_fdd)
summary(lm_agg_snow)

plot_model(lm_agg_gdd,type="pred",terms=c("Heat_load_index_25"))
plot_model(lm_agg_gdd,type="pred",terms=c("","mnt_25_vosges"))
plot_model(lm_agg_gdd,type="pred",terms=c("ipv_1000_25"))

plot_model(lm_agg_fdd,type="pred",terms=c("mnt_25_vosges"))

plot_model(lm_agg_snow,type="pred",terms=c("mnt_25_vosges"))
plot_model(lm_agg_snow,type="pred",terms=c("Heat_load_index_25"))
plot_model(lm_agg_snow,type="pred",terms=c("tree_density_2018"))

plot_model(lm_agg_mean,type="pred",terms=c("mnt_25_vosges"))
plot_model(lm_agg_mean,type="pred",terms=c("mnt_25_vosges","ipv_1000_25"))

plot_model(lm_agg_mean,type="pred",terms=c("ipv_1000_25"))
plot_model(lm_agg_mean,type="pred",terms=c("Heat_load_index_25"))
plot_model(lm_agg_mean,type="pred",terms=c("tree_density_2018"))

### validation agg model


ggplot(all_tms_data_agg[!is.na(GDD5_sum),],aes(x=x,y=y,color=residuals(lm_agg_gdd)))+geom_point()

all_tms_data_day_ext<-all_tms_data_day[! locality_id %in% c(suspicious_loggers,new_loggers,broken_loggers),]
all_tms_data_day_ext<-all_tms_data_day_ext[date%in% Ventron_station_day_sub[extreme_day_abs==1,date],]
all_tms_data_day_ext<-all_tms_data_day_ext[date%in% Ventron_station_day_sub[extreme_cold_day==1,date],]


all_tms_data_day_ext<-all_tms_data_day_ext[,.(T3_mean=mean(T3_mean),T3_max=mean(T3_max,na.rm=T),T3_min=mean(T3_min),var_T3_max=sd(T3_max)),by=locality_id]



all_tms_data_agg[,canopy_cover_glama:=str_replace_all(canopy_cover_glama,",",".")]
all_tms_data_agg[,canopy_cover_glama:=as.numeric(canopy_cover_glama)]



all_tms_data_day_ext<-merge(all_tms_data_day_ext,meta_data,by.x="locality_id",by.y="site_ID")

all_tms_data_day_ext[,canopy_cover_glama:=str_replace_all(canopy_cover_glama,",",".")]
all_tms_data_day_ext[,canopy_cover_glama:=as.numeric(canopy_cover_glama)]


summary(all_tms_data_day_ext)

#### extreme vents 

lm_agg_mean_ext<-lm(T3_max~mnt_25_vosges+(I(Heat_load_index_25>0.70)+tree_density_2018+ipv_1000_25) ,
                data=all_tms_data_day_ext)

lm_agg_mean_ext<-lm(T3_mean~mnt_25_vosges+(Heat_load_index_25+tree_density_2018+ipv_1000_25) ,
                    data=all_tms_data_day_ext)

lm_agg_mean_ext<-lm(T3_max~mnt_25_vosges+(Heat_load_index_25*canopy_cover_glama+ipv_1000_25) ,
                    data=all_tms_data_day_ext)

lm_agg_mean_ext<-lm(T3_mean~mnt_25_vosges+(Heat_load_index_25*canopy_cover_glama+ipv_1000_25) ,
                    data=all_tms_data_day_ext)

lm_agg_mean_ext<-lm(T3_mean~mnt_25_vosges+(Heat_load_index_25+ipv_1000_25)+tree_density_2018 ,
                    data=all_tms_data_day_ext)

lm_agg_mean_glm<-glm(T3_mean~mnt_25_vosges+Heat_load_index_25+ipv_1000_25 +tree_density_2018, data=all_tms_data_agg)
lm_agg_mean_glm<-glm(T3_mean~mnt_25_vosges+Heat_load_index_25+ipv_1000_25 , data=all_tms_data_agg)



grp<-data.frame(var=c("mnt_25_vosges","Heat_load_index_25","ipv_1000_25"),group=c("Elevation","HL","Topo"))
grp<-data.frame(var=c("mnt_25_vosges","Heat_load_index_25","ipv_1000_25","tree_density_2018"),group=c("Elevation","Topo","Topo","canopy"))
grp<-data.frame(var=c("mnt_25_vosges","Heat_load_index_25","ipv_1000_25","tree_density_2018"),group=c("Elevation","Topo","Topo","canopy"))


grp<-data.frame(var=c("mnt_25_vosges","Heat_load_index_25","ipv_1000_25","canopy_cover_glama"),group=c("Elevation","Topo","Topo","canopy"))

library(modEvA)
varPart(model=lm_agg_mean_glm,groups = grp)

summary(lm_agg_mean_glm)

lm_agg_mean_ext<-lm(T3_min~mnt_25_vosges+(ipv_1000_25) ,
                    data=all_tms_data_day_ext)

table(all_tms_data_day_ext$locality_id)

summary(lm_agg_mean_ext)
summary(lm_agg_mean)

plot(lm_agg_mean_ext)
hist(residuals(lm_agg_mean_ext),nc=15)
hist(all_tms_data_day_ext$Heat_load_index_25,nc=15)
all_tms_data_day_ext[,mean(Heat_load_index_25),by=aspect_NS_25_vosges>0]
all_tms_data_day_ext[,table(Heat_load_index_25,aspect_NS_25_vosges>0)]

plot_model(lm_agg_mean_ext,type="pred",terms=c("Heat_load_index_25","canopy_cover_glama"))

ggplot(all_tms_data_day_ext,aes(x=mnt_25_vosges,y=T3_mean))+geom_point()+geom_smooth()
ggplot(all_tms_data_day_ext,aes(x=Heat_load_index_25,y=T3_max + -0.0079476 *mnt_25_vosges ))+geom_point()+geom_smooth()
ggplot(all_tms_data_day_ext,aes(x=ipv_1000_25,y=T3_max + -0.0079476 *mnt_25_vosges ))+geom_point()+geom_smooth()

ggplot(all_tms_data_day_ext,aes(x=ipv_1000_25,y=residuals(lm_agg_mean_ext)))+geom_point()+geom_smooth(method = "lm")+geom_smooth(color="orange")

#### prediction and validation ####

validation_data<-all_tms_data_day
validation_data[,fitted_micro:=predict(daily_model,newdata=validation_data)]
validation_data[,fitted_micro_no_random:=predict(daily_model,newdata=validation_data,re.form=NA)]
validation_data[,residual_micro:=T3_min-fitted_micro]
validation_data[,residual_micro_no_random:=T3_min-fitted_micro_no_random]

sqrt(mean(validation_data$residual_micro_no_random^2,na.rm=T))


ggplot(validation_data,aes(x=fitted_micro_no_random,y=T3_min))+theme_bw()+geom_point(size=0.5,alpha=0.7)+geom_smooth()
summary(lm(T3_mean~fitted_micro_no_random,data=validation_data))
ggplot(validation_data,aes(x=fitted_micro,y=T3_mean))+theme_bw()+geom_point(size=0.5,alpha=0.7)+geom_smooth()

ggplot(validation_data,aes(x=residual_micro_no_random))+theme_bw()+geom_histogram(color="grey",fill="lightgreen",bins = 30,boundary =0)

ggplot(validation_data,aes(x=mnt_25_vosges,y=residual_micro_no_random))+theme_bw()+geom_point(size=0.5,alpha=0.7)+geom_smooth()


validation_data[!is.na(T3_mean),.(mean(T3_mean),mean(fitted_micro_no_random)),by=locality_id]

ggplot(validation_data,aes(x=date,y=residual_micro))+theme_bw()+geom_point()+geom_smooth()
ggplot(validation_data,aes(x=date,y=residual_micro_no_random))+theme_bw()+geom_point()+geom_smooth()

ggplot(validation_data,aes(x=date,y=fitted_micro_no_random,color=mnt_25_vosges,group=locality_id))+theme_bw()+geom_line()
ggplot(validation_data,aes(x=date,y=T3_mean-fitted_micro_no_random,color=mnt_25_vosges,group=locality_id))+theme_bw()+geom_line()
ggplot(validation_data,aes(x=date,y=T3_mean-fitted_micro,color=mnt_25_vosges,group=locality_id))+theme_bw()+geom_line()


ggplot(validation_data[locality_id=="H_N_3",],aes(x=date,y=T3_mean,color="true"))+theme_bw()+geom_line()+geom_line(mapping=aes(y=fitted_micro_no_random,color="fitted"))
ggplot(validation_data[locality_id==sample(unique(validation_data$locality_id),1),],aes(x=date,y=T3_max,color="true"))+theme_bw()+geom_line()+geom_line(mapping=aes(y=fitted_micro_no_random,color="fitted"))

#### coef table ####

sum_mean<-summary(daily_model)
export_model<-signif(sum_mean$coefficients,3)

write.table(export_model,file.path("figure_result","tmean_model_fiexd.csv"),sep=";",row.names = T,dec=",")

#### validate that out of the soil loggers did not affect the results ####

str(daily_model_max)
daily_model_max

  
  ## end out of soil flagging
  
  ## workflow to test out of soil sensor validity
periode_1<-ymd("2021/07/29","2021/08/25")
periode_2<-ymd("2022/04/29","2022/05/26")
#Ventron_station_day[,n_temp:=NULL]
all_tms_data_day_backup_2<-all_tms_data_day_backup

all_tms_data_day_backup_2[date %between% study_period,weight:=1/sum(!is.na(T3_mean)),by=locality_id]

all_tms_data_day_backup_2<-merge(all_tms_data_day_backup_2,meta_data[variable=="sensor_ID",],by.x="locality_id",by.y="site_ID")
all_tms_data_day_backup_2<-merge(all_tms_data_day_backup_2,Ventron_station_day,by="date",all.x=T)
all_tms_data_day_backup_2<-all_tms_data_day_backup_2[order(locality_id)]

data_periode_1<-all_tms_data_day_backup_2[date%between% periode_1,]
data_periode_1<-data_periode_1[! locality_id %in% c(suspicious_loggers,new_loggers,broken_loggers)]

data_periode_2<-all_tms_data_day_backup[date%between% periode_2,]
data_periode_2<-data_periode_2[! locality_id %in% c(suspicious_loggers,new_loggers,broken_loggers)]

data_periode_2_unamed<-all_tms_data_day_backup_2[date%between% periode_2,]
data_periode_2_unamed<-data_periode_2_unamed[! locality_id %in% c(suspicious_loggers,new_loggers,broken_loggers)]


categ_loggers<-merge(data_periode_1[,.(in_soil_1=sum(in_soil)),by=locality_id],data_periode_2[,.(in_soil_2=sum(in_soil)),by=locality_id],by="locality_id")
categ_loggers[,categ_num:= in_soil_1- in_soil_2]
categ_loggers[,categ:=ifelse(in_soil_1==0,"never_in",ifelse(categ_num==0,"always_in","transi"))]

data_periode_1<-merge(data_periode_1,categ_loggers[,c(1,5)],by="locality_id")
data_periode_2<-merge(data_periode_2,categ_loggers[,c(1,5)],by="locality_id")
data_periode_2_unamed<-merge(data_periode_2_unamed,categ_loggers[,c(1,5)],by="locality_id")

colnames(data_periode_2)<-paste0(colnames(data_periode_2),"_p2")


Validation_model<-lme4::lmer(formula(daily_model),data=data_periode_2_unamed[categ!="always_in",])
Validation_model_max<-lme4::lmer(formula(daily_model_max),data=data_periode_2_unamed[categ!="always_in",])
Validation_model_min<-lme4::lmer(formula(daily_model_min),data=data_periode_2_unamed[categ!="always_in",])


data_periode_1$T3_mean_pred<-predict(Validation_model,newdata=data_periode_1,re.form=NA)
data_periode_1$T3_max_pred<-predict(Validation_model_max,newdata=data_periode_1,re.form=NA)
data_periode_1$T3_min_pred<-predict(Validation_model_min,newdata=data_periode_1,re.form=NA)

tab_model(Validation_model)

for ( i in c("T3_mean","T3_max","T3_min")){
  print(i)
  print(summary(lm(get(i)~get(paste0(i,"_pred")),data=data_periode_1[categ=="transi",])))
  print(summary(lm(get(i)~get(paste0(i,"_pred")),data=data_periode_1[categ=="never_in",])))
  print("- - - -")
  
  
}



summary(lm(get("T3_mean")~T3_mean_pred,data=data_periode_1[categ=="transi",]))
summary(lm(T3_mean~T3_mean_pred,data=data_periode_1[categ=="never_in",]))

# ls4 hs4 problematic for tmax ?
summary(lm(T3_mean~T3_mean_pred*categ,data=data_periode_1[categ!="always_in",]))
summary(lm(T3_max~T3_max_pred*categ,data=data_periode_1[categ!="always_in",]))
summary(lm(T3_min~T3_min_pred*categ,data=data_periode_1[categ!="always_in",]))

ggplot(data_periode_1[categ!="always_in",],aes(x=T3_mean_pred,y=T3_mean,color=locality_id))+theme_bw()+geom_point()+facet_wrap(~categ)
ggplot(data_periode_1[categ!="always_in",],aes(x=T3_max_pred,y=T3_max,color=locality_id))+theme_bw()+geom_point()+facet_wrap(~categ)
ggplot(data_periode_1[categ!="always_in",],aes(x=T3_min_pred,y=T3_min,color=locality_id))+theme_bw()+geom_point()+facet_wrap(~categ)

#### extreme events #####
#compute effect size ?

ggplot(Ventron_station_day,aes(x=date,y=temp_ventr_station))+geom_line()+geom_smooth()
ggplot(Ventron_station_day,aes(x=date,y=temp_ventr_station))+geom_line()+geom_smooth()+geom_line(data=all_tms_data_day[locality_id=="H_N_8_2",],aes(x=date,y=T3_mean),color="darkred")

library(mgcv)
Ventron_station_day[,date_num:=as.numeric(date)]
gam_ventron<-gam(temp_ventr_station~s(date_num),data=Ventron_station_day)
plot(gam_ventron)
Ventron_station_day[,resid_gam:=residuals(gam_ventron)]

ggplot(Ventron_station_day_sub,aes(x=date,y=temp_ventr_station))+
  geom_line()+
  geom_smooth()+
  geom_point(data=Ventron_station_day_sub[extreme_cold_day==1,])+
  geom_line(data=all_tms_data_day[locality_id=="H_N_8_2",],aes(x=date,y=T3_mean),color="darkred")
Ventron_station_day_sub<-Ventron_station_day[date%between%study_period,]

hist(Ventron_station_day_sub$resid_gam,nc=20)

Ventron_station_day_sub[,extreme_day:=resid_gam>4]
Ventron_station_day_sub[,extreme_cold_day:=resid_gam< (-4) & date<dmy("15/05/2022")]
Ventron_station_day_sub[,extreme_day_abs:=temp_ventr_station>20]


 #### extreme events 2 ####

strange<-c("L_N_4")

all_tms_data_day[,buffering:=T3_mean - -6.783e-03 * mnt_25_vosges - 1* temp_ventr_station]
all_tms_data_day[,buffering:=T3_mean-6.783e-03 *( 1204-mnt_25_vosges ) - 1* temp_ventr_station]
all_tms_data_day[,buffering:=T3_mean - 1* temp_ventr_station]


ggplot(all_tms_data_day,aes(x = temp_ventr_station , y = T3_mean  -6.783e-03 *( 1204-mnt_25_vosges ) ))+geom_point()+theme_bw()+geom_smooth()+geom_abline(slope=1,intercept = 0,color="darkgreen")
ggplot(all_tms_data_day,aes(x = temp_ventr_station , y = T3_mean  -6.783e-03 * mnt_25_vosges - 1* temp_ventr_station))+geom_point()+theme_bw()+geom_smooth()+geom_abline(slope=0,intercept = 10,color="darkgreen")
ggplot(all_tms_data_day[temp_ventr_station>5,],aes(x = temp_ventr_station , y = T3_mean  -6.783e-03 *( 1204-mnt_25_vosges ) - 1* temp_ventr_station))+geom_point()+theme_bw()+geom_smooth()+geom_smooth(method='lm',color="orange")+geom_abline(slope=0,intercept = 0,color="darkgreen")

ggplot(all_tms_data_day[temp_ventr_station>0&!locality_id%in% strange ,],aes(color=Heat_load_index_25,x = temp_ventr_station , y = T3_mean - -6.783e-03 * mnt_25_vosges - 1* temp_ventr_station))+geom_point()+theme_bw()+geom_smooth(method = "lm")+geom_abline(slope=0,intercept = 10,color="darkgreen")
ggplot(all_tms_data_day[temp_ventr_station>5&!locality_id%in% strange,],aes(color=Heat_load_index_25>0.7,x = temp_ventr_station , y = buffering))+geom_point()+theme_bw()+geom_smooth(method='lm')+geom_abline(slope=0,intercept = 0 ,color="darkgreen")
ggplot(all_tms_data_day[temp_ventr_station>5&!locality_id%in% strange,],aes(color=mnt_25_vosges>650,x = temp_ventr_station , y = buffering))+geom_point()+theme_bw()+geom_smooth(method='lm')+geom_abline(slope=0,intercept = 0 ,color="darkgreen")

lm_buff_ext<-lm(buffering~temp_ventr_station*(Heat_load_index_25+ipv_1000_25+tree_density_2018),data=all_tms_data_day[temp_ventr_station>0 & !locality_id%in% strange,])


ggplot(all_tms_data_day[temp_ventr_station>5&!locality_id%in% strange,],aes(color=Heat_load_index_25>0.7,x = mnt_25_vosges , y = T3_max))+geom_point()+theme_bw()+geom_smooth(method='lm')+geom_abline(slope=-6.783e-03,intercept = 22)
ggplot(all_tms_data_day[temp_ventr_station>5&!locality_id%in% strange,],aes(color=Heat_load_index_25>0.7,x = mnt_25_vosges , y = buffering))+geom_point()+theme_bw()+geom_smooth(method='lm')+geom_abline(slope=-6.783e-03,intercept = 7.75)
ggplot(all_tms_data_day[temp_ventr_station>5&!locality_id%in% strange,],aes(color=Heat_load_index_25>0.7,x = Heat_load_index_25 , y = T3_mean))+geom_point()+theme_bw()+geom_smooth(method='lm')+geom_abline(slope=-6.783e-03,intercept = 7.75)




lm_buff_ext<-lme4::lmer(buffering~Heat_load_index_25+ipv_1000_25+mnt_25_vosges+temp_ventr_station+ Heat_load_index_25:temp_ventr_station +

# lm_buff_ext<-lme4::lmer(buffering~Heat_load_index_25+ipv_1000_25+tree_density_2018+poly(mnt_25_vosges,2)+temp_ventr_station+

                                                    
#lm_buff_ext<-lme4::lmer(buffering~temp_ventr_station+(Heat_load_index_25+ipv_1000_25+tree_density_2018+mnt_25_vosges)+Heat_load_index_25:temp_ventr_station+mnt_25_vosges:temp_ventr_station + 

#lm_buff_ext<-lme4::lmer(buffering~mnt_25_vosges+temp_ventr_station*(Heat_load_index_25+ipv_1000_25+tree_density_2018)+Heat_load_index_25:temp_ventr_station+
                                                    
                          (1|locality_id)+(1|date),
                        data=all_tms_data_day[temp_ventr_station>5 & !locality_id%in% strange ,],
                        weights = all_tms_data_day[temp_ventr_station>5 & !locality_id%in% strange,]$weight,
)


summary(lm_buff_ext)
tab_model(lm_buff_ext)


plot_model(lm_buff_ext,type="pred",terms = c("temp_ventr_station","Heat_load_index_25"))
plot_model(lm_buff_ext,type="pred",terms = c("temp_ventr_station","ipv_1000_25"))
plot_model(lm_buff_ext,type="pred",terms = c("temp_ventr_station","tree_density_2018"))
plot_model(lm_buff_ext,type="pred",terms = c("temp_ventr_station","mnt_25_vosges"))
plot_model(lm_buff_ext,type="pred",terms = c("Heat_load_index_25"))

#### spatial prediction #### 

summary(daily_model)
dt_raster_2<-dt_raster
dt_raster[,temp_ventr_station:=20]
dt_raster_2[,predict_microclim_mean:=predict(daily_model,dt_raster,re.form=NA)]
dt_raster_2[,med:= median()]

median(dt_raster_2$predict_microclim_mean,na.rm=T)

dt_raster_2[,mean(predict_microclim_mean,na.rm=T)]
summary(dt_raster_2[!is.na(predict_microclim_mean),]$predict_microclim_mean)

plot_comp_2<-plot_comp
plot_comp_2[,temp_ventr_station:=20]
plot_comp_2[,tree_density_2018:=tree_cover_density_thur_2018_copernicus]

plot_comp_2[,predict_microclim_mean:=predict(daily_model,plot_comp_2,re.form=NA)]
plot_comp_2[,mean(cit_climplant,na.rm=T),by=.(period,ifelse(predict_microclim_mean<median(predict_microclim_mean),"cold","warm"))]
plot_comp_2[,mean(predict_microclim_mean,na.rm=T),by=.(ifelse(predict_microclim_mean<median(predict_microclim_mean),"cold","warm"))]


dt_raster<-data.table(as.data.frame(stack_all_variable_vallee))
dt_raster<-cbind(dt_raster,xyFromCell(stack_all_variable_vallee$mnt_25_vosges,1:length(stack_all_variable_vallee$mnt_25_vosges)))

dt_raster[,tree_density_2018:=tree_cover_density_thur_2018_copernicus]


dt_raster[,temp_ventr_station:=15]
dt_raster[,max_ventron:=20]
dt_raster[,min_ventron:=10]
dt_raster[,range_ventron:=10]
dt_raster[,temp_ventr_station_corrected:=temp_ventr_station + ( get_coef["mnt_25_vosges"]*(mnt_25_vosges-1210))]



#dt_raster[,temp_ventr_max_corrected:= max_ventron + ( coef_grad_elevation_max*(mnt_25_vosges-1210))]


get_coef<-summary(daily_model)$coefficients[,1]
dt_raster[,pred_macro_mean:=get_coef["(Intercept)"] + temp_ventr_station*(get_coef["temp_ventr_station"] + mnt_25_vosges*get_coef["temp_ventr_station:mnt_25_vosges"]) + mnt_25_vosges*get_coef["mnt_25_vosges"]   ]
#dt_raster[,pred_effect_ipv:= + ipv_1000_25 * get_coef["ipv_1000_25"] + (get_coef["temp_ventr_station"] + mnt_25_vosges*get_coef["temp_ventr_station:mnt_25_vosges"]) + mnt_25_vosges*get_coef["mnt_25_vosges"]   ]

get_coef<-summary(daily_model_max)$coefficients[,1]
dt_raster[,pred_macro_max:=get_coef["(Intercept)"] + temp_ventr_station*(get_coef["max_ventron"] + mnt_25_vosges*get_coef["max_ventron:mnt_25_vosges"]) + mnt_25_vosges*get_coef["mnt_25_vosges"]   ]
get_coef<-summary(daily_model_min)$coefficients[,1]
dt_raster[,pred_macro_min:=get_coef["(Intercept)"] + temp_ventr_station*(get_coef["min_ventron"] + mnt_25_vosges*get_coef["min_ventron:mnt_25_vosges"]) + mnt_25_vosges*get_coef["mnt_25_vosges"]   ]

#dt_raster[,bdv2_25_vosges:=as.character(bdv2_25_vosges)]
#dt_raster[,bdv2_25_vosges:=ifelse(bdv2_25_vosges=="4","3",bdv2_25_vosges)]
dt_raster[,predict_microclim_mean:=predict(daily_model,dt_raster,re.form=NA)]
dt_raster[,predict_microclim_max:=predict(daily_model_max,dt_raster,re.form=NA)]
dt_raster[,predict_microclim_min:=predict(daily_model_min,dt_raster,re.form=NA)]
dt_raster[,predict_microclim_range:=predict(daily_model_range,dt_raster,re.form=NA)]

dt_raster[,predict_microclim_gdd:=predict(lm_agg_gdd,dt_raster)]
dt_raster[,predict_microclim_fdd:=predict(lm_agg_fdd,dt_raster)]
dt_raster[,predict_microclim_snow:=predict(lm_agg_snow,dt_raster)]

dt_raster[,fitted_microclim:=NA]
dt_raster[,residual_microclim:=NA]

dt_raster[,fitted_microclim:=ifelse( !is.na(predict_microclim_mean),fitted(lapse_rate),NA  )]
dt_raster[,fitted_microclim:=predict(lapse_rate,dt_raster,re.form=NA)]


dt_raster[,residual_microclim:=ifelse( !is.na(predict_microclim_mean),residuals(lapse_rate),NA  )]
dt_raster[,residual_microclim:=predict_microclim_mean-fitted_microclim]


plot_comp$fitted_microclim<-plot_comp$mnt_25_vosges* summary(daily_model)$coefficients[,1]["mnt_25_vosges"] + summary(daily_model)$coefficients[,1]["(Intercept)"] +15*summary(daily_model)$coefficients[,1]["temp_ventr_station"]
plot_comp$resid_microclim<-  plot_comp$predict_microclim_mean - plot_comp$fitted_microclim


dt_raster[,fitted_microclim:=ifelse( !is.na(predict_microclim_mean),mnt_25_vosges* summary(daily_model)$coefficients[,1]["mnt_25_vosges"] + summary(daily_model)$coefficients[,1]["(Intercept)"] +15*summary(daily_model)$coefficients[,1]["temp_ventr_station"]  ,NA  )]
dt_raster[,resid_microclim:=predict_microclim_mean-fitted_microclim]



lapse_rate<-lm(predict_microclim_mean~mnt_25_vosges,data=dt_raster[!is.na(predict_microclim_mean),])

dt_raster[,predict_snow:=predict(glm_snow,dt_raster)]

# dt_raster[,predict_mesoclim:=predict(daily_meso_model,dt_raster,re.form=NA)]
# dt_raster[,predict_offset:=predict(daily_offset,dt_raster,re.form=NA)]
# 

ggplot(dt_raster,aes(x,y,fill=predict_microclim_mean<17))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_d(na.value ="transparent")+geom_point(data=vege_plot_meta_data[Y_L93>675858,],aes(x=X_L93,y=Y_L93),fill="white",shape=3)


ggplot(dt_raster,aes(x,y,fill=predict_microclim_mean-pred_macro_mean))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )+labs(fill="offset")
ggplot(dt_raster,aes(x,y,fill=predict_microclim_mean))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")+geom_point(data=vege_plot_meta_data[Y_L93>675858,],aes(x=X_L93,y=Y_L93),fill="white",shape=3)
ggplot(dt_raster,aes(x,y,fill=predict_microclim_max))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=predict_microclim_min))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=predict_microclim_range))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")

ggplot(dt_raster,aes(x,y,fill=predict_microclim_gdd))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=predict_microclim_fdd))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=predict_microclim_snow))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")

ggplot(dt_raster,aes(x,y,fill=fitted_microclim))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=residual_microclim))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent",limits=c(-1.5,1.5),oob=scales::squish)
ggplot(dt_raster,aes(x,y,fill=residual_microclim*0.22241 ))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent",oob=scales::squish)

ggplot(dt_raster,aes(x,y,fill=ipv_1000_25* 0.0006515 ))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent",limits=c(0,2))+labs(fill="var")
ggplot(dt_raster,aes(x,y,fill=Heat_load_index_25*  1.5292014 ))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent",limits=c(0,2),oob=scales::squish)+labs(fill="var")

ggplot(dt_raster,aes(x,y,fill=tree_density_2018 ))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent",limits=c(0,40),oob=scales::squish)+labs(fill="var")

ggplot(dt_raster,aes(x,y,fill=Heat_load_index_25*  1.5292014+ipv_1000_25* 0.0006515   ))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent",limits=c(0,2),oob=scales::squish)+labs(fill="var")


ggplot(dt_raster,aes(x,y,fill=mnt_25_vosges))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )+labs(fill="offset")

ggplot(dt_raster,aes(x,y,fill=mnt_25_vosges))+theme_classic()+
  geom_raster()+
  coord_fixed()+scale_fill_gradientn(colors=terrain.colors(100),na.value ="transparent")+labs(fill="offset")+
  new_scale_fill() +
  geom_raster(mapping = aes(fill=hillshade),alpha=0.3)+
   scale_fill_gradient(low = "black", high = "white",na.value="transparent") 



#### Maps figure ####

library(ggspatial)
library(ggnewscale)
extended_mnt<-raster(file.path("data","extented_dem.tif"))
extended_mnt<-crop(big_mnt,extent(stack_all_variable_vallee$bdv2_25_vosges)+2000)
slope<-slopeAspect (extended_mnt)
hs<-hillShade(slope$slope,slope$aspect)
extended_mnt<-stack(extended_mnt,hs)

dt_raster[,extended_dem:=getValues(big_mnt)]
dt_raster[,hillshade:=getValues(hs)]


valley_border<-st_read("data",layer="valley_shp")

coords_sf_map<-coords_sf[coords_sf$site_ID%in%unique(all_tms_data_day$locality_id),]
coords_vege_map<-vege_plot_meta_data_sf[vege_plot_meta_data_sf$plot_ID%in%plot_comp$plot_ID,]

map_fig_1_dt<-data.table(as.data.frame(extended_mnt))
map_fig_1_dt<-cbind(map_fig_1_dt,xyFromCell(extended_mnt$mnt_25_vosges,1:length(extended_mnt$mnt_25_vosges)))


sampling_map<-ggplot(map_fig_1_dt,aes(x,y,fill=mnt_25_vosges))+
  theme_classic()+
  geom_raster()+
  coord_fixed()+
  labs(x="",y="",fill="Elevation (m)",shape="Study sites")+
  scale_fill_gradientn(colors=terrain.colors(100),na.value ="transparent",breaks=c(300,600,900,1200,1400))+
  new_scale_fill() +
  geom_raster(mapping = aes(fill=layer),alpha=0.35,show.legend = F)+
  scale_fill_gradient(low = "black", high = "white",na.value="transparent") +
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.5,fill=NA,alpha=0,color="grey15") + 
  geom_sf(data=coords_vege_map,aes(shape="Flora surveys"),inherit.aes=F,size=1.75)+
  geom_sf(data=coords_sf_map,aes(shape="Temperature \nloggers"),inherit.aes=F,size=2,fill="white")+
  scale_shape_manual(values=c(3,21))+
  coord_sf(expand=F)+
  annotation_north_arrow( height = unit(1, "cm"),width = unit(1, "cm"),location="tr")+
  annotation_scale()+
  theme(panel.border= element_rect(fill=NA,color="grey5",linewidth=1))

color_three_fact<-c("#0080C2FF", "#EFC000FF", "#CD534CFF")

dt_raster[,hillshade_2:=ifelse(is.na(predict_microclim_mean),NA,hillshade)]
# Heat_load_index_25*  1.5292014+ipv_1000_25* 0.0006515 
# fitted_microclim
# residual_microclim
pred_macro_mean
temp_ventr_station_corrected

annotate("label",645,20,label="Lapse rate",fill="#0080C2FF",alpha=0.3)
  annotate("label",790,17.9,label="Microclimate",fill="#EFC000FF",alpha=0.3)

  viridis::viridis(4)
  viridis::viridis(5)
  viridis::viridis(6)
  viridis::viridis(12)
  viridis::viridis(7)
  c( "#440154FF" ,"#3B528BFF" ,"#21908CFF", "#5DC863FF" ,"#FDE725FF")
  
  c("#440154FF", "#414487FF" ,"#2A788EFF", "#22A884FF" ,"#7AD151FF" ,"#FDE725FF")
  
  "#440154FF" "#482173FF" "#433E85FF" "#38598CFF" "#2D708EFF" "#25858EFF" "#1E9B8AFF" "#2BB07FFF" "#51C56AFF"
   "#85D54AFF" "#C2DF23FF" "#FDE725FF"
  
  
  dt_raster[,cut_micr:= cut(resid_microclim,breaks=c(3,-0.33,-0.66,-0.83,-1,-1.5,-3))]
  
  plot_comp$cut_microclim
  
map_microclimate <- ggplot(dt_raster,aes(x,y,fill= cut_micr))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  # scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,limits=c(-2,0))+
  # scale_fill_gradientn(colors= c("#440154FF", "#414487FF" ,"#2A788EFF", "#22A884FF" ,"#7AD151FF" ,"#FDE725FF","#FDE725FF"),values=scales::rescale(  c(-2,-1.4,-0.8,-0.3,0)) ,breaks=(  c(-2,-1.8,-1,-0.6,0)) ,na.value ="transparent",oob=scales::squish,limits=c(-2,0))+
  scale_fill_manual(values= c("#482173FF", "#38598CFF" ,"#1E9B8AFF", "#2BB07FFF" ,"#C2DF23FF" ,"#FDE725FF"),na.value ="transparent")+
  
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  coord_sf(expand=F)+
  annotation_scale()+ 
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))+
  labs(x="",y="",fill="Buffered \ntemperature °C")+annotate("label",1001325 ,6774050 ,label="Microclimate",fill="#EFC000FF",alpha=0.3,size=6)



map_elevation_climate<-ggplot(dt_raster,aes(x,y,fill=fitted_microclim))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster$fitted_microclim,probs=c(0.025,0.975),na.rm=T))+
  #scale_fill_gradientn(colors = brewer.pal(7, "PuOr")[order(7:1)],na.value ="transparent",oob=scales::squish)+
  # scale_fill_gradient2(low="lightskyblue3",mid="grey99",high="#0080C2FF",midpoint = 22,na.value ="transparent",oob=scales::squish)+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Temperature °C")+ 
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))+
  annotate("label",1001325 ,6774050 ,label="Elevation",fill="#0080C2FF",alpha=0.3,size=6)



map_whole_pred<-ggplot(dt_raster,aes(x,y,fill=predict_microclim_mean))+
  theme_classic()+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster$predict_microclim_mean,probs=c(0.025,0.975),na.rm=T))+
  #scale_fill_gradientn(colors = brewer.pal(7, "PuOr")[order(7:1)],na.value ="transparent",oob=scales::squish)+
  # scale_fill_gradient2(low="lightskyblue3",mid="grey99",high="#0080C2FF",midpoint = 22,na.value ="transparent",oob=scales::squish)+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Temperature °C")



display.brewer.all()

#  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="grey15")

  new_scale_fill() +
  geom_raster(mapping = aes(fill=hillshade_2),alpha=0.15,show.legend = F)  +
  scale_fill_gradient(low = "black", high = "white",na.value="transparent")


  
ggsave(file.path("Figure_result","sampling_map.jpg"),sampling_map,dpi=400,unit="mm",width = 88,height = 75,scale=2)

ggsave(file.path("Figure_result","climate_map.jpg"),ggarrange(map_elevation_climate,map_microclimate,align = "hv",labels=c("a)","b)")),dpi=400,unit="mm",width = 180,height = 100,scale=2)


ggplot(dt_raster,aes(x,y,fill=tree_density_2018))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")

dt_raster[,tree_density_2018:=ifelse(is.na(predict_microclim_gdd),NA,tree_density_2018)]
dt_raster[,cor(predict_microclim_mean,tmoy_digitalis_8610_25_vosges,use="complete.obs")]


ggplot(dt_raster,aes(x,y,fill=Heat_load_index_25))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=Heat_load_index_25))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")

ggplot(dt_raster,aes(x,y,fill=pred_macro_mean))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=temp_ventr_station_corrected))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")
ggplot(dt_raster,aes(x,y,fill=tmoy_digitalis_8610_25_vosges))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")

ggplot(dt_raster,aes(x=tmoy_digitalis_8610_25_vosges,y=pred_macro_mean))+theme_bw()+geom_point(alpha=0.01)+geom_smooth(method="lm")
ggplot(dt_raster,aes(x=tmoy_digitalis_8610_25_vosges,y=predict_microclim_mean))+theme_bw()+geom_point(alpha=0.01)+geom_smooth(method="lm")


tmoy_digitalis_8610_25_vosges
predict_microclim_mean
ggplot(dt_raster,aes(x,y,fill=predict_snow))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_viridis_c(na.value ="transparent")

var_to_pred<-c("predict_microclim_mean","predict_microclim_max","predict_microclim_min","predict_microclim_range","predict_microclim_gdd","predict_microclim_fdd","predict_microclim_snow")

raster_of_prediction<-foreach(variable =  var_to_pred,.combine = raster::stack)%do%{
  tmp_rast<-stack(stack_all_variable_vallee$mnt_25_vosges)
  tmp_rast<-setValues(tmp_rast,dt_raster[,get(variable)])
  return(tmp_rast)
  
  
}

names(raster_of_prediction)<-var_to_pred


vege_pred<-data.table(extract(raster_of_prediction,vege_plot_data_sf))


ggplot(all_tms_data_day,aes(x=date,y=T3_mean,color=locality_id))+geom_line()+theme_bw()

#### terrain ?####

ggplot(all_tms_data_day_backup[date>ymd("2022-08-25"),],aes(x=date,y=T3_mean,color=locality_id))+geom_line()

all_tms_data_day_backup[locality_id=="new1",][1,]
