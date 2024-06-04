source("Packages.R")

#### load some predictors ####

#if(Sys.info()["nodename"]=="CALCULUS") setwd("//tsclient/C/Users/borderieux/Desktop/data_microclim") else 
  

stack_all_variable_vallee<-stack(file.path("data","processed_raster",list.files(file.path("data","processed_raster"))))


# get copernicus data

tree_density<-raster(file.path("data","tree_density_raster","tree_cover_density_thur_2018_copernicus.tif"))
tree_density<-projectRaster(tree_density,crs=crs(stack_all_variable_vallee))

tree_density_reprojected<-projectRaster(tree_density,stack_all_variable_vallee)
stack_all_variable_vallee$tree_density_projected<-tree_density_reprojected


## precise coordinates with gnss trimble
coords<-fread(file.path("data","metadata_loggers","Logger_gps_data_trimble.csv"),dec=",")
coords[,sensor_ID:=NULL]
coords_sf<-st_as_sf(coords,coords =c("x","y") ,crs=st_crs(4326))
coords_sf_laea<-st_transform(coords_sf,st_crs(3035))
coords_sf<-st_transform(coords_sf,st_crs(2154))

coords[,xl93:=st_coordinates(coords_sf)[,1]]
coords[,yl93:=st_coordinates(coords_sf)[,2]]

env_var<- extract(stack_all_variable_vallee,coords_sf)
tree_density_2018<- extract(tree_density,coords_sf_laea)
env_var<-cbind(env_var,tree_density_2018);rm(tree_density_2018)

coords<-cbind(coords,env_var)

#### load and clean loggers data ####
setwd(file.path("~","AdapFor_analysis"))

## meta data such as date, local canopy, ste ID
meta_data<-fread("Raw_data/loggers_microclim/meta_data_2022_points.csv",dec = ",")
meta_data<-meta_data[site_ID!="",]
colnames(meta_data)[4]<-"date_setup"
meta_data[,date_setup:=dmy(date_setup)]
meta_data[,date_revisite_1:=ifelse(date_revisite_1=="",NA,date_revisite_1)]
meta_data[,date_revisite_1:=dmy(date_revisite_1)]
meta_data<-melt(meta_data,measure.vars = c("sensor_ID","sensor_ID_2"),value.name = "sensor_ID")
meta_data<-meta_data[order(sensor_ID),]
meta_data<-meta_data[!is.na(sensor_ID),]

## calibration curves for T3, obtained with a thermocouple
calibration_curves_T3<-fread("Raw_data/loggers_microclim/Calibration_T3.csv")
sum(meta_data$sensor_ID%in%calibration_curves_T3$sensor_ID)
mean(calibration_curves_T3$offset);mean(calibration_curves_T3$slope)
## data from a meteological statino of grand ventron
Ventron_station<-fread("Raw_data/loggers_microclim/station_ventron_21_22.csv")
Ventron_station[,time:=ymd_hms(time)]

## canopy cover in a 25m radius data
cover_25m_sites<-fread("Raw_data/loggers_microclim/cover_25m_sites.csv")
cover_25m_sites_aggr<-cover_25m_sites[,.(canopy_cover_25m=sum(free_cover),under_canopy_cover=sum(under_cano_cover),broadleaved_prop=(sum(ifelse(Broadleaved,free_cover,0))/sum(free_cover)  )),by=site_ID]
meta_data<-merge(meta_data,cover_25m_sites_aggr,by="site_ID",all.x=T)

## tms 4 data
folder<-file.path("C:","Program Files (x86)","Lolly","data","data_visit_2022")
folder_2<-file.path("C:","Program Files (x86)","Lolly","data","data_visit_2022_2")

data_file<-grep("data_9",list.files(folder_2,full.names = F),value=T)

get_tms_number<-function(stringfile){
  res<-str_remove_all(stringfile,"data_")
  res<-str_remove_all(res,"_0.csv")
  as.integer(res)
}
tms_number<-get_tms_number(data_file)

all_tms_data<-foreach(i=data_file,.combine = rbind)%do%{
  number_i<-get_tms_number(i)
  reading<-fread(file.path(folder_2,i))
  reading$tms_number<-number_i
  reading
}

colnames(all_tms_data)<-c("index","date","time_zone","T1","T2","T3","soil_moisture","shake","err_flag","sensor_ID")
all_tms_data<-merge(all_tms_data,meta_data[,c("site_ID" ,"sensor_ID","date_setup")],by.x="sensor_ID",by.y="sensor_ID",all.x=T)


meta_data<-merge(meta_data,coords,by.x="site_ID",by.y="site_ID")

all_tms_data[,date:=ymd_hm(date)]
all_tms_data[,date:=with_tz(date,"Europe/Paris")]
all_tms_data[,day:=floor_date(date,"day")]
all_tms_data[,month:=floor_date(date,"month")]

all_tms_data[,t1_later:=c(NA,T1[1:(.N-1)]),by=site_ID]
all_tms_data[,var_t1:=abs(T1-t1_later),by=site_ID]

all_tms_data[,t2_later:=c(NA,T2[1:(.N-1)]),by=site_ID]
all_tms_data[,var_t2:=abs(T2-t2_later),by=site_ID]

### apply the calibration procedure for T3

all_tms_data<-merge(all_tms_data,calibration_curves_T3,by="sensor_ID",all.x=T)
all_tms_data[is.na(offset),offset:=mean(calibration_curves_T3$offset)]
all_tms_data[is.na(slope),offset:=mean(calibration_curves_T3$slope)]
all_tms_data[,raw_T3:=T3]
all_tms_data[,T3:=offset+T3*slope]

hist(all_tms_data$T3,nc=50);hist(all_tms_data$raw_T3,nc=50)
all_tms_data<-all_tms_data[T1>-25&T2>-25&T3>-25,]



### remove the days where the loggers wasn't in the field

all_tms_data[,not_in_field:=floor_date(date,unit="day")<=date_setup]
all_tms_data<-all_tms_data[not_in_field==0,]
all_tms_data[,last_recorded_day:=max(day),by=sensor_ID]
all_tms_data<-all_tms_data[day!=last_recorded_day,] ## the last day of a sensor is by definition incomplete

### some loggers were replaced with others, removing the erroneous reading
ymd("2022/07/06")
meta_data[site_ID=="cad_L_8",]
all_tms_data<-all_tms_data[!(sensor_ID=="94195464" & day<=ymd("2022/07/06")),]

### remove measurements of intensive loggers ( more than the usual four 00 15 30 45 measurements)
all_tms_data[,minutes:=minute(date)]
all_tms_data<-all_tms_data[minutes%in%c(0,15,30,45),]
all_tms_data[,minutes:=NULL]

### compute daily data
Ventron_station[,time:=with_tz(time=time,"Europe/Paris")]
Ventron_station[,day:=floor_date(time,"day")]
Ventron_station[,month:=floor_date(time,"month")]
Ventron_station_day<-Ventron_station[,.(temp_ventr_station=mean(temperature),max_ventron=max(temperature)),by=day]

all_tms_data_day<-all_tms_data[,.(T1=mean(T1),T2=mean(T2),T3=mean(T3),soil_moisture=mean(soil_moisture),
                                  ecart=mean(T3)-mean(T1),amplitude=max(T2)-min(T2),
                                  max_T3=quantile(T3,probs=0.95),T3_2=mean(quantile(T3,probs=c(0.05,0.95))),
                                  var_extreme_errflag_t1=mean(var_t1,na.rm=T),var_extreme_errflag_t2=max(var_t2,na.rm=T)),by=.(day,site_ID,sensor_ID)]

all_tms_data_day<-merge(all_tms_data_day,meta_data,by=c("site_ID","sensor_ID"))

### correcte the station temperature by elevation


all_tms_data_day<-merge(all_tms_data_day,Ventron_station_day,by="day")

#### flagging the day the tms is out of the soil ####

all_tms_data_day[,max(day),by=site_ID]
nice_date<-list(scale_x_date(name   = "Date",
                             breaks = function(date) seq.Date(from = lubridate::dmy("15/07/2021"),  to = lubridate::today(),  by = "1 month"),
                             limits = c(lubridate::dmy("15/07/2021"),lubridate::today()),expand  = c(0,0),labels = scales::date_format("%b %Y")))

ggplot(all_tms_data_day[day>dmy("01/03/2021"),],aes(x=(day),y=T3,color=site_ID))+theme_bw()+geom_line(lwd=1)+geom_line(color="black",data=Ventron_station,mapping=aes(x=time,y=temperature))
ggplot(all_tms_data_day[day>dmy("01/07/2021"),],aes(x=day,y=T3_2,color=site_ID))+theme_bw()+geom_line(lwd=1)
ggplot(all_tms_data_day[day>dmy("01/07/2021"),],aes(x=day,y=offset_T3,color=site_ID))+theme_bw()+geom_line(lwd=1)

ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1%in%c(0),],aes(x=day,y=rolling_mini_err,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)
to_check_trsh<-ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1%in%c(0,1),],aes(x=day,y=var_extreme_errflag_t1,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.3,lty=2)
plotly::ggplotly(to_check_trsh)
plotly::ggplotly(ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1%in%c(0),],aes(x=day,y=rolling_mean_err,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3))
ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1==1  ,],aes(x=day,y=rolling_mean_err,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2)+geom_hline(yintercept = 0.2,lty=2)+coord_cartesian(ylim=c(-1,0.42))
plotly::ggplotly(ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1==0 ,],aes(x=day,y=delta_extreme,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2)+geom_hline(yintercept = 0.2,lty=2)+geom_hline(yintercept = -0.1,lty=2)#+coord_cartesian(ylim=c(-1,0.42))
)
ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1==0 &site_ID=="H_S_1" ,],aes(x=day,y=rolling_mean_err,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2)
ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1==0 &site_ID!="H_S_1" ,],aes(x=day,y=rolling_mean_err,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2)

ggplot(all_tms_data_day[day>dmy("01/07/2021"),],aes(x=day,y=sum_rolling,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2)

ggplot(all_tms_data_day[day>dmy("01/07/2021")& site_ID%in%c("H_S_1","L_S_2","L_S_3","H_S_2","L_S_4","cad_L_4","L_S_8","H_N_1","H_N_2") ,],aes(x=day,y=var_extreme_errflag_t1,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2)
ggplot(all_tms_data_day[day>dmy("01/07/2021")& site_ID%in%c("H_S_1","L_S_2","L_S_3","H_S_2","L_S_4","cad_L_4","L_S_8","H_N_1") ,],aes(x=day,y=sum_rolling,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2)

ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1==1& site_ID!="H_S_1" |site_ID%in%c("H_S_2","L_S_4","cad_L_4","L_S_8","H_N_1") ,],aes(x=day,y=rolling_mean_err,color=as.character(in_soil_1),group=site_ID))+theme_bw()+geom_line(lwd=1)+geom_hline(yintercept = 0.4,lty=2)#+facet_wrap(~in_soil_1,nrow = 3)
ggplot(all_tms_data_day[day>dmy("01/07/2021")& in_soil_1==1& site_ID!="H_S_1" |site_ID%in%c("H_S_2","L_S_4","cad_L_4","L_S_8","H_N_1") ,],aes(x=day,y=rolling_mini_err_l,color=as.character(in_soil_1),group=site_ID))+theme_bw()+geom_line(lwd=1)+geom_hline(yintercept = 0.4,lty=2)#+facet_wrap(~in_soil_1,nrow = 3)


all_tms_data_day[in_soil==1 & site_ID != "H_S_1",quantile(var_extreme_errflag_t1,probs=0.75)+1.5*(quantile(var_extreme_errflag_t1,probs=0.75)-quantile(var_extreme_errflag_t1,probs=0.25))]


meta_data[,back_in_soil_date:=date_revisite_1]


meta_data[site_ID=="H_S_1",back_in_soil_date:=ymd("2022/06/04")]
meta_data[site_ID=="H_N_2",back_in_soil_date:=ymd("2022/01/19")]


all_tms_data_day[,rolling_mini_err_r:=frollapply(var_extreme_errflag_t1,15,FUN = min,align = "right"),by=site_ID]
all_tms_data_day[,rolling_mini_err_l:=frollapply(var_extreme_errflag_t1,15,FUN = min,align = "left"),by=site_ID]

all_tms_data_day[,rolling_mean_err:=frollapply(var_extreme_errflag_t1,10,FUN = mean,align = "left"),by=site_ID]
all_tms_data_day[,rolling_mean_err_lag:=data.table::shift(rolling_mean_err,10,type="lead"),by=site_ID]
all_tms_data_day[,delta_extreme:=(rolling_mean_err-rolling_mean_err_lag),by=site_ID]


all_tms_data_day[,sum_rolling:=abs(rolling_mini_err_r-rolling_mini_err_l)]

all_tms_data_day[,var_extreme_errflag_t1_shift:=data.table::shift(var_extreme_errflag_t1,1),by=site_ID]
all_tms_data_day[,delta_extreme:=(var_extreme_errflag_t1-var_extreme_errflag_t1_shift),by=site_ID]
all_tms_data_day[,delta_extreme_rolling:=frollapply(delta_extreme,15,FUN = max,align = "right"),by=site_ID]



treshold_out<-0.35


meta_data[,put_out_date_1:=NA]

all_tms_data_day[,treshold_crossed:=rolling_mean_err>treshold_out]
all_tms_data_day[,in_soil:=1]
for (sensor in unique(all_tms_data_day[,sensor_ID])){
  day_put_in<-meta_data[sensor_ID==sensor,date_setup]
  
  day_crossing<-all_tms_data_day[treshold_crossed==T& sensor_ID==sensor& day>day_put_in,day][1]
  print(paste0("Sensor n° ",sensor," was put out of soil on ",day_crossing))
  all_tms_data_day[,in_soil:=ifelse(sensor_ID==sensor,ifelse(day>=day_crossing,0,1),in_soil)]
  meta_data[,put_out_date_1:=ifelse(sensor_ID==sensor,as.character(day_crossing),put_out_date_1)]
  
  
}


meta_data[,put_out_date_1:=ymd(put_out_date_1)]


table(meta_data$put_out_date_1)
table(meta_data$in_soil_1,meta_data$put_out_date_1)
table(meta_data$in_soil_1,is.na(meta_data$put_out_date_1))

sum(is.na(meta_data$put_out_date_1))

put_out_logg<-all_tms_data_day[!is.na(put_out_date_1),unique(site_ID)]
all_tms_data_day[,in_soil:=ifelse(is.na(in_soil),1,in_soil)]

all_tms_data_day[in_soil_1==2,table(in_soil)]


### two special case : 94207915 and  94207913
### 94207915 was in a building for several days, then it was planted again

ggplotly( ggplot(all_tms_data_day[site_ID=="H_N_2" ,],aes(x=day,y=var_extreme_errflag_t1,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2))
date_out<-ymd(c("2021/08/27","2022/07/13"))
date_int<-ymd("2022/01/19")
all_tms_data_day[site_ID=="H_N_2"& day>date_out[1] ,in_soil:=0]
all_tms_data_day[site_ID=="H_N_2"& day>date_int[1] ,in_soil:=1]
all_tms_data_day[site_ID=="H_N_2"& day>date_out[2],in_soil:=0]
ggplotly( ggplot(all_tms_data_day[site_ID=="H_S_1" ,],aes(x=day,y=var_extreme_errflag_t1,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2))
date_out<-ymd(c("2021/06/03"))
date_int<-ymd("2022/06/04")
all_tms_data_day[site_ID=="H_S_1" ,in_soil:=1]
all_tms_data_day[site_ID=="H_S_1"& day>date_out[1] ,in_soil:=0]
all_tms_data_day[site_ID=="H_S_1"& day>date_int[1],in_soil:=1]

## end special case


#### checking if the results are influenced by out of soil loggers ####

prop_soil_day<-all_tms_data_day_backup[,.(prop=mean(in_soil),prop_2=sum(in_soil)/49,n_log=.N),by=date]

Ventron_station_day[,n_temp:=(temp_ventr_station-min(temp_ventr_station))/(max(temp_ventr_station)-min(temp_ventr_station))]

ggplot(prop_soil_day,aes(x=date,y=prop))+geom_line()+theme_bw()+geom_line(aes(y=as.numeric(n_log==49)),color="darkred")#+ geom_line(color="black",lwd=1,mapping =aes(x=day,y=n_temp) ,data =Ventron_station_day)

ggplotly(ggplot(prop_soil_day,aes(x=day,y=prop))+geom_line()+theme_bw()+geom_line(aes(y=as.numeric(n_log==49)),color="darkred")+ geom_line(color="black",lwd=1,mapping =aes(x=day,y=n_temp) ,data =Ventron_station_day))

### methodology: after 07/07/2021 take the 18 first thermometer to be put out of soil, and compare them with 18 "good" loggers

good_reading<-meta_data[is.na(put_out_date_1),site_ID]
out_to_test_reading<-meta_data[order(put_out_date_1),][1:30,][order(date_revisite_1)][,.(site_ID,put_out_date_1,date_revisite_1)][10:28,site_ID]
all_tms_data_day_test<-all_tms_data[site_ID%in%c(good_reading,out_to_test_reading),]
all_tms_data_day_test<-all_tms_data_day[site_ID%in%c(good_reading,out_to_test_reading),]
all_tms_data_day_test[,type_of_reading:=ifelse(site_ID%in%good_reading,"Good","out_of_soil")]

meta_data[site_ID%in%out_to_test_reading,sort(put_out_date_1)]

dates_first_last_out<-meta_data[site_ID%in%out_to_test_reading,.(min=min(put_out_date_1),max=max(put_out_date_1))]

plot_sense<-ggplot(all_tms_data_day_test[day>ymd("2022/08/05"),],aes(x=date,y=T3,color=as.character(sensor_ID)))+geom_line()+theme_bw()+
  geom_vline(lty=2,xintercept =unlist(dates_first_last_out[,]) )+facet_wrap(~type_of_reading,nrow=2)+geom_hline(yintercept = 30,lty=2,lwd=0.5)+geom_hline(yintercept = 20,lty=2,lwd=0.5)

plot_sense<-ggplot(all_tms_data_day_test[day>ymd("2021/07/08")& day<ymd("2021/08/01"),],aes(x=date,y=T3,color=as.character(sensor_ID)))+geom_line()+theme_bw()+
  geom_vline(lty=2,xintercept =unlist(dates_first_last_out[,]) )+facet_wrap(~type_of_reading,nrow=2)+geom_hline(yintercept = 30,lty=2,lwd=0.5)+geom_hline(yintercept = 20,lty=2,lwd=0.5)

plot_sense<-ggplot(all_tms_data_day_test[day>ymd("2022/03/15")& day<ymd("2022/06/15"),],aes(x=date,y=T3,color=as.character(sensor_ID)))+geom_line()+theme_bw()+
  geom_vline(lty=2,xintercept =unlist(dates_first_last_out[,]) )+facet_wrap(~type_of_reading,nrow=2)+geom_hline(yintercept = 30,lty=2,lwd=0.5)+geom_hline(yintercept = 20,lty=2,lwd=0.5)


plot_sense<-ggplot(all_tms_data_day_test[day%between% dmy(c("17-08-2021","05-07-2022")),],aes(x=day,y=T3,color=as.character(sensor_ID)))+geom_line()+theme_bw()+
  geom_vline(lty=2,xintercept =unlist(dates_first_last_out[,]) )#+facet_wrap(~type_of_reading,nrow=2)



all_tms_data_day_test<-all_tms_data[site_ID%in%c(put_out_logg),]
all_tms_data_day_test[,gap:=grepl("^L_",site_ID)]

plot_sense<-ggplot(all_tms_data_day_test[day>ymd("2021/07/08")& day<ymd("2021/08/01"),],aes(x=date,y=T3,color=as.character(sensor_ID)))+geom_line()+theme_bw()+
  geom_vline(lty=2,xintercept =unlist(dates_first_last_out[,]) )+geom_hline(yintercept = 30,lty=2,lwd=0.5)+geom_hline(yintercept = 20,lty=2,lwd=0.5)#+facet_wrap(~gap,nrow=2)

plot_sense<-ggplot(all_tms_data_day_test[,],aes(x=date,y=T3,color=as.character(sensor_ID)))+geom_line()+theme_bw()+
  geom_vline(lty=2,xintercept =unlist(dates_first_last_out[,]) )+geom_hline(yintercept = 30,lty=2,lwd=0.5)+geom_hline(yintercept = 20,lty=2,lwd=0.5)+facet_wrap(~gap,nrow=2)





library(plotly)
ggplotly(plot_sense)
htmlwidgets::saveWidget(ggplotly(plot_sense),"widget_communication/All_tms_15min.html")

## it seems that H_S_1 has been moved out the soil, replaced again

ggplot()


### second methodologie, study two time period and 3 categ of loggers
## constant out of soil, constant in, and transition

periode_1<-ymd("2021/07/29","2021/08/25")
periode_2<-ymd("2022/04/29","2022/05/26")
Ventron_station_day[,n_temp:=NULL]


data_periode_1<-all_tms_data_day[day%between% periode_1,]
data_periode_2<-all_tms_data_day[day%between% periode_2,]
data_periode_2_unamed<-all_tms_data_day[day%between% periode_2,]


categ_loggers<-cbind(data_periode_1[,.(in_soil_1=sum(in_soil)),by=site_ID],data_periode_2[,.(in_soil_2=sum(in_soil)),by=site_ID])
categ_loggers[,categ_num:= in_soil_1- in_soil_2]
categ_loggers[,categ:=ifelse(in_soil_1==0,"never_in",ifelse(categ_num==0,"always_in","transi"))]

data_periode_1<-merge(data_periode_1,categ_loggers[,c(1,6)],by="site_ID")
data_periode_2<-merge(data_periode_2,categ_loggers[,c(1,6)],by="site_ID")
data_periode_2_unamed<-merge(data_periode_2_unamed,categ_loggers[,c(1,6)],by="site_ID")

colnames(data_periode_2)<-paste0(colnames(data_periode_2),"_p2")

#merge_periode<-cbind(data_periode_1,data_periode_2)
merge_periode<-merge(data_periode_1,data_periode_2,by.x=c("site_ID","day"),by.y=c("site_ID_p2","day_p2"))
merge_periode[,delta_t3:=T3-T3_p2]
merge_periode[,true_delta:= (T3_p2)/(T3)]
merge_periode[,delta_ref:=temp_ventr_station_p2 /temp_ventr_station]
#merge_periode<-merge(data_periode_1,data_periode_2,by=c("site_ID"),suffixes=c("_1","_2"))

merge_periode[site_ID%in% balanced_plot ,lapply(.SD,mean),by=categ,.SDcols=c("T3","T3_p2","delta_t3","true_delta","temp_ventr_station","temp_ventr_station_p2","delta_ref","mnt_25_vosges") ]

lm_find_delta<-lm(true_delta~poly(delta_ref,2)+categ,data=merge_periode)
lm_find_delta<-lm(true_delta~delta_ref*categ+mnt_25_vosges,data=merge_periode[site_ID%in% balanced_plot])
lm_find_delta<-lm(true_delta~poly(delta_ref,2)*categ+mnt_25_vosges,data=merge_periode[site_ID%in% balanced_plot])

plot_model(lm_find_delta,type="pred",terms=c("delta_ref","categ"))
summary(lm_find_delta)

ggplot(merge_periode[],aes(x=delta_ref,y=true_delta,color=categ))+theme_bw()+geom_point()+geom_smooth()

merge_periode[,delta_t3:=T3-T3_p2]

tmp_comput<-merge_periode[,.(T3=mean(T3),T3_p2=mean(T3_p2),delta_t3=mean(delta_t3),temp_ref=mean(temp_ventr_station),temp_ref_p2=mean(temp_ventr_station_p2),true_delta=mean(true_delta),delta_ref=mean(delta_ref)),by=.(site_ID,categ)]
tmp_comput[,true_delta:= (T3_p2)/(T3)]
tmp_comput[,delta_ref:=temp_ref_p2/temp_ref]
tmp_comput[site_ID%in% balanced_plot,lapply(.SD,mean),by=categ,.SDcols=c("T3","T3_p2","delta_t3","true_delta","temp_ref","temp_ref_p2","delta_ref") ]


table(categ_loggers$site_ID,categ_loggers$categ)

balanced_plot<-c("slope_H_1","slope_H_2","slope_H_4", 
                 "L_S_5" ,"L_S_6" ,"L_S_7",
                 "H_S_1", "H_S_1_2" ,"H_S_2",
                 "H_S_3" ,"H_S_4", "H_S_5",
                 "H_N_8_2")


## perform models
suspicious_loggers<-c("L_S_1","L_S_3","slope_H_1","L_S_7")

daily_model<-lme4::lmer(T3~temp_ventr_station+mnt_25_vosges+aspect_NS_25_vosges+tree_density_2018 +(1|site_ID),
                        data=data_periode_1[categ!="never_in",])

daily_model<-lme4::lmer(T3~mnt_25_vosges+temp_ventr_station_corrected*tree_density_projected+aspect_NS_25_vosges+I(as.factor(bdv2_25_vosges))+ipv_1000_25 +(1|site_ID),
                        data=data_periode_2_unamed[categ!="always_in",][!site_ID%in%suspicious_loggers,])

daily_model<-lme4::lmer(max_T3~max_ventron+mnt_25_vosges+aspect_NS_25_vosges+tree_density_2018 +(1|site_ID),
                        data=data_periode_2_unamed[categ!="always_in",][!site_ID%in%suspicious_loggers,])

[!site_ID%in%suspicious_loggers,]

tab_model(daily_model)

summary(daily_model)
merge_periode_2<-merge_periode[categ!="always_in",][order(site_ID,day),][!site_ID%in%suspicious_loggers,]
merge_periode_2$prediction<- predict(daily_model,newdata=merge_periode_2[order(site_ID,day),][!site_ID%in%suspicious_loggers,],re.form=NA)
merge_periode_2[,se:=(T3-prediction)^2]
merge_periode_2[,error:=T3-prediction]
merge_periode_2[,ratio:=T3/prediction]

ggplot(merge_periode_2,aes(prediction,T3,color=site_ID))+theme_bw()+geom_point()+facet_wrap(~categ)
ggplot(merge_periode_2,aes(prediction,T3))+theme_bw()+geom_point()+facet_wrap(~categ)+geom_smooth()

ggplot(merge_periode_2,aes(mnt_25_vosges,error))+theme_bw()+geom_point()+facet_wrap(~categ)+geom_smooth()+geom_smooth(method = "lm",col="orange")


merge_periode_2[,lapply(.SD,mean),by=categ,.SDcols=c("se","error","T3","ratio","prediction","aspect_NS_25_vosges","mnt_25_vosges","temp_ventr_station","temp_ventr_station_p2") ]

summary(lm(T3~prediction,data=merge_periode_2[categ=="transi",]))
summary(lm(T3~prediction,data=merge_periode_2[categ=="always_in",]))
summary(lm(T3~prediction,data=merge_periode_2[categ=="never_in",]))


#### performe analysis and modelling ####
ggplot(all_tms_data_day_backup[locality_id=="H_N_2",],aes(x=date,y=T3_mean))+geom_line()
ggplot(all_tms_data[site_ID=="H_N_2",],aes(x=date,y=T3))+geom_line()


suspicious_loggers<-c("L_S_1","L_S_3","slope_H_1","L_S_7")
new_loggers<-c("new1")
broken_loggers<-c("cad_L_8")

all_tms_data_day_backup<-all_tms_data_day
all_tms_data_day<-all_tms_data_day[! site_ID %in% c(suspicious_loggers,new_loggers)]

study_period<-c(dmy("01/05/2022"),dmy("21/08/2022"))
#all_tms_data_day<-all_tms_data_day[day%between%study_period,]
all_tms_data_day[day %between% study_period,weight:=1/.N,by=site_ID]


#### correct the weather station with elevation data
daily_meso_model_get_coef<-lme4::lmer(T3~mnt_25_vosges+temp_ventr_station+(1|site_ID),
                                      data=all_tms_data_day[day %between% study_period,],
                                      weights = all_tms_data_day[day %between% study_period,weight])
daily_meso_model_get_coef<-lme4::lmer(max_T3~mnt_25_vosges+max_ventron+(1|site_ID),
                                      data=all_tms_data_day[day %between% study_period,],
                                      weights = all_tms_data_day[day %between% study_period,weight])

summary(daily_meso_model_get_coef)
tab_model(daily_meso_model_get_coef)
coef_grad_elevation<-    -0.0063960  
coef_grad_elevation_max<-   -0.0071134  
all_tms_data_day[,temp_ventr_station_corrected:= temp_ventr_station + ( coef_grad_elevation*(mnt_25_vosges-1210))]
all_tms_data_day[,temp_ventr_max_corrected:= max_ventron + ( coef_grad_elevation_max*(mnt_25_vosges-1210))]



ggplot(all_tms_data_day[day%between% study_period,],aes(x=temp_ventr_station_corrected,y=T1))+theme_bw()+geom_point()+geom_smooth(method = 'lm')

ggplot(all_tms_data_day[day%between% study_period,],aes(x=mnt_25_vosges,y=T3,color=as.factor(day)))+theme_bw()+geom_point(show.legend=F)+geom_smooth(show.legend=F,method = 'lm')
ggplot(all_tms_data_day[day%between% study_period,],aes(x=temp_ventr_station,y=T3))+theme_bw()+geom_point(show.legend=F)+geom_smooth(show.legend=F,method = 'lm')
ggplot(all_tms_data_day[day%between% study_period,],aes(x=aspect_NS_25_vosges,y=T3))+theme_bw()+geom_point(show.legend=F)+geom_smooth(show.legend=F,method = 'lm')

ggplot(aggr_daily,aes(x=mnt_25_vosges,y=T3  -0.0056772*(1300-mnt_25_vosges)  ))+theme_bw()+geom_point(show.legend=F)+geom_smooth(show.legend=F,method = 'lm')
ggplot(aggr_daily,aes(x=ipv_1000_25,y=T3  -0.0056772*(1300-mnt_25_vosges)  ))+theme_bw()+geom_point(show.legend=F)+geom_smooth(show.legend=F,method = 'lm')


plot_ventron<-ggplot(all_tms_data_day[day%between% study_period,],aes(x=day,y=T3,color=site_ID))+theme_bw()+geom_line()+
  geom_line(color="black",lwd=1,mapping =aes(x=day,y=temp_ventr_station) ,data =all_tms_data_day[site_ID=="L_S_6" & day%between% study_period,] )

plot_ventron<-function(sites,study_period,lwd=1){
  res<-ggplot(all_tms_data_day[site_ID%in%sites &day%between% study_period,],aes(x=day,y=T3,color=site_ID))+theme_bw()+geom_line(lwd=lwd,lty=1,mapping =aes(x=day,y=temp_ventr_station_corrected,color=paste0(site_ID,"_air")) ,data =all_tms_data_day[site_ID%in%sites & day%between% study_period,] )+
    geom_line(lwd=lwd,lty=1)
  res<-res+scale_color_manual(breaks=c(sites,paste0(sites,"_air")) ,
                         values=c("slategray3","tomato","slateblue4","indianred3"))
  res<-res
  res
  }
  

plot_ventron(c("H_S_7","cad_L_1"),dmy(c("15/07/21","15/07/22")))
res<-plot_ventron(c("H_S_7","cad_L_1"),dmy(c("01/03/22","15/08/22")),0.5)+labs(x="Date 2022",y="Temperature C°",title="Temperature à l'exterieur et l'interieur de la forêt",subtitle="Avril - juin 2022")
res_2<-plot_ventron(c("H_S_7","cad_L_1"),dmy(c("01/08/21","15/03/22")),0.5)+labs(x="Date 2022",y="Temperature C°",title="Temperature à l'exterieur et l'interieur de la forêt",subtitle="Aout 2021 - mars 2022")

plot_ventron_ly<-ggplotly(res)

ggsave("widget_communication/two_curves.pdf",res)
#ggsave("widget_communication/two_curves.tif",res,device="tiff")
#ggsave("widget_communication/two_curves_winter.tif",res_2,device="tiff")

ventron_plotly<-ggplotly(plot_ventron)

htmlwidgets::saveWidget(plot_ventron_ly,"widget_communication/two_curves_summer.html")

aggr_daily<-all_tms_data_day[day%between% study_period,lapply(.SD,quantile,probs=0.95),by=site_ID, .SDcols=c("T3","max_T3","offset_T3","temp_ventr_station","mnt_25_vosges","aspect_NS_25_vosges","bdv2_25_vosges","ipv_1000_25","tree_density_projected","canopy_cover_25m")]
aggr_daily<-all_tms_data_day[day%between% study_period,lapply(.SD,function(x)if(class(x)=="character")unique( x) else mean(x)),by=site_ID, .SDcols=c("T3","max_T3","offset_T3","temp_ventr_station","temp_ventr_station_corrected","mnt_25_vosges","aspect_NS_25_vosges","bdv2_25_vosges","ipv_1000_25","tree_density_projected","canopy_cover_25m")]


aggr_daily_model<-lm(T3~temp_ventr_station_corrected*tree_density_projected+aspect_NS_25_vosges+I(as.factor(bdv2_25_vosges))+ipv_1000_25+tree_density_projected,
                data=aggr_daily)
summary(aggr_daily_model)



mean(all_tms_data_day[day %between% study_period,]$temp_ventr_station)
weighted.mean(all_tms_data_day[day %between% study_period,]$temp_ventr_station,all_tms_data_day[day %between% study_period,]$weight)
mean(all_tms_data_day[day %between% study_period,]$max_ventron)
weighted.mean(all_tms_data_day[day %between% study_period,]$temp_ventr_max_corrected,all_tms_data_day[day %between% study_period,]$weight)

mean(all_tms_data_day[day %between% study_period,]$temp_ventr_station_corrected)

all_tms_data_day



library(lme4)
daily_model<-lme4::lmer(max_T3~temp_ventr_max_corrected*tree_density_projected+aspect_NS_25_vosges:temp_ventr_station_corrected+temp_ventr_station_corrected:ipv_1000_25+aspect_NS_25_vosges+I(as.factor(bdv2_25_vosges))+ipv_1000_25 +
                          #(1|site_ID)+(1|day),
                          (1|site_ID),
                        data=all_tms_data_day[day %between% study_period,],
                        weights = all_tms_data_day[day %between% study_period,weight])
daily_meso_model<-lme4::lmer(max_T3~temp_ventr_max_corrected+(1|site_ID),
                     data=all_tms_data_day[day %between% study_period,])

daily_model<-lme4::lmer(T3~temp_ventr_station_corrected*tree_density_projected+aspect_NS_25_vosges+I(as.factor(bdv2_25_vosges))+ipv_1000_25 +(1|site_ID),
                        data=all_tms_data_day[day %between% study_period,])

all_tms_data_day[day %between% study_period,weight:=1/.N,by=site_ID  ]
all_tms_data_day[,temp_ventr_station_corrected_2:=temp_ventr_station_corrected^2 ]
all_tms_data_day[,temp_ventr_station_corrected_day_before:=c(NA,temp_ventr_station_corrected[1:(.N-1)]),by=site_ID ]

all_tms_data_day[,bdv2_25_vosges:=as.character(bdv2_25_vosges)]

daily_model<-lme4::lmer(T3~temp_ventr_station_corrected*tree_density_projected+aspect_NS_25_vosges:temp_ventr_station_corrected+temp_ventr_station_corrected:ipv_1000_25+aspect_NS_25_vosges+aspect_NS_25_vosges:temp_ventr_station_corrected+bdv2_25_vosges+ipv_1000_25 +
                          #(1|site_ID)+(1|day),
                          (1|site_ID),
                        data=all_tms_data_day[day %between% study_period,],
                        weights = all_tms_data_day[day %between% study_period,weight])


daily_model_autocor<-lme(T3~temp_ventr_station_corrected*tree_density_projected+aspect_NS_25_vosges:temp_ventr_station_corrected+temp_ventr_station_corrected:ipv_1000_25+aspect_NS_25_vosges+aspect_NS_25_vosges:temp_ventr_station_corrected+bdv2_25_vosges+ipv_1000_25 ,
                         random= ~day|site_ID ,
                         correlation=corCAR1(form=~day|site_ID),
                         data=all_tms_data_day[day %between% study_period,],
                         weights = ~weight)


daily_meso_model<-lme4::lmer(T3~temp_ventr_station_corrected+(1|site_ID),
                             data=all_tms_data_day[day %between% study_period,])

tab_model(daily_model)
tab_model(daily_model_autocor)

library(nlme)

nl_daily_model<-lme4::nlmer(T3~resu(temp_ventr_station,coef_grad_elevation,mnt_25_vosges)  ~tree_density_projected+aspect_NS_25_vosges+bdv2_25_vosges+ipv_1000_25+(coef_grad_elevation|site_ID) ,
                        data=all_tms_data_day[day %between% study_period,],
                        start=c(coef_grad_elevation=-0.0062926))


nl_daily_model<-gnls(T3~ (coef_cover*tree_density_projected+ coef_clim)*(temp_ventr_station + ( coef_elev*(mnt_25_vosges-1210))) +intercept+coef_aspect*aspect_NS_25_vosges + coef_ipv*ipv_1000_25 ,
     all_tms_data_day[day %between% study_period,],
     start=c(coef_clim=1,coef_elev=-0.0066365,intercept=10 ,coef_aspect= -0.3360430,coef_cover=0 ,coef_ipv=0.0006061 ),na.action="na.pass",
     control=gnlsControl(nlsTol=0.01, nlsMaxIter=20))




nl_daily_model<-nlme(T3~ coef_clim*(temp_ventr_station + ( coef_elev*(mnt_25_vosges-1210))) +intercept+coef_aspect*aspect_NS_25_vosges + coef_cover*tree_density_projected+ coef_ipv*ipv_1000_25 ,
                     all_tms_data_day[day %between% study_period,],
                     fixed=coef_elev+intercept+coef_aspect+coef_cover+coef_ipv~1,
                     random=site_ID~1,
                     start=c(site_ID=1,coef_elev=-0.0066365,intercept=5,coef_aspect= -0.3360430,coef_cover=-0.0341302 ,coef_ipv=0.0006061 ),na.action="na.pass",
                     )


resu<-deriv(~temp_ventr_station + ( coef_grad_elevation*(mnt_25_vosges-1210)),
            namevec="coef_grad_elevation",
            function.arg = c("temp_ventr_station","coef_grad_elevation","mnt_25_vosges"))
selfStart(resu)
daily_meso_model<-lme4::lmer(T3~temp_ventr_station_corrected+(1|site_ID),
                             data=all_tms_data_day[day %between% study_period,])


car::vif(daily_model)
summary(daily_model)
tab_model(daily_model)
tab_model(daily_meso_model)
summary(nl_daily_model)

### 1st try to model offests instead of absolute value: 


ll_tms_data_day[,offset_T3:=T3 -temp_ventr_station_corrected ]

ggplot(all_tms_data_day,aes(max_ventron,max_T3))+theme_bw()+geom_point()+geom_smooth()

plot_model(daily_model,type="pred",terms=c("mnt_25_vosges","temp_ventr_station"))

all_tms_data_day[,bdv2_25_vosges:=as.character(bdv2_25_vosges)]

daily_offset<-lme4::lmer(offset_T3~temp_ventr_station+mnt_25_vosges+aspect_NS_25_vosges+slope_25_vosges+bdv2_25_vosges+ipv_1000_25+tree_density_projected +(1|site_ID),
                        data=all_tms_data_day[day %between% study_period,])

tab_model(daily_offset)
summary(daily_offset)
plot(daily_offset)
ipv<-plot_model(daily_model,type="pred",terms=c("ipv_1000_25","temp_ventr_station_corrected"))
plot_model(daily_model,type="pred",terms=c("tree_density_2018","temp_ventr_station_corrected"))
aspect<-plot_model(daily_model,type="pred",terms=c("aspect_NS_25_vosges","temp_ventr_station_corrected"))

custom<-list(theme_bw(),labs(title=" ",y="Temperature inside forest",color="Outside forest T°"),scale_color_manual(values = c("cyan4","gold3","coral3")),scale_fill_manual(values = c("cyan4","gold3","coral3")))

plot_factor<-ggarrange(ipv+custom+labs(x="Bottom - top of the valley position"),aspect+custom+labs(x="South - North exposure"),widths = 2,common.legend = T,legend="bottom")
ggsave('figures/microclim/factors_micro.png',plot_factor,units = "in",height =4 ,width = 6,dpi=340)

#### validation ? ####
daily_model<-daily_model_autocor

validation_data<-all_tms_data_day[day%between% study_period,]
validation_data[,fitted_micro:=predict(daily_model,newdata=validation_data)]
validation_data[,fitted_micro_no_random:=predict(daily_model,newdata=validation_data,re.form=NA)]
validation_data[,residual_micro:=residuals(daily_model)]
validation_data[,residual_micro_no_random:=T3-fitted_micro_no_random]

ggplot(validation_data,aes(x=fitted_micro_no_random,y=T3))+theme_bw()+geom_point(size=0.5,alpha=0.7)+geom_smooth()
ggplot(validation_data,aes(x=fitted_micro,y=T3))+theme_bw()+geom_point(size=0.5,alpha=0.7)+geom_smooth()

plot_valid<-ggplot(validation_data,aes(y=fitted_micro_no_random,x=T3))+theme_bw()+geom_point(size=0.5,alpha=0.7)+geom_smooth()+labs(x="Measurements",y="Model prediction")+geom_abline(slope = 1,intercept = 0,lty=2)
ggsave('figures/microclim/valid_plot.png',plot_valid,units = "in",height =4 ,width = 6,dpi=340)

ggplot(validation_data,aes(x=T3,y=residual_micro_no_random))+theme_bw()+geom_point()+geom_smooth()
ggplot(validation_data,aes(x=temp_ventr_station_corrected,y=residual_micro_no_random))+theme_bw()+geom_point()+geom_smooth()
ggplot(validation_data,aes(x=tree_density_projected,y=T3 ))+theme_bw()+geom_point()+geom_smooth()
ggplot(validation_data,aes(x=residual_micro ))+theme_bw()+geom_histogram(color="grey",fill="lightgreen")

ggplot(validation_data,aes(x=day,y=residual_micro))+theme_bw()+geom_point()+geom_smooth()
ggplot(validation_data,aes(x=day,y=residual_micro_no_random))+theme_bw()+geom_point()+geom_smooth()

ggplot(validation_data,aes(x=day,y=fitted_micro_no_random,color=mnt_25_vosges,group=site_ID))+theme_bw()+geom_line()
ggplot(validation_data,aes(x=day,y=T3-fitted_micro_no_random,color=mnt_25_vosges,group=site_ID))+theme_bw()+geom_line()
ggplot(validation_data,aes(x=day,y=T3-fitted_micro,color=mnt_25_vosges,group=site_ID))+theme_bw()+geom_line()

ggplot(validation_data[site_ID=="H_N_3",],aes(x=day,y=T3,color="true"))+theme_bw()+geom_line()+geom_line(mapping=aes(y=fitted_micro,color="fitted"))
ggplot(validation_data[site_ID==sample(unique(validation_data$site_ID),1),],aes(x=day,y=T3,color="true"))+theme_bw()+geom_line()+geom_line(mapping=aes(y=fitted_micro,color="fitted"))


validation_site<-validation_data[,lapply(.SD,mean),by=site_ID,.SDcols =c ("T3","fitted_micro","fitted_micro_no_random","residual_micro","residual_micro_no_random")]
mean(validation_site$residual_micro_no_random);mean(abs(validation_site$residual_micro_no_random))


ggplot(validation_site,aes(x=residual_micro_no_random))+geom_histogram()


#### vegetation surveys dataset ####

vege_plot_data<-fread("C:/Users/borderieux/Desktop/docs/reconciliation_ecoplant_sophy_ancien_ifn/bd_vosges/plot_data_vosges.csv")
vege_plot_data[is.na(Y_L93),]
vege_plot_data<-vege_plot_data[!is.na(Y_L93),]
vege_plot_data<-vege_plot_data[Y_L93>6464402 ,]
vege_plot_data<-vege_plot_data[SOURCE!="",]
vege_plot_data<-vege_plot_data[ANNEE>=2000,]
vege_plot_data<-vege_plot_data[Vallee=="Thur",]

vege_plot_data<-vege_plot_data[SOURCE=="APT",]

vege_plot_data[,.(.N,mean(`altitude MNT`)),by=ANNEE]

table(vege_plot_data$ANNEE,vege_plot_data$SOURCE)

table(vege_plot_data$ANNEE,vege_plot_data$Accompagnateur)


flora_survey<-fread("C:/Users/borderieux/Desktop/docs/reconciliation_ecoplant_sophy_ancien_ifn/bd_vosges/flora_vosges.csv")

sp_it_climplant<-fread(file.path("Raw_data","CIT_data","climplant_names_trait_V1.2.csv"))
flora_survey<-merge(flora_survey,sp_it_climplant,by.x="Nom initial",by.y="Nom initial",all.x=T)
sum(is.na(vege_plot_data$YearMeanMean))
nrow(vege_plot_data)



table(vege_plot_data$Vallee)
vege_plot_data[ANNEE>2009,table(Vallee,`Qualité des reelvés: 1 bon, 2 mauvaise`)]
vege_plot_data[ANNEE>2009,table(Vallee,`precision localisation`)]
vege_plot_data[ANNEE>2009,table(Vallee,SOURCE)]


cti_climplant_vosges<-flora_survey[! lb_nom_final %in% c("Picea abies","Abies alba","Fagus sylvatica","Acer campestre","Pseudotsuga menziesii","Acer pseudoplatanus","Acer platanoides"),## remove of the tree species and the main woody species
                            .(cit_climplant=mean(YearMeanMean,na.rm=T),
                              min_optimum=min(YearMeanMean,na.rm=T),
                              median_climplant=mean(YearMeanMedian,na.rm=T),
                              cit_tmax_climplant=mean(YearMaxMean,na.rm=T),
                              cit_climplant_05=mean(YearMean05,na.rm=T),
                              cit_climplant_95=mean(YearMean95,na.rm=T),
                              cit_ecoplant=mean(topt,na.rm=T),
                              cit_ecoplant_picq=mean(topt_picq,na.rm=T),
                              n_sp_climplant=sum(!is.na(YearMeanMean)),
                              mean_area=mean(Area,na.rm=T),
                              freq_for_mean=mean(indFor_Freq,na.rm=T),
                              indfor_Chytry=mean(indFor_Chytry,na.rm=T),
                              mean_azote=mean(azote,na.rm=T),
                              mean_N=mean(Nopt,na.rm=T),
                              mean_R=mean(R_ellenberg,na.rm=T),
                              mean_pH=mean(pHopt,na.rm=T),
                              mean_CN=mean(vi_CN,na.rm=T),
                              mean_L=mean(Li,na.rm=T)),by=ID]

vege_plot_data<-merge(vege_plot_data,cti_climplant_vosges,by="ID",all.x=T)
vege_plot_data_sf<-st_as_sf(vege_plot_data,coords = c("X_L93","Y_L93"),crs=st_crs(2154))


mapview(vege_plot_data_sf[vege_plot_data_sf$Vallee=="Thur",],zcol="ANNEE")


ggplot(vege_plot_data,aes(x=ANNEE,y=cit_climplant))+theme_bw()+geom_point()+geom_smooth(method = "lm")
ggplot(vege_plot_data,aes(x=ANNEE,y=`altitude MNT`))+theme_bw()+geom_point()+geom_smooth(method = "lm")


summary(lm(cit_climplant~ANNEE+`altitude MNT`,data=vege_plot_data))

#### resurveyd plot

revisite_flora_survey<-fread("C:/Users/borderieux/Desktop/docs/reconciliation_ecoplant_sophy_ancien_ifn/bd_vosges/revisite_2022_flora.csv")
revisite_flora_survey<-revisite_flora_survey[Releve!="",]
revisite_flora_survey<-merge(revisite_flora_survey,sp_it_climplant,by.x="Espece",by.y="Nom initial",all.x=T)

cti_climplant_revisite<-revisite_flora_survey[! lb_nom_final %in% c("Picea abies","Abies alba","Fagus sylvatica","Acer campestre","Pseudotsuga menziesii","Acer pseudoplatanus","Acer platanoides"),## remove of the tree species and the main woody species
                                   .(cit_climplant=mean(YearMeanMean,na.rm=T),
                                     min_optimum=min(YearMeanMean,na.rm=T),
                                     median_climplant=mean(YearMeanMedian,na.rm=T),
                                     cit_tmax_climplant=mean(YearMaxMean,na.rm=T),
                                     cit_climplant_05=mean(YearMean05,na.rm=T),
                                     cit_climplant_95=mean(YearMean95,na.rm=T),
                                     cit_ecoplant=mean(topt,na.rm=T),
                                     cit_ecoplant_picq=mean(topt_picq,na.rm=T),
                                     n_sp_climplant=sum(!is.na(YearMeanMean)),
                                     mean_area=mean(Area,na.rm=T),
                                     freq_for_mean=mean(indFor_Freq,na.rm=T),
                                     indfor_Chytry=mean(indFor_Chytry,na.rm=T),
                                     mean_azote=mean(azote,na.rm=T),
                                     mean_N=mean(Nopt,na.rm=T),
                                     mean_R=mean(R_ellenberg,na.rm=T),
                                     mean_pH=mean(pHopt,na.rm=T),
                                     mean_CN=mean(vi_CN,na.rm=T),
                                     mean_L=mean(Li,na.rm=T)),by=Releve]


revisite_plot_data<-fread("C:/Users/borderieux/Desktop/docs/reconciliation_ecoplant_sophy_ancien_ifn/bd_vosges/revisite_2022_plot_data.csv")
revisite_plot_data<-revisite_plot_data[,c(1,2,7,8,9,10,11,12,13)]
revisite_plot_data<-merge(revisite_plot_data,cti_climplant_revisite,by.x="Id releve 2022",by.y="Releve",all.x=T)

revisite_plot_data<-merge(revisite_plot_data,cti_climplant_vosges,by.x="Id Releve ancien",by.y="ID",all.x=T,suffixes =c("_recent","_past"))
revisite_plot_data[,delta_cti_climplant:=cit_climplant_recent- cit_climplant_past]
revisite_plot_data[,delta_cti_ecoplant:=cit_ecoplant_recent- cit_ecoplant_past]
revisite_plot_data[,delta_cti_ecoplant_picq:=cit_ecoplant_picq_recent- cit_ecoplant_picq_past]


#### mapping for representation ####

dt_raster<-data.table(as.data.frame(stack_all_variable_vallee))
dt_raster<-cbind(dt_raster,xyFromCell(stack_all_variable_vallee$mnt_25_vosges,1:length(stack_all_variable_vallee$mnt_25_vosges)))
dt_raster[,temp_ventr_station:=15.06]
dt_raster[,temp_ventr_station_corrected:=temp_ventr_station + ( coef_grad_elevation*(mnt_25_vosges-1210))]

dt_raster[,max_ventron:=21.20475]
dt_raster[,temp_ventr_max_corrected:= max_ventron + ( coef_grad_elevation_max*(mnt_25_vosges-1210))]

dt_raster[,bdv2_25_vosges:=as.character(bdv2_25_vosges)]
dt_raster[,bdv2_25_vosges:=ifelse(bdv2_25_vosges=="4","3",bdv2_25_vosges)]
dt_raster[,predict_microclim:=predict(daily_model,dt_raster,re.form=NA)]
dt_raster[,predict_mesoclim:=predict(daily_meso_model,dt_raster,re.form=NA)]
dt_raster[,predict_offset:=predict(daily_offset,dt_raster,re.form=NA)]

ggplot(dt_raster,aes(x,y,fill=bdv2_25_vosges))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_discrete(na.value ="transparent" )

ggplot(dt_raster,aes(x,y,fill=temp_ventr_station_corrected))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient(na.value ="transparent" )

ggplot(dt_raster,aes(x,y,fill=predict_mesoclim))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient(na.value ="transparent" )
ggplot(dt_raster[predict_microclim<20,],aes(x,y,fill=predict_microclim))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient(na.value ="transparent" )

ggplot(dt_raster,aes(x,y,fill=predict_microclim-predict_mesoclim))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )
ggplot(dt_raster,aes(x,y,fill=predict_microclim-temp_ventr_station_corrected))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )+labs(fill="")
ggplot(dt_raster,aes(x,y,fill=predict_microclim-temp_ventr_max_corrected))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )
ggplot(dt_raster,aes(x,y,fill=ipv_1000_25))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )

ggplot(dt_raster,aes(x,y,fill=predict_offset))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )

ggplot(dt_raster,aes(x=predict_microclim-temp_ventr_station_corrected,y=predict_offset))+theme_bw()+coord_fixed()+geom_point(size=0.5)+geom_smooth()

                                                                           
ggplot(vege_plot_data_sf)+geom_raster(data=dt_raster,aes(x,y,fill=predict_microclim))+geom_sf(pch=3)+theme_bw()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent")
ggplot(vege_plot_data_sf)+geom_raster(data=dt_raster,aes(x,y,fill=predict_microclim-predict_mesoclim))+geom_sf(pch=3)+theme_bw()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent")

### data visu
dt_raster[,difference_clim:=predict_microclim-temp_ventr_station_corrected]
offset_map<-ggplot(dt_raster[difference_clim< (-min(difference_clim,na.rm=T)),],aes(x,y,fill=difference_clim))+theme_bw()+geom_raster()+coord_fixed()+scale_fill_gradient2(low = "darkblue",mid="white",high="firebrick",na.value ="transparent" )+theme(panel.background = element_rect(fill = "gray98"))+labs(title="Fraicheur du microclimat   avril - juin 2022",fill="t°")
micro_map<-ggplot(dt_raster[predict_microclim<quantile(predict_microclim,probs=0.95,na.rm=T),],aes(x,y,fill=predict_microclim))+theme_bw()+geom_raster()+coord_fixed()+theme(panel.background = element_rect(fill = "gray98"))+scale_fill_viridis_c(na.value = "transparent")+labs(title="Microclimat (sous la forêt)  mai - aout 2022",fill="t°")
mmeso_map<-ggplot(dt_raster[temp_ventr_station_corrected<quantile(temp_ventr_station_corrected,probs=0.95,na.rm=T),],aes(x,y,fill=temp_ventr_station_corrected))+theme_bw()+geom_raster()+coord_fixed()+theme(panel.background = element_rect(fill = "gray98"))+scale_fill_viridis_c(na.value = "transparent")+labs(title="Climat (T) plein air  avril - juin 2022",fill="t°")
ggsave('figures/microclim/map_micro.png',micro_map,units = "in",height =5 ,width = 5,dpi=300)

htmlwidgets::saveWidget(ggplotly(micro_map),"widget_communication/microclimat.html")
htmlwidgets::saveWidget(ggplotly(mmeso_map),"widget_communication/mesoclimat.html")
htmlwidgets::saveWidget(ggplotly(offset_map),"widget_communication/offset_map.html")

###



# raster_pred_meso<-raster(stack_all_variable_vallee$mnt_25_vosges)
# raster_pred_meso<-setValues(raster_pred_meso,dt_raster$predict_mesoclim)

raster_pred_micro<-raster(stack_all_variable_vallee$mnt_25_vosges)
raster_pred_micro<-setValues(raster_pred_meso,dt_raster$predict_microclim)

# raster_pred_offset<-raster(stack_all_variable_vallee$mnt_25_vosges)
# raster_pred_offset<-setValues(raster_pred_meso,dt_raster$predict_offset)
# 

# stack_predictions<-stack(list(pred_micro=raster_pred_micro,pred_meso=raster_pred_meso,
#                               dif_pred=raster_pred_micro-raster_pred_meso,pred_offset=raster_pred_offset))

stack_predictions<-stack(list(pred_micro=raster_pred_micro))


vege_pred<-data.table(extract(stack_predictions,vege_plot_data_sf))

vege_pred<-cbind(vege_pred,extract(stack_all_variable_vallee,vege_plot_data_sf))
vege_pred[,temp_ventr_station:= 20]
vege_pred[,temp_ventr_station_corrected:= temp_ventr_station + ( coef_grad_elevation*(mnt_25_vosges-1210))]
#### link vege and pred of microclim ####
vege_plot_data_pred<-cbind(vege_plot_data,vege_pred)
vege_plot_data_pred<-vege_plot_data_pred[!`precision localisation`%in%c(3),]
vege_plot_data_pred<-vege_plot_data_pred[pred_micro<30,]
vege_plot_data_pred[,dif_pred:=pred_micro - temp_ventr_station_corrected]

ggplot(vege_plot_data_pred,aes(x=pred_micro,y=cit_climplant))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(vege_plot_data_pred,aes(x=temp_ventr_station_corrected,y=cit_climplant))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(vege_plot_data_pred,aes(x=dif_pred,y=cit_climplant))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(vege_plot_data_pred[pred_offset<1,],aes(x=pred_offset,y=cit_climplant))+theme_bw()+geom_point()+geom_smooth(method="lm")
# cit_climplant

raster_pred_offset
ggplot(vege_plot_data_pred,aes(x=aspect_NS_25_vosges,y=cit_climplant))+theme_bw()+geom_point()+geom_smooth(method="lm")


summary(lm(cit_climplant~temp_ventr_station_corrected,data=vege_plot_data_pred))
summary(lm(cit_climplant~pred_micro,data=vege_plot_data_pred))
summary(lm(cit_climplant~mnt_25_vosges+aspect_NS_25_vosges,data=vege_plot_data_pred))

summary(lm(cit_ecoplant_picq~I(pred_micro-temp_ventr_station_corrected),data=vege_plot_data_pred))

summary(lm(cit_ecoplant_picq~dif_pred,data=vege_plot_data_pred))


hist(vege_pred$pred_micro,nc=20)
hist(vege_pred$dif_pred,nc=20)





#### get daily climate ####

library(easyclimate)
coords_sf_84<-st_transform(coords_sf,crs=st_crs(4326))

era5_test<-get_daily_climate(coords_sf_84,c("Tmax","Tmin"),period="2015-01-01:2019-08-15")

era5_test_rast<-get_daily_climate(coords_sf_84,c("Tmax","Tmin"),period="2015-01-01:2019-08-15","raster")






#### save and load ####
save.image("microclim_vosges.RData")


#### export for grand Ventron ####

reserve_thermo<-c("H_N_8_2","H_S_6","H_N_7","slope_L_6","H_N_6","slope_L_5","L_N_5","cad_L_5")

all_tms_data_ventron<-all_tms_data[site_ID%in%reserve_thermo,c( "site_ID" , "sensor_ID" ,"date","day", "T1","T2", "T3" ,"soil_moisture")]

all_tms_data_day_ventron<-all_tms_data_day[site_ID%in%reserve_thermo,c( "site_ID" , "sensor_ID" ,"day", "T1","T2", "T3" ,"soil_moisture","amplitude","max_T3","in_soil","temp_ventr_station","temp_ventr_station_corrected")]
all_tms_data_day_ventron[,effet_tampon:=T3-temp_ventr_station_corrected]

all_tms_data_day_ventron[,mean(effet_tampon,na.rm=T),by=site_ID]
effet_tempon_ete_22<-all_tms_data_day_ventron[day%between% study_period,.(effet_tempon_ete_22=mean(effet_tampon,na.rm=T)),by=site_ID]


ggplot(all_tms_data_day_ventron,aes(x=day,y=effet_tampon,color=site_ID))+
  theme_bw()+
  geom_line()+
  geom_hline(yintercept = 0,lty=2)+
  coord_cartesian(ylim=c(-5,5))#+geom_line(data=Ventron_station_day,aes(y=temp_ventr_station,color="Ventron_station"))

Ventron_station_day

coords_ventron<-coords[site_ID%in%reserve_thermo,c("site_ID" ,"x","y","xl93","yl93","mnt_25_vosges")]
coords_ventron<-merge(coords_ventron,effet_tempon_ete_22,all.x=T)
coords_ventron[,fonctionel:=as.numeric(!is.na(effet_tempon_ete_22))]

write.table(coords_ventron,"thermometre_metadonnee.csv",sep = ";")
write.table(Ventron_station_day,"temperature_station_ventron.csv",sep = ";")
write.table(all_tms_data_ventron,"temperature_brute.csv",sep = ";")
write.table(all_tms_data_day_ventron,"temperature_jour.csv",sep = ";")


library(mapview)

coords_ventron_sf<-st_as_sf(coords_ventron,coords = c("xl93","yl93"),crs=st_crs(2154))

mapview(coords_ventron_sf,zcol="effet_tempon_ete_22")

mapshot(mapview(coords_ventron_sf,zcol="effet_tempon_ete_22"),url="carte_thermometre.html")
