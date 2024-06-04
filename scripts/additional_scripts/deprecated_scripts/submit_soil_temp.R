#### soiltemp formating ####

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
all_tms_data<-merge(all_tms_data,meta_data[,c("site_ID" ,"sensor_ID","date_setup","put_out_date_1")],by.x="sensor_ID",by.y="sensor_ID",all.x=T)


meta_data<-merge(meta_data,coords,by.x="site_ID",by.y="site_ID")

all_tms_data[,date:=ymd_hm(date)]
all_tms_data[,date:=with_tz(date,"Europe/Paris")]

all_tms_data[,day:=floor_date(date,"day")]
all_tms_data[,month:=floor_date(date,"month")]
all_tms_data[,year:=floor_date(date,"year")]

all_tms_data[,t1_later:=c(NA,T1[1:(.N-1)]),by=site_ID]
all_tms_data[,var_t1:=abs(T1-t1_later),by=site_ID]

all_tms_data[,t2_later:=c(NA,T2[1:(.N-1)]),by=site_ID]
all_tms_data[,var_t2:=abs(T2-t2_later),by=site_ID]


all_tms_data[,not_in_field:=floor_date(date,unit="day")<=date_setup]
all_tms_data[,out_of_soil:=ifelse(is.na(put_out_date_1),FALSE,floor_date(date,unit="day")>put_out_date_1)]


all_tms_data<-all_tms_data[not_in_field==0,]

all_tms_data<-all_tms_data[out_of_soil==FALSE,]



all_tms_data[,last_recorded_day:=max(day),by=sensor_ID]
all_tms_data<-all_tms_data[day!=last_recorded_day,] ## the last day of a sensor is by definition incomplete

table(all_tms_data$sensor_ID)/96
all_tms_data$sensor_ID
meta_data
View(all_tms_data[sensor_ID==94207944,])
all_tms_data[,.(min(date),max(date)),by=sensor_ID]


### remove measurements of intensive loggers ( more than the usual four 00 15 30 45 measurements)
all_tms_data[,minutes:=minute(date)]
all_tms_data<-all_tms_data[minutes%in%c(0,15,30,45),]
all_tms_data[,minutes:=NULL]




too_short<-c(94207928 ,94207929,94207936)

### some loggers were replaced with others, removing the erroneous reading
ymd("2022/07/06")
meta_data[site_ID=="cad_L_8",]
all_tms_data<-all_tms_data[!(sensor_ID=="94195464" & day<=ymd("2022/07/06")),]


### two special case : 94207915 and  94207913
### 94207915 was in a building for several days, then it was planted again

suspicious_loggers<-c("L_S_1","L_S_3","slope_H_1","L_S_7")
new_loggers<-c("new1")
broken_loggers<-c("cad_L_8")

all_tms_data[site_ID=="cad_L_8",]

all_tms_data<-all_tms_data[!site_ID%in%suspicious_loggers,]


all_tms_data[,in_soil:=1]


plotly::ggplotly( ggplot(all_tms_data_day[site_ID=="H_N_2" ,],aes(x=day,y=var_extreme_errflag_t1,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2))
date_out<-ymd(c("2021/08/27","2022/07/13"))
date_int<-ymd("2022/01/19")
all_tms_data[site_ID=="H_N_2"& day>date_out[1] ,in_soil:=0]
all_tms_data[site_ID=="H_N_2"& day>date_int[1] ,in_soil:=1]
all_tms_data[site_ID=="H_N_2"& day>date_out[2],in_soil:=0]
plotly::ggplotly( ggplot(all_tms_data_day[site_ID=="H_S_1" ,],aes(x=day,y=var_extreme_errflag_t1,color=site_ID))+theme_bw()+geom_line(lwd=1)+facet_wrap(~in_soil_1,nrow = 3)+geom_hline(yintercept = 0.4,lty=2))
date_out<-ymd(c("2021/06/03"))
date_int<-ymd("2022/06/04")
all_tms_data[site_ID=="H_S_1" ,in_soil:=1]
all_tms_data[site_ID=="H_S_1"& day>date_out[1] ,in_soil:=0]
all_tms_data[site_ID=="H_S_1"& day>date_int[1],in_soil:=1]

all_tms_data[in_soil==0,"T1"]<-NA
all_tms_data[in_soil==0,"T2"]<-NA
all_tms_data[in_soil==0,"T3"]<-NA
all_tms_data[in_soil==0,"T1"]<-NA

all_tms_data[,time_h:=hour(date)]
all_tms_data[,time_m:=minute(date)]

all_tms_data[,time_hm:=substr(as.character(date),12,16)]

nchar("2022-06-19 09:00:00")

all_tms_data_export<-all_tms_data[,c("sensor_ID","year","month","day","time_hm","T1","T2","T3","soil_moisture","date")]

colnames(all_tms_data_export)<-c("Raw_data_identifier","Year","Month","Day","Time (24h)","T1","T2","T3","Soil_moisture_raw","date")

all_tms_data_export[,Year:=year((date))]
all_tms_data_export[,Month:=month((date))]
all_tms_data_export[,Day:=day((date))]

all_tms_data_export[,date:=NULL]

write.table(all_tms_data_export,"RAW-TIME-SERIES-DATA.csv",row.names = F,col.names = T,sep=";")

#### meta data ######

meta_data_export<-meta_data[sensor_ID%in%unique(all_tms_data_export$Raw_data_identifier),]

meta_data_export<-meta_data_export[,c("sensor_ID","site_ID","x.x","y.x","canopy_cover_glama","mnt_25_vosges.x","hat_t3_1","put_out_date_1")]

meta_data_export<-rbind(meta_data_export,meta_data_export,meta_data_export,meta_data_export)
meta_data_export<-meta_data_export[order(sensor_ID),]
meta_data_export[,Sensor_code:=rep(c("T1","T2","T3","Soil_moisture_raw"),50)]

meta_data_export[,comment:=ifelse(is.na(put_out_date_1),"","This time series is truncated, the logger has been removed out of the soil by wild animals, we removed the invalid readings with the sudden change of daily T1 variation")]
meta_data_export[,comment2:=ifelse(hat_t3_1==0,paste0("This logger was found without T3 sunshield"),comment)]

meta_data_export[,comment2:=ifelse(hat_t3_1==0  &  Sensor_code =="T3",paste0("This logger was found without T3 sunshield"),"")]

table((meta_data_export$hat_t3_1))

sum(is.na(meta_data_export$hat_t3_1))

write.table(meta_data_export,"meta_data_soil_temp.csv",row.names = F,col.names = T,sep=";")

meta_data_export<-merge(meta_data_export,
                        all_tms_data[,.(min=min(date),max=max(date)),by=sensor_ID],by="sensor_ID")


meta_data_export[,min_y:=year(min)]
meta_data_export[,min_m:=month(min)]
meta_data_export[,min_d:=day(min)]


meta_data_export[,max_y:=year(max)]
meta_data_export[,max_m:=month(max)]
meta_data_export[,max_d:=day(max)]


meta_soil_temp<-fread("PEOPLE.csv")

write.table(meta_soil_temp,"PEOPLE.csv",row.names = F,col.names = T,sep=";")



meta_soil_temp<-fread("METADATA.csv")

write.table(meta_soil_temp,"METADATA.csv",row.names = F,col.names = T,sep=";")
