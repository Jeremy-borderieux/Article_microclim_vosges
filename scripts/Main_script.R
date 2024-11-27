#### Packages ####
library(data.table)
library(sf)
library(stringr)
library(raster)
library(ggspatial)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(lubridate)
library(car)
library(modEvA)

#### Spatial covariable and logger meta data reading ####
meta_data<-fread(file.path("data","metadata_loggers","meta_data_2022_points.csv"),dec = ",")

## covariables
stack_all_variable_vallee<-stack(file.path("data","processed_raster",list.files(file.path("data","processed_raster"))))

## get copernicus data
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

## extraction of the covariables
env_var<- extract(stack_all_variable_vallee,coords_sf)
tree_density_2018<- extract(tree_density,coords_sf_laea)
env_var<-cbind(env_var,tree_density_2018);rm(tree_density_2018)

coords<-cbind(coords,env_var)
meta_data<-merge(meta_data,coords,by="site_ID")

## canopy cover in a 25m radius data
cover_25m_sites<-fread(file.path("data","metadata_loggers","cover_25m_sites.csv"))
cover_25m_sites_aggr<-cover_25m_sites[,.(canopy_cover_25m=sum(free_cover),under_canopy_cover=sum(under_cano_cover),broadleaved_prop=(sum(ifelse(Broadleaved,free_cover,0))/sum(free_cover)  )),by=site_ID]
cover_25m_sites_aggr[,absolute_cover:=canopy_cover_25m+under_canopy_cover]

meta_data<-merge(meta_data,cover_25m_sites_aggr,by="site_ID",all.x=T)


#### Temperature Logger reading and cleaning with myClim ####
#install.packages("http://labgis.ibot.cas.cz/myclim/myClim_latest.tar.gz", repos=NULL, build_vignettes=TRUE)
library(myClim)


## preparing the locality table (aka metadata) for myclim, and load the path in the correct order
folder<-file.path("data","microclimate_loggers")

data_file<-grep("data_9",list.files(folder,full.names = T),value=T)
local_table<- data.table(locality_id = coords$site_ID,elevation = coords$mnt_25_vosges,lon_wgs84= coords$x ,lat_wgs84=coords$y ,tz_offset= 120)
local_table_total<-merge(local_table,meta_data[,c("site_ID","sensor_ID")],by.x="locality_id",by.y="site_ID")
local_table_total[,path_to_file:=file.path(folder,paste0("data_",sensor_ID,"_0.csv"))]
local_table_working_logs<-local_table_total[path_to_file%in%data_file,]

##reading with myclim
tms_data<-mc_read_data(files_table = data.table(path= local_table_working_logs$path_to_file,locality_id = local_table_working_logs$locality_id, data_format ="TOMST" ),localities_table = local_table_working_logs)

infos<-mc_info(tms_data)

## calibration og the t3 sensors with curves obtained from lab calibration (T3 : sensor at 12cm )
calibration_curves_T3<-fread(file.path("data","metadata_loggers","Calibration_T3.csv"))
infos<-merge(infos,calibration_curves_T3,by.x="serial_number",by.y="sensor_ID",all.x=T)

prep_for_calib<-data.frame(serial_number = infos$serial_number,
                           sensor_id = infos$sensor_id,
                           datetime = as.POSIXct("2021-05-15",tz="UTC"),
                           cor_factor = infos$offset,
                           cor_slope = infos$slope-1)

tms_load <- mc_prep_calib_load(tms_data, prep_for_calib)
tms_data <- mc_prep_calib(tms_load, sensors = "TMS_T3")

rm(tms_load)
mc_info_count(tms_data) #which returns the number of localities, loggers and sensors in myClim object
mc_info(tms_data)# returning data frame with summary per sensor
mc_info_meta(tms_data)# returning the data frame with locality metadata
mc_info_clean(tms_data) #returning the data frame with cleaning log

# min max mean temp ok, 
tms_error_flag <- mc_agg(tms_data, fun = list(TMS_T1 = "error_flag"), period = "day",
                         custom_functions = list(error_flag = function(x){mean(na.rm=T,abs(x-c(NA,x[1:(length(x-1))])))}))

tms_error_flag <- data.table(mc_reshape_long(tms_error_flag))
# errflag ok


temp_env <- mc_agg(tms_data, 
                   fun=list(TMS_T3=c("mean_t","percentile","range_t")),
                   period = "day", min_coverage = 0.95,percentiles = c(5,95),
                   custom_functions = list(mean_t=function(x) mean(x[x%between%  quantile(x,na.rm=T,probs=c(0.05,0.95))],na.rm=T) ,
                                           range_t=function(x) diff(range(x[x%between%  quantile(x,na.rm=T,probs=c(0.05,0.95))],na.rm=T) )))


## The temperature are then aggregated to the growing season of 2022, from 1st april to mid august
agg_growing_season<-mc_agg(temp_env,
                           fun=list(TMS_T3_mean_t="mean",TMS_T3_percentile95="mean",TMS_T3_percentile5="mean"),
                           period="custom",
                           custom_start ="04-01" ,
                           custom_end = "08-15")

## we transform the output from myclim to a data.table
all_tms_data_agg<-data.table(mc_reshape_long(agg_growing_season))
all_tms_data_agg[,height:=NULL]
all_tms_data_agg[,serial_number:=NULL]
all_tms_data_agg[,datetime:=ymd(datetime)]
all_tms_data_agg<-all_tms_data_agg[datetime==ymd("2022-04-01"),]
all_tms_data_agg[,datetime:=NULL]
all_tms_data_agg[,time_to:=NULL]

## renaming
all_tms_data_agg[sensor_name=="TMS_T3_mean_t_mean","sensor_name"]<- "T3_mean"
all_tms_data_agg[sensor_name=="TMS_T3_percentile5_mean","sensor_name"]<- "T3_min"
all_tms_data_agg[sensor_name=="TMS_T3_percentile95_mean","sensor_name"]<- "T3_max"

all_tms_data_agg<-data.table::dcast(all_tms_data_agg, locality_id~sensor_name)

#### Microclimate model fitting ####
loggers_to_remove<-c("L_S_1","L_S_3","slope_H_1","L_S_7", # these loggers were found out of soil in plain sun, their time series reveals an abnormal heat, safer to remove them
                     "cad_L_8", # this logger was malfunction, the new replacing it did not record temperature during the study period 
                     "new1") # this logger was not part of the original sampling, it don't have the whole growing season recorded


all_tms_data_agg_fit<-all_tms_data_agg[ !locality_id %in% loggers_to_remove,]

all_tms_data_agg_fit<-merge(all_tms_data_agg_fit,meta_data,by.x="locality_id",by.y="site_ID")
## correction of the canopy_cover_glama variable (was character)
## some data are missing, filling out with the mean
all_tms_data_agg_fit[,canopy_cover_glama:=as.numeric(str_replace_all(canopy_cover_glama,",","."))]
all_tms_data_agg_fit[,canopy_cover_glama:=ifelse(is.na(canopy_cover_glama),mean(all_tms_data_agg_fit$canopy_cover_glama,na.rm=T),canopy_cover_glama)]
all_tms_data_agg_fit[,canopy_cover_25m:=ifelse(is.na(canopy_cover_25m),mean(all_tms_data_agg_fit$canopy_cover_25m,na.rm=T),canopy_cover_25m)]

## the aggregated temperature are modeled with a simple linear regression
lm_agg_mean<-lm(T3_mean~mnt_25_vosges+Heat_load_index_25+ipv_25+tree_density_2018 , data=all_tms_data_agg_fit)
lm_agg_mean_canopy_25m<-lm(T3_mean~mnt_25_vosges+Heat_load_index_25+canopy_cover_25m+ipv_25 , data=all_tms_data_agg_fit,model=T)
lm_agg_mean_canopy_glama<-lm(T3_mean~mnt_25_vosges+Heat_load_index_25*canopy_cover_glama+ipv_25 , data=all_tms_data_agg_fit,model=T)

summary(lm_agg_mean)
summary(lm_agg_mean_canopy_25m)
summary(lm_agg_mean_canopy_glama)

lm_agg_max<-lm(T3_max~mnt_25_vosges+Heat_load_index_25+ipv_25+tree_density_2018 , data=all_tms_data_agg_fit)
lm_agg_min<-lm(T3_min~mnt_25_vosges+Heat_load_index_25+ipv_25+tree_density_2018 , data=all_tms_data_agg_fit)


##cor between max and mean temperature 
all_tms_data_agg_fit[,cor(T3_max,T3_mean,method = "pearson")]
all_tms_data_agg_fit[,plot(T3_max,T3_mean)]
## finale sampling scheme
table(str_remove(all_tms_data_agg_fit$locality_id,"_[:digit:]"))

## Small test of elevation/topography/canopy spatial autocorrelation
library(foreach)
library(terra)
stack_all_variable_vallee <- rast(stack_all_variable_vallee)

pred_topo <- stack_all_variable_vallee$Heat_load_index_25 * get_coef_mean["Heat_load_index_25",1] + stack_all_variable_vallee$ipv_25 * get_coef_mean["ipv_25",1]
names(pred_topo) <- "pred_topo"

pred_cano <- stack_all_variable_vallee$tree_density_projected * get_coef_mean["tree_density_2018",1] 
names(pred_cano) <- "pred_cano"

pred_elev <- stack_all_variable_vallee$mnt_25_vosges * get_coef_mean["mnt_25_vosges",1] 
names(pred_elev) <- "pred_elev"


stack_all_variable_vallee <- stack(stack_all_variable_vallee,pred_topo,pred_cano,pred_elev)
stack_all_variable_vallee$pred_cano <- mask(stack_all_variable_vallee$pred_cano,stack_all_variable_vallee$pred_elev)


library(usdm)
variogram_topo <- Variogram(stack_all_variable_vallee$pred_topo,lag = 125,cutoff = 4160*2)
variogram_vege <- Variogram(stack_all_variable_vallee$pred_cano,lag = 25,cutoff =2000)
variogram_vege_2 <- Variogram(stack_all_variable_vallee$pred_cano,lag = 25,cutoff =2000)

variogram_topo_2 <- Variogram(stack_all_variable_vallee$pred_topo,lag = 25,cutoff = 2000)

variogram_elev <- Variogram(stack_all_variable_vallee$pred_elev,lag = 125,cutoff = 4160*2)
variogram_elev_2 <- Variogram(stack_all_variable_vallee$pred_elev,lag = 25,cutoff =2000)

plot(variogram_topo)
plot(variogram_topo_2)
plot(variogram_vege)
plot(variogram_vege_2)

plot(variogram_elev)
plot(variogram_elev_2)

variogram_data <- rbind(variogram_vege@variogram,variogram_topo_2@variogram,variogram_elev_2@variogram)
variogram_data$class <- rep(c("Canopy cooling","Topographic effect","Lapse rate"),each=80)

variogram_data <- data.table(variogram_data)
variogram_data[,scaled_gamma := scale(gamma), by = class]
scale(variogram_data$gamma)

semivariogram <- ggplot(variogram_data,aes(x = distance, y = scaled_gamma , color = class))+
#ggplot(variogram_data,aes(x = distance, y = ifelse(class=="Lapse rate" ,gamma/5 ,gamma) , color = class))+
    theme_classic()+
  theme(legend.position = c(0.7,0.25))+
  geom_point(alpha=0.75,size=1)+
  geom_path()+
  scale_color_manual(values=c("#4EB8C2", "#5BBA58", "#6B6B6B")[c(2,3,1)])+
  labs(y = "Scaled semivariance", x= "Distance (m)", color ="Predictors")
  
ggsave(file.path("figure_result","semivariogram.png"),semivariogram,width = 180,height = 110,unit="mm",dpi=300)

# 
# 
# autocor_raster <- foreach( window = c(seq(from = 3 , to = 31, by=2),41,51),.combine = rbind)%do%{
#   f <- matrix(1, nrow=window,ncol=window)
#   cat(paste0(window," / "))
#   a1 <- autocor(stack_all_variable_vallee$pred_topo,f, method="moran")
#   a2 <- autocor(stack_all_variable_vallee$tree_density_projected,f, method="moran")
#   a3 <- autocor(stack_all_variable_vallee$mnt_25_vosges,f, method="moran")
#   
#   return(data.table(window=window, moran_pred_topo = a1,moran_canopy =a2 , moran_lapse = a3))
#   
# }
# autocor_raster[,window_meter := window * 25]
# autocor_plot <- ggplot(autocor_raster[window< 52,])+
#   theme_classic()+
#   geom_point(mapping=aes(x = window_meter, y=moran_pred_topo , color= "Topographic effect"),size = 2.5)+
#   geom_point(mapping=aes(x = window_meter, y=moran_canopy , color= "Canopy cooling"),size = 2.5)+
#   geom_point(mapping=aes(x = window_meter, y=moran_lapse , color= "Lapse rate"),size = 2.5)+
#   
#   geom_line(mapping=aes(x = window_meter, y=moran_pred_topo , color= "Topographic effect"),alpha = 0.5)+
#   geom_line(mapping=aes(x = window_meter, y=moran_canopy , color= "Canopy cooling"),alpha = 0.5)+
#   geom_line(mapping=aes(x = window_meter, y=moran_lapse , color= "Lapse rate"),alpha = 0.5)+
# 
#   geom_hline(yintercept = 1,col="grey",lty=2)+
#   theme(legend.position = c(0.2,0.75))+
#   labs(x="Focal distance (m)",y="Global Moran I",color="Predictions")+
#   scale_color_manual(values=c("#4EB8C2", "#5BBA58", "#6B6B6B")[c(2,3,1)])
# 
# ggsave(file.path("figure_result","global_moranI.png"),autocor_plot,width = 180,height = 130,unit="mm",dpi=300)

## a function to format the model estimates for publication
create_table_model<-function(model,dovarpart=T,interaction=F,canopy_variable="tree_density_2018"){

  get_coef<-summary(model)$coefficients
  if(canopy_variable=="canopy_cover_25m")get_coef<-get_coef[c(1,2,3,5,4),]
  if(canopy_variable=="canopy_cover_glama")get_coef<-get_coef[c(1,2,3,5,4,6),]
get_range<-function(x){
  paste0(round(range(x),3),collapse=" : ")
  
}  
all_tms_data_agg_fit<-model$model

## ompute the range of predictor
range_predictor<-c(NA,apply(all_tms_data_agg_fit[,c("mnt_25_vosges","Heat_load_index_25","ipv_25",canopy_variable)],2,function(x)diff(range(x))))
range_predictor_char<-c(NA,apply(all_tms_data_agg_fit[,c("mnt_25_vosges","Heat_load_index_25","ipv_25",canopy_variable)],2,get_range))

sd_predictor <- c(NA,apply(all_tms_data_agg_fit[,c("mnt_25_vosges","Heat_load_index_25","ipv_25",canopy_variable)],2,sd))

export_table_model<-data.table(get_coef,keep.rownames=T)
export_table_model[,range:=if (interaction) c(range_predictor,NA) else range_predictor]
export_table_model[,sd:=if (interaction) c(sd_predictor,NA) else sd_predictor]

export_table_model[,effect_size:=sd*Estimate]
export_table_model[,sd:=NULL]



if(dovarpart){
## variance partitionning
grp<-data.frame(var=c("mnt_25_vosges","Heat_load_index_25","ipv_25",canopy_variable),group=c("Elevation","Topo","Topo","Canopy"))

coef_varpart<-as.vector(t(varPart(model=glm(model),groups = grp)))
names(coef_varpart)<-colnames(t(varPart(model=glm(model),groups = grp)))

var_elev<-coef_varpart["Elevation"] + coef_varpart["Elevation_Topo"]/2 + coef_varpart["Elevation_Canopy"]/2  +   coef_varpart["Elevation_Topo_Canopy"]/3
var_topo<- coef_varpart["Topo"] + coef_varpart["Elevation_Topo"]/2 + coef_varpart["Topo_Canopy"]/2  +   coef_varpart["Elevation_Topo_Canopy"]/3
var_canopy<- coef_varpart["Canopy"] + coef_varpart["Topo_Canopy"]/2 + coef_varpart["Elevation_Canopy"]/2  +   coef_varpart["Elevation_Topo_Canopy"]/3

export_table_model[,variance:= c(NA,var_elev,var_topo,var_topo,var_canopy)*100]

}
export_table_model<-export_table_model[,-c("t value")]

export_table_model<-cbind(export_table_model[,"rn"],signif(export_table_model[,-"rn"],3))
export_table_model[,`Pr(>|t|)`:=ifelse(`Pr(>|t|)`<10e-4,"<10-4",`Pr(>|t|)`)]
export_table_model[,type:=if(interaction)c(NA,"Elevation","Topography","Topography","Canopy",NA) else c(NA,"Elevation","Topography","Topography","Canopy")]

export_table_model[,range:= if(interaction) c(range_predictor_char,NA) else range_predictor_char ]

if(dovarpart) export_table_model<-export_table_model[,c(1,8,2,3,5,6,7,4)] else export_table_model<-export_table_model[,c(1,7,2,3,5,6,4)]
export_table_model

}


export_table_model<-create_table_model(lm_agg_mean)
export_table_model_max<-create_table_model(lm_agg_max)

export_table_model_cano_1<-create_table_model(lm_agg_mean_canopy_25m,F,F,"canopy_cover_25m")
export_table_model_cano_2<-create_table_model(lm_agg_mean_canopy_glama,F,T,"canopy_cover_glama")

write.table(export_table_model,file.path("figure_result","Model_Tmean.csv"),row.names = F,sep=";",dec=".")
write.table(export_table_model_max,file.path("figure_result","Model_Tmax.csv"),row.names = F,sep=";",dec=",")

write.table(export_table_model_cano_1,file.path("figure_result","Model_canopy_25m.csv"),row.names = F,sep=";",dec=",")
write.table(export_table_model_cano_2,file.path("figure_result","Model_canopy_glama.csv"),row.names = F,sep=";",dec=",")

## Please provide plots to better justify the expected linearity between temperature and topography or canopy closure.
## visualization of the linearity of the temperature predictors 

all_tms_data_agg_fit_melt <- all_tms_data_agg_fit[,c("mnt_25_vosges","Heat_load_index_25","ipv_25","tree_density_2018","T3_mean","T3_max")]
colnames(all_tms_data_agg_fit_melt)[1:4] <- c("Elevation (m.a.s.l)","Heat load index","Topographic position","Tree density (%)")
all_tms_data_agg_fit_melt <- melt(all_tms_data_agg_fit_melt,measure.vars = c("Elevation (m.a.s.l)","Heat load index","Topographic position","Tree density (%)"))

plot_linear_pred_mean <- ggplot(all_tms_data_agg_fit_melt,aes(x=value,y=T3_mean))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "lm",color='orange',fill="orange")+
  geom_smooth()+
  facet_wrap(~variable,scales ="free_x")+
  labs(y="Mean understory temperature °C",x="Predictor")

plot_linear_pred_max<- ggplot(all_tms_data_agg_fit_melt,aes(x=value,y=T3_max))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "lm",color='orange',fill="orange")+
  geom_smooth()+
  facet_wrap(~variable,scales ="free_x")+
  labs(y="Max understory temperature °C",x="Predictor")
## visual check of the residuals

ggsave(file.path("Figure_result","linear_predictor_mean_temp.jpg"),
       ggarrange(plot_linear_pred_mean,plot_linear_pred_max,labels = c("a)","b)"),align="hv",nrow=2,common.legend = T),
       dpi=250,width = 180,height=210,unit="mm")


#### vegetation data reading ####
vege_plot_meta_data<-fread(file.path("data","flora_data_metadata","plot_data_2022_final.csv"))
vege_plot_meta_data_sf<-st_as_sf(vege_plot_meta_data,coords=c("X_L93","Y_L93"),crs=st_crs(2154))

## extract from the same used for microclimate various env variables
vege_covariable<-data.table(extract(stack_all_variable_vallee,vege_plot_meta_data_sf))
vege_covariable$plot_ID<-vege_plot_meta_data$plot_ID

vege_covariable$tree_density_projected<-extract(tree_density,vege_plot_meta_data_sf)

vege_covariable<-vege_covariable[!is.na(mnt_25_vosges),]

vege_plot_meta_data<-merge(vege_plot_meta_data,vege_covariable,by="plot_ID")

## database of species attributes (thermal otpimum, pH optimum etc...)
indicator_value<-fread(file.path("data","climplant_names_trait_V1.2.csv"),encoding = "Latin-1")

## eurforplant reading (list of habitat preference of plant species)
EurForPlant<-fread(file.path("data","EurForPlant_melt.csv")) 
EurForPlant<-merge(EurForPlant[biogeo_region=="France_mountains",],indicator_value[,c("species","lb_nom_final")],by="species",all.x=T)
EurForPlant[,nom_espece:=lb_nom_final]
EurForPlant[,lb_nom_final:=NULL]
EurForPlant<-EurForPlant[!is.na(nom_espece),]

more_than_one_sp<-EurForPlant[,.N,by=nom_espece]

# some species have synonymes in the list because of subspecies 
selected_sub_sp<-c("Aconitum lycoctonum L. subsp. lycoctonum","Adenostyles alliariae (Gouan) Kern.","Arabidopsis arenosa (L.) Lawalrée subsp. arenosa",
                   "Asplenium adiantum-nigrum L. subsp. adiantum-nigrum","Carex divulsa Stokes subsp. divulsa","Carex muricata L. subsp. muricata","Centaurea phrygia subsp. pseudophrygia (C. A. Mey.) Gugler",
                   "Cerastium fontanum Baumg. subsp. fontanum","Dactylis glomerata L. subsp. glomerata","Juniperus communis L. subsp. communis",
                   "Lamium galeobdolon (L.) Crantz subsp. galeobdolon","Luzula sylvatica (Huds.) Gaudin subsp. sylvatica","Pinus nigra J. F. Arnold subsp. nigra",
                   "Pyrus communis L. subsp. communis","Salix cinerea L.","Stellaria nemorum L. subsp. nemorum",
                   "Viola canina L.","Viscum album L. subsp. album","Dactylis glomerata L. subsp. glomerata",
                   "Adenostyles alliariae (Gouan) Kern.","Asplenium adiantum-nigrum L. subsp. adiantum-nigrum",
                   "Carex divulsa Stokes subsp. divulsa","Stellaria nemorum L. subsp. nemorum")
EurForPlant[,more_than_one_subsp:=nom_espece%in%more_than_one_sp[N>=2]$nom_espece]
EurForPlant<-EurForPlant[more_than_one_subsp==FALSE| Scientific_name%in%selected_sub_sp,]
EurForPlant[,forest_species:= ifelse(habitat_categ%in%c("1,1","1,2"),"Forest_species",ifelse(habitat_categ%in%c("2,1","2,2","O"),"Open_species","NC"))]
EurForPlant
# 1.1: Taxa which can be found mainly in the closed forest
# 1.2: Taxa which are mainly typical along forest edges and in forest openings. This includes species of forest edges, species which mainly occur in windthrow, burned or clear-cut areas, or during the regeneration phase after such events, species which mainly occur on exploitation roads and unpaved forest paths, species which are restricted to open forests due to extreme site conditions 
# 2.1: Taxa which can be found in forest as well as open vegetation
# 2.2: Taxa which can be found partly in forest, mainly in open vegetation
# O: Taxa which can be ± only found in open vegetation
# nk: Habitat preference not known
# ext: Taxon extinct
# / : Taxon does not occur

## we append the database of otpima with euforplant
indicator_value<-merge(indicator_value,EurForPlant[,c("nom_espece","Family","Layer","forest_species","habitat_categ")],all.x=T,by.x="lb_nom_final",by.y="nom_espece")

## reading of the actual surveys
flora_survey<-fread(file.path("data","flora_data_metadata","flora_vosges_2022_tax_hom.csv"))
flora_survey<-flora_survey[plot_ID%in%vege_plot_meta_data$plot_ID,]
## adding the species attributes
flora_survey<-merge(flora_survey,indicator_value,all.x=T,by.x="species_name",by.y="lb_nom_final")
flora_survey<-flora_survey[tree!=1,]# removing the trees

## function to create a [site_id,species] matrix of absence-presence from the flora survey
create_table_sp<-function(survey,id_names="idp"){
  
  table_survey<-table(survey[,get(id_names)],survey$species_name)
  table_survey<-as.data.frame.matrix(table_survey)
  
  di<-dim(table_survey)
  sp<-colnames(table_survey)
  id<-rownames(table_survey)
  
  table_survey<-as.matrix(table_survey)
  table_survey<-as.numeric(table_survey)
  table_survey<-matrix(table_survey,nrow=di[1],ncol =di[2] )
  table_survey<-as.data.frame(table_survey)
  
  rownames(table_survey)<-id
  colnames(table_survey)<-sp
  
  return(table_survey)
}

flora_survey<-flora_survey[!grepl("subsp.",species_name),]
table_sp_vosges<-create_table_sp(flora_survey,"plot_ID")
table_sp_vosges[1:6,1:6]
dim(table_sp_vosges)

## create plot scale averages  # cit = Community thermal index , RS = species richness
cti_indicator<-flora_survey[,
                             .(cit_climplant=mean(YearMeanMean,na.rm=T),
                               sd_cit_climplant=sd(YearMeanMean,na.rm=T),
                               median_climplant=mean(YearMeanMedian,na.rm=T),
                               cit_tmax_climplant=mean(YearMaxMean,na.rm=T),
                               cit_climplant_05=mean(YearMean05,na.rm=T),
                               cit_climplant_95=mean(YearMean95,na.rm=T),
                               cit_ecoplant=mean(topt,na.rm=T),
                               cit_ecoplant_picq=mean(topt_picq,na.rm=T),
                               sd_cit_ecoplant_picq=sd(topt_picq,na.rm=T),
                               n_sp_climplant=sum(!is.na(YearMeanMean)),
                               mean_azote=mean(azote,na.rm=T),
                               mean_N=mean(Nopt,na.rm=T),
                               mean_R=mean(R_ellenberg,na.rm=T),
                               mean_pH=mean(pHopt,na.rm=T),
                               mean_CN=mean(vi_CN,na.rm=T),
                               mean_L=mean(Li,na.rm=T),
                               RS=.N,
                               RS_forest=sum(forest_species=="Forest_species",na.rm=T),
                               RS_generalist=sum(forest_species=="Open_species",na.rm=T)),
                             by=plot_ID]

vege_plot_meta_data<-merge(vege_plot_meta_data,cti_indicator,by="plot_ID")

## use the microclimate model to predict the mean understory temperature
## and the contribution of each estimated parameter
vege_plot_meta_data[,tree_density_2018:=tree_density_projected,]
vege_plot_meta_data[,pred_T_mean:= predict(lm_agg_mean,newdata=vege_plot_meta_data)]
vege_plot_meta_data[,pred_T_max:= predict(lm_agg_max,newdata=vege_plot_meta_data)]

vege_plot_meta_data[,pred_elev:= get_coef_mean[1,1] + get_coef_mean["mnt_25_vosges",1]*mnt_25_vosges ]
vege_plot_meta_data[,pred_micro:=  pred_T_mean - pred_elev  ]

vege_plot_meta_data[,pred_elev_max:= get_coef_max[1,1] + get_coef_max["mnt_25_vosges",1]*mnt_25_vosges ]
vege_plot_meta_data[,pred_micro_max:=  pred_T_max - pred_elev_max  ]

vege_plot_meta_data[,pred_heat_load:=  get_coef_mean["Heat_load_index_25",1]*Heat_load_index_25 ]
vege_plot_meta_data[,pred_ipv:=  get_coef_mean["ipv_25",1]*ipv_25 ]
vege_plot_meta_data[,pred_canopy:=  get_coef_mean["tree_density_2018",1]*tree_density_2018 ]
vege_plot_meta_data[,pred_topo:=pred_heat_load + pred_ipv ]

vege_plot_meta_data[,pred_heat_load_max:=  get_coef_max["Heat_load_index_25",1]*Heat_load_index_25 ]
vege_plot_meta_data[,pred_ipv_max:=  get_coef_max["ipv_25",1]*ipv_25 ]
vege_plot_meta_data[,pred_canopy_max:=  get_coef_max["tree_density_2018",1]*tree_density_2018 ]
vege_plot_meta_data[,pred_topo_max:=pred_heat_load_max + pred_ipv_max ]

## pred_topo is relative to the maximum temperature created by topogrpahy: exposed ridges (ipv25 = 1, heatload = 1)
#maximum_pred_topo =   get_coef_mean["Heat_load_index_25",1]* mean(vege_plot_meta_data$ipv_25) + get_coef_mean["ipv_25",1]*mean(vege_plot_meta_data$Heat_load_index_25)
mean_pred_topo =   get_coef_mean["Heat_load_index_25",1]* 0.66 + get_coef_mean["ipv_25",1]* 0.5

summary(vege_plot_meta_data)

vege_plot_meta_data[,pred_topo:=pred_topo -mean_pred_topo ]

## we remove the few plots with positive prediction of microclimate, they were plots that fall into an open canopy
## and are unreliable in light of what the microclimate model was trained with
vege_plot_meta_data<-vege_plot_meta_data[pred_micro<0,]

table_sp_vosges<-table_sp_vosges[vege_plot_meta_data$plot_ID,]

## these model are fitted with only bioindicatded soil pH, so that the effect of acidity can be latter removed
## these corrected values are not used anymore as suggested by the reviewers
model_out_ph<-lm(RS~ mean_pH,data=vege_plot_meta_data)
summary(model_out_ph)
vege_plot_meta_data$RS_pH_corrected<-residuals(model_out_ph)-16.0019 +  6.7577* 5.163606

model_out_ph<-lm(cit_climplant~ mean_pH,data=vege_plot_meta_data[!is.na(cit_climplant),])
summary(model_out_ph)
vege_plot_meta_data[!is.na(cit_climplant),cit_climplant_pH_corrected:=residuals(model_out_ph)+ 6.12442 + 0.32454 *5.1707]

mean(vege_plot_meta_data$ANNEE)

### we create 3 (of topo or micro climate) categories with equal number of plots (102), to study i na discrete manner CTI and Specific richness
vege_plot_meta_data[,cut_canopy:=cut_number(pred_canopy,3)]
vege_plot_meta_data[,cut_topo:=cut_number(pred_topo,3)]

unique(vege_plot_meta_data$cut_canopy)
levels(vege_plot_meta_data$cut_canopy)<-c("[-2.70 : -2.50 °]\nCold microclimate","[-2.50 : -2.30 °]\nModerate microclimate","[-2.30 : -1.30°]\nWarm microclimate")


unique(vege_plot_meta_data$cut_topo)
levels(vege_plot_meta_data$cut_topo)<-c("[-0.70 : -0.15°]\nCold topoclimate","[-0.15: 0.20°]\nModerate topoclimate","[0.20 : 0.70 °]\nWarm topoclimate")

### Vegetation analysis ####
# small ordination for supplementary
library(ade4)
influencial_plots<-c("521341_15","3F4MO7_19","4052605_21")# these 3 plots distord the ordination

table_sp_vosges_tmp<-as.data.table(table_sp_vosges[vege_plot_meta_data[!plot_ID %in% influencial_plots]$plot_ID,])

ca<-dudi.coa(table_sp_vosges_tmp,nf = 4, scannf = FALSE) 

vege_plot_meta_data_ordi<-cbind(vege_plot_meta_data[!plot_ID %in% influencial_plots],data.table(cbind(ca$li,ca$l1)))

plot_ordination<-ggplot(vege_plot_meta_data_ordi,aes(x=Axis1,y=Axis2,color=cut_topo,fill=cut_topo))+
  theme_bw()+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.spacing =   unit(5, 'cm'))+
  geom_point(alpha=0.75,pch=21,color="grey20")+
  stat_ellipse(type="t",geom = "polygon",alpha=0.15,segments =100)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+
  scale_color_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+
  labs(fill="Topoclimate \ncategories",color="Topoclimate \ncategories")

ggsave(file.path("Figure_result","plot_ordination.jpg"),plot_ordination,dpi=175,width = 180,height=110,unit="mm")

rm(table_sp_vosges_tmp)

#### vegetation model fitting ####

## table to export the results of linear model prediction specific richness (RS) or community thermla index (CTI)
create_table_model_flora<-function(model,do_varpart=F,ismax=F){
  
  get_coef<-summary(model)$coefficients
  colnames(get_coef)[3] <- "testval"
  colnames(get_coef)[4] <- "pval"
  
  get_range<-function(x){
    paste0(round(range(x),3),collapse=" : ")
    
  }  
  
  names_predictor <- c("pred_elev","pred_topo","pred_canopy","mean_pH")
  if(ismax) names_predictor <- c("pred_elev_max","pred_topo_max","pred_canopy_max","mean_pH")
  
  range_predictor<-c(NA,apply(vege_plot_meta_data[,..names_predictor],2,function(x)diff(range(x))))
  range_predictor_char<-c(NA,apply(vege_plot_meta_data[,..names_predictor],2,get_range))
  
  mean_predictor <- apply(vege_plot_meta_data[,..names_predictor],2,mean)
  min_pred <- as.matrix(rbind(mean_predictor,mean_predictor,mean_predictor,mean_predictor))
  max_pred <- as.matrix(rbind(mean_predictor,mean_predictor,mean_predictor,mean_predictor))
  
  min_sd <- function(x) mean(x)
  max_sd <- function(x) mean(x) + sd(x)
  
  diag(min_pred) <- apply(vege_plot_meta_data[,..names_predictor],2,min_sd)
  diag(max_pred) <- apply(vege_plot_meta_data[,..names_predictor],2,max_sd)
  
  min_range <- predict(model,data.table(min_pred),type="response")
  max_range <- predict(model,data.table(max_pred),type="response")
  max_range - min_range
  
  
  export_table_model<-data.table(get_coef,keep.rownames=T)
  export_table_model[,range:=range_predictor]
  export_table_model[,effect_size:=c(NA,max_range - min_range)]
  
  # ## variance partitionning
  # if(do_varpart){
  #   grp<-data.frame(var=c("pred_elev","pred_topo","pred_canopy","mean_pH"),group=c("Elevation","Topo","Canopy","Canopy"))
  #   
  #   coef_varpart<-as.vector(t(varPart(model=glm(model),groups = grp)))
  #   names(coef_varpart)<-colnames(t(varPart(model=glm(model),groups = grp)))
  #   
  #   var_elev<-coef_varpart["Elevation"] + coef_varpart["Elevation_Topo"]/2 + coef_varpart["Elevation_Canopy"]/2  +   coef_varpart["Elevation_Topo_Canopy"]/3
  #   var_topo<- coef_varpart["Topo"] + coef_varpart["Elevation_Topo"]/2 + coef_varpart["Topo_Canopy"]/2  +   coef_varpart["Elevation_Topo_Canopy"]/3
  #   var_canopy<- coef_varpart["Canopy"] + coef_varpart["Topo_Canopy"]/2 + coef_varpart["Elevation_Canopy"]/2  +   coef_varpart["Elevation_Topo_Canopy"]/3
  #   
  #   export_table_model[,variance:= c(NA,var_elev,var_topo,var_canopy,var_canopy)*100]
  # }
  # 
  export_table_model<-export_table_model[,-c("testval")]
  export_table_model<-cbind(export_table_model[,"rn"],signif(export_table_model[,-"rn"],3))
  export_table_model[,pval:=ifelse(pval<10e-4,"<10-4",pval)]
  export_table_model[,range:=range_predictor_char]
  
  if(!do_varpart) export_table_model<-export_table_model[,c(1,2,3,5,6,4)] else export_table_model<-export_table_model[,c(1,2,3,5,6,4,7)]
  export_table_model[,r_squared:=summary(model)$r.squared]
  
  return(export_table_model)
  
  
}

## linear regression
lm_cit_no_ph_correction<-lm(cit_climplant~pred_elev+pred_topo+pred_canopy+mean_pH,data=vege_plot_meta_data)
lm_cit_no_ph_correction_max<-lm(cit_climplant~pred_elev_max+pred_topo_max+pred_canopy_max+mean_pH,
                            data=vege_plot_meta_data)

summary(lm_cit_no_ph_correction)
library(MASS)
glm_rs_no_ph_correction<-glm.nb(RS~pred_elev+pred_topo+pred_canopy+mean_pH,data=vege_plot_meta_data)
glm_rs_no_ph_correction_max<-glm.nb(RS~pred_elev_max+pred_topo_max+pred_canopy_max+mean_pH,
                                data=vege_plot_meta_data)

summary(lm_cit_no_ph_correction)
summary(glm_rs_no_ph_correction)

vege_plot_meta_data[,cor(pred_elev,mean_pH)]

car::vif(lm_cit_no_ph_correction)
car::vif(glm_rs_no_ph_correction)

library(DHARMa)

testResiduals(lm_cit_no_ph_correction)
testResiduals(glm_rs_no_ph_correction)
testResiduals(lm_agg_mean)

all_tms_data_agg_fit[,resid_T3_mean:=residuals(lm_agg_mean)]
ggplot(all_tms_data_agg_fit,aes(x=resid_T3_mean))+theme_classic()+geom_histogram(fill="grey80",color="grey40",boundary=0,binwidth =0.125)
ggplot(all_tms_data_agg_fit,aes(x=xl93,y=resid_T3_mean))+theme_classic()+geom_point()+geom_smooth()
ggplot(all_tms_data_agg_fit,aes(x=yl93,y=resid_T3_mean))+theme_classic()+geom_point()+geom_smooth()
ggplot(all_tms_data_agg_fit,aes(x=Heat_load_index_25,y=resid_T3_mean))+theme_classic()+geom_point()+geom_smooth()


create_plot_res <- function(model){
  p_a <- ggplot(mappin = aes(x = fitted(model), y= residuals(model)))+
    geom_point(size=1)+
    geom_smooth(method="lm")+
    theme_bw()+
    labs(x="Fitted values",y="Residuals")
  
  p_b <- ggplot(mappin = aes(x = residuals(model)))+
    theme_bw()+
    geom_histogram(fill="#8ABFFC",color="grey10",bins =16,origin = 0)+
    labs(x="Residuals")
 
  return( ggarrange(p_a,p_b,nrow = 1,widths = c(1,0.7)))
  
}

to_export_residuals <- ggarrange(create_plot_res(lm_agg_mean),
          create_plot_res(lm_cit_no_ph_correction),
          create_plot_res(glm_rs_no_ph_correction),
          nrow=3,
          labels= c("a)","b)","c)"))

ggsave(file.path("figure_result","residuals_of_models.png"),to_export_residuals,width = 180,height = 130,units = "mm")

## T mean predicting flora
table_cit<-create_table_model_flora(lm_cit_no_ph_correction,F,ismax = F)
table_rs<-create_table_model_flora(glm_rs_no_ph_correction,F)
table_rs[,r_squared:=NA]

export_flora_model<-rbind(table_rs,table_cit,fill=T)

write.table(export_flora_model,file.path("figure_result","Model_flora.csv"),row.names = F,sep=";",dec=".")

## T max predicting flora
table_cit<-create_table_model_flora(lm_cit_no_ph_correction_max,F,T)
table_rs<-create_table_model_flora(glm_rs_no_ph_correction_max,F,T)
table_rs[,r_squared:=NA]
export_flora_model_max<-rbind(table_rs,table_cit,fill=T)

write.table(export_flora_model_max,file.path("figure_result","Model_flora_max.csv"),row.names = F,sep=";",dec=".")


#### descritpives statistics ####

## microclimate prediction
summary(all_tms_data_agg_fit$T3_mean)
summary(all_tms_data_agg_fit$T3_min)
summary(all_tms_data_agg_fit$T3_max)

## specific richness
summary(vege_plot_meta_data$RS)
sd(vege_plot_meta_data$RS)
## CTI
summary(vege_plot_meta_data$cit_climplant)
sd(vege_plot_meta_data$cit_climplant,na.rm=T)
## per class of topoclimate
vege_plot_meta_data[,mean(mnt_25_vosges),by=cut_topo]
vege_plot_meta_data[,mean(RS),by=cut_topo]
vege_plot_meta_data[,mean(cit_climplant,na.rm=T),by=cut_topo]
vege_plot_meta_data[,mean(RS_pH_corrected),by=cut_topo]
vege_plot_meta_data[,mean(cit_climplant_pH_corrected,na.rm=T),by=cut_topo]

all_tms_data_agg_fit[,canopy_cover_glama:= as.numeric(str_replace_all(canopy_cover_glama,",","."))]
## correlation between field and remote sensed canopy density
cor.test(all_tms_data_agg_fit$canopy_cover_glama,all_tms_data_agg_fit$tree_density_projected,use="complete.obs")
cor.test(all_tms_data_agg_fit$canopy_cover_25m,all_tms_data_agg_fit$tree_density_2018,use="complete.obs")
cor.test(all_tms_data_agg_fit$canopy_cover_25m,all_tms_data_agg_fit$canopy_cover_glama,use="complete.obs")

## descritpive number for species
flora_survey<-flora_survey[plot_ID%in% vege_plot_meta_data$plot_ID]

## unique species
flora_survey[,length(unique(species_name))]
## unique species with a know Topt
flora_survey[! is.na(YearMeanMean) ,length(unique(species_name))]
## unique species with a know habitat preference
flora_survey[! habitat_categ%in%c(NA,"/") ,length(unique(species_name))]
## proportions
1-flora_survey[plot_ID%in% vege_plot_meta_data$plot_ID ,sum(is.na(YearMeanMean))  ]/nrow(flora_survey)
1-flora_survey[plot_ID%in% vege_plot_meta_data$plot_ID ,sum( habitat_categ%in%c(NA,"/"))  ]/nrow(flora_survey)

## beta div descriptors
table_sp_vosges_2<-data.table(table_sp_vosges)[,-"plot_ID"]

plot_cold<-vege_plot_meta_data[cut_topo==levels(vege_plot_meta_data$cut_topo)[1],plot_ID]
plot_moder<-vege_plot_meta_data[cut_topo==levels(vege_plot_meta_data$cut_topo)[2],plot_ID]
plot_warm<-vege_plot_meta_data[cut_topo==levels(vege_plot_meta_data$cut_topo)[3],plot_ID]

for(i in list(plot_cold,plot_moder,plot_warm)){
  subset_table<-table_sp_vosges_2[rownames(table_sp_vosges)%in% i,]
  atleastone<-apply(subset_table,2,sum)
  atleastone<-atleastone[atleastone!=0]
  gamma<-length(atleastone)
 
  alpha<-mean(apply(subset_table,1,sum))
  
  print(gamma)
  print(alpha)
  print(gamma/alpha)
  print("- - - ")
  
}

## all the following code is to detect sepcies unique to one topoclimate class
table_sp_vosges_2$plot_ID<-rownames(table_sp_vosges)
table_sp_vosges_cut<-as.data.table(merge(table_sp_vosges_2,vege_plot_meta_data[,c("plot_ID","cut_topo")]))
table_sp_vosges_cut[,plot_ID:=NULL]

table_sp_vosges_cut<-melt(table_sp_vosges_cut)
table_sp_vosges_cut<-table_sp_vosges_cut[,.(occ=sum(value)),by=.(cut_topo,variable)]
table_sp_vosges_cut[,occ_tot:=sum(occ),by=variable]
table_sp_vosges_cut<-table_sp_vosges_cut[occ!=0,]
table_sp_vosges_cut


table_sp_vosges_cut[,unique_to_cold:=ifelse(occ == occ_tot & cut_topo=="[-1.55 : -1.0°]\nCold topoclimate",1,0)]
table_sp_vosges_cut[,unique_to_moder:=ifelse(occ == occ_tot & cut_topo=="[-1.0 : -0.65 °]\nModerate topoclimate",1,0)]
table_sp_vosges_cut[,unique_to_warm:=ifelse(occ == occ_tot & cut_topo=="[-0.65 : -0.15 °]\nWarm topoclimate",1,0)]

table_sp_vosges_cut[unique_to_cold==1,]
table_sp_vosges_cut[unique_to_moder==1,]
table_sp_vosges_cut[unique_to_warm==1,]

nrow(table_sp_vosges_cut[unique_to_cold==1,])
nrow(table_sp_vosges_cut[unique_to_moder==1,])
nrow(table_sp_vosges_cut[unique_to_warm==1,])

#### Figures ####

### boxplots of CTI and specific richness per classes

## list used to display significance of wilcox.test
my_comparisons <- list( levels(vege_plot_meta_data$cut_topo)[1:2], 
                        levels(vege_plot_meta_data$cut_topo)[2:3], 
                        levels(vege_plot_meta_data$cut_topo)[c(1,3)] )


RS_plot<-ggplot(vege_plot_meta_data[,],aes(x=cut_topo,y=RS,fill=cut_topo))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="grey20",size=1.5,show.legend = F)+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F)+
  #geom_violin(alpha=0.5,show.legend = F,draw_quantiles = 0.5)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  coord_cartesian(ylim=c(-2,68))+
  labs(y="Species richness",x="Topoclimate class")+
  stat_compare_means(comparisons = my_comparisons, label.y = c(55,58,61),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns")))


CTI_plot<-ggplot(vege_plot_meta_data[,],aes(x=cut_topo,y=cit_climplant,fill=cut_topo))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="grey20",size=1.5,show.legend = F)+
  scale_y_continuous(limits=c(5.75,9.65))+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F,draw_quantiles = 0.5)+
  #geom_violin(alpha=0.5, outlier.alpha = 0,show.legend = F,draw_quantiles = 0.5)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  labs(y="Community thermal index (°C)",x="Topoclimate class")+
  stat_compare_means(comparisons = my_comparisons ,label.y = c(9,9.25,9.45),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns")))



Export_boxplot<-ggarrange(plotlist=list(RS_plot,CTI_plot),ncol=2,labels=c("a)","b)"))

ggsave(plot=Export_boxplot,file.path("figure_result","Figure_boxplot.jpg"),width=180,height =80 ,scale=1.25,unit="mm",dpi=650)


my_comparisons <- list( levels(vege_plot_meta_data$cut_canopy)[1:2], 
                        levels(vege_plot_meta_data$cut_canopy)[2:3], 
                        levels(vege_plot_meta_data$cut_canopy)[c(1,3)] )

RS_plot_canopy<-ggplot(vege_plot_meta_data,aes(x=cut_canopy,y=RS,fill=cut_canopy))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="white",size=1.5,show.legend = F)+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  coord_cartesian(ylim=c(-2,68))+
  labs(y="Species richness",x="Microclimate categories")+ #\ncorrected for pH
  stat_compare_means(comparisons = my_comparisons, label.y = c(55,58,61),
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")))


CTI_plot_canopy<-ggplot(vege_plot_meta_data[,],aes(x=cut_canopy,y=cit_climplant,fill=cut_canopy))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="white",size=1.5,show.legend = F)+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  labs(y="Community thermal index (°C)",x="Microclimate categories")+ # \ncorrected for pH
  stat_compare_means(comparisons = my_comparisons ,label.y = c(9.25,9.5,9.75)-0.1,
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")))



Export_boxplot_canopy<-ggarrange(plotlist=list(RS_plot_canopy,CTI_plot_canopy),ncol=2,labels=c("a)","b)"))


ggsave(plot=Export_boxplot_canopy,file.path("figure_result","Figure_boxplot_suppl_canopy.jpg"),width=180,height =80 ,scale=1.25,unit="mm",dpi=400)

### continuous version( no boxplots)

RS_plot_continuous<-ggplot(vege_plot_meta_data,aes(x=pred_topo ,y=RS,fill=cut_topo))+
  theme_bw()+
  geom_point(alpha=0.5,pch=21,color="grey10",size=3,show.legend = T)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  geom_smooth(method="lm",fill="grey75",color="grey10")+
  #coord_cartesian(ylim=c(0,72))+
  labs(y="Species richness",x="Topoclimate effect",fill="Topoclimatic classes")


CTI_plot_continuous<-ggplot(vege_plot_meta_data,aes(x=pred_topo ,y=cit_climplant,fill=cut_topo))+
  theme_bw()+
  geom_point(alpha=0.5,pch=21,color="grey10",size=3,show.legend = T)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  geom_smooth(method="lm",fill="grey75",color="grey10")+
  #coord_cartesian(ylim=c(0,72))+
  labs(y="Community thermal index (°C)",x="Topoclimate effect",fill="Topoclimatic classes")

export_continuous<-ggarrange(plotlist=list(RS_plot_continuous,CTI_plot_continuous),nrow=2,labels=c("a)","b)"),common.legend = T,legend = "bottom")

ggsave(plot=export_continuous,file.path("figure_result","Figure_continuous.jpg"),width=180,height =180 ,scale=1,unit="mm",dpi=400)

### histogram of species Topt and habitat preferences per classes

vege_plot_meta_data_melt<-melt(vege_plot_meta_data,measure.vars = c("RS_forest","RS_generalist"))
vege_plot_meta_data_melt[,variable:=ifelse(variable=="RS_forest","Forest specialist","Generalist")]

# RS_plot_2<-ggplot(vege_plot_meta_data_melt[,],aes(x=cut_topo,y=value,fill=cut_topo,lty=variable))+
#   theme_bw()+
#   geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.35,pch=21,size=1.2,show.legend = F)+
#   
#   geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = T)+
#   scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"),guide="none")+ 
#   coord_cartesian(ylim=c(0,35))+
#   labs(y="Specific richness",x="Topoclimate categories",lty=" ")+
#   theme(legend.position = c(0.9,0.9))
# 
# RS_plot_2
# 


col<-c("plot_ID",grep("cut_",colnames(vege_plot_meta_data),value = T))
flora_survey_sub<-merge(flora_survey,vege_plot_meta_data[,..col],by="plot_ID")
flora_survey_sub[,cut_YearMeanMean:=cut(YearMeanMean,10)]



climplant_histogram<-ggplot(flora_survey_sub,aes(x=YearMeanMean ,fill=cut_topo))+theme_bw()+
  geom_histogram(binwidth =1,alpha=0.5,origin=0,color="grey40",position ="dodge")+
  scale_fill_manual(values=c("#440154FF","#22A884FF","#FDE725FF"))+
  scale_x_continuous(breaks = seq(from=3,to=13,by =1))+
  scale_y_continuous(sec.axis = sec_axis( transform=~./(nrow(vege_plot_meta_data)/3), name="  "))+
  coord_cartesian(expand=0.1)+
  #geom_vline(xintercept = seq(from=3,to=12.5,by =1))+
  theme( panel.grid.major.x = element_line(color="grey60"),panel.grid.minor.x = element_line(color="white"),legend.spacing.y = unit(0.15, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))+
  labs(x="Thermal optimum (°C)",y="Occurences",fill="Topoclimate\ncategories")#+facet_wrap(~forest_species)


habitat_histogram<-ggplot(flora_survey_sub[habitat_categ%in%c("1,1","1,2","2,1","2,2","O"),],aes(x=habitat_categ ,fill=cut_topo))+theme_bw()+
  geom_histogram(binwidth =1,alpha=0.5,origin=0,color="grey40",position ="dodge",stat="count")+
  scale_fill_manual(values=c("#440154FF","#22A884FF","#FDE725FF"))+
  scale_y_continuous(breaks=c(0,200,400,600,800),sec.axis = sec_axis( trans=~./(nrow(vege_plot_meta_data)/3), name="Plot scale"))+
  #geom_vline(xintercept = seq(from=3,to=12.5,by =1))+
  theme( panel.grid.minor.x = element_line(color="grey60"),panel.grid.major.x = element_line(color="white"),legend.spacing.y = unit(0.15, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))+
  labs(x="Species habitat affinity",y="Occurences",fill="Topoclimate\ncategories")#+facet_wrap(~forest_species)


export_histogram<-ggarrange(plotlist=list(climplant_histogram,habitat_histogram),ncol=2,labels=c("a)","b)"),widths = c(1.8,1),common.legend = T,legend = "bottom")

ggsave(plot=export_histogram,file.path("figure_result","Figure_histogram.jpg"),width=180,height =90 ,scale=1,unit="mm",dpi=550)


####  Creation of microclimate maps ####
## dt_raster is a data?tble object that will contain all the spatialized variable and the predicted understory temperature
dt_raster<-data.table(as.data.frame(stack_all_variable_vallee))
dt_raster<-cbind(dt_raster,xyFromCell(stack_all_variable_vallee$mnt_25_vosges,1:length(stack_all_variable_vallee$mnt_25_vosges)))
dt_raster[,tree_density_2018:=tree_density_projected]

dt_raster[,pred_T_mean:= predict(lm_agg_mean,newdata=dt_raster)]
dt_raster[,pred_elev:= get_coef_mean[1,1] + get_coef_mean["mnt_25_vosges",1]*mnt_25_vosges ]
dt_raster[,pred_micro:=  pred_T_mean - pred_elev  ]

dt_raster[,pred_heat_load:=   get_coef_mean["Heat_load_index_25",1]*Heat_load_index_25 ]
dt_raster[,pred_ipv:=   get_coef_mean["ipv_25",1]*ipv_25  ]
dt_raster[,pred_topo:=   pred_heat_load+pred_ipv ]
dt_raster[,pred_topo:=   pred_topo - mean_pred_topo]

dt_raster[,pred_canopy:=   get_coef_mean["tree_density_2018",1]*tree_density_2018  ]

## proportion of public forest computation (!is.na(pred_elev  delineates the study region)
dt_raster[,foret_pub:= !is.na(pred_elev) & is.na(foret_pub_25_vosges)]
dt_raster[!is.na(pred_elev),sum(!foret_pub)/length(foret_pub)]


#### Maps #### 
## usefull for later to show where the inset is
square_inset<-st_buffer(st_as_sf(data.table(x=1000125,y=6764575),coords=c("x","y"),crs=st_crs(2154)),1000,endCapStyle="SQUARE")
valley_border<-st_read("data",layer="valley_shp")

map_elevation_climate<-ggplot(dt_raster,aes(x,y,fill=pred_elev+ get_coef_mean["tree_density_2018",1]*90))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster$pred_elev+ get_coef_mean["tree_density_2018",1]*90,probs=c(0.025,0.975),na.rm=T))+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Mean \ntemperature °C")+ 
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))


map_pred_topo<-ggplot(dt_raster,aes(x,y,fill=pred_topo ))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=c(-0.7,0.7),breaks=c(-0.7,-0.35,-0,0.35,0.7),labels = c("-0.70","-0.35","0","0.35","0.70"))+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Topoclimatic \neffect °C")+
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))

## We bound the extreme value for a easier display 
dt_raster[,pred_canopy_2:=ifelse(is.na(pred_T_mean),NA,pred_canopy)]
dt_raster[,pred_canopy_2:=ifelse(pred_canopy_2> -1.5,-1.5,pred_canopy_2)]
## used for the mini maps below
subset_map_dt<-dt_raster[x%between% c(999125-100,999625+1500)& y %between%c(6765075-1500,6766075-500) ,]

map_pred_canopy<-ggplot(dt_raster[tree_density_projected  == 0],aes(x,y,fill=pred_canopy_2 ))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster$pred_canopy_2,probs=c(0.025,0.975),na.rm=T))+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Canopy closure\neffect °C")+
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))

ggplot(dt_raster[ ],aes( x= mnt_25_vosges, y= tree_density_projected))+
  geom_point(alpha=0.015)+
  theme_classic()+
  geom_point(data = vege_plot_meta_data,aes(fill="Vegetation surveys"),color="black",pch=21,size=1.75)+
  geom_point(data = coords[site_ID %in% all_tms_data_agg_fit$locality_id,],aes(fill = "Loggers"),color="black",pch=21,size=2.25)+
  scale_fill_manual(values = c("white","#55C463"))+
  labs(x = "Elevation (M a.s.l.)",y = "canopy closure",fill="Plots")


## create a zoomed inset of one raster
create_mini_map<-function(what){

mini_map_1<-ggplot(subset_map_dt,aes(x,y,fill=get(what) ))+
  theme_void()+
  geom_raster(show.legend = F)+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster[,..what],probs=c(0.025,0.975),na.rm=T))+
  labs(fill="var")+
  coord_sf(expand=F)+
  theme(panel.border = element_rect(fill=NA,color="grey20",linewidth =3) , plot.margin = margin(10,10,10,10))+
  labs(x="",y="",fill="Temperature °C")
  return(mini_map_1)

}

minimaps<-ggarrange(plotlist = list( create_mini_map("pred_elev"),create_mini_map("pred_topo"),create_mini_map("pred_canopy_2")),ncol=2,nrow=2,align = "hv")

out_4_panel<-ggarrange(plotlist=list(map_elevation_climate+ggtitle("    "),
                                     map_pred_topo,
                                     map_pred_canopy,minimaps),
                       ncol=2,nrow=2,labels = c("a) Lapse rate","b) Topograpy ","c) Canopy    ","d)        "),align = "hv",
                       font.label = list(size = 18, color = "grey5", face = "bold", family = NULL),hjust = 0)

ggsave(file.path("Figure_result","climate_map_2.jpg"),
       out_4_panel,
       dpi=400,unit="mm",width = 180,height = 180,scale=1.5)


#### code to create the sampling map ####
extended_mnt<-raster(file.path("data","extented_dem.tif"))
extended_mnt<-crop(extended_mnt,extent(stack_all_variable_vallee$bdv2_25_vosges)+2000)
slope<-slopeAspect (extended_mnt)
hs<-hillShade(slope$slope,slope$aspect)
extended_mnt<-stack(extended_mnt,hs)

dt_raster[,extended_dem:=getValues(extended_mnt$extented_dem)]
dt_raster[,hillshade:=getValues(hs)]


valley_border<-st_read("data",layer="valley_shp")
river_sf<-read_sf("data",layer="river_valley")
river_sf<-st_intersection(river_sf,valley_border)

coords_sf_map<-coords_sf[coords_sf$site_ID%in%unique(all_tms_data_agg$locality_id),]
coords_vege_map<-vege_plot_meta_data_sf[vege_plot_meta_data_sf$plot_ID%in%vege_plot_meta_data$plot_ID,]

map_fig_1_dt<-data.table(as.data.frame(extended_mnt))
map_fig_1_dt<-cbind(map_fig_1_dt,xyFromCell(extended_mnt$extented_dem,1:length(extended_mnt$extented_dem)))


sampling_map<-ggplot(map_fig_1_dt,aes(x,y,fill=extented_dem))+
  theme_classic()+
  geom_raster()+
  coord_fixed()+
  labs(x="",y="",fill="Elevation (m)",shape="Study sites")+
  scale_fill_gradientn(colors=terrain.colors(100),na.value ="transparent",breaks=c(300,600,900,1200,1400))+
  new_scale_fill() +
  geom_raster(mapping = aes(fill=layer),alpha=0.35,show.legend = F)+
  scale_fill_gradient(low = "black", high = "white",na.value="transparent") +
  geom_sf(data=river_sf,color="dodgerblue3",inherit.aes=F)+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.5,fill=NA,alpha=0,color="grey15") + 
  geom_sf(data=coords_vege_map,aes(shape="Floristic surveys"),inherit.aes=F,size=1.75)+
  geom_sf(data=coords_sf_map,aes(shape="Temperature \nloggers"),inherit.aes=F,size=2,fill="white")+
  scale_shape_manual(values=c(3,21))+
  coord_sf(expand=F)+
  annotation_north_arrow( height = unit(0.9, "cm"),width = unit(0.9, "cm"),location="tl")+
  annotation_scale()+
  theme(panel.border= element_rect(fill=NA,color="grey5",linewidth=1))


library(cowplot)
rnaturalearth<-rnaturalearth::ne_countries(
  scale = 10,
  type = "countries",
  continent = "Europe",
  country = NULL,
  geounit = NULL,
  sovereignty = NULL,
  returnclass = c( "sf")
)


vosges_sf<-read_sf("data",layer="vosges_border")
valley_border_point<-st_centroid(valley_border)

europe_country<-(rnaturalearth[rnaturalearth$name%in%c("France","Spain","poland","United Kingdom","Andorra","Portugal","Belgium","Switzerland","Italy","Germany","Austria","Liechtenstein","Netherlands","Luxembourg"),])
europe_country<-rnaturalearth
europe_country<-st_transform(europe_country,st_crs(2154))

inset_france_vosges<-ggplot(europe_country)+
  geom_sf(fill="grey98",color="grey60")+
  geom_sf(data=vosges_sf,fill="grey60",color="grey60")+
  geom_sf(data=valley_border_point,color="grey5",size=0.75)+
  coord_sf(xlim= c(-100000+300*1000,1000000+500*1000 +300*1000 ),ylim =c(5850000+200*1000,7350000+200*1000))+
  theme_void()+
  theme(panel.border =  element_rect(color="grey30",fill=NA,linewidth = 1),panel.background = element_rect(color="grey10",fill="white"))

sampling_map_inset<-ggdraw(sampling_map)+draw_plot(inset_france_vosges,0.45, 0.647, 0.32, 0.32)

ggsave(file.path("Figure_result","sampling_map.jpg"),sampling_map_inset,dpi=250,unit="mm",width = 88,height = 75,scale=2)

#### supplementary figures: remote and field canopy ####

## we gather the locality of the sensors we used, excluding the one with missing field measurements 
all_tms_data_supp_cano_plot<-all_tms_data_agg_fit[locality_id%in% meta_data[!is.na(canopy_cover_25m),site_ID],]

## A simple linear model to gather Pvalue, estimates and R2
lm_test<-lm(tree_density_2018~canopy_cover_25m,data=all_tms_data_supp_cano_plot)

r2<-signif(summary(lm_test)$r.squared,4)
pval<-summary(lm_test)$coefficients[2,4]
a<-signif(summary(lm_test)$coefficients[1,1],3)
b<-signif(summary(lm_test)$coefficients[2,1],3)
labs<-paste0("Y = ",a," + ",b," * Remote tree density \nR² = ",r2,ifelse(pval<0.05," ***","   ")," ")

cover_25m_supp<-ggplot(all_tms_data_supp_cano_plot,aes(x=tree_density_2018,y=canopy_cover_25m))+
  geom_point(color="black",fill="grey25",shape=21)+
  geom_smooth(method='lm')+
  theme_bw()+
  ylim(c(50,100))+
  labs(x="Remote sensed tree density (%)",y="Canopy cover smartphone\nglama photography (%)")+
  geom_text(aes(x=78,y=100,label=labs),inherit.aes = F,hjust   = 0,vjust   = 1,size=4)


## same process
lm_test<-lm(tree_density_2018~canopy_cover_glama,data=all_tms_data_supp_cano_plot)
r2<-signif(summary(lm_test)$r.squared,4)
pval<-summary(lm_test)$coefficients[2,4]
a<-signif(summary(lm_test)$coefficients[1,1],3)
b<-signif(summary(lm_test)$coefficients[2,1],3)
labs_glama<-paste0("Y = ",a," + ",b," * Remote tree density \nR² = ",r2,ifelse(pval<0.05," ***"," .")," ")


glama_supp<-ggplot(all_tms_data_supp_cano_plot,aes(x=tree_density_2018,y=canopy_cover_glama))+
  geom_point(color="black",fill="grey25",shape=21)+
  geom_smooth(method='lm')+
  theme_bw()+
  ylim(c(30,100))+
  labs(x="Remote sensed tree density (%)",y="Visually estimated \ncanopy cover (%)")+
  geom_text(aes(x=78,y=100,label=labs_glama),inherit.aes = F,hjust   = 0,vjust   = 1,size=4)

Export_suppl_canopy<-ggarrange(plotlist=list(cover_25m_supp,glama_supp),nrow=2,labels=c("a)","b)"))

ggsave(plot=Export_suppl_canopy,file.path("figure_result","Canopy_methods_suppl.jpg"),width=180,height =155 ,scale=1,unit="mm",dpi=300)

#### Appendice: List of species ####
Table_sum<-flora_survey_sub[,.(.N,topt = unique(YearMeanMean),
                               pH=unique(pHopt),
                               habitat_categ=unique(habitat_categ)),
                            by=.(species_name,cut_topo)]
Table_sum<-Table_sum[!grepl("subsp.",species_name),,]

Table_sum[is.na(topt),]
Table_sum[!is.na(topt),]

Table_sum<-Table_sum[order(species_name),]

Table_sum<-dcast(Table_sum,species_name+topt+pH+habitat_categ~cut_topo,value.var = "N")
colnames(Table_sum)[c(5,6,7)]<-c("cold_topo","inter_topo","warm_topo")

Table_sum[,cold_topo:=ifelse(is.na(cold_topo),0,cold_topo)]
Table_sum[,inter_topo:=ifelse(is.na(inter_topo),0,inter_topo)]
Table_sum[,warm_topo:=ifelse(is.na(warm_topo),0,warm_topo)]

Table_sum[,count:=cold_topo+inter_topo+warm_topo]
Table_sum<-Table_sum[,c(1,8,5,6,7,2,3,4)]

Table_sum[,topt:=round(topt,2)]

write.table(Table_sum,file.path("figure_result","liste_species.csv"),row.names = F,sep=";",dec=",")
