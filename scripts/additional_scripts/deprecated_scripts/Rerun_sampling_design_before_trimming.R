#### format of the vosges workflow  ####
setwd(file.path("~","AdapFor_analysis"))
source("Packages.R")

if(Sys.info()["nodename"]=="CALCULUS") setwd("//tsclient/C/Users/borderieux/Desktop/data_microclim") else setwd("C:/Users/borderieux/Desktop/data_microclim")



canopy_gedi<-raster("canopy/Forest_height_2019_NAFR_fance_l93.tif")
cover_2010<-raster("canopy/treecover2010_50N_000E.tif")
mnt_25_vosges<-raster("topo/mnt_25_vosges.tif")

tmoy_an_50m_86_10<-raster("mean_clim/digitalis_tmoy_50m_8610_L93.tif")
pp_an_50m_86_10<-raster("mean_clim/digitalis_pp_an_50m_8610_L93.tif")

border_vosges_shp<-read_sf("shp",layer="Vosges_border_shp")
vallee_vosge<-read_sf("shp",layer="Vallee_plaine")
ser_d11_vosges<-read_sf("shp",layer="ser_D11")

ipv_1000_50<-raster("D:/quercus/BD_SIG/topo/bdalti50/indices_topo/ipv/ipv_x_1000.tif")
extent_l2<-extent(projectRaster(mnt_25_vosges_crop,crs=crs(ipv_1000_50)) )
ipv_1000_50<-crop(ipv_1000_50,extent_l2)
ipv_1000_50<-projectRaster(ipv_1000_50,to=mnt_25_vosges_crop)
extent<-st_bbox(border_vosges_shp)
#extent_lamb2<-st_bbox(st_transform(vallee_vosge,crs=st_crs(  27572 )))
extent_vosges<-extent(x=extent[1],xmax=extent[3],ymin=extent[2],ymax=extent[4])
#extent_vosges_lamb2<-extent(x=extent_lamb2[1],xmax=extent_lamb2[3],ymin=extent_lamb2[2],ymax=extent_lamb2[4])

extent<-st_bbox(vallee_vosge)
#extent_lamb2<-st_bbox(st_transform(vallee_vosge,crs=st_crs(  27572 )))
extent_vallee<-extent(x=extent[1],xmax=extent[3],ymin=extent[2],ymax=extent[4])
#extent_vosges_lamb2<-extent(x=extent_lamb2[1],xmax=extent_lamb2[3],ymin=extent_lamb2[2],ymax=extent_lamb2[4])

#○extents are in l 93



mnt_25_vosges_crop<-crop(mnt_25_vosges,extent_vosges)
canopy_gedi_crop<-crop(canopy_gedi,extent_vosges)
canopy_gedi_crop<-projectRaster(canopy_gedi_crop,mnt_25_vosges_crop)


tmoy_an_50m_86_10<-projectRaster(tmoy_an_50m_86_10,mnt_25_vosges_crop)
pp_an_50m_86_10<-projectRaster(pp_an_50m_86_10,mnt_25_vosges_crop)





mnt_25_tpi<-terrain(mnt_25_vosges_crop,opt="tpi")

## other tpi value for cad and cap calculation
mnt_25_tpi_250<-terrain(mnt_25_vosges_crop,opt="tpi")
weight_tpi<-matrix(1,3,3) 
ceiling(length(weight_tpi)/2)
weight_tpi_125<-matrix(1,13,13)
weight_tpi_250<-matrix(1,21,21)
weight_tpi_500<-matrix(1,41,41)
weight_tpi_1000<-matrix(1,81,81)
weight_tpi_1500<-matrix(1,121,121)
#weight_tpi_2500<-matrix(1,201,201)

mnt_25_tpi_1000<-focal(mnt_25_vosges_crop,w=weight_tpi_1000,fun=function(x, ...) x[ceiling(length(weight_tpi_1000)/2)] - mean(x[-(ceiling(length(weight_tpi_1000)/2))]), pad=TRUE, padValue=NA)


#cap_25_tpi_150<-focal(mnt_25_vosges_crop,w=weight_tpi_125,fun=function(x, ...) log(x[5] - min(x[-5])), pad=TRUE, padValue=NA)
#cap_25_tpi_250<-focal(mnt_25_vosges_crop,w=weight_tpi_250,fun=function(x, ...) log(x[5] - min(x[-5])), pad=TRUE, padValue=NA)

cap_25_tpi_500<-focal(mnt_25_vosges_crop,w=weight_tpi_500,fun=function(x, ...) log(x[ceiling(length(weight_tpi_500)/2)] - min(x[-(ceiling(length(weight_tpi_500)/2))])), pad=TRUE, padValue=NA)
cap_25_tpi_1000<-focal(mnt_25_vosges_crop,w=weight_tpi_1000,fun=function(x, ...) log(x[ceiling(length(weight_tpi_1000)/2)] - min(x[-(ceiling(length(weight_tpi_1000)/2))])), pad=TRUE, padValue=NA)
cap_25_tpi_1500<-focal(mnt_25_vosges_crop,w=weight_tpi_1500,fun=function(x, ...) log(x[ceiling(length(weight_tpi_1500)/2)] - min(x[-(ceiling(length(weight_tpi_1500)/2))])), pad=TRUE, padValue=NA)
#cap_25_tpi_2500<-focal(mnt_25_vosges_crop,w=weight_tpi_2500,fun=function(x, ...) log(x[ceiling(length(weight_tpi_2500)/2)] - min(x[-(ceiling(length(weight_tpi_2500)/2))])), pad=TRUE, padValue=NA)


plot(cap_25_tpi_500)

plot( cap_25_tpi_500)

mnt_25_slope<-terrain(mnt_25_vosges_crop,opt="slope",unit="degrees")
mnt_25_aspect<-terrain(mnt_25_vosges_crop,opt="aspect")

mnt_25_aspect<-cos(mnt_25_aspect)



if(Sys.info()["nodename"]=="CALCULUS") setwd("D:/quercus/BD_SIG/occup_sol/IFN/_BD_FORET") else setwd("S:/BD_SIG/occup_sol/IFN/_BD_FORET")

big_shp_dep<-readRDS("France_entiere/France_all_forest_bdv2.RData")


if(Sys.info()["nodename"]=="CALCULUS") setwd("D:/quercus/BD_SIG/occup_sol/Forets_publiques_domaniales_communales") else setwd("S:/BD_SIG/occup_sol/Forets_publiques_domaniales_communales")

foret_pub<-read_sf(".",layer="Foret_pub")

if(Sys.info()["nodename"]=="CALCULUS") setwd("//tsclient/C/Users/borderieux/Desktop/data_microclim") else setwd("C:/Users/borderieux/Desktop/data_microclim")




library(fasterize)

raster_bdv2_25m<-fasterize(big_shp_dep,mnt_25_vosges_crop,field = "integer_forest")
reclass_matric_all_forest<-matrix(c(1:32,1:30,NA,NA),ncol = 2)
raster_bdv2_25m<-reclassify(raster_bdv2_25m,reclass_matric_all_forest)
## remvoe lande and prairie ?


foret_pub$domaniale<-ifelse(foret_pub$cdom_frt=="OUI",1,0)
raster_foret_pub_25m<-fasterize(foret_pub,mnt_25_vosges_crop,field = "domaniale")


plot(raster_distance_25m)



reclass_matric_distance_to_edge<-matrix(c(1:32,NA,rep(NA,30),0,0,0),ncol = 2)

raster_distance_25m<-reclassify(raster_bdv2_25m,reclass_matric_distance_to_edge)

#raster_distance_25m<-distance(raster_distance_25m)


#longer calculation

mnt_25_tri<-terrain(mnt_25_vosges_crop,opt="tri")


raster_distance_25m<-gridDistance(raster_distance_25m, 0) 


plot(canopy_gedi_crop)

plot(mask(mnt_25_vosges_crop,as_Spatial(border_vosges_shp)))




plot(raster_distance_25m)




writeRaster(mnt_25_vosges_crop,file="export_final_crop_l93_25m/mnt_25_vosges.tif",overwrite=TRUE)
writeRaster(mnt_25_tpi,file="export_final_crop_l93_25m/tpi_25_vosges.tif",overwrite=TRUE)
writeRaster(mnt_25_tri,file="export_final_crop_l93_25m/tri_25_vosges.tif",overwrite=TRUE)
writeRaster(mnt_25_slope,file="export_final_crop_l93_25m/slope_25_vosges.tif",overwrite=TRUE)
writeRaster(mnt_25_aspect,file="export_final_crop_l93_25m/aspect_NS_25_vosges.tif",overwrite=TRUE)
writeRaster(canopy_gedi_crop,file="export_final_crop_l93_25m/canopy_height_25_vosges.tif",overwrite=TRUE)
writeRaster(tmoy_an_50m_86_10,file="export_final_crop_l93_25m/tmoy_digitalis_8610_25_vosges.tif",overwrite=TRUE)
writeRaster(pp_an_50m_86_10,file="export_final_crop_l93_25m/pp_digitalis_8610_25_vosges.tif",overwrite=TRUE)
writeRaster(raster_bdv2_25m,file="export_final_crop_l93_25m/bdv2_25_vosges.tif",overwrite=TRUE)
writeRaster(raster_foret_pub_25m,file="export_final_crop_l93_25m/foret_pub_25_vosges.tif",overwrite=TRUE)
writeRaster(raster_distance_25m,file="export_final_crop_l93_25m/distance_edge_25_vosges.tif",overwrite=TRUE)
writeRaster(cap_25_tpi_1500,file="export_final_crop_l93_25m/cap_index_25_1500.tif",overwrite=TRUE)
writeRaster(cap_25_tpi_1000,file="export_final_crop_l93_25m/cap_index_25_1000.tif",overwrite=TRUE)
writeRaster(ipv_1000_50,file="export_final_crop_l93_25m/ipv_1000_25.tif",overwrite=TRUE)

setwd("C:/Users/borderieux/OneDrive - agroparistech.fr/Docs/data_microclim")


mnt_25_vosges_crop<-raster( "export_final_crop_l93_25m/mnt_25_vosges.tif" )
#mnt_25_tpi<-raster( "export_final_crop_l93_25m/tpi_25_vosges.tif" )
#mnt_25_tri<-raster( "export_final_crop_l93_25m/tri_25_vosges.tif" )
mnt_25_slope<-raster( "export_final_crop_l93_25m/slope_25_vosges.tif" )
mnt_25_aspect<-raster( "export_final_crop_l93_25m/aspect_NS_25_vosges.tif" )
canopy_gedi_crop<-raster( "export_final_crop_l93_25m/canopy_height_25_vosges.tif" )
tmoy_an_50m_86_10<-raster( "export_final_crop_l93_25m/tmoy_digitalis_8610_25_vosges.tif" )
pp_an_50m_86_10<-10*raster( "export_final_crop_l93_25m/pp_digitalis_8610_25_vosges.tif" )
raster_bdv2_25m<-raster( "export_final_crop_l93_25m/bdv2_25_vosges.tif" )
raster_foret_pub_25m<-raster( "export_final_crop_l93_25m/foret_pub_25_vosges.tif" )
raster_distance_25m<-raster( "export_final_crop_l93_25m/distance_edge_25_vosges.tif" )
ipv_1000_25m<-raster( "export_final_crop_l93_25m/ipv_1000_25.tif" )
cap_25_tpi_1500<-raster("export_final_crop_l93_25m/cap_index_25_1500.tif")
cover_2010_crop_25m<-raster("export_final_crop_l93_25m/canopy_cover_2010_25m_vosges.tif")

#cover_2010_crop_25m<-projectRaster(cover_2010,canopy_gedi_crop)
#writeRaster(cover_2010_crop_25m,file="export_final_crop_l93_25m/canopy_cover_2010_25m_vosges.tif",overwrite=TRUE)

mapview(cover_2010_crop_25m)
plot(ipv_1000_25m)

mask_forest_vosges<-trim(mask(raster_bdv2_25m,ser_d11_vosges))
ext_mask_foret_vosges<-extent(mask_forest_vosges)

mask_forest_vallee<-trim(mask(mask_forest_vosges,vallee_vosge[vallee_vosge$Vallee%in%c("Thann"),]))
ext_mask_foret_valle<-extent(mask_forest_vallee)

stack_all_variable<-stack(list(mnt_25_vosges_crop,mnt_25_aspect,ipv_1000_25m,cap_25_tpi_1500,mnt_25_slope,
                               canopy_gedi_crop,cover_2010_crop_25m,tmoy_an_50m_86_10,
                               pp_an_50m_86_10,raster_bdv2_25m,
                               raster_foret_pub_25m,raster_distance_25m))



stack_all_variable_vosges<-crop(stack_all_variable,ext_mask_foret_vosges)
stack_all_variable_vosges<-mask(stack_all_variable_vosges,mask_forest_vosges)


reclass_matrix_microclim_forest<-matrix(c(1:30,4,rep(1,8),rep(2,14),rep(3,2),rep(4,4),1),ncol = 2)
stack_all_variable_vosges$bdv2_25_vosges<-reclassify(stack_all_variable_vosges$bdv2_25_vosges,reclass_matrix_microclim_forest)


stack_all_variable_vallee<-stack_all_variable
stack_all_variable_vallee<-crop(stack_all_variable_vallee,ext_mask_foret_valle)
stack_all_variable_vallee<-mask(stack_all_variable_vallee,mask_forest_vallee)

stack_all_variable_vallee$bdv2_25_vosges<-reclassify(stack_all_variable_vallee$bdv2_25_vosges,reclass_matrix_microclim_forest)

plot(stack_all_variable_vallee$ipv_1000_25)
plot(stack_all_variable_vallee$mnt_25_vosges)
# mnt_25_tpi_250_vallee<-crop(mnt_25_tpi_250,ext_mask_foret_valle)
# mnt_25_tpi_250_vallee<-mask(mnt_25_tpi_250_vallee,mask_forest_vallee)
# 
# mnt_25_cap_150_vallee<-crop(cap_25_tpi_150,ext_mask_foret_valle)
# mnt_25_cap_150_vallee<-mask(mnt_25_cap_150_vallee,mask_forest_vallee)
# 
# mnt_25_cap_250_vallee<-crop(cap_25_tpi_250,ext_mask_foret_valle)
# mnt_25_cap_250_vallee<-mask(mnt_25_cap_250_vallee,mask_forest_vallee)
# 
# mnt_25_cap_500_vallee<-crop(cap_25_tpi_500,ext_mask_foret_valle)
# mnt_25_cap_500_vallee<-mask(mnt_25_cap_500_vallee,mask_forest_vallee)
# 
# stack_all_variable_vallee$mnt_25_tpi_250_vallee<-mnt_25_tpi_250_vallee
# stack_all_variable_vallee$mnt_25_cap_150_vallee<-mnt_25_cap_150_vallee
# stack_all_variable_vallee$mnt_25_cap_250_vallee<-mnt_25_cap_250_vallee
# stack_all_variable_vallee$mnt_25_cap_500_vallee<-mnt_25_cap_500_vallee


plot(stack_all_variable_vallee)
plot(stack_all_variable_vallee$slope_25_vosges)


summary(stack_all_variable_vosges)
summary(stack_all_variable_vallee)

for(i in 1:dim(stack_all_variable_vallee)[3]){
  
  
writeRaster(stack_all_variable_vallee[[i]],file=paste0("export_final_mask_valley/",names(stack_all_variable_vallee[[i]]),".tif"),overwrite=TRUE)

}


library(terra)
#### classification for natural sampling
stack_all_variable_vallee<-stack(paste0("export_final_mask_valley/",list.files("export_final_mask_valley")))

# aspect_NS_25_vosges<-cos(terrain(mnt_25_vosges_crop,opt="aspect"))
# aspect_EW_25_vosges<-sin(terrain(mnt_25_vosges_crop,opt="aspect"))
# 
# aspect_NS_25_vosges<-mask(crop(aspect_NS_25_vosges,stack_all_variable_vallee$bdv2_25_vosges),stack_all_variable_vallee$bdv2_25_vosges)
# aspect_EW_25_vosges<-mask(crop(aspect_EW_25_vosges,stack_all_variable_vallee$bdv2_25_vosges),stack_all_variable_vallee$bdv2_25_vosges)
# 
# writeRaster(aspect_NS_25_vosges,file=paste0("export_final_mask_valley/","aspect_NS_25_vosges",".tif"),overwrite=TRUE)
# writeRaster(aspect_EW_25_vosges,file=paste0("export_final_mask_valley/","aspect_EW_25_vosges",".tif"),overwrite=TRUE)
# 


plot(aspect_EW_25_vosges)

min_elev<-cellStats(stack_all_variable_vallee$mnt_25_vosges,stat="min")
max_elev<-cellStats(stack_all_variable_vallee$mnt_25_vosges,stat="max")

min_elev<-quantile(stack_all_variable_vallee$mnt_25_vosges,probs=0.025)
max_elev<-quantile(stack_all_variable_vallee$mnt_25_vosges,probs=0.975)


from_elev<-seq(from=300,to=1350,by=80)
to_elev<-seq(from=350,to=1400,by=80)



from_elev<-seq(from=min_elev,to=max_elev,length.out = 8)
by_elev<-from_elev[3]-from_elev[2]
to_elev<-seq(from=from_elev[2],to=from_elev[length(from_elev)]+by_elev,by=by_elev)

mid_value_elev<-from_elev+by_elev/2


label_elev<-paste(floor(from_elev),"-",floor(to_elev),sep="")
id_elev<-1:(length(from_elev))

ref_elev<-data.table(id_elev,label_elev,mid_value_elev)

reclass_matrix_elev_quali<-matrix(c(from_elev,to_elev,id_elev),ncol=3)
reclass_matrix_elev_quali<-rbind(reclass_matrix_elev_quali,matrix(c(0,min(from_elev)+1,NA),ncol=3))
reclass_matrix_elev_quali_2<-matrix(c(15,2000,NA),ncol=3)

stack_all_variable_vallee$elev_classification<-reclassify(stack_all_variable_vallee$mnt_25_vosges,reclass_matrix_elev_quali)
stack_all_variable_vallee$elev_classification<-reclassify(stack_all_variable_vallee$elev_classification,reclass_matrix_elev_quali_2)

plot(stack_all_variable_vallee$mnt_25_cap_500_vallee)


#reclass_matrix_elev_foret_pub<-matrix(c(0,1,NA,NA,1,1),ncol=2)

#stack_all_variable_vallee$foret_pub_25_vosges<-reclassify(stack_all_variable_vallee$foret_pub_25_vosges,reclass_matrix_elev_foret_pub)

getwd()# "\\\\tsclient/C/Users/borderieux/Desktop/data_microclim"


plot(stack_all_variable_vallee$elev_classification)

plot(stack_all_variable_vallee$ipv_1000_25)


stack_all_variable_vallee<-stack(paste0("export_final_mask_valley/",list.files("export_final_mask_valley")))




all_vallee_cell<-data.table(xyFromCell(stack_all_variable_vallee,1:ncell(stack_all_variable_vallee)))
all_vallee_cell[,ident:=1:nrow(all_vallee_cell)]



all_vallee_cell_value<-cbind(all_vallee_cell,data.table(extract(stack_all_variable_vallee,1:ncell(stack_all_variable_vallee))))
#all_vallee_cell_value<-all_vallee_cell_value[!is.na(mnt_25_vosges),]
#all_vallee_cell_value<-all_vallee_cell_value[!is.na(mnt_25_cap_500_vallee),]
#all_vallee_cell_value<-all_vallee_cell_value[!is.infinite(mnt_25_cap_500_vallee),]

all_vallee_cell_value<-all_vallee_cell_value[canopy_height_25_vosges>4,]
all_vallee_cell_value<-all_vallee_cell_value[canopy_height_25_vosges<35 ,]

all_vallee_cell_value$bdv2_25_vosges<-as.factor(all_vallee_cell_value$bdv2_25_vosges)

all_vallee_cell_value<-all_vallee_cell_value[bdv2_25_vosges!=4,]

all_vallee_cell_value_dt<-all_vallee_cell_value


all_vallee_cell_value<-as.data.frame(all_vallee_cell_value)
rownames(all_vallee_cell_value)<-all_vallee_cell_value$ident


plot(all_vallee_cell_value$mnt_25_vosges[1:100000],all_vallee_cell_value$tmoy_digitalis_8610_25_vosges[1:100000])
plot(all_vallee_cell_value$mnt_25_vosges[1:100000],all_vallee_cell_value$pp_digitalis_8610_25_vosges[1:100000])
plot(all_vallee_cell_value$mnt_25_vosges,all_vallee_cell_value$ipv_1000_25,cex=0.25,pch=16)
plot(all_vallee_cell_value$canopy_cover_2010_25m_vosges,all_vallee_cell_value$canopy_height_25_vosges,cex=0.25,pch=16)

cor(all_vallee_cell_value$mnt_25_vosges,all_vallee_cell_value$ipv_1000_25,use = "complete.obs",method="spearman")
cor(all_vallee_cell_value$canopy_cover_2010_25m_vosges,all_vallee_cell_value$canopy_height_25_vosges,use = "complete.obs",method="spearman")


summary(lm(all_vallee_cell_value$tmoy_digitalis_8610_25_vosges[1:100000]~all_vallee_cell_value$aspect_NS_25_vosges[1:100000]))


summary(all_vallee_cell_value)


all_vallee_cell_value_dt[,.N,by=.(elev_classification,foret_pub_25_vosges)]
all_vallee_cell_value_dt[,mean(mnt_25_vosges),by=.(elev_classification,foret_pub_25_vosges)]


quantile(all_vallee_cell_value_dt$canopy_height_25_vosges,probs=c(0,0.025,0.05,0.25,0.5,0.75,0.95,0.975))

hist(all_vallee_cell_value_dt$canopy_height_25_vosges)
hist(all_vallee_cell_value_dt$aspect_NS_25_vosges)
hist(all_vallee_cell_value_dt$slope_25_vosges)
hist(all_vallee_cell_value_dt$ipv_1000_25)

all_vallee_cell_value_dt<-data.table(all_vallee_cell_value)

all_vallee_cell_value_dt[,high_cano:=ifelse(canopy_height_25_vosges>25,1,0)]
all_vallee_cell_value_dt[,expo_south:=ifelse(aspect_NS_25_vosges<= -0.71,"S",ifelse(aspect_NS_25_vosges>= 0.71,"N",NA))]
all_vallee_cell_value_dt[,expo_east:=ifelse(aspect_EW_25_vosges< -0.71,"W",ifelse(aspect_EW_25_vosges> 0.71,"E",NA))]

#all_vallee_cell_value_dt[,expo_east:=ifelse(aspect_NS_25_vosges<= -0.7,"S",ifelse(aspect_NS_25_vosges>= 0.7,"N",NA))]

all_vallee_cell_value_dt[,cap_pot:=ifelse(ipv_1000_25<= 200,1,ifelse(ipv_1000_25>=800,2,0))]
all_vallee_cell_value_dt[,slope_quanti:=ifelse(slope_25_vosges<= 10,"low",ifelse(slope_25_vosges>=25,"hight","moderate"))]


all_vallee_cell_value_dt<-merge(all_vallee_cell_value_dt,ref_elev,by.x="elev_classification",by.y="id_elev",all.x=T)

all_vallee_cell_value_dt[foret_pub_25_vosges==1 &!is.na(expo_south)& !is.na(elev_classification),.N,by=.(elev_classification,high_cano,expo_south)][order(elev_classification,high_cano,expo_south)]
all_vallee_cell_value_dt[foret_pub_25_vosges==1 &!is.na(expo_south)& !is.na(elev_classification),.N,by=.(elev_classification,cap_pot)][order(elev_classification)][order(elev_classification,cap_pot)]

summary(all_vallee_cell_value_dt_sample)

all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt[!is.na(foret_pub_25_vosges) & !is.na(elev_classification)]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[distance_edge_25_vosges>51,]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[slope_25_vosges<35,]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[canopy_height_25_vosges>10,]


all_vallee_cell_value_dt_sample[,foret_pub_25_vosges:=1]


all_vallee_cell_value_dt_sample[,dif_mid:=abs(mid_value_elev-mnt_25_vosges)]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[order(dif_mid),]

all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[dif_mid<10,]

merge(all_vallee_cell_value_dt_sample[,ifelse(.N==0,0,.N),by=.(elev_classification,high_cano,expo_south)][order(elev_classification,high_cano,expo_south)],ref_elev[,c("id_elev","label_elev")],by.x="elev_classification",by.y="id_elev",all.x=T)
merge(all_vallee_cell_value_dt_sample[cap_pot!=0,.N,by=.(elev_classification,cap_pot)][order(elev_classification)][order(elev_classification,cap_pot)],ref_elev[,c("id_elev","label_elev")],by.x="elev_classification",by.y="id_elev",all.x=T)
all_vallee_cell_value_dt_sample[,expo_all:=ifelse(!is.na(expo_south),expo_south,ifelse(!is.na(expo_east),expo_east,NA))]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[!is.na(expo_all),]

all_vallee_cell_value_dt_sample[,.N,by=.(elev_classification,label_elev,slope_quanti)][order(elev_classification),]
all_vallee_cell_value_dt_sample[cap_pot==0 & slope_quanti!="moderate",.N,by=.(elev_classification,label_elev,slope_quanti)][order(elev_classification),]
all_vallee_cell_value_dt_sample[slope_quanti=="moderate"&cap_pot!=0,.N,by=.(elev_classification,label_elev,cap_pot)][order(elev_classification,cap_pot),]
all_vallee_cell_value_dt_sample[slope_quanti=="moderate"&cap_pot!=0,.N,by=.(elev_classification,label_elev,cap_pot,expo_south,expo_east)][order(elev_classification,cap_pot,expo_south,expo_east),]
all_vallee_cell_value_dt_sample[slope_quanti=="moderate"&cap_pot%in%c(0,2),.N,by=.(elev_classification,label_elev,expo_south,high_cano)][order(elev_classification,expo_south,high_cano),]
all_vallee_cell_value_dt_sample[slope_quanti=="moderate"&cap_pot%in%c(0,2),.N,by=.(elev_classification,label_elev,expo_south,high_cano)][order(elev_classification,expo_south,high_cano),]


sample_point<-function(avcv_dts){
  selected=NULL
  
  for(i in 1:8){
    if(i%in%c(1,5,8))desired_expo<-"S"## to have the extreme
    if(i%in%c(2,6,7))desired_expo<-"N"
    if(i%in%c(3))desired_expo<-"E"
    if(i%in%c(4))desired_expo<-"W"
    slope_test<-avcv_dts[elev_classification==i&expo_all==desired_expo & cap_pot==0 & slope_quanti!="moderate",sample(ident,1),by=.(elev_classification,label_elev,slope_quanti)]$V1
    selected<-c(selected,slope_test)
  }
 # return(avcv_dts[ident%in%selected,])
  # slope_test<-avcv_dts[cap_pot==0 & slope_quanti!="moderate",sample(ident,1),by=.(elev_classification,label_elev,slope_quanti)]$V1
  # selected<-c(selected,slope_test)
  # 
  
  for(i in 1:6){
  if(i%in%c(1))desired_expo<-"S"## to have the extreme
  if(i%in%c(2))desired_expo<-"N"
  if(i%in%c(3,5))desired_expo<-"E"
  if(i%in%c(4,6))desired_expo<-"W"
  cap_ipv_test_all_but_7_8<- avcv_dts[expo_all==desired_expo& elev_classification==i & !ident%in%selected & slope_quanti=="moderate"&cap_pot!=0,sample(ident,1),by=.(elev_classification,label_elev,cap_pot)]$V1
  selected<-c(selected,cap_ipv_test_all_but_7_8)
  }
  
  cap_ipv_test_7<- avcv_dts[ expo_all=="S" & elev_classification%in%c(7) & cap_pot==2 & !ident%in%selected & slope_quanti=="moderate",sample(ident,1),by=.(elev_classification,label_elev,cap_pot)]$V1
  selected<-c(selected,cap_ipv_test_7)
  cap_ipv_test_8<- avcv_dts[ expo_all=="N" & elev_classification%in%c(8) & cap_pot==2 & !ident%in%selected & slope_quanti=="moderate",sample(ident,1),by=.(elev_classification,label_elev,cap_pot)]$V1
  selected<-c(selected,cap_ipv_test_8)
  
  
  
  canop_expo_test_all_but_8<- avcv_dts[!is.na(high_cano) &!is.na(expo_south)& elev_classification!=8 & !ident%in%selected & slope_quanti=="moderate"&cap_pot%in%c(0),sample(ident,1),by=.(elev_classification,label_elev,expo_south,high_cano)]$V1
  selected<-c(selected,canop_expo_test_all_but_8)
  canop_expo_8<- avcv_dts[!is.na(high_cano) & !is.na(expo_south)& elev_classification==8 & !ident%in%selected & slope_quanti=="moderate"&cap_pot%in%c(0,2),sample(ident,1),by=.(elev_classification,label_elev,expo_south,high_cano)]$V1
  selected<-c(selected,canop_expo_8)
  
  res<-avcv_dts[ident%in%selected,]
   return(res)
  
  
}


tst<-all_vallee_cell_value_dt_sample[expo_all=="E"&cap_pot==0 & slope_quanti!="moderate",.N,by=.(elev_classification,label_elev,slope_quanti,expo_all)][order(elev_classification,slope_quanti),]

sample_point(all_vallee_cell_value_dt_sample)

points_to_sample<-all_vallee_cell_value_dt_sample[all_vallee_cell_value_dt_sample[,.I[1],by=.(elev_classification,high_cano,expo_south)]$V1]
all_vallee_cell_value_dt_sample_reduced<-all_vallee_cell_value_dt_sample[!ident%in%points_to_sample$ident &cap_pot!=0,]
points_to_sample<-rbind(points_to_sample,all_vallee_cell_value_dt_sample_reduced[all_vallee_cell_value_dt_sample_reduced[,.I[1],by=.(elev_classification,cap_pot)]$V1])

summary(points_to_sample)
hist(points_to_sample$slope_25_vosges)
plot(stack_all_variable_vallee)

set.seed(1)
# points_to_sample<-all_vallee_cell_value_dt_sample[all_vallee_cell_value_dt_sample[,sample(ident,1),by=.(elev_classification,high_cano,expo_south)]$V1]
# all_vallee_cell_value_dt_sample_reduced<-all_vallee_cell_value_dt_sample[!ident%in%points_to_sample$ident &cap_pot!=0,]
# points_to_sample<-rbind(points_to_sample,all_vallee_cell_value_dt_sample_reduced[all_vallee_cell_value_dt_sample_reduced[,sample(.I,1),by=.(elev_classification,cap_pot)]$V1])
#
points_to_sample<-sample_point(all_vallee_cell_value_dt_sample)
points_to_sample[order(elev_classification)]
points_to_sample_sf<-st_as_sf(points_to_sample,coords=c("x","y"),crs=st_crs(2154))


dist_point<-(st_distance(points_to_sample_sf))
minmax<-mean(apply(dist_point,1,function(x)min(x[as.numeric(x)!=0])))
min(dist_point[as.numeric(dist_point)!=0])

points_to_sample_final<-points_to_sample
minmax_final<-mean(apply(dist_point,1,function(x)min(x[x!=0])))
dist_point_final<-dist_point

dist_point<-as.numeric(st_distance(points_to_sample_sf))
dist_point<-dist_point[dist_point!=0]
### optimize distance between

set.seed(1)

for(i in 1:10000){
 # points_to_sample<-all_vallee_cell_value_dt_sample[all_vallee_cell_value_dt_sample[,sample(.I,1),by=.(elev_classification,high_cano,expo_south)]$V1]
 # all_vallee_cell_value_dt_sample_reduced<-all_vallee_cell_value_dt_sample[!ident%in%points_to_sample$ident &cap_pot!=0,]
 # points_to_sample<-rbind(points_to_sample,all_vallee_cell_value_dt_sample_reduced[all_vallee_cell_value_dt_sample_reduced[,sample(.I,1),by=.(elev_classification,cap_pot)]$V1])
  points_to_sample<-sample_point(all_vallee_cell_value_dt_sample)
  
  
  
  points_to_sample_sf<-st_as_sf(points_to_sample,coords=c("x","y"),crs=st_crs(2154))
  dist_point<-(st_distance(points_to_sample_sf))
  minmax<-mean(apply(dist_point,1,function(x)min(x[x!=0])))
  
  #dist_point<-dist_point[dist_point!=0]
  
  #if(mean(dist_point)<mean(dist_point_final) &min(dist_point)>min(dist_point_final) ){
    if(minmax > minmax_final & as.numeric( min(dist_point[as.numeric(dist_point)!=0])) > 250){
    points_to_sample_final<-points_to_sample
    dist_point_final<-dist_point
    minmax_final<-minmax
  }
  
  cat(paste0(i,"-"))
  
}


mean(dist_point_final)
min(dist_point_final)
min(dist_point_final[as.numeric(dist_point_final)!=0])
summary(points_to_sample_final)
points_to_sample_final_sf<-st_as_sf(points_to_sample_final,coords=c("x","y"),crs=st_crs(2154))
points_to_sample_final_sf$foret_communale<-st_intersection(points_to_sample_final_sf,foret_pub)$llib_frt
points_to_sample_final$foret_communale<-st_intersection(points_to_sample_final_sf,foret_pub)$llib_frt

foret_pub_points<-st_intersection(points_to_sample_final_sf,foret_pub)
commu<-names(summary(as.factor(foret_pub_points$llib_frt)))

mapview::mapview(points_to_sample_final_sf,zcol="mnt_25_vosges")
mapview(list(foret_pub[foret_pub$llib_frt%in%commu,],points_to_sample_final_sf),label=c("llib_frt","ident"),zcol=c("llib_frt","mnt_25_vosges"))

points_to_sample_final[,c("ident_class","foret_communale")]
summary_village<-points_to_sample_final[,.(N_sensors=.N,names_site=paste(ident_class,collapse = " "),elev=paste(elev_classification,collapse = " ")),by=foret_communale][order(-N_sensors),]

points_to_sample_final$ident_class
kept_var_sf<-c("ident","elev_classification","label_elev","high_cano","cap_pot","slope_quanti","expo_all",
               "foret_communale","bdv2_25_vosges" ,"mnt_25_vosges","ipv_1000_25","slope_25_vosges","tmoy_digitalis_8610_25_vosges","pp_digitalis_8610_25_vosges", 
               "canopy_cover_2010_25m_vosges" ,"canopy_height_25_vosges","distance_edge_25_vosges")
kep_var<-c("ident","x","y","elev_classification","label_elev","high_cano","cap_pot","slope_quanti","expo_all",
            "foret_communale","bdv2_25_vosges" ,"mnt_25_vosges","ipv_1000_25","slope_25_vosges","tmoy_digitalis_8610_25_vosges","pp_digitalis_8610_25_vosges", 
            "canopy_cover_2010_25m_vosges" ,"canopy_height_25_vosges","distance_edge_25_vosges")


kep_var_2<-c("ident","x","y","elev_classification","label_elev","high_cano","cap_pot","slope_quanti","expo_all",
           "bdv2_25_vosges" ,"mnt_25_vosges","ipv_1000_25","slope_25_vosges","tmoy_digitalis_8610_25_vosges","pp_digitalis_8610_25_vosges", 
           "canopy_cover_2010_25m_vosges" ,"canopy_height_25_vosges","distance_edge_25_vosges")
points_to_sample_final_sf[,kept_var_sf]

write_sf(points_to_sample_final_sf[,kept_var_sf],"export_shp_sampling_final","sampling_points_microclim_vosges",driver="ESRI Shapefile")
write_sf(st_buffer(points_to_sample_final_sf[,kept_var_sf],12.5,endCapStyle = "SQUARE"),"export_shp_sampling_final","sampling_sites_microclim_vosges",driver="ESRI Shapefile")

backup_site<-st_as_sf(all_vallee_cell_value_dt_sample[,..kep_var_2],coords=c("x","y"),crs=st_crs(2154))
write_sf(backup_site,"export_shp_sampling_final","all_potential_sites",driver="ESRI Shapefile")
write_sf(st_buffer(backup_site,12.5,endCapStyle = "SQUARE"),"export_shp_sampling_final","all_potential_sites_pixels",driver="ESRI Shapefile")



write.table(as.data.frame(points_to_sample_final)[,kep_var],row.names = F,sep=";",file = "export_shp_sampling_final/sample_points_.csv")

write.table(as.data.frame(summary_village),row.names = F,sep=";",file = "export_shp_sampling_final/communale_forests.csv")


mapview(points_to_sample_final_sf,zcol="mnt_25_vosges",label="ident_class")

mapview(st_as_sf(all_vallee_cell_value_dt_sample,coords=c("x","y"),crs=st_crs(2154)),zcol="mnt_25_vosges")
mapview(st_buffer(points_to_sample_final_sf,30,endCapStyle = "SQUARE"),zcol="mnt_25_vosges")
mapview(points_to_sample_sf,zcol="mnt_25_vosges")
mapview(selection_sf,zcol="mnt_25_vosges")

hist(points_to_sample_final$slope_25_vosges,nc=20)
hist(selection_finale$slope_25_vosges,nc=20)
hist(points_to_sample_final$ipv_1000_25,nc=20)
hist(selection_finale$ipv_1000_25,nc=20)

sort(points_to_sample_final$mnt_25_vosges)
sort(selection_finale$mnt_25_vosges)

sort(points_to_sample_final$slope_25_vosges)
sort(selection_finale$slope_25_vosges)

foret_pub_points<-st_intersection(points_to_sample_final_sf,foret_pub)
commu<-names(summary(as.factor(foret_pub_points$llib_frt)))

mapview(list(foret_pub[foret_pub$llib_frt%in%commu,],points_to_sample_final_sf),label=c("llib_frt","ident"),zcol=c("llib_frt","mnt_25_vosges"))
mapview(list(foret_pub[foret_pub$llib_frt%in%commu,]),label=c("llib_frt"),zcol=c("llib_frt"))
mapview(list(foret_pub[foret_pub$llib_frt%in%commu,]),label=c("llib_frt"),zcol=c("llib_frt"),col.regions=rainbow(18))


save_map<-mapview(list(forets_communales=foret_pub[foret_pub$llib_frt%in%commu,],Points=points_to_sample_final_sf),label=c("llib_frt","ident_class"),zcol=c("llib_frt","mnt_25_vosges"),col.regions=list(rainbow(18),grey(seq(from=0.9,to=0.1,length.out = 8))),alpha.regions=list(0.8,1))
save_map

mapview(list(forets_communales=foret_pub[foret_pub$llib_frt%in%commu,],Points=points_to_sample_final_sf),label=c("llib_frt","ident_class"),zcol=c("llib_frt","mnt_25_vosges"),col.regions=list(rainbow(18),"white"),alpha.regions=list(0.8,1))

mapshot(save_map ,url="communale_forest_map.html",selfcontained = TRUE)
remotes::install_github("r-spatial/mapview")
webshot::install_phantomjs()
library(mapview)
mapview



hist(points_to_sample_final$canopy_height_25_vosges)
hist(selection_finale$canopy_height_25_vosges)


hist(selection_finale$slope_25_vosges)
summary(points_to_sample_final)
summary(selection_finale)
summary(points_to_sample)
nrow(points_to_sample)

plot(stack_all_variable_vallee$mnt_25_cap_500_vallee)

library(mapview)

library(ade4)

library(factoextra)

colnames(all_vallee_cell_value)
pca_variable<-c("mnt_25_vosges","ipv_1000_25","slope_25_vosges","canopy_height_25_vosges","aspect_NS_25_vosges")

pca_valle_1st<-dudi.pca(df = all_vallee_cell_value[,pca_variable], scannf = FALSE, nf = 3)

mix_variable<-c(pca_variable,"bdv2_25_vosges")


pca_valle_1st<-dudi.mix(df = all_vallee_cell_value[,mix_variable], scannf = FALSE, nf = 3)

cor(all_vallee_cell_value[,c(pca_variable,"slope_25_vosges")])

inertia.dudi(pca_valle_1st)


s.arrow(pca_valle_1st$co,1,2)
s.arrow(pca_valle_1st$co,2,3)
fviz_pca_var(pca_valle_1st,axes=c(1,2),col.var = "contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE     )
fviz_pca_var(pca_valle_1st,axes=c(1,3),col.var = "contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE     )
fviz_pca_var(pca_valle_1st,axes=c(2,3),col.var = "contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE     )

coord_pca<-pca_valle_1st$li

ggplot(coord_pca,mapping=aes(x=Axis1,y= Axis2))+theme_bw()+geom_point(size=0.45,alpha=0.1)+geom_point(data=coord_pca[all_vallee_cell_value$ident%in%points_to_sample_final$ident,],mapping=aes(x=Axis1,y= Axis2),size=1.6,col="red")
ggplot(coord_pca,mapping=aes(x=Axis1,y= Axis2))+theme_bw()+geom_point(size=0.45,alpha=0.1)+geom_point(data=coord_pca[all_vallee_cell_value$ident%in%selection_finale$ident,],mapping=aes(x=Axis1,y= Axis2),size=1.6,col="red")

ggplot(coord_pca,mapping=aes(x=Axis2,y= Axis3))+theme_bw()+geom_point(size=0.45,alpha=0.1)+geom_point(data=coord_pca[all_vallee_cell_value$ident%in%points_to_sample_final$ident,],mapping=aes(x=Axis2,y= Axis3),size=1.6,col="red")
ggplot(coord_pca,mapping=aes(x=Axis2,y= Axis3))+theme_bw()+geom_point(size=0.45,alpha=0.1)+geom_point(data=coord_pca[all_vallee_cell_value$ident%in%selection_finale$ident,],mapping=aes(x=Axis2,y= Axis3),size=1.6,col="red")

cor(all_vallee_cell_value)

summary(all_vallee_cell_value)
ggplot(melt(all_vallee_cell_value),aes(x=value))+theme_bw()+geom_histogram()+facet_wrap(~variable,scales = "free")


sqrt(pca_valle_1st$li)






#### sqrt root correction of the pc coordinates
all_vallee_cell_value_end<-cbind(all_vallee_cell_value,pca_valle_1st$li)
dim(all_vallee_cell_value)

Min<-all_vallee_cell_value_end[all_vallee_cell_value_end[,"Axis1"]<=0,]
Plus<-all_vallee_cell_value_end[all_vallee_cell_value_end[,"Axis1"]>0,]
Plus[,"Axis1"]<-sqrt(Plus[,"Axis1"])
Min[,"Axis1"]<-(-sqrt(-Min[,"Axis1"]+0.0001))
all_vallee_cell_value_end<-rbind(Min,Plus)


Min<-all_vallee_cell_value_end[all_vallee_cell_value_end[,"Axis2"]<=0,]
Plus<-all_vallee_cell_value_end[all_vallee_cell_value_end[,"Axis2"]>0,]
Plus[,"Axis2"]<-sqrt(Plus[,"Axis2"])
Min[,"Axis2"]<-(-sqrt(-Min[,"Axis2"]+0.0001))
all_vallee_cell_value_end<-rbind(Min,Plus)


Min<-all_vallee_cell_value_end[all_vallee_cell_value_end[,"Axis3"]<=0,]
Plus<-all_vallee_cell_value_end[all_vallee_cell_value_end[,"Axis3"]>0,]
Plus[,"Axis1"]<-sqrt(Plus[,"Axis3"])
Min[,"Axis3"]<-(-sqrt(-Min[,"Axis3"]+0.0001))
all_vallee_cell_value_end<-rbind(Min,Plus)




#### selection along the pca axis with cube 


Selection<-c()

Points<-64 #number of locations to select. 
#Note: algorithm does NOT cycle until all these points have been achieved, 
#although it has two optional clauses to select additional locations till this number is reached (see below)

Groups<-ceiling(Points^0.33) #0.33 is used, as we are using environmental space as a cube, which we want to cut in little cubes,
#with selection of points in each cube. 0.33 gives a reasonable selection, but number can be changed to have a coarser or more strict selection

#divide each PCA axis in the same number of bins with equal sizes along each axis
Steps_Axis1<-(max(all_vallee_cell_value_end$Axis1)-min(all_vallee_cell_value_end$Axis1))/Groups
Steps_Axis2<-(max(all_vallee_cell_value_end$Axis2)-min(all_vallee_cell_value_end$Axis2))/Groups
Steps_Axis3<-(max(all_vallee_cell_value_end$Axis3)-min(all_vallee_cell_value_end$Axis3))/Groups


set.seed(40)


for (i in 1:Groups) {
  #hierarchical selection procedure, first based on dimension 1, 
  #then on dim 2 and 3, each time within the groups of Dim. 1
  
  #select all points within the i'st bin of axis 1
  Selection1 <- all_vallee_cell_value_end[all_vallee_cell_value_end$Axis1>(min(all_vallee_cell_value_end$Axis1)+(i-1)*Steps_Axis1) & 
                             all_vallee_cell_value_end$Axis1<(min(all_vallee_cell_value_end$Axis1)+i*Steps_Axis1),]
  for (j in 1:Groups) {
    #select all points from the i'st bin that fall in the j'st bin of axis 2
    Selection2 <- Selection1[Selection1$Axis2>(min(Selection1$Axis2)+(j-1)*Steps_Axis2) & 
                               Selection1$Axis2<(min(Selection1$Axis2)+j*Steps_Axis2),]
    for (k in 1:Groups) {
      #select all points in that j'st bin that fall within the k'st bin of axis 3
      Selection3 <- Selection2[Selection2$Axis3>(min(Selection2$Axis3)+(k-1)*Steps_Axis3) & 
                                 Selection2$Axis3<(min(Selection2$Axis3)+k*Steps_Axis3),]
      #randomly select 3 points within that cube with dimensions i,j,k
      #depending on distribution of locations in the PCA space (and size of the bins suggested earlier), 
      #there will be more or less cubes with dimension i,j,k that have no points, where thus no points will be selected
      Selection <- rbind(Selection,Selection3[sample(1:nrow(Selection3), 1),])
    }
  }
  
  if(nrow(Selection[complete.cases(Selection[,1:5]),])<(i*Points/Groups)){ 
    #if not enough potential sampling points in the cubes, 
    #we can fill out with random selections from the main category, as many as are in there
    #this is optional and not done in the paper, as it only matters if a fixed amount of locations is required
    if((i*Points/Groups)-nrow(Selection[complete.cases(Selection[,1:5]),])< nrow(Selection1)){
      Selection4 <- Selection1[sample(1:nrow(Selection1), (i*Points/Groups)-nrow(Selection[complete.cases(Selection[,1:5]),])),]
      Selection <- rbind(Selection, Selection4)
    }
    else {
      Selection4<-Selection1
      Selection <- rbind(Selection, Selection4)
    }
  }
}


Selection<-Selection[complete.cases(Selection[,1:5]),]
if(nrow(Selection[complete.cases(Selection[,1:5]),])<Points){ 
  #if  still not enough points to go to the final requested number, just fill up with random points
  Selection4 <- all_vallee_cell_value_end[sample(1:nrow(all_vallee_cell_value_end), Points-nrow(Selection[complete.cases(Selection[,1:5]),])),]
  Selection <- rbind(Selection, Selection4)
}




dim(Selection[complete.cases(Selection),])


#♥Selection$ident<-as.numeric(rownames(Selection))
summary(Selection)
selection_finale<-merge(Selection,all_vallee_cell,by="ident")


Selection<-Selection[complete.cases(Selection),]
selection_finale<-Selection
selection_sf<-st_as_sf(selection_finale,coords=c("x","y"),crs=st_crs(2154))

library(mapview)


mapview(selection_sf,zcol="tmoy_digitalis_8610_25_vosges")
mapview(selection_sf,zcol="mnt_25_vosges")
mapview(selection_sf,zcol="canopy_height_25_vosges")
mapview(selection_sf,zcol="aspect_NS_25_vosges")
mapview(selection_sf,zcol="bdv2_25_vosges")

mapview(points_to_sample_final_sf,zcol="mnt_25_vosges")
n_step<-10

n_step*2*2+2*n_step

n_step*2*2+2*n_step

summary(stack_all_variable_vallee)

library(terra)

xyFromCell(mask_forest_vallee,1:ncell(mask_forest_vallee))


plot(mask_forest_vallee)


mapview(mask_forest_vallee)

summary(as.factor(vallee_vosge$Vallee))

#### manipulaion  ####

if(Sys.info()["nodename"]=="CALCULUS") setwd(file.path("D:","quercus","BD_SIG","climat","france")) else setwd(file.path("S:","BD_SIG","climat","france"))
current_climate<-raster(file.path("DIGITALIS_v2","1999_2014","1km","tmoy_v2_1","tmoy_9914_an","w001001.adf"))
current_climate_v2<-raster(file.path("DIGITALIS_v2b","tmoy_v2b","tmoy_0919_13.tif"))

setwd(file.path("~","AdapFor_analysis"))

crs_wgs84<-crs(canopy_gedi)
current_climate_v2_newcrs<-projectRaster(current_climate_v2,crs=crs_wgs84)
extent_france_wgs84<-extent(current_climate_v2_newcrs)
crs_l93<-crs(current_climate_v2)
crs_l2<-crs(current_climate)

canopy_gedi_2<-crop(canopy_gedi,extent_france_wgs84)

writeRaster(canopy_gedi_2,"C:/Users/borderieux/Desktop/chesla_micro_macroclim_data/canopy/Forest_height_2019_NAFR_fance_wgs84.tif",overwrite=T)

plot(canopy_gedi_2)

canopy_gedi<-projectRaster(canopy_gedi,crs=crs(past_climate))



vallee_vosge<-read_sf(file.path("Raw_data","GIS","Field_work_planning"),layer="Vallee_plaine")
extent<-st_bbox(vallee_vosge)
extent_lamb2<-st_bbox(st_transform(vallee_vosge,crs=st_crs(  27572 )))
extent_vosges<-extent(x=extent[1],xmax=extent[3],ymin=extent[2],ymax=extent[4])

if(Sys.info()["nodename"]=="CALCULUS") setwd("D:/quercus/BD_SIG/topo") else setwd("S:/BD_SIG/topo")


mnt_25<-raster(file.path("BDalti25","BDALTIV2_2-0_25M_ASC_LAMB93-IGN69_D068_2017-08-29","BDALTIV2_merge_068.tif"))
mnt_25<-projectRaster(mnt_25,crs=crs(past_climate))
mnt_25<-crop(mnt_25,extent_vosges_lamb2)
plot(mnt_25)
setwd(file.path("~","AdapFor_analysis"))

canopy_gedi_vosges<-crop(canopy_gedi,extent_vosges_lamb2)
canopy_gedi_france<-crop(canopy_gedi,extent(past_climate))
plot(canopy_gedi_vosges)


#### field work ####


library(mapview)

tst<-points_to_sample_final[cap_pot%in%c(0,2) & slope_quanti== "moderate"& elev_classification%in%8,]

points_to_sample_final$ident_class<-NA
points_to_sample_final[,ident_class:=ifelse(slope_quanti!="moderate",ifelse(slope_quanti=="hight",paste0("slope_H_",elev_classification),paste0("slope_L_",elev_classification)),ident_class)]
points_to_sample_final[,ident_class:=ifelse(slope_quanti=="moderate"& cap_pot==0 & elev_classification%in%1:7,ifelse(high_cano==1,paste0("H_",expo_all,"_",elev_classification),paste0("L_",expo_all,"_",elev_classification)),ident_class)]

## special case for the elev class8, because it was lacking plots with cap_pot==0, so they are cap_pot==2
points_to_sample_final[,ident_class:=ifelse(slope_quanti=="moderate"& cap_pot%in%c(0,2) & elev_classification==8,ifelse(high_cano==1,paste0("H_",expo_all,"_",elev_classification),paste0("L_",expo_all,"_",elev_classification)),ident_class)]
points_to_sample_final[,ident_class:=ifelse(is.na(ident_class),ifelse(cap_pot==1,paste0("CAD_H_",elev_classification),paste0("CAD_L_",elev_classification)),ident_class)]

# du to the lack of low cad points ih the highest elevation class (8),two sites have the same caracteristics,
# that is aspect north, low canopy, elev_class= 8 and cad ==2   (low cad potential)
# they will be renamed to avoid confusions later on

points_to_sample_final[ident_class=="L_N_8","ident_class"]<-c("L_N_8_1","L_N_8_2")


sort(points_to_sample_final$ident_class)


points_to_sample_final_sf<-st_as_sf(points_to_sample_final,coords=c("x","y"),crs=st_crs(2154))
mapview(points_to_sample_final_sf,zcol="mnt_25_vosges")
mapview(points_to_sample_final_sf[points_to_sample_final_sf$ident_class=="L_N_8",],zcol="mnt_25_vosges")
mapview(st_buffer(points_to_sample_final_sf,12.5,endCapStyle = "SQUARE"),zcol="mnt_25_vosges")

export_gpx<-points_to_sample_final_sf[,"geometry"]
export_gpx<-st_transform(export_gpx,crs=st_crs(4326 ))
write_sf(export_gpx,"export_shp_sampling_final/sampling_points_no_name_f.gpx",driver="GPX",overwrite=T)
file.remove("export_shp_sampling_final/sampling_points.gpx")

#### visualize the new forest temps map ####

map_forest_temp<-raster("C:/Users/borderieux/Downloads/ForestBIO1.tif")

border_vosges_shp<-read_sf("shp",layer="Vosges_border_shp")
vallee_vosge<-read_sf("shp",layer="Vallee_plaine")
border_vosges_shp<-st_transform(border_vosges_shp,st_crs(crs(map_forest_temp)))
vallee_vosge<-st_transform(vallee_vosge,st_crs(crs(map_forest_temp)))

extent<-st_bbox(border_vosges_shp)
extent_vosges<-extent(x=extent[1],xmax=extent[3],ymin=extent[2],ymax=extent[4])

extent<-st_bbox(vallee_vosge)
extent_vallee<-extent(x=extent[1],xmax=extent[3],ymin=extent[2],ymax=extent[4])
extent_vallee<-extent(x=4090000,xmax=4110000,ymin=2740000,ymax=2770000)


map_forest_temp_vosges<-crop(map_forest_temp,extent_vosges)
map_forest_temp_vallee<-crop(map_forest_temp,extent_vallee)

plot(map_forest_temp_vallee)
mapview::mapview(map_forest_temp_vallee)

#### heat load and twi computation  ####
setwd("C:/Users/borderieux/OneDrive - agroparistech.fr/Docs/data_microclim")
mnt_25_vosges_crop<-raster( "export_final_crop_l93_25m/mnt_25_vosges.tif" )
mnt_25_vosges_crop<-crop(mnt_25_vosges_crop,stack_all_variable_vallee)
library(terra)
library(spatialEco)

mnt_25_vosges_crop<-rast(mnt_25_vosges_crop)

HL_25_vosges<-hli(mnt_25_vosges_crop)
plot(HL_25_vosges)

HL_25_vosges<-raster(HL_25_vosges)

HL_25_vosges<-mask(HL_25_vosges,stack_all_variable_vallee$bdv2_25_vosges)

writeRaster(HL_25_vosges,file="export_final_crop_l93_25m/Heat_load_index_25.tif",overwrite=TRUE)



library(dynatopmodel)
a.atb <- upslope.area(mnt_25_vosges_crop, atb=TRUE)
a.atb <- upslope.areatst(mnt_25_vosges_crop, atb=TRUE)
sp::plot(a.atb, main=c("Upslope area (log(m^2/m))", "TWI log(m^2/m)"))

TWI_25_vosges<-mask(a.atb$atb,stack_all_variable_vallee$bdv2_25_vosges)
plot(TWI_25_vosges)
hist(TWI_25_vosges)
writeRaster(TWI_25_vosges,file="export_final_crop_l93_25m/Topographic_wetness_index_25.tif",overwrite=TRUE)

#### twi function####
upslope.areatst <- function(dem, log=TRUE, atb=FALSE, deg=0.1, fill.sinks=TRUE)
{
  # check
  if(xres(dem)!=yres(dem))
  {
    stop("Raster has differing x and y cell resolutions. Check that it is in a projected coordinate system (e.g. UTM) and use raster::projectRaster to reproject to one if not. Otherwise consider using raster::resample")
  }
  # any sinks still present may give strange results
  #  sink(file="e:/junk/sink.txt")
  #  on.exit(sink(NULL))
  if(fill.sinks)
  {
    # use capture.output to supress the function console output
    capture.output(dem <- invisible(raster::setValues(dem, topmodel::sinkfill(raster::as.matrix(dem), res=xres(dem), degree=deg))))
  }
  topidx <- topmodel::topidx(raster::as.matrix(dem), res=xres(dem))
  
  a <- raster::setValues(dem, topidx$area)
  if(log)
  {
    a <- log(a)
  }
  if(atb)
  {
    atb <- raster::setValues(dem, topidx$atb)
    # add the topographic index ln(a/tanB)
    a <- addLayer(a, atb)
    names(a)<-c("a", "atb")
  }
  return(a)
}

# use function from TOPMODEL
flow.lens <- function(dem,
                      src=NULL,  #  starting cells, defaults to all cells i dem
                      agg=1, # initial aggregation factor
                      max.agg=4,
                      
                      outlet=NULL)  # A vector containing the row and column indices of the pixel representing the catchment outlet. if a single value, treated as the index of a DEM cell
{
  lens <- raster::setValues(dem, NA)
  if(length(outlet)>0)
  {
    outlet.sp <- xyFromCell(dem, outlet, spatial=TRUE)
  }
  
  dem.agg <- dem
  while(agg <= max.agg & max(c(0,lens[]), na.rm=TRUE)==0)
  {
    if(agg>1)
    {
      message("Trying a aggregated dem to determine flow lengths")
      # try a coarser
      dem.agg <- raster::aggregate(dem, agg)
      #	reaches <- aggregate(reaches, )
    }
    
    if(length(outlet)>0)
    {
      outlet <- extract.cells(dem.agg, outlet.sp)
      iout<- rowColFromCell(dem.agg, outlet)
    }
    else{iout <- NA}
    
    dem.agg <- fill.sinks(dem.agg, deg=0.1)
    lens <- raster::setValues(dem.agg, flowlength(as.matrix(dem.agg), outlet=iout))
    if(!is.null(src))
    {
      lens[setdiff(1:ncell(dem), src)]<- NA
    }
    agg <- agg+1
  }
  agg <- agg-1
  # disaggregate
  if(agg>1)
  {
    message("Disaggregating back to original resolution...")
    lens <- raster::disaggregate(lens, agg, method="bilinear")
  }
  return(raster::xres(dem)*lens)
}

extract.cells <- function(dem, drn,...)
{
  target <- extract(dem, drn, cellnumbers=TRUE,...)
  if(is.list(target))
  {
    return(do.call(rbind, target)[,1])
  }
  else
  {
    return(target[,1])
  }
}


fill.sinks <- function(dem, deg=0.01,
                       silent=TRUE,
                       ipass=1,   # perform sinkfill a maximum of this times or until all sinks filled
                       fail.if.not.complete=FALSE)
{
  DEM <- as.matrix(dem)
  res <- xres(dem)
  # stopifnot(is(DEM, "matrix"))
  # if (min(as.vector(DEM[!is.na(DEM)])) < -9000)
  #    stop("DEM contains unrealistic values (< -9000)")
  #   DEM[is.na(DEM)] <- -9999
  #  nrow <- dim(DEM)[1]
  #  ncol <- dim(DEM)[2]
  
  i <- 1
  sinks.remain <- TRUE
  while(i <= ipass & sinks.remain)
  {
    prev <- DEM
    #  DEM[is.na(DEM)] <- -9999
    capture.output(DEM <- topmodel::sinkfill(DEM, res, deg))
    diff <- sum(DEM[]-prev[], na.rm=TRUE)
    if(diff==0)
    {
      # there are definitely no sinks left now
      sinks.remain <- FALSE
    }
    #
    #       .C("sinkfill", PACKAGE = "topmodel", as.double(DEM),
    #                  result = double(nrow * ncol + 2), as.integer(nrow), as.integer(ncol),
    #                  as.double(res), as.double(deg))$result
    #    result[result > 999998] <- NA
    #     DEM <- matrix(result[3:(nrow * ncol + 2)], nrow = nrow)
    #     # 100 is max number of iterations, so if reached thsi then have to run the sinkfill again
    #     if(result[1]< 100 & result[1]>0)
    #     {
    #       sinks.remain <- FALSE
    #
    #     }
    #     else if(!silent)
    #     {
    #       cat("Sinkfill pass #", i, " No. of sinks remaining = ", result[2], "\n")
    #
    #     }
    
    i <- i + 1
  }
  #   if (result[1] == -1)
  #   {
  #     warning("incomplete sink removal")
  #   }
  #   if (result[1] == 100)
  #   {
  #     msg <- paste("Maximum no. iterations reached (100). No. sinks remaining=", result[2])
  #     message(msg)
  #   }
  
  #  mat <- matrix(result[3:(nrow * ncol + 2)], nrow = nrow)
  
  return(setValues(dem, DEM))
}


#### plot a PCA of the sampled plot in the whole space ####
library(ade4)
library(factoextra)

colnames(all_vallee_cell_value)
dt_raster
dt_raster[,tree_density_2018:=tree_density_projected]
dt_raster[,tree_density_2018:=ifelse(is.na(mnt_25_vosges),NA,tree_density_2018)]


pca_variable<-c("mnt_25_vosges","ipv_25","slope_25_vosges","tree_density_2018","Heat_load_index_25")

pca_variable_col<-c(pca_variable,"site_ID")


pca_data<-na.omit(dt_raster[,..pca_variable])
pca_data[,site_ID:=NA]


coords_pca<-coords[site_ID!="new1",]
pca_data<-rbind(pca_data,coords_pca[,..pca_variable_col])


pca_valle_1st<-dudi.pca(df = pca_data[,..pca_variable], scannf = FALSE, nf = 3)

pca_valle_1st$co

score(pca_valle_1st)

res.var <- get_pca_var(pca_valle_1st)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2    


inertia(pca_valle_1st)

coord_pca<-data.table(pca_valle_1st$li)

unique(all_tms_data_agg_fit$locality_id)
coord_pca[,site:=(pca_data$site_ID)]

coord_pca[,state:=ifelse(is.na(site),"no",ifelse(site%in% unique(all_tms_data_agg_fit$locality_id),"kept","Malfunction"))]


a12<-ggplot(coord_pca[is.na(site)],aes(x=Axis1,y=Axis2))+
  geom_point(alpha=0.025, size=0.4)+
  geom_point(data=coord_pca[!is.na(site)],aes(color=state),size=1.5)+
  scale_color_manual(breaks=c("kept","Malfunction"),values=c("deepskyblue1","firebrick2"),label=c("Kept loggers","lost or discarded\nloggers"))+
  theme_bw()+labs(color="Logger state")

a3<-ggplot(coord_pca[is.na(site)],aes(x=Axis1,y=Axis3))+
  geom_point(alpha=0.025, size=0.4)+
  geom_point(data=coord_pca[!is.na(site)],aes(color=state),size=1.5)+
  scale_color_manual(breaks=c("kept","Malfunction"),values=c("deepskyblue1","firebrick2"),label=c("Kept loggers","lost or discarded\nloggers"))+
  theme_bw()+labs(color="Logger state")


exportpca<-ggarrange(plotlist = list(a12,a3),nrow = 2,common.legend = T,legend = "bottom",labels = c("a)","b)"))


ggsave("pca.png",exportpca+theme(plot.background = element_rect(color=NA,fill = "white")),unit="mm",width=180,height = 180,dpi=180)


