#### pckg ####
library(raster)
library(data.table)
library(sf)

#### start ####

###lexicon:
## mnt_25_vosges elevation (masl) from DEM 
## bdv2 is the description of type of forest
## IPV: metric of topographic position ranging from 0 to 1
## canopy_height_25_vosges : canopy height from GEDI
## aspect_NS_25_vosges north - south aspect, from -1 to 1
## aspect_EW_25_vosges east - west aspect, from -1 to 1
## slope_25_vosges slope in degree
## 25 stands for the 25 meter resolutino of the raster cells


##import the raster
stack_all_variable_vallee<-stack(file.path("data","processed_raster",list.files(file.path("data","processed_raster"))))

## get copernicus data
tree_density<-raster(file.path("data","tree_density_raster","tree_cover_density_thur_2018_copernicus.tif"))
tree_density<-projectRaster(tree_density,crs=crs(stack_all_variable_vallee))
tree_density_reprojected<-projectRaster(tree_density,stack_all_variable_vallee)
stack_all_variable_vallee$tree_density_projected<-tree_density_reprojected

## create the elevation classes
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
## end of creation of elevation classes

## we create a data.table  all_vallee_cell_value_dt that will hold the relevant variable 
all_vallee_cell<-data.table(xyFromCell(stack_all_variable_vallee,1:ncell(stack_all_variable_vallee)))
all_vallee_cell[,ident:=1:nrow(all_vallee_cell)]
all_vallee_cell_value_dt<-cbind(all_vallee_cell,data.table(extract(stack_all_variable_vallee,1:ncell(stack_all_variable_vallee))))

rm(all_vallee_cell)

## filter extreme values of canopy height
all_vallee_cell_value_dt<-all_vallee_cell_value_dt[canopy_height_25_vosges>4,]
all_vallee_cell_value_dt<-all_vallee_cell_value_dt[canopy_height_25_vosges<35 ,]

all_vallee_cell_value_dt$bdv2_25_vosges<-as.factor(all_vallee_cell_value_dt$bdv2_25_vosges)

## we remove open forests
all_vallee_cell_value_dt<-all_vallee_cell_value_dt[bdv2_25_vosges!=4,]

## specification of the classes used for the stratified sampling
all_vallee_cell_value_dt[,high_cano:=ifelse(canopy_height_25_vosges>25,1,0)]# breakpoint translate to 90% tree density
all_vallee_cell_value_dt[,expo_south:=ifelse(aspect_NS_25_vosges<= -0.71,"S",ifelse(aspect_NS_25_vosges>= 0.71,"N",NA))]
all_vallee_cell_value_dt[,expo_east:=ifelse(aspect_EW_25_vosges< -0.71,"W",ifelse(aspect_EW_25_vosges> 0.71,"E",NA))]

## ipv: metric of topographic position (see manuscript), cap mean cold air pooling potential, expo = aspect
all_vallee_cell_value_dt[,cap_pot:=ifelse(ipv_25<= 0.2,1,ifelse(ipv_25>=0.8,2,0))]
all_vallee_cell_value_dt[,slope_quanti:=ifelse(slope_25_vosges<= 10,"low",ifelse(slope_25_vosges>=25,"hight","moderate"))]

all_vallee_cell_value_dt<-merge(all_vallee_cell_value_dt,ref_elev,by.x="elev_classification",by.y="id_elev",all.x=T)

## we keep plot within public forest, without an excessive slope, without a small canopy height, which translate to open forest
## and we remove pixel too close to the forest edge to avoid edge effects
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt[!is.na(foret_pub_25_vosges) & !is.na(elev_classification)]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[distance_edge_25_vosges>51,]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[slope_25_vosges<35,]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[canopy_height_25_vosges>10,]

all_vallee_cell_value_dt_sample[,foret_pub_25_vosges:=1]


all_vallee_cell_value_dt_sample[,dif_mid:=abs(mid_value_elev-mnt_25_vosges)]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[order(dif_mid),]
## only plots 10 m appart from the mean of an elevatino classes are kept 
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[dif_mid<10,]
##we remove plot with undefined aspect
all_vallee_cell_value_dt_sample[,expo_all:=ifelse(!is.na(expo_south),expo_south,ifelse(!is.na(expo_east),expo_east,NA))]
all_vallee_cell_value_dt_sample<-all_vallee_cell_value_dt_sample[!is.na(expo_all),]

## this function generate a stratified sampling of 60 points across the elevation classes
## special conditions are 
sample_point<-function(avcv_dts){
  selected=NULL
  
  ## loop that create the sampling accross low and high slope (2 points per elevation class)
  for(i in 1:8){
    if(i%in%c(1,5,8))desired_expo<-"S"## to have the extreme
    if(i%in%c(2,6,7))desired_expo<-"N"
    if(i%in%c(3))desired_expo<-"E"
    if(i%in%c(4))desired_expo<-"W"
    slope_test<-avcv_dts[elev_classification==i&expo_all==desired_expo & cap_pot==0 & slope_quanti!="moderate",sample(ident,1),by=.(elev_classification,label_elev,slope_quanti)]$V1
    selected<-c(selected,slope_test)
  }
  
  ## for the 6 fist elevation classes, loop that sample low and high cold air pooling potential regardless of canopy, with preset aspect
  ## (2 points per elevation class)
  for(i in 1:6){
  if(i%in%c(1))desired_expo<-"S"## to have the extreme
  if(i%in%c(2))desired_expo<-"N"
  if(i%in%c(3,5))desired_expo<-"E"
  if(i%in%c(4,6))desired_expo<-"W"
  cap_ipv_test_all_but_7_8<- avcv_dts[expo_all==desired_expo& elev_classification==i & !ident%in%selected & slope_quanti=="moderate"&cap_pot!=0,sample(ident,1),by=.(elev_classification,label_elev,cap_pot)]$V1
  selected<-c(selected,cap_ipv_test_all_but_7_8)
  }
  
  ## special selection is performed for the class 7 and 8 because of the lack of valley bottom in these classes
  ## (1 ponit instead of 2)
  cap_ipv_test_7<- avcv_dts[ expo_all=="S" & elev_classification%in% 7 & cap_pot==2 & !ident%in%selected & slope_quanti=="moderate",sample(ident,1),by=.(elev_classification,label_elev,cap_pot)]$V1
  selected<-c(selected,cap_ipv_test_7)
  cap_ipv_test_8<- avcv_dts[ expo_all=="N" & elev_classification%in% 8 & cap_pot==2 & !ident%in%selected & slope_quanti=="moderate",sample(ident,1),by=.(elev_classification,label_elev,cap_pot)]$V1
  selected<-c(selected,cap_ipv_test_8)
  
  ## sampling of low and high canopy height, with a south and a north facing slope (4 points per elevation class)
  ## in moderated cold air pooling potential
  canop_expo_test_all_but_8<- avcv_dts[!is.na(high_cano) &!is.na(expo_south)& elev_classification!=8 & !ident%in%selected & slope_quanti=="moderate"&cap_pot%in%c(0),sample(ident,1),by=.(elev_classification,label_elev,expo_south,high_cano)]$V1
  selected<-c(selected,canop_expo_test_all_but_8)
  
  ## exception for the highest elevation to also include ridges (aka low cold air pooling potential), otherwise not points are found
  canop_expo_8<- avcv_dts[!is.na(high_cano) & !is.na(expo_south)& elev_classification==8 & !ident%in%selected & slope_quanti=="moderate"&cap_pot%in%c(0,2),sample(ident,1),by=.(elev_classification,label_elev,expo_south,high_cano)]$V1
  selected<-c(selected,canop_expo_8)
  
  res<-avcv_dts[ident%in%selected,]
   return(res)
  
}



set.seed(1)

## a test with one sampling
points_to_sample<-sample_point(all_vallee_cell_value_dt_sample)
points_to_sample[order(elev_classification)]
points_to_sample_sf<-st_as_sf(points_to_sample,coords=c("x","y"),crs=st_crs(2154))
plot(points_to_sample_sf[,c("mnt_25_vosges")])
table(points_to_sample$expo_all)
table(points_to_sample$slope_quanti) # 2 low slope points are missing, this condition was not met for every class
table(points_to_sample$cap_pot) 


## This loop run 10 000 times and compute the mean minimum distance between points and try to maximize it with each iteration
## We keep the set of points with the larger mean minimum distance between points
## points should also be at least 250 meters appart
minmax_final<-230

set.seed(1)

for(i in 1:10000){
  points_to_sample<-sample_point(all_vallee_cell_value_dt_sample)
  points_to_sample_sf<-st_as_sf(points_to_sample,coords=c("x","y"),crs=st_crs(2154))
  dist_point<-st_distance(points_to_sample_sf)
  minmax<-mean(apply(dist_point,1,function(x)min(x[x!=0])))
  cat(paste0(i,"-"))
  
  if(minmax > minmax_final & as.numeric(min(dist_point[as.numeric(dist_point)!=0])) > 250){
    points_to_sample_final<-points_to_sample
    dist_point_final<-dist_point
    minmax_final<-minmax
    print(paste0("New maximum found : ",round(minmax_final,1)," m"))
  }
  
}


mean(dist_point_final[as.numeric(dist_point_final)!=0])
min(dist_point_final[as.numeric(dist_point_final)!=0])

points_to_sample_final_sf<-st_as_sf(points_to_sample_final,coords=c("x","y"),crs=st_crs(2154))

library(mapview)
mapview(points_to_sample_final_sf,zcol="mnt_25_vosges")


#### plot a PCA of the sampled plot in the whole space ####
library(ade4)
library(factoextra)

## for this part of the script we use dt_raster from the main script, DT_raster summarize only the pixel studied wit hthe microclimate model

dt_raster[,tree_density_2018:=tree_density_projected]
dt_raster[,tree_density_2018:=ifelse(is.na(mnt_25_vosges),NA,tree_density_2018)]


pca_variable<-c("mnt_25_vosges","ipv_25","slope_25_vosges","tree_density_2018","Heat_load_index_25")
pca_variable_col<-c(pca_variable,"site_ID")

pca_data<-na.omit(dt_raster[,..pca_variable])
pca_data[,site_ID:=NA]


coords_pca<-coords[site_ID!="new1",] # object from main script
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
coord_pca <- cbind(coord_pca,pca_valle_1st)
unique(all_tms_data_agg_fit$locality_id)
coord_pca[,site:=(pca_data$site_ID)]

coord_pca[,state:=ifelse(is.na(site),"no",ifelse(site%in% unique(all_tms_data_agg_fit$locality_id),"kept","Malfunction"))]

a12<-ggplot(coord_pca[is.na(site)],aes(x=Axis1,y=Axis2))+
  geom_point(alpha=0.025, size=0.4,aes(color = ""))+
  geom_point(data=coord_pca[!is.na(site)],aes(color=state),size=1.5)+
  scale_color_manual(breaks=c("kept","Malfunction"),values=c("deepskyblue1","firebrick2"),label=c("Kept loggers","lost or discarded\nloggers"))+
  theme_bw()+labs(color="Logger state")

a3<-ggplot(coord_pca[is.na(site)],aes(x=Axis1,y=Axis3))+
  geom_point(alpha=0.025, size=0.4)+
  geom_point(data=coord_pca[!is.na(site)],aes(color=state),size=1.5)+
  scale_color_manual(breaks=c("kept","Malfunction"),values=c("deepskyblue1","firebrick2"),label=c("Kept loggers","lost or discarded\nloggers"))+
  theme_bw()+labs(color="Logger state")


exportpca<-ggarrange(plotlist = list(a12,a3),nrow = 2,common.legend = T,legend = "bottom",labels = c("a)","b)"))


ggsave(file.path("figure_result","pca.png"),
       exportpca+theme(plot.background = element_rect(color=NA,fill = "white")),
       unit="mm",width=180,height = 180,dpi=180)
