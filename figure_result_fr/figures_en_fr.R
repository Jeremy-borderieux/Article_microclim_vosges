## fr

map_elevation_climate<-ggplot(dt_raster,aes(x,y,fill=pred_elev+ get_coef_mean["tree_density_2018",1]*90))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster$pred_elev+ get_coef_mean["tree_density_2018",1]*90,probs=c(0.025,0.975),na.rm=T))+
  #scale_fill_gradientn(colors = brewer.pal(7, "PuOr")[order(7:1)],na.value ="transparent",oob=scales::squish)+
  # scale_fill_gradient2(low="lightskyblue3",mid="grey99",high="#0080C2FF",midpoint = 22,na.value ="transparent",oob=scales::squish)+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Temperature\nmoyenne °C")+ 
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))



map_pred_topo<-ggplot(dt_raster,aes(x,y,fill=pred_topo ))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=c(-1.5,-0.15),breaks=c(-0.15,-0.5,-1,-1.5))+
  #scale_fill_gradientn(colors = brewer.pal(7, "PuOr")[order(7:1)],na.value ="transparent",oob=scales::squish)+
  #scale_fill_gradient2(low="lightskyblue3",mid="grey99",high="#0080C2FF",midpoint = 22,na.value ="transparent",oob=scales::squish)+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Refroidissement \ntopoclimatique °C")+
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))


dt_raster[,pred_canopy_2:=ifelse(is.na(pred_T_mean),NA,pred_canopy)]
dt_raster[,pred_canopy_2:=ifelse(pred_canopy_2> -1.5,-1.5,pred_canopy_2)]


map_pred_canopy<-ggplot(dt_raster,aes(x,y,fill=pred_canopy_2 ))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster$pred_canopy_2,probs=c(0.025,0.975),na.rm=T))+
  #scale_fill_gradientn(colors = brewer.pal(7, "PuOr")[order(7:1)],na.value ="transparent",oob=scales::squish)+
  #scale_fill_gradient2(low="lightskyblue3",mid="grey99",high="#0080C2FF",midpoint = 22,na.value ="transparent",oob=scales::squish)+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Refroidissement \nmicroclimatique °C")+
  theme(legend.position = c(0.88,0.75),legend.text = element_text(size=11),legend.title = element_text(size=12.5))



square_inset<-st_buffer(st_as_sf(data.table(x=1000125,y=6764575),coords=c("x","y"),crs=st_crs(2154)),1000,endCapStyle="SQUARE")
subset_map_dt<-dt_raster[x%between% c(999125-100,999625+1500)& y %between%c(6765075-1500,6766075-500) ,]


create_mini_map<-function(what){
  
  mini_map_1<-ggplot(subset_map_dt,aes(x,y,fill=get(what) ))+
    theme_void()+
    geom_raster(show.legend = F)+
    coord_fixed()+
    scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                         limits=quantile(dt_raster[,..what],probs=c(0.025,0.975),na.rm=T))+
    labs(fill="var")+
    # geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
    coord_sf(expand=F)+
    #annotation_scale()+ 
    theme(panel.border = element_rect(fill=NA,color="grey20",linewidth =3) , plot.margin = margin(10,10,10,10))+
    labs(x="",y="",fill="Temperature °C")
  return(mini_map_1)
  
}


minimaps<-ggarrange(plotlist = list( create_mini_map("pred_elev"),create_mini_map("pred_topo"),create_mini_map("pred_canopy_2")),ncol=2,nrow=2,align = "hv")


out_4_panel<-ggarrange(plotlist=list(map_elevation_climate,map_pred_topo,map_pred_canopy,minimaps),
                       ncol=2,nrow=2,labels = c("a)","b)","c)","d)"),align = "hv",
                       font.label = list(size = 18, color = "grey5", face = "bold", family = NULL))



ggsave(file.path("Figure_result_fr","climate_map_2.jpg"),
       out_4_panel,
       dpi=200,unit="mm",width = 180,height = 160,scale=1.5)


#### boxplot  ####


unique(vege_plot_meta_data$cut_topo)
levels(vege_plot_meta_data$cut_topo)<-c("[-1.80 : -1.25°]\nTopoclimat froid","[-1.25 : -0.90 °]\nTopoclimat modéré","[-0.90 : -0.40 °]\nTopoclimat chaud")

my_comparisons <- list( levels(vege_plot_meta_data$cut_topo)[1:2], 
                        levels(vege_plot_meta_data$cut_topo)[2:3], 
                        levels(vege_plot_meta_data$cut_topo)[c(1,3)] )



RS_plot<-ggplot(vege_plot_meta_data,aes(x=cut_topo,y=RS_pH_corrected,fill=cut_topo))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="grey20",size=1.5,show.legend = F)+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  coord_cartesian(ylim=c(-2,68))+
  labs(y="Richesse spécifique",x="Catégories de topoclimat")+
  stat_compare_means(comparisons = my_comparisons, label.y = c(55,58,61),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns")))


CTI_plot<-ggplot(vege_plot_meta_data[n_sp_climplant>5,],aes(x=cut_topo,y=cit_climplant_pH_corrected,fill=cut_topo))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="grey20",size=1.5,show.legend = F)+
  scale_y_continuous(limits=c(6.5,9.65))+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  labs(y="Affinité climatique \nde la communauté (°C)",x="Catégories de topoclimat")+
  stat_compare_means(comparisons = my_comparisons ,label.y = c(9,9.25,9.5),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns")))



Export_boxplot<-ggarrange(plotlist=list(RS_plot,CTI_plot),ncol=2,labels=c("a)","b)"))

ggsave(plot=Export_boxplot,file.path("figure_result_fr","Figure_boxplot.jpg"),width=180,height =80 ,scale=1.25,unit="mm",dpi=400)




unique(vege_plot_meta_data$cut_topo)
levels(vege_plot_meta_data$cut_topo)<-c("[-1.55 : -1.0°]\nTopoclimat froid","[-1.0 : -0.65 °]\nTopoclimat modéré","[-0.65 : -0.15 °]\nTopoclimat chaud")

my_comparisons <- list( levels(vege_plot_meta_data$cut_topo)[1:2], 
                        levels(vege_plot_meta_data$cut_topo)[2:3], 
                        levels(vege_plot_meta_data$cut_topo)[c(1,3)] )

RS_plot<-ggplot(vege_plot_meta_data[,],aes(x=cut_topo,y=RS_pH_corrected,fill=cut_topo))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="grey20",size=1.5,show.legend = F)+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F)+
  #geom_violin(alpha=0.5,show.legend = F,draw_quantiles = 0.5)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  coord_cartesian(ylim=c(-2,68))+
  labs(y="Richesse spécifique",x="Catégories de topoclimat")+
  stat_compare_means(comparisons = my_comparisons, label.y = c(55,58,61),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns")))


CTI_plot<-ggplot(vege_plot_meta_data[,],aes(x=cut_topo,y=cit_climplant_pH_corrected,fill=cut_topo))+
  theme_bw()+
  geom_jitter(alpha=0.35,pch=21,color="grey20",size=1.5,show.legend = F)+
  scale_y_continuous(limits=c(6.49,9.65))+
  geom_boxplot(alpha=0.5, outlier.alpha = 0,show.legend = F,draw_quantiles = 0.5)+
  #geom_violin(alpha=0.5, outlier.alpha = 0,show.legend = F,draw_quantiles = 0.5)+
  scale_fill_manual(values=c("#482173FF","#22A884FF","#FDE725FF"))+ 
  labs(y="Affinité climatique \nde la communauté (°C)",x="Catégories de topoclimat")+
  stat_compare_means(comparisons = my_comparisons ,label.y = c(9,9.25,9.45),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns")))




Export_boxplot<-ggarrange(plotlist=list(RS_plot,CTI_plot),ncol=2,labels=c("a)","b)"))

ggsave(plot=Export_boxplot,file.path("figure_result_fr","Figure_boxplot.jpg"),width=180,height =80 ,scale=1.25,unit="mm",dpi=400)
