

map_pred_topo<-ggplot(dt_raster,aes(x,y,fill=pred_topo-2.3+0.4))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,
                       limits=quantile(dt_raster$pred_topo-2.3+0.4,probs=c(0.05,0.95),na.rm=T))+
  #scale_fill_gradientn(colors = brewer.pal(7, "PuOr")[order(7:1)],na.value ="transparent",oob=scales::squish)+
  #scale_fill_gradient2(low="lightskyblue3",mid="grey99",high="#0080C2FF",midpoint = 22,na.value ="transparent",oob=scales::squish)+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Refroidissement par\nla topographie (°C)")+
  theme(legend.position = c(0.88,0.78),legend.text = element_text(size=11),legend.title = element_text(size=12.5))



map_pred_canopy<-ggplot(dt_raster,aes(x,y,fill=pred_canopy_2+1.84 ))+
  theme_classic(base_size =14)+
  geom_raster()+
  coord_fixed()+
  
 scale_fill_gradientn(colors=c("#5E4FA2","#5E4FA2","#3288BD","#3288BD","#66C2A5","#ABDDA4", "#E6F598","#FFFFBF",
                                        "#FEE08B","#FDAE61","#F46D43" ,"#D53E4F"),na.value ="transparent",
                                        oob=scales::squish,limits=quantile(dt_raster$pred_canopy_2+1.84,probs=c(0.075,0.925),na.rm=T))+
  #scale_fill_gradient2(low="#5E4FA2",mid="#FFFFBF",high="#9E0142",midpoint = -0.9,na.value ="transparent",oob=scales::squish,limits=quantile(dt_raster$pred_canopy_2+1.1,probs=c(0.075,0.925),na.rm=T))+
  

  #scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,limits=quantile(dt_raster$pred_canopy_2+1.1,probs=c(0.075,0.925),na.rm=T))+
  labs(fill="var")+
  geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
  geom_sf(data=square_inset,fill=NA,inherit.aes=F,color="firebrick3",linewidth= 1)+
  coord_sf(expand=F)+
  annotation_scale()+ 
  labs(x="",y="",fill="Refroidissement\npar la forêt (°C)")+
  theme(legend.position = c(0.88,0.78),legend.text = element_text(size=11),legend.title = element_text(size=12.5))

map_pred_canopy


create_mini_map<-function(what){
  
  mini_map_1<-ggplot(subset_map_dt,aes(x,y,fill=get(what) ))+
    theme_void()+
    geom_raster(show.legend = F)+
    coord_fixed()+
   # scale_fill_viridis_c(na.value ="transparent",oob=scales::squish,limits=quantile(dt_raster[,..what],probs=c(0.075,0.925),na.rm=T))+
  
    scale_fill_gradientn(colors=c("#5E4FA2","#5E4FA2","#3288BD","#3288BD","#66C2A5","#ABDDA4", "#E6F598","#FFFFBF",
                                           "#FEE08B","#FDAE61","#F46D43" ,"#D53E4F"),na.value ="transparent",oob=scales::squish,limits=quantile(dt_raster[,..what],probs=c(0.075,0.925),na.rm=T))+
                                             
      labs(fill="var")+
    # geom_sf(data=valley_border,inherit.aes=F,linewidth=0.4,fill=NA,alpha=0,color="white")+
    coord_sf(expand=F)+
    #annotation_scale()+ 
    theme(panel.border = element_rect(fill=NA,color="grey20",linewidth =3) , plot.margin = margin(10,10,10,10))+
    labs(x="",y="",fill="Temperature °C")
  return(mini_map_1)
  
}


create_mini_map("pred_canopy_2")

out_2_panel<-ggarrange(plotlist=list(map_pred_topo,map_pred_canopy),
                       ncol=2,nrow=1,labels = c("a)","b)","c)","d)"),align = "hv",
                       font.label = list(size = 18, color = "grey5", face = "bold", family = NULL))




ggsave(file.path("Figure_result_fr","Microclimat_salon_agri.jpg"),
       out_2_panel,
       dpi=600,unit="mm",width = 180,height = 80,scale=1.55)

ggsave( file.path("Figure_result_fr","mini_topo.jpg")   ,
           )




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



create_mini_map("pred_topo")
