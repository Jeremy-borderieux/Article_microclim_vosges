#### read flora ####

vege_plot_meta_data<-fread(file.path("data","flora_data_metadata","plot_data_2022_final.csv"))
vege_plot_meta_data_sf<-st_as_sf(vege_plot_meta_data,coords=c("X_L93","Y_L93"),crs=st_crs(2154))

vege_pred<-data.table(extract(raster_of_prediction,vege_plot_meta_data_sf))
vege_pred$plot_ID<-vege_plot_meta_data$plot_ID
vege_pred<-vege_pred[!is.na(predict_microclim_mean),]

vege_covariable<-data.table(extract(stack_all_variable_vallee,vege_plot_meta_data_sf[vege_plot_meta_data_sf$plot_ID%in%vege_pred$plot_ID,]))
vege_covariable$plot_ID<-vege_pred$plot_ID


vege_plot_meta_data<-merge(vege_plot_meta_data,vege_pred,by="plot_ID")
vege_plot_meta_data<-merge(vege_plot_meta_data,vege_covariable,by="plot_ID")


table(vege_plot_meta_data$`precision localisation`)

indicator_value<-fread(file.path("data","climplant_names_trait_V1.2.csv"),encoding = "Latin-1")

flora_survey<-fread(file.path("data","flora_data_metadata","flora_vosges_2022_tax_hom.csv"))
flora_survey<-flora_survey[plot_ID%in%vege_plot_meta_data$plot_ID,]
flora_survey<-merge(flora_survey,indicator_value,all.x=T,by.x="species_name",by.y="lb_nom_final")

flora_survey<-flora_survey[tree!=1,]

##  function to create a [site_id,species] matrix of absence-presence from the flora survey
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

table_sp_vosges<-create_table_sp(flora_survey,"plot_ID")
table_sp_vosges[1:6,1:6]

count_sp<-apply(table_sp_vosges,2,sum)
sort(count_sp)
names(count_sp[count_sp>=3])

table_sp_vosges<-table_sp_vosges[,names(count_sp[count_sp>=2])]
table_sp_vosges<-table_sp_vosges[,names(count_sp[count_sp>=1])]


influencial_species<-c("Galium uliginosum","Gentiana lutea","Lychnis flos-cuculi","Petasites hybridus","Lysimachia nemorum")
influencial_plots<-c("521341_15","3F4MO7_19","4052605_21")

table_sp_vosges<-table_sp_vosges[!rownames(table_sp_vosges)%in%influencial_plots,!colnames(table_sp_vosges)%in%influencial_species]

RS<-apply(data.table(table_sp_vosges)[,-("plot_ID")],1,sum)
RS<-data.table(rs=RS,plot_ID=table_sp_vosges$plot_ID)

##create plot scale averages

cti_indicator<-flora_survey[ ,
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
                                  mean_area=mean(Area,na.rm=T),
                                  freq_for_mean=mean(indFor_Freq,na.rm=T),
                                  indfor_Chytry=mean(indFor_Chytry,na.rm=T),
                                  mean_azote=mean(azote,na.rm=T),
                                  mean_N=mean(Nopt,na.rm=T),
                                  mean_R=mean(R_ellenberg,na.rm=T),
                                  mean_pH=mean(pHopt,na.rm=T),
                                  mean_CN=mean(vi_CN,na.rm=T),
                                  mean_L=mean(Li,na.rm=T)),
                                by=plot_ID]

vege_plot_meta_data<-merge(vege_plot_meta_data,cti_indicator,by="plot_ID")

#### flora exploratory analysis ####
library(vegan)
data_table<-plot_comp[plot_ID%in%rownames(table_sp_vosges),]
table_sp_vosges_2<-table_sp_vosges[rownames(table_sp_vosges)%in%data_table$plot_ID,]
var_part<-varpart(table_sp_vosges_2,~fitted_microclim_scale,~resid_microclim_scale,~mean_pH,data=data_table)
var_part<-varpart(table_sp_vosges_2,~fitted_microclim_max_scale,~resid_microclim_max_scale,~mean_pH,data=data_table)

rda_vosges<-rda(table_sp_vosges_2~ resid_microclim+fitted_microclim+mean_pH,data=data_table)
rda_vosges<-cca(table_sp_vosges_2~ resid_microclim,data=data_table)
rda_vosges<-cca(table_sp_vosges_2~ Heat_load_index_25+ ipv_1000_25+mnt_25_vosges+tree_cover_density_thur_2018_copernicus+mean_pH,data=data_table)


anova.cca(rda_vosges)
summary(rda_vosges)

cca_vosges<-cca(table_sp_vosges_2~ (fitted_microclim+resid_microclim)+mean_pH ,data=data_table)

constraint_vosges<-capscale(table_sp_vosges_2~ (fitted_microclim+resid_microclim)+mean_pH ,data=data_table)

rda_vosges<-rda(table_sp_vosges_2~ fitted_microclim ,data=data_table)
rda_vosges<-rda(table_sp_vosges_2~ resid_microclim ,data=data_table)
rda_vosges<-rda(table_sp_vosges_2~ mean_pH ,data=data_table)
rda_vosges<-rda(table_sp_vosges_2~ period ,data=data_table)

rda_vosges$CCA$u
rda_vosges$CCA$v

summary(var_part)
anova(var_part)
plot(rda_vosges,choices = c(1, 2),scaling=3)
plot(rda_vosges,choices = c(1, 2),scaling=3)
plot(constraint_vosges)

rda_vosges$CCA$wa

rda_plot_dt<-data.table(scores(rda_vosges,choices = c(1,2,3))$sites)
rda_plot_dt[,plot_ID:=rownames(scores(rda_vosges)$sites)]

rda_sp_dt<-data.table(scores(rda_vosges)$species)
rda_sp_dt[,species_name:=rownames(scores(rda_vosges)$species)]

rda_plot_dt<-merge(rda_plot_dt,vege_plot_meta_data,by.x="plot_ID",by.y="plot_ID")



plot(rda_plot_dt[,1],rda_plot_dt[,2])
plot(rda_sp_dt[,1],rda_sp_dt[,2])

plot(rda_vosges$CCA$wa[,1],rda_vosges$CCA$wa[,2])
plot(rda_vosges$CA$u[,1],rda_vosges$CA$u[,2])

library(ade4)

ca<-dudi.coa(table_sp_vosges,nf = 4, scannf = FALSE) 
inertia.dudi(ca)

plot_comp<-data.table(cbind(ca$li,ca$l1))
plot_comp[,plot_ID:=rownames(ca$li)]
plot_comp<-merge(plot_comp,vege_plot_meta_data,by.x="plot_ID",by.y="plot_ID")

plot_comp<-merge(plot_comp,RS)

plot_comp[order(-RS2),][1:30,]
plot_comp[order(RS2),][1:30,]

sp_comp<-data.table(ca$co)
sp_comp[,species_name:=rownames(ca$co)]
sp_comp[,species_name:=str_replace_all(species_name,"[.]"," ")]
sp_comp$count<-colSums(table_sp_vosges[,sp_comp$species_name])

sp_comp[,abbrev_name:=paste0(str_extract_all(species_name,"[:upper:][:alpha:]{3}"),str_extract_all(species_name," [:alpha:]{4}"))]


sp_comp[order(-Comp1),][1:30,]



sp_comp<-merge(sp_comp,indicator_value,by.x="species_name",by.y="lb_nom_final")
rda_sp_dt<-merge(rda_sp_dt,indicator_value,by.x="species_name",by.y="lb_nom_final")
rda_sp_dt[,abbrev_name:=paste0(str_extract_all(species_name,"[:upper:][:alpha:]{3}"),str_extract_all(species_name," [:alpha:]{4}"))]
rda_sp_dt$count<-colSums(table_sp_vosges_2[,rda_sp_dt$species_name])

ggplot(sp_comp[,],aes(x=Comp1,y=Comp2,color=YearMeanMean,size=vi_pH))+geom_point(alpha=0.5)

ggplot(sp_comp[,],aes(x=Comp1,y=Comp2,color=YearMeanMean,size=vi_pH))+geom_point(alpha=0.75)

ggplot(sp_comp[,],aes(x=Comp1,y=Comp2))+geom_point(alpha=0.5)

ggplot(rda_sp_dt[,],aes(x=RDA1,y=PC1,label=species_name))+geom_point(alpha=0.5)+geom_text(data=rda_sp_dt[,])
ggplot(rda_sp_dt[,],aes(x=CCA1,y=CA1,label=species_name))+geom_point(alpha=0.5)+geom_text(data=rda_sp_dt[,])

ggplot(rda_sp_dt[,],aes(x=RDA1,y= RDA2,label=species_name))+geom_point(alpha=0.5)+geom_text(data=rda_sp_dt[,])


ggplot(sp_comp[,],aes(x=Comp1,y=Comp2,color=vi_pH,size=count))+geom_point()
plotly::ggplotly()

ca_plot<-ggplot(sp_comp[,],aes(x=Comp1,y=Comp2,color=vi_pH,label=abbrev_name))+
  geom_point(size=0.25)+
  geom_text(data=sp_comp[count>0,])+
  theme_bw()


cca_plot<-ggplot(rda_sp_dt[,],aes(x=CCA1,y=CCA2,color=vi_pH,label=abbrev_name))+
  geom_point(size=0.25)+
  geom_text(data=rda_sp_dt[count>10,])+
  theme_bw()


plotly::ggplotly(ca_plot)
plotly::ggplotly(cca_plot)


ggplot(sp_comp[,],aes(x=Comp2,y=Comp3,color=YearMeanMean))+geom_point()

ggplot(sp_comp[,],aes(x=Comp1,y=vi_pH))+geom_point()+geom_smooth(method="lm")
ggplot(sp_comp[,],aes(x=YearMean95,y=Comp2))+geom_point()+geom_smooth(method="lm")
ggplot(sp_comp[Comp3>( -5),],aes(x=YearMean95,y=Comp3))+geom_point()+geom_smooth(method="lm")
ggplot(sp_comp[,],aes(x=YearMeanMean,y=Comp4))+geom_point()+geom_smooth(method="lm")





ggplot(rda_plot_dt[,],aes(x= CCA1,y= CCA2,color=predict_microclim_mean,size=mean_pH))+geom_point(alpha=0.5)+theme_bw()


ggplot(rda_plot_dt[predict_microclim_mean<21,],aes(x=predict_microclim_mean,y= CCA2))+geom_point()+geom_smooth()
ggplot(rda_plot_dt[predict_microclim_mean<21,],aes(x=mean_pH,y= CCA1))+geom_point()+geom_smooth()


ggplot(plot_comp[,],aes(x=Axis1,y=Axis2,color=predict_microclim_mean,size=mean_pH))+geom_point(alpha=0.5)

ggplot(plot_comp[,],aes(x=Axis1,y=Axis2,color=bool_microclim))+geom_point(alpha=0.75)


ggplot(plot_comp,aes(x=Axis2,fill=bool_microclim))+geom_density(alpha=0.5)+scale_fill_manual(values=c("lightblue","tomato"))+theme_bw()


ggplot(plot_comp[predict_microclim_mean<21,],aes(x=predict_microclim_mean,y=Axis2))+geom_point()+geom_smooth()

ggplot(plot_comp[predict_microclim_mean<21,],aes(x=predict_microclim_mean,y=cit_climplant))+geom_point()+geom_smooth()

ggplot(plot_comp[,],aes(x=mnt_25_vosges,y=RS2))+geom_point()+geom_smooth(method="lm")

ggplot(plot_comp[,],aes(x=mean_pH,y=RS1))+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp[RS3 > (-20),],aes(x=cit_climplant,y=RS3))+geom_point()+geom_smooth(method="lm")

plot_comp[,cv_cit_climplant:=sd_cit_climplant/cit_climplant]
plot_comp[,cv_cit_ecoplant_picq:=sd_cit_ecoplant_picq/cit_ecoplant_picq]

plot_comp[,.(.N,mean(n_sp_climplant)),by=quali_plot]
plot_comp$quali_plot<-plot_comp$`Qualité des reelvés: 1 bon, 2 A voir`
plot_comp[is.na(quali_plot),"quali_plot"]<-2
plot_comp[,quali_plot:=as.character(quali_plot)]

ggplot(plot_comp,aes(x=quali_plot,y=resid_microclim))+theme_bw()+geom_boxplot()

ggplot(plot_comp,aes(x=tree_density_2018,y=n_sp_climplant))+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp[ rs < 40 & ANNEE!=2014 ],aes(x=resid_microclim,y=rs))+geom_point()+geom_smooth(method="gam")
ggplot(plot_comp[ abs(resid_rs)<18],aes(x=resid_microclim,y=resid_rs+2))+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp[cit_climplant>7.2 & rs < 50,],aes(x=resid_microclim,y=rs))+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp[cit_climplant>7.2 & rs < 50,],aes(x=resid_microclim,y=rs,size=cit_ecoplant_picq))+geom_point()+geom_smooth(method="lm")


ggplot(plot_comp[,],aes(x=Heat_load_index_25,y=cit_ecoplant_picq))+geom_point()+geom_smooth(method="lm")



plot_comp[,bool_microclim:=ifelse(resid_microclim<median(resid_microclim),"cold_microclim","warm_microclim")]
plot_comp[,cut_microclim:=cut_number(resid_microclim,3)]
plot_comp[,cut_ipv:=cut_number(ipv_effect,3)]
plot_comp[,cut_canopy:=cut_number(canopy_effect,3)]
plot_comp[,cut_hl:=cut_number(heat_load_effect,3)]


plot_comp[,.(mean(cit_ecoplant_picq),mean(cit_climplant,na.rm=T),mean(mean_pH),median(resid_rs),mean(rs)),by=cut_microclim]

flora_survey_sub[,.(.N,mean(topt_picq,na.rm=T),mean(YearMeanMean ,na.rm=T),mean(vi_pH,na.rm=T)),by=bool_microclim]
flora_survey_sub[,.(quantile(topt_picq,na.rm=T,probs=c(0.90)),quantile(YearMeanMean ,na.rm=T,probs=c(0.90)),quantile(vi_pH,na.rm=T,probs=c(0.90))),by=bool_microclim]
flora_survey_sub[,.(quantile(topt_picq,na.rm=T,probs=c(0.1)),quantile(YearMeanMean ,na.rm=T,probs=c(0.1)),quantile(vi_pH,na.rm=T,probs=c(0.1))),by=bool_microclim]


flora_survey_sub<-merge(flora_survey,plot_comp[,c("plot_ID","bool_microclim","cut_microclim")],by="plot_ID")


ggplot(flora_survey_sub,aes(x=topt_picq ,fill=as.factor(bool_microclim)))+theme_bw()+geom_density(alpha=0.5)+scale_fill_manual(values=c("lightblue","tomato"))
ggplot(flora_survey_sub,aes(x=YearMeanMean ,fill=as.factor(bool_microclim)))+theme_bw()+geom_histogram(binwidth =0.5,alpha=0.5,color="grey40",position ="identity")+scale_fill_manual(values=c("lightblue","tomato"))

ggplot(flora_survey_sub,aes(x=YearMeanMean ,fill=cut_microclim))+theme_bw()+geom_density(alpha=0.5)+scale_fill_manual(values=c("#440154FF","#22A884FF","#FDE725FF"))
ggplot(flora_survey_sub,aes(x=YearMeanMean ,fill=cut_microclim))+theme_bw()+geom_histogram(binwidth =0.5,alpha=0.5,color="grey40",position ="identity")+scale_fill_manual(values=c("#440154FF","#22A884FF","#FDE725FF"))


plot_comp$resid_rs<-residuals(rs_lm_mod)

ggplot(plot_comp,aes(x=bool_microclim,y=resid_rs,fill=bool_microclim))+theme_bw()+geom_boxplot(alpha=0.5)+scale_fill_manual(values=c("#440154FF","#22A884FF","#FDE725FF"))
ggplot(plot_comp,aes(x=cut_microclim,y=resid_rs+2,fill=cut_microclim))+theme_bw()+geom_boxplot(alpha=0.5)+scale_fill_manual(values=c("#440154FF","#22A884FF","#FDE725FF"))
levels(plot_comp$cut_microclim)
my_comparisons <- list( levels(plot_comp$cut_microclim)[1:2],  levels(plot_comp$cut_microclim)[2:3],  levels(plot_comp$cut_microclim)[c(1,3)] )
ggplot(plot_comp,aes(x=cut_microclim,y=mnt_25_vosges,fill=cut_microclim))+theme_bw()+geom_boxplot(alpha=0.5)+scale_fill_manual(values=c("#440154FF","#22A884FF","#FDE725FF"))+ stat_compare_means(comparisons = my_comparisons)

ggplot(plot_comp,aes(x=cut_microclim,y=cit_climplant,fill=cut_microclim))+theme_bw()+geom_violin(draw_quantiles = c(0.25,0.5,0.75),alpha=0.5)+scale_fill_manual(values=c("lightblue","grey80","tomato"))
ggplot(plot_comp,aes(x=cut_microclim,y=cit_ecoplant_picq,fill=cut_microclim))+theme_bw()+geom_violin(draw_quantiles = c(0.25,0.5,0.75),alpha=0.5)+scale_fill_manual(values=c("lightblue","grey80","tomato"))
ggplot(plot_comp,aes(x=cut_microclim,y=Axis2,fill=cut_microclim))+theme_bw()+geom_violin(draw_quantiles = c(0.25,0.5,0.75),alpha=0.5)+scale_fill_manual(values=c("lightblue","grey80","tomato"))

ggplot(plot_comp,aes(x=bool_microclim,y=resid_rs,fill=bool_microclim))+theme_bw()+geom_violin(alpha=0.5)+scale_fill_manual(values=c("lightblue","grey80","tomato"))
plot_comp$cit_climplant

ggplot(flora_survey_sub,aes(x=cut_microclim,y=topt_picq,fill=cut_microclim))+theme_bw()+geom_boxplot(draw_quantiles = c(0.25,0.5,0.75),alpha=0.5)+scale_fill_manual(values=c("lightblue","grey80","tomato"))

ggplot(flora_survey_sub,aes(x=Li ,fill=as.factor(bool_microclim)))+theme_bw()+geom_histogram(binwidth =0.5,alpha=0.5,color="grey40",position ="identity")+scale_fill_manual(values=c("lightblue","tomato"))

summary(lm(rs~fitted_microclim+resid_microclim+quali_plot+mean_pH,data=plot_comp[cit_climplant>7.2 & rs < 50,]))


summary(lm(rs~predict_microclim_mean+quali_plot+mean_pH,data=plot_comp[cit_climplant>7.2 & rs < 50,]))
summary(lm(n_sp_climplant~resid_microclim+mean_pH,data=plot_comp))
summary(lm(n_sp_climplant~fitted_microclim+tree_density_2018+Heat_load_index_25+ipv_1000_25+quali_plot+mean_pH,data=plot_comp))
summary(lm(Axis2~fitted_microclim+resid_microclim+quali_plot+mean_pH,data=plot_comp))

rs_lm_mod<-lm(rs~fitted_microclim+resid_microclim+quali_plot+mean_pH,data=plot_comp[ rs < 50,])
rs_lm_mod<-lm(rs~fitted_microclim+quali_plot+mean_pH,data=plot_comp)
rs_lm_mod<-lm(cit_ecoplant_picq~fitted_microclim+mean_pH,data=plot_comp[ ,])

summary(rs_lm_mod)

ggplot(plot_comp,aes(x=resid_microclim,y=residuals(rs_lm_mod)))+geom_point()+geom_smooth(method="lm")


summary(lm(residuals(rs_lm_mod)~resid_microclim,data=plot_comp))


summary(lm(n_sp_climplant~pred_extrem+micro_extrem+mean_pH,data=plot_comp))

ggplot(plot_comp[,],aes(x=predict_microclim_mean,y=mnt_25_vosges,color=Heat_load_index_25))+geom_point()+geom_smooth(method="lm")

summary(lm(RS1~mean_pH,data=plot_comp))

summary(lm(Axis2~mnt_25_vosges+Heat_load_index_25+ipv_1000_25,data=plot_comp))

summary(lm(Axis2~cit_climplant,data=plot_comp))
summary(lm(Axis2~poly(cit_climplant,2),data=plot_comp[!is.na(cit_climplant),]))
summary(lm(Axis2~poly(predict_microclim_mean,1)+predict_microclim_min,data=plot_comp))
summary(lm(Axis2~poly(predict_microclim_max,1),data=plot_comp))
summary(lm(Axis2~poly(predict_microclim_range,1),data=plot_comp))



summary(lm(Axis2~fitted_microclim+resid_microclim,data=plot_comp[predict_microclim_mean<21,]))

summary(lm(CCA2~fitted_microclim+resid_microclim,data=rda_plot_dt[predict_microclim_mean<21,]))

summary(lm(Axis2~poly(mnt_25_vosges,2)+poly(Heat_load_index_25,2),data=plot_comp))

summary(lm(cit_climplant~poly(predict_microclim_mean,2),data=plot_comp[predict_microclim_mean<21,]))
summary(lm(cit_climplant~predict_microclim_mean,data=plot_comp[predict_microclim_mean<21,]))
summary(lm(cit_climplant~mnt_25_vosges,data=plot_comp[predict_microclim_mean<21,]))


summary(lm(RS4~predict_microclim_snow,data=plot_comp))


summary(lm(Axis2~tmoy_digitalis_8610_25_vosges+pp_digitalis_8610_25_vosges+aspect_NS_25_vosges,data=plot_comp))

scatter(ca)
plot(ca)

library(FactoMineR)
ca<-CA(table_sp_vosges)

plot(ca)


#### plot_comp maps ####
plot_comp_sf<-st_as_sf(plot_comp,coords=c("X_L93","Y_L93"),crs=st_crs(2154))


ggplot(plot_comp_sf,aes(color=period))+theme_bw()+geom_sf()
ggplot(plot_comp_sf,aes(color=RS1))+theme_bw()+geom_sf()+scale_color_gradient2()

ggplot(plot_comp_sf,aes(color=ipv_1000_25))+theme_bw()+geom_sf()+scale_color_gradient2()

ggplot(plot_comp_sf,aes(color=resid_rs))+theme_bw()+geom_sf()+scale_color_gradient2()

ggplot(plot_comp_sf,aes(color=Heat_load_index_25))+theme_bw()+geom_sf()+scale_color_gradient2()

ggplot(plot_comp_sf,aes(color=resid_rs))+theme_bw()+geom_sf()+scale_color_gradient2(limits=c(-20,20),oob=scales::squish)

ggplot(plot_comp_sf,aes(color=cut_microclim))+theme_bw()+geom_sf()


ggplot(plot_comp_sf,aes(color=RS3))+theme_bw()+geom_sf()+scale_color_gradient2(limits=c(-0.5,0.5),oob=scales::squish)

ggplot(plot_comp_sf,aes(color=resid_microclim_scale))+theme_bw()+geom_sf()+scale_color_gradient2(limits=c(-2,2),oob=scales::squish,low="cadetblue",high="tomato")


#### temporal analysis ? ####

cumsum( table(plot_comp$ANNEE))
plot_comp[,period:=ifelse(ANNEE<= 2015,"Past","Recent")]

plot_comp<-plot_comp[predict_microclim_mean<21,]


ggplot(plot_comp,aes(x=as.character(ANNEE),y=predict_microclim_mean))+theme_bw()+geom_boxplot()
ggplot(plot_comp,aes(x=as.character(ANNEE),y=mean_pH))+theme_bw()+geom_boxplot()
ggplot(plot_comp[n_sp_climplant>=5,],aes(x=as.character(period),y=predict_microclim_mean,fill=period))+theme_bw()+geom_violin(alpha=0.5)
ggplot(plot_comp[n_sp_climplant>=5,],aes(x=as.character(period),y=mean_pH,fill=period))+theme_bw()+geom_violin(alpha=0.5)

ggplot(plot_comp[n_sp_climplant>=5,],aes(x=as.character(period),y=mnt_25_vosges,fill=period))+theme_bw()+geom_violin(alpha=0.5)


ggplot(plot_comp,aes(x=as.character(period),y=Heat_load_index_25))+theme_bw()+geom_boxplot()+geom_violin()

ggplot(plot_comp,aes(x=mean_pH,y=cit_climplant,color=period))+theme_bw()+geom_point()+geom_smooth(method="lm")

library(corrplot)
nnnn<-grep("predict",colnames(plot_comp),value=T)
corrplot(cor(plot_comp[,..nnnn]))


ggplot(plot_comp,aes(x=predict_microclim_mean,y=RS2))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp,aes(x=mnt_25_vosges,y=cit_climplant,color=period))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp,aes(x=mnt_25_vosges,y=predict_microclim_mean))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp,aes(x=resid_microclim,y=RS2))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp[n_sp_climplant>=5,],aes(x=predict_microclim_max,y=cit_climplant,color=period))+theme_bw()+geom_point()+geom_smooth(method="lm")

ggplot(plot_comp[,],aes(x=mnt_25_vosges,y=predict_microclim_max))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp[n_sp_climplant>=5,],aes(x=resid_microclim,y=cit_climplant,color=period))+theme_bw()+geom_point()+geom_smooth(method="lm")



ggplot(plot_comp,aes(x=fitted_microclim,y=Axis2))+theme_bw()+geom_point()+geom_smooth(method="lm")
ggplot(plot_comp,aes(x=resid_microclim,y=Axis2))+theme_bw()+geom_point()+geom_smooth(method="lm")


lm_resid<-lm(predict_microclim_mean~mnt_25_vosges,data=plot_comp)

summary(lm_resid)
summary(daily_model)$coefficients[,1]$mnt_25_vosges

plot_comp$fitted_microclim<-plot_comp$mnt_25_vosges* summary(daily_model)$coefficients[,1]["mnt_25_vosges"] + summary(daily_model)$coefficients[,1]["(Intercept)"] +15*summary(daily_model)$coefficients[,1]["temp_ventr_station"]
plot_comp$resid_microclim<-  plot_comp$predict_microclim_mean - plot_comp$fitted_microclim

plot_comp$fitted_microclim_max<-plot_comp$mnt_25_vosges* summary(daily_model_max)$coefficients[,1]["mnt_25_vosges"]
plot_comp$resid_microclim_max<-  plot_comp$predict_microclim_max - plot_comp$fitted_microclim


hist(plot_comp$fitted_microclim)
hist(plot_comp$resid_microclim)

plot_comp<-plot_comp[resid_microclim<0,]


plot_comp$resid_microclim_scale<-scale(plot_comp$resid_microclim)
plot_comp$fitted_microclim_scale<-scale(plot_comp$fitted_microclim)

plot_comp$resid_microclim_max_scale<-scale(plot_comp$fitted_microclim_max)
plot_comp$fitted_microclim_max_scale<-scale(plot_comp$resid_microclim_max)

plot_comp$mean_pH_scale<-scale(plot_comp$mean_pH)

plot_comp[,tree_density_2018:=tree_cover_density_thur_2018_copernicus]
plot_comp[,pred_extrem:=predict(lm_agg_mean_ext,newdata=plot_comp)]
plot_comp[,elev_extrem:=mnt_25_vosges*coef(lm_agg_mean_ext)["mnt_25_vosges"] + coef(lm_agg_mean_ext)[1] ]
plot_comp[,micro_extrem:=pred_extrem-elev_extrem]



rda_plot_dt$resid_microclim_scale<-scale(rda_plot_dt$resid_microclim)
rda_plot_dt$fitted_microclim_scale<-scale(rda_plot_dt$fitted_microclim)
rda_plot_dt$mean_pH_scale<-scale(rda_plot_dt$mean_pH)






plot_comp[,resid_microclim_bool:=as.character(resid_microclim>0)]
lm_temp<-(lm(cit_climplant~mnt_25_vosges*period+mean_pH,data=plot_comp[n_sp_climplant>=5,]))

lm_temp<-(lm(cit_climplant~predict_microclim_mean*period+mean_pH,data=plot_comp[n_sp_climplant>=5,]))

lm_temp<-(lm(cit_climplant~fitted_microclim*period+mean_pH+resid_microclim*period,data=plot_comp[n_sp_climplant>=5,]))

lm_temp<-(lm(cit_climplant~fitted_microclim*period*resid_microclim+mean_pH,data=plot_comp[n_sp_climplant>=5 & resid_microclim<2 ,]))

lm_temp<-(lm(cit_climplant~fitted_microclim*period+mean_pH+quant_residual*period,data=plot_comp[n_sp_climplant>=5,]))

lm_temp<-(lm(cit_climplant~fitted_microclim*period+quant_residual*period+mean_pH,data=plot_comp[n_sp_climplant>=5,]))

lm_temp<-(lm(cit_climplant~mnt_25_vosges*period*resid_microclim_bool+mean_pH,data=plot_comp[n_sp_climplant>=5,]))


lm_temp<-(lm(cit_climplant~mnt_25_vosges*period + Heat_load_index_25*period+mean_pH,data=plot_comp[n_sp_climplant>=5,]))

summary(lm_temp)
car::vif(lm_temp,type = 'predictor')

summary(lm(RS2~fitted_microclim*resid_microclim,data=plot_comp))
summary(lm(cit_climplant~fitted_microclim*resid_microclim+mean_pH,data=plot_comp[n_sp_climplant>=5,]))
summary(lm(cit_climplant~fitted_microclim+resid_microclim+mean_pH,data=plot_comp[n_sp_climplant>=5,]))

summary(lm(RS2~fitted_microclim,data=plot_comp))
summary(lm(RS2~resid_microclim,data=plot_comp))

lm_temp<-(lm(Axis2~mean_pH_scale+(fitted_microclim_scale*resid_microclim_scale),data=plot_comp[`precision localisation`!=5,]))

lm_temp<-(lm(cit_climplant~mean_pH_scale+(fitted_microclim_scale*resid_microclim_scale),data=plot_comp))
lm_temp<-(lm(cit_ecoplant_picq~mean_pH_scale+(fitted_microclim_scale*resid_microclim_scale),data=plot_comp))

lm_temp<-(lm(Axis2~fitted_microclim_scale*resid_microclim_scale,data=plot_comp))
lm_temp<-(lm(Axis2~fitted_microclim_max_scale*resid_microclim_max_scale,data=plot_comp)) # redondant info
lm_temp<-(lm(Axis2~predict_microclim_mean,data=plot_comp))
lm_temp<-(lm(Axis2~predict_microclim_max,data=plot_comp))

lm_temp<-(lm(cit_climplant~fitted_microclim_scale*resid_microclim_scale+mean_pH,data=plot_comp))
lm_temp<-(lm(cit_climplant~predict_microclim_mean+mean_pH,data=plot_comp))
lm_temp<-(lm(cit_ecoplant_picq~fitted_microclim_scale*resid_microclim_scale+mean_pH,data=plot_comp))
lm_temp<-(lm(cit_ecoplant_picq~predict_microclim_mean+mean_pH,data=plot_comp))

summary(lm(RS2~period*(fitted_microclim+resid_microclim),data=plot_comp))

summary(lm(RS2~predict_microclim_mean,data=plot_comp))


summary(lm_with_resid)

summary(lm_temp)

#### Figure effect size and variation partitioning ####

library(modEvA)
library(ggvenn)
library(ggeffects)
create_eff_size_table<-function(model,x1="fitted_microclim_scale",x2="resid_microclim_scale"){
  tmp1<-ggpredict(model,terms=x1)
  tmp1$factor<-"Elevation"
  tmp2<-ggpredict(model,terms=x2)
  tmp2$factor<-"Microclimate"
  tmp2$predicted<-tmp2$predicted-0.15
  tmp2$conf.low<-tmp2$conf.low-0.15
  tmp2$conf.high<-tmp2$conf.high-0.15
  tmp3<-ggpredict(model,terms="mean_pH_scale")
  tmp3$factor<-"Soil pH"
  tmp3$predicted<-tmp3$predicted-0.3
  tmp3$conf.low<-tmp3$conf.low-0.3
  tmp3$conf.high<-tmp3$conf.high-0.3
  effect_sizes<-rbind(tmp1,tmp2,tmp3)
  effect_sizes
}

create_variation_part<-function(model,groups=c("fitted_microclim_scale","resid_microclim_scale","mean_pH_scale")){
  
  grp<-data.frame(var=groups,group=c("Elevation","Microclimate","soil"))
  
  resvarpart<-varPart(model=model,groups=grp)
  print(resvarpart)
  resvarpart<-round(resvarpart$Proportion*100,1)
  resvarpart<-paste0(resvarpart,"%")
  
  list_venn<-list(Elevation=c((resvarpart[1]),(resvarpart[4]),(resvarpart[6]),resvarpart[7]),
                  Microclimate=c((resvarpart[2]),(resvarpart[4]),(resvarpart[5]),resvarpart[7]),
                  `Soil pH`=c((resvarpart[3]),(resvarpart[5]),(resvarpart[6]),resvarpart[7]))
  
  list(resvarpart,list_venn)
  
}

plot_comp_climplant<-plot_comp[n_sp_climplant>=5,]
lm_ca<-lm(Axis2~mean_pH_scale+fitted_microclim_scale+resid_microclim_scale,data=plot_comp)
lm_ca<-lm(Axis2~mean_pH_scale+mean_L+fitted_microclim_scale+resid_microclim_scale,data=plot_comp)

lm_ca<-glm(Axis2~mean_pH_scale+fitted_microclim_scale+resid_microclim_scale,data=plot_comp)
car::vif(lm_ca)
lm_climplant<-glm(cit_climplant~mean_pH_scale+fitted_microclim_scale+resid_microclim_scale,data=plot_comp_climplant)

lm_cca<-glm(CCA2~fitted_microclim_scale+resid_microclim_scale+mean_pH_scale,data=rda_plot_dt[predict_microclim_mean<21,])

lm_extreme<-glm(Axis2~elev_extrem + micro_extrem+mean_pH_scale,data=plot_comp[ ,])
lm_extreme_2<-glm(cit_climplant~elev_extrem + micro_extrem+mean_pH_scale,data=plot_comp_climplant[ ,])

lm_rs<-glm(rs~fitted_microclim_scale + resid_microclim_scale+mean_pH_scale,data=plot_comp)

color_three_fact<-c("#0080C2FF", "#EFC000FF", "#CD534CFF")

summary(lm_ca)

list_venn<-create_variation_part(lm_ca)
list_venn_climplant<-create_variation_part(lm_climplant)
list_venn_cca2<-create_variation_part(lm_cca)
list_venn_extreme<-create_variation_part(lm_extreme,c("elev_extrem","micro_extrem","mean_pH_scale"))
list_venn_extreme<-create_variation_part(lm_extreme_2,c("elev_extrem","micro_extrem","mean_pH_scale"))

list_venn_rs<-create_variation_part(lm_rs)


plot_varpart_ca<-ggvenn(list_venn[[2]],show_elements=T,fill_color =color_three_fact,text_size = 3.7,fill_alpha = 0.3,stroke_size = 0.5,set_name_size =4.5)+
  annotate("text",x = -0.8, y = -2.2, label = paste0("Unexplained = ",list_venn[[1]][8]),size=3.5)+theme(plot.margin = margin(t=2,b=4,unit = "mm"),plot.background = element_rect(fill=NA,color=NA))

plot_varpart_climplant<-ggvenn(list_venn_climplant[[2]],show_elements=T,fill_color =color_three_fact,text_size = 3.7,fill_alpha = 0.3,stroke_size = 0.5,set_name_size =4.5)+
  annotate("text",x = -0.8, y = -2.2, label = paste0("Unexplained = ",list_venn_climplant[[1]][8]),size=3.5)+theme(plot.margin = margin(t=2,b=4,unit = "mm"),plot.background = element_rect(fill=NA,color=NA))

plot_varpart_cca2<-ggvenn(list_venn_extreme[[2]],show_elements=T,fill_color =color_three_fact,text_size = 3.7,fill_alpha = 0.3,stroke_size = 0.5,set_name_size =4.5)+
  annotate("text",x = -0.8, y = -2.2, label = paste0("Unexplained = ",list_venn_extreme[[1]][8]),size=3.5)+theme(plot.margin = margin(t=2,b=4,unit = "mm"),plot.background = element_rect(fill=NA,color=NA))


plot(stack_all_variable_vallee$Heat_load_index_25)


effect_sizes<-as.data.table(create_eff_size_table(lm_ca))
effect_sizes_climplant<-as.data.table(create_eff_size_table(lm_climplant))
effect_sizes_ext<-as.data.table(create_eff_size_table(lm_extreme,"elev_extrem","micro_extrem"))


table_effect_size_plot<-data.table(color=color_three_fact,
                                   name=c("Elevation","Microclimate","Soil pH"),
                                   mean=c(mean(plot_comp$fitted_microclim),mean(plot_comp$resid_microclim),mean(plot_comp$mean_pH)),
                                   sd=c(sd(plot_comp$fitted_microclim),sd(plot_comp$resid_microclim),sd(plot_comp$mean_pH)))

table_effect_size_plot_climplant<-data.table(color=color_three_fact,
                                   name=c("Elevation","Microclimate","Soil pH"),
                                   mean=c(mean(plot_comp_climplant$fitted_microclim),mean(plot_comp_climplant$resid_microclim),mean(plot_comp_climplant$mean_pH)),
                                   sd=c(sd(plot_comp_climplant$fitted_microclim),sd(plot_comp_climplant$resid_microclim),sd(plot_comp_climplant$mean_pH)))

table_effect_size_plot_ext<-data.table(color=color_three_fact,
                                   name=c("Elevation","Microclimate","Soil pH"),
                                   mean=c(mean(plot_comp$elev_extrem),mean(plot_comp$micro_extrem),mean(plot_comp$mean_pH)),
                                   sd=c(sd(plot_comp$elev_extrem),sd(plot_comp$micro_extrem),sd(plot_comp$mean_pH)))


table_effect_size_plot[,label_ca:=paste0(name," µ = " ,round(mean,1),", σ = ",round(sd,2))]
table_effect_size_plot_climplant[,label_ca:=paste0(name," µ = " ,round(mean,1),", σ = ",round(sd,2))]
table_effect_size_plot_ext[,label_ca:=paste0(name," µ = " ,round(mean,1),", σ = ",round(sd,2))]


plot_ca_effect_size<-ggplot(effect_sizes[abs(x) <= 2 ,],aes(x=x, y=predicted,color=factor,fill=factor))+
  geom_line(linewidth=0.9)+
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2,color="grey40",linewidth=0.2,lty=2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3,color="grey40",linewidth=0.2)+
  theme_bw()+
  scale_color_manual(values=color_three_fact,label=table_effect_size_plot$label_ca)+
  scale_fill_manual(values=color_three_fact,label=table_effect_size_plot$label_ca)+
  labs(x="Scaled predictor",y="Flora composition (CA Axis 2)",color="Predictor",fill="Predictor")+ 
  theme(legend.position = c(0.35, 0.2),legend.text=element_text(size=9),legend.title = element_text(size=10),legend.background = element_rect(fill=NA))

plot_ca_effect_size_climplant<-ggplot(effect_sizes_climplant[abs(x) <= 2 ,],aes(x=x, y=predicted,color=factor,fill=factor))+
  geom_line(linewidth=0.9)+
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2,color="grey40",linewidth=0.2,lty=2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3,color="grey40",linewidth=0.2)+
  theme_bw()+
  scale_color_manual(values=color_three_fact,label=table_effect_size_plot_climplant$label_ca)+
  scale_fill_manual(values=color_three_fact,label=table_effect_size_plot_climplant$label_ca)+
  labs(x="Scaled predictor",y="Mean thermal optimum °C",color="Predictor",fill="Predictor")+ 
  theme(legend.position = c(0.72, 0.2),legend.text=element_text(size=9),legend.title = element_text(size=10),legend.background = element_rect(fill=NA))


plot_ca_effect_size_ext<-ggplot(effect_sizes_climplant[abs(x) <= 2 ,],aes(x=x, y=predicted,color=factor,fill=factor))+
  geom_line(linewidth=0.9)+
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2,color="grey40",linewidth=0.2,lty=2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3,color="grey40",linewidth=0.2)+
  theme_bw()+
  scale_color_manual(values=color_three_fact,label=table_effect_size_plot_ext$label_ca)+
  scale_fill_manual(values=color_three_fact,label=table_effect_size_plot_ext$label_ca)+
  labs(x="Scaled predictor",y="Mean thermal optimum °C",color="Predictor",fill="Predictor")+ 
  theme(legend.position = c(0.72, 0.2),legend.text=element_text(size=9),legend.title = element_text(size=10),legend.background = element_rect(fill=NA))




export_figure_eff_size_varpart<-ggarrange(plot_ca_effect_size,plot_varpart_ca,plot_ca_effect_size_climplant,plot_varpart_climplant,
                                          labels=c("a)","","b)",""),ncol=2,nrow=2,widths =c(1.4,1),align = "v")
export_figure_eff_size_varpart<-export_figure_eff_size_varpart+theme(plot.background = element_rect(fill="white",color=NA))

ggsave(plot=export_figure_eff_size_varpart,file.path("figure_result","Figure_eff_varpart.jpg"),width=180,height =160 ,scale=1.1,unit="mm",dpi=400)

get_coef<-summary(daily_model)$coefficients[,1]


summary(daily_model)

plot_comp[,ipv_effect:= ipv_1000_25*get_coef["ipv_1000_25"] + 15 *(get_coef["temp_ventr_station"]+(ipv_1000_25*get_coef["temp_ventr_station:ipv_1000_25"])) ]
plot_comp[,canopy_effect:= tree_cover_density_thur_2018_copernicus*get_coef["tree_density_2018"] + 15 *(get_coef["temp_ventr_station"]+(tree_cover_density_thur_2018_copernicus*get_coef["temp_ventr_station:tree_density_2018"])) ]
plot_comp[,heat_load_effect:= Heat_load_index_25*get_coef["Heat_load_index_25"] + 15 *(get_coef["temp_ventr_station"]+(Heat_load_index_25*get_coef["temp_ventr_station:Heat_load_index_25"])) ]

get_coef<-summary(lm_agg_mean_glm)$coefficients[,1]
plot_comp[,ipv_effect:= ipv_1000_25*get_coef["ipv_1000_25"] ]
plot_comp[,canopy_effect:= tree_cover_density_thur_2018_copernicus*get_coef["tree_density_2018"]  ]
plot_comp[,heat_load_effect:= Heat_load_index_25*get_coef["Heat_load_index_25"] ]



plot_comp[,ipv_effect:= scale(ipv_effect) ]
plot_comp[,canopy_effect:= scale(canopy_effect)]
plot_comp[,heat_load_effect:= scale(heat_load_effect)]


plot_comp[,sum_comp:=ipv_effect+canopy_effect+heat_load_effect]

lm_microclim<-lm(Axis2~resid_microclim_scale,data=plot_comp)
lm_microclim<-glm(Axis2~ipv_effect+canopy_effect+heat_load_effect,data=plot_comp)
lm_microclim<-glm(rs~ipv_effect+canopy_effect+heat_load_effect,data=plot_comp)
lm_microclim<-glm(cit_climplant~ipv_effect+canopy_effect+heat_load_effect,data=plot_comp)
lm_microclim<-glm(cit_ecoplant_picq~ipv_effect+canopy_effect+heat_load_effect,data=plot_comp)


create_variation_part_2<-function(model){
  
  grp<-data.frame(var=c("heat_load_effect","ipv_effect","canopy_effect"),group=c("Heat load","Topographic position","Canopy cover"))
  
  resvarpart<-varPart(model=model,groups=grp)
  print(resvarpart)
  resvarpart<-round(resvarpart$Proportion*100,1)
  resvarpart<-paste0(resvarpart,"%")
  
  list_venn<-list(`Heat load`=c((resvarpart[1]),(resvarpart[4]),(resvarpart[6]),resvarpart[7]),
                  `Topographic position`=c((resvarpart[2]),(resvarpart[4]),(resvarpart[5]),resvarpart[7]),
                  `Canopy cover`=c((resvarpart[3]),(resvarpart[5]),(resvarpart[6]),resvarpart[7]))
  
  list(resvarpart,list_venn)
  
}

create_variation_part_2(lm_microclim)
car::vif(lm_microclim)



lm_climplant<-lm(cit_climplant~mean_pH_scale+fitted_microclim_scale+resid_microclim_scale,data=plot_comp_climplant)
lm_climplant<-glm(cit_climplant_05~mean_pH_scale+fitted_microclim_scale+resid_microclim_scale,data=plot_comp_climplant)
lm_climplant<-glm(cit_climplant_95~mean_pH_scale+fitted_microclim_scale+resid_microclim_scale,data=plot_comp_climplant)

hist(plot_comp_climplant$cit_climplant_05)

list_venn_climplant<-create_variation_part(lm_climplant)


plot_varpart_ca<-ggvenn(list_venn[[2]],show_elements=T,fill_color =color_three_fact,text_size = 3.7,fill_alpha = 0.3,stroke_size = 0.5,set_name_size =4.5)+
  annotate("text",x = -0.8, y = -2.2, label = paste0("Unexplained = ",list_venn[[1]][8]),size=3.5)+theme(plot.margin = margin(t=2,b=4,unit = "mm"),plot.background = element_rect(fill=NA,color=NA))

plot_varpart_climplant<-ggvenn(list_venn_climplant[[2]],show_elements=T,fill_color =color_three_fact,text_size = 3.7,fill_alpha = 0.3,stroke_size = 0.5,set_name_size =4.5)+
  annotate("text",x = -0.8, y = -2.2, label = paste0("Unexplained = ",list_venn_climplant[[1]][8]),size=3.5)+theme(plot.margin = margin(t=2,b=4,unit = "mm"),plot.background = element_rect(fill=NA,color=NA))





summary(lm_climplant)
summary(lm_microclim)

ggplot(plot_comp,aes(x=fitted_microclim_scale,y=cit_climplant))+theme_bw()+geom_point(alpha=0.5)+geom_smooth(method="lm")
ggplot(plot_comp,aes(x=fitted_microclim_scale,y=cit_climplant_95))+theme_bw()+geom_point(alpha=0.5)+geom_smooth(method="lm")

ggplot(plot_comp,aes(x=mnt_25_vosges,y=predict_microclim_mean))+theme_bw()+geom_point(alpha=0.75)+geom_line(mapping = aes(y=fitted_microclim),linewidth=1,color="#0080C2FF")
ggplot(plot_comp,aes(x=resid_microclim_scale,y=cit_climplant_95))+theme_bw()+geom_point(alpha=0.5)+geom_smooth(method="lm")



#### explanatory figure ####

pred_elev_micro_figure_2<-ggplot(plot_comp,aes(x=mnt_25_vosges,y=predict_microclim_mean,fill=cit_climplant))+
  theme_bw()+
  geom_point(alpha=1,size=2,shape=21,color="grey60",fill="grey80")+
  geom_line(mapping = aes(y=fitted_microclim),linewidth=1,color="#0080C2FF")+
  geom_segment(data=plot_comp[plot_ID%in%plot_with_arrows,],aes(xend=mnt_25_vosges,y=predict_microclim_mean+0.05,yend=fitted_microclim-0.05),linewidth=1,color="#EFC000FF",arrow = arrow(length = unit(0.03, "npc"), ends = "both"))+
  #scale_color_viridis_c(direction = -1)
  #scale_color_gradientn(colors = brewer.pal(7, "PiYG")[order(7:1)])
  scale_fill_gradient2(midpoint = 7.5,low="#4D9221",high="#C51B7D",mid="white")+
  labs(x="Elevation (m)",y="Predicted understory temperature °C",fill="Thermal optimum")+
  annotate("label",645,20,label="Lapse rate",fill="#0080C2FF",alpha=0.3)+
  annotate("label",790,17.9,label="Microclimate",fill="#EFC000FF",alpha=0.3)


plotly::ggplotly(pred_elev_micro_figure)


ggsave(file.path("figure_result","Figure_prediction_elevation.jpg"),pred_elev_micro_figure,scale=1.25,width=180,unit="mm",height=100)
ggsave(file.path("figure_result","Figure_prediction_elevation_2.jpg"),pred_elev_micro_figure_2,scale=1.25,width=180,unit="mm",height=100)

plot_with_arrows<-c("620911_16","4112501_21","30901_09","510981_15","4113402_21")

plot_model(lm_temp,type="pred",terms=c("predict_microclim_mean","period"))
plot_model(lm_temp,type="pred",terms=c("predict_microclim_max","period"))
plot_model(lm_temp,type="pred",terms=c("mnt_25_vosges","period"))
plot_model(lm_temp,type="pred",terms=c("fitted_microclim_scale","resid_microclim_scale"))
plot_model(lm_temp,type="pred",terms=c("resid_microclim","period"))
plot_model(lm_temp,type="pred",terms=c("fitted_microclim","period","resid_microclim"))
plot_model(lm_temp,type="pred",terms=c("resid_microclim","period","fitted_microclim"))

plot_model(lm_temp,type="pred",terms=c("resid_microclim_scale","mean_pH_scale"))

plot_model(lm_temp,type="pred",terms=c("fitted_microclim_scale"))
plot_model(lm_temp,type="pred",terms=c("resid_microclim_scale"))
plot_model(lm_temp,type="pred",terms=c("mean_pH_scale"))



plot_model(lm_with_resid,type="pred",terms=c("ipv_1000_25"))






library(POV)
POV(Axis2~fitted_microclim_scale*resid_microclim_scale,Data=plot_comp, Complete = TRUE)
POV(Axis2~fitted_microclim_max_scale*resid_microclim_max_scale,Data=plot_comp, Complete = TRUE)


POV(cit_climplant~fitted_microclim_scale*resid_microclim_scale+mean_pH,Data=plot_comp[n_sp_climplant>=5,], Complete = TRUE)


ggplot(plot_comp,aes(x=as.character(period),y=cit_climplant))+theme_bw()+geom_boxplot()+geom_violin()

ggplot(plot_comp,aes(x=as.character(period),y=RS2))+theme_bw()+geom_boxplot()+geom_violin()

plot_comp[,mean(RS2,na.rm=T),by=period]
plot_comp[,mean(cit_climplant,na.rm=T),by=period]

plot_comp[,mean(cit_climplant,na.rm=T),by=.(period,ifelse(predict_microclim_mean<17,"cold","warm"))]

plot_comp[,mean(cit_climplant,na.rm=T),by=.(period,ifelse(resid_microclim <0,"cold","warm"))]


plot_comp[,mean(cit_ecoplant_picq,na.rm=T),by=period]
plot_comp[n_sp_climplant>=5,mean(cit_climplant,na.rm=T),by=period]
plot_comp[n_sp_climplant>=5,mean(ANNEE,na.rm=T),by=period]
plot_comp[n_sp_climplant>=5,.N,by=period]

wilcox.test(plot_comp[period=="Past",mnt_25_vosges],plot_comp[period=="Recent",mnt_25_vosges])

mean(plot_comp[period=="Past",mnt_25_vosges])
mean(plot_comp[period=="Recent",mnt_25_vosges][1:150])
ggplot(plot_comp_sf,aes(color=cit_climplant,shape=period))+theme_bw()+geom_sf(size=2)+scale_color_viridis_c()
ggplot(plot_comp_sf,aes(color=resid_microclim,shape=period))+theme_bw()+geom_sf(size=2)+scale_color_gradient2(limits=c(-1,1),oob=scales::squish,low="cadetblue",high = "tomato")

ggplot(plot_comp_sf,aes(color=period,shape=period))+theme_bw()+geom_sf(size=2,alpha=0.5)



ggplot(plot_comp,aes(x=mean_pH,y=mnt_25_vosges,color=period))+theme_bw()+geom_point(alpha=0.5)
ggplot(plot_comp,aes(x=mean_pH,y=mnt_25_vosges))+theme_bw()+geom_density()

ggplot(plot_comp,aes(x=mean_pH,y=cit_climplant))()+geom_point(alpha=0.5)



sp_indicator_value<-indicator_value

plot_comp[,table(period,ifelse(predict_microclim_mean<18,"cold","warm"))]

id_past<-plot_comp[period=="Past" & predict_microclim_mean<18 ,plot_ID]
id_recent<-plot_comp[period=="Recent" & predict_microclim_mean<18,plot_ID][1:84]
  
id_past<-plot_comp[period=="Past" & predict_microclim_mean>18 ,plot_ID]
id_recent<-plot_comp[period=="Recent" & predict_microclim_mean>18,plot_ID][1:66]



table_sp_vosges_past<-table_sp_vosges[plot_comp[period=="Past",plot_ID],]
table_sp_vosges_recent<-table_sp_vosges[plot_comp[period=="Recent",plot_ID][1:150],]


table_compare<-create_comparaison_table(table_sp_vosges_past,table_sp_vosges_recent)

contrib_thermo_vosges<-get_contrib_one_ser(id_past,id_recent,table_sp_vosges,topt_name = "YearMeanMean")

View(contrib_thermo_vosges[[2]])



#### species response curves ####

frequent_sp<-sort(apply(table_sp_vosges,2,sum))
table_sp_vosges$plot_ID<-rownames(table_sp_vosges)

sp_of_interest<-table_sp_vosges[,c("plot_ID","Dryopteris filix-mas")]
colnames(sp_of_interest)<-c("plot_ID","presence_sp")
plot_comp<-merge(plot_comp[,-c("presence_sp")],sp_of_interest,by="plot_ID",all.x=T)
plot_comp_model_sp<-plot_comp[!is.na(presence_sp),]
plot_comp_model_sp<-plot_comp_model_sp[resid_microclim< (-0.5),]

# fitted_microclim resid_microclim    predict_microclim_mean
ggplot(plot_comp_model_sp,aes(x=predict_microclim_max,y=presence_sp))+theme_bw()+geom_point()+geom_smooth()
ggplot(plot_comp_model_sp,aes(x=fitted_microclim,y=presence_sp))+theme_bw()+geom_point()+geom_smooth()
ggplot(plot_comp_model_sp,aes(x=resid_microclim,y=presence_sp))+theme_bw()+geom_point()+geom_smooth()

glm_sp<-glm(presence_sp~ poly(fitted_microclim,2) + poly(resid_microclim,2) , family = "binomial",data=plot_comp_model_sp)
glm_sp<-glm(presence_sp~ poly(predict_microclim_mean,2) , family = "binomial",data=plot_comp_model_sp)
glm_sp<-glm(presence_sp~ poly(mnt_25_vosges,2), family = "binomial",data=plot_comp_model_sp)
glm_sp<-glm(presence_sp~ poly(resid_microclim,2) , family = "binomial",data=plot_comp_model_sp)
# +  poly(Heat_load_index_25,2)
summary(glm_sp)
plot_comp_model_sp[,pred_sp:=fitted(glm_sp)]
roc(plot_comp_model_sp$presence_sp,plot_comp_model_sp$pred_sp)


ggplot(plot_comp_model_sp,aes(x=predict_microclim_mean,y=pred_sp))+theme_bw()+geom_point()
ggplot(plot_comp_model_sp,aes(x=resid_microclim,y=pred_sp))+theme_bw()+geom_point()
ggplot(plot_comp_model_sp,aes(x=fitted_microclim,y=pred_sp,color=resid_microclim))+theme_bw()+geom_point()

library(pROC)

roc(plot_comp_model_sp$presence_sp,plot_comp_model_sp$pred_sp)
