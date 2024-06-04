data_meteo_france<-fread(file.path("data","weather_station","MENSQ_68_previous-1950-2022.csv"))


data_meteo_france<-data_meteo_france[NOM_USUEL=="MARKSTEIN CRETE",]

data_meteo_france[,AAAAMM_chr:=as.character(AAAAMM)]
data_meteo_france[,annee:=as.numeric(substr(AAAAMM_chr,1,4))]
data_meteo_france[,mois:=as.numeric(substr(AAAAMM_chr,5,6))]
data_meteo_france[,moy_pres:= !is.na(TM)]


data_meteo_france[annee%in%c(2005:2020) & mois %in% c(4:8) ,mean(TM)]
data_meteo_france[annee%in%c(2022) &  mois %in% c(4:8) ,mean(TM)]

