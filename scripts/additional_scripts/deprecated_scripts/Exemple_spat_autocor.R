###########################################################################################################
##################                   R Code accompanying De Frenne et al.                ##################
##################        Nutrient fertilization by dogs in peri-urban ecosystems        ##################
###########################################################################################################

setwd("D:/Users/pdfrenne/Documents/R/Dog-fertilization")
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)
library(readxl)
library(nlme)


mean.na <- function(x) mean(x, na.rm=T)
max.na <- function(x) max(x, na.rm=T)
se.na <- function(x) sd(x, na.rm=T)/sqrt(length(x))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



Data<-read_excel("Countdata-V7.xlsx",sheet="Countdata")
Data$Date<-as.Date(Data$Date, "%d/%m/%y")
Data$Day <- factor(Data$Day,levels = c("Monday","Tuesday","Wednesday","Thursday","Friday", "Saturday", "Sunday"))
Data$Site<-as.factor(Data$Site)
Data$Site.name.full<-as.factor(Data$Site.name.full)
Data$SiteID<-as.factor(Data$SiteID)
Data$Area.paths<-as.numeric(as.character(Data$Area.paths))


Data<-data.frame(Data)
str(Data)
head(Data)

###########################################################################################################
###############           To obtain dog densities                                   ######################
###########################################################################################################
#Sum of dogs: use the raw counts, not taking into account the hours
Data$Sum<-Data$Dog.off.leash+Data$Dog.on.leash  #Sum of off and on leash dogs
#Sum of dogs PER HOUR: taking into account the hours and search efforts:
Data$Sum.hour<-Data$Dog.off.hour+Data$Dog.on.hour  #Sum of off and on leash  per hour

daylength = 12 #(hours) average daylength 

#Add standardized dog count to Data file, per unit area (ha) and unit time (day)
Data$Off.ha.day<-Data$Dog.off.hour*daylength/(Data$Area) #off leash
Data$On.ha.day<-Data$Dog.on.hour*daylength/(Data$Area)   #on leash
Data$Sum.ha.day<-Data$Sum.hour*daylength/(Data$Area)    #Sum
Data$ryear = Data$Sum.ha.day*365.25 #(counts) average number of dogs per ha per year




############## Scenario analysis - ONLY run this if you want to perform the scenario analysis of no off leash dogs
#Assuming no off leash dogs and all deposits are within an area 2m left and 2m right of each path
# Data$Off.ha.day<-Data$Dog.off.hour*daylength/(Data$Area.paths) #off leash
# Data$On.ha.day<-Data$Dog.on.hour*daylength/(Data$Area.paths)   #on leash
# Data$Sum.ha.day<-Data$Sum.hour*daylength/(Data$Area.paths)    #Sum
# Data$ryear = Data$Sum.ha.day*365.25 #(counts) average number of dogs per ha per year




###########################################################################################################
###############           To obtain fertilization rates                            ######################
###########################################################################################################
#1) deposition from urine
vol.urine = 0.736 #(L) volume of urine produced per dog per day following Paradeis et al. (2013)
Data$urine.dep = vol.urine*Data$ryear/4 #(L) volume of urine deposited in the area per year for all dogs, divide by four assuming dogs deposit one quarter of their urine during this walk
Pconc.urine = 484.61 #(mg P/L) mean P concentration in dog urine
Nconc.urine = 18681.41 #(mg N/L) mean N concentration in dog urine

Data$Pdep.urine = Pconc.urine*Data$urine.dep/1000/1000  #(kg P) total amount of P deposited in the area #(kg P / ha / yr)
Data$Ndep.urine = Nconc.urine*Data$urine.dep/1000/1000  #(kg N) total amount of N deposited in the area #(kg N / ha / yr)




#2) deposition from faeces
mass.faeces = 0.1 #(kg) 100g of dry faeces mass per dog, assuming one deposit per walk
Data$total.mass.faeces = mass.faeces*Data$ryear #(kg) total mass of faeces deposited in the area

Pconc.faeces = 31.95 #(mg P / g dry mass dung) mean concentration in dog faeces
Nconc.faeces = 44.27 #(mg N / g dry mass dung) mean concentration in dog faeces


Data$Pdep.faeces = Pconc.faeces*Data$total.mass.faeces/1000  #(kg P) total amount of P deposited in the area via faeces 
#(kg P / ha / yr)
Data$Ndep.faeces = Nconc.faeces*Data$total.mass.faeces/1000  #(kg N) total amount of N deposited in the area from faeces
#(kg N / ha / yr)




#4) summed deposition
Data$totalP = Data$Pdep.faeces + Data$Pdep.urine ##(kg P / ha / yr)
Data$totalN = Data$Ndep.faeces + Data$Ndep.urine ##(kg N / ha / yr)

#Raw means
summary(Data$Pdep.faeces)
summary(Data$Ndep.faeces)

summary(Data$Pdep.urine)
summary(Data$Ndep.urine)

summary(Data$totalP)
summary(Data$totalN)


###########################################################################################################
######     To obtain results with mixed models including temporal autocorrelation       #################
###########################################################################################################
#Faeces
summary(lme(Pdep.faeces~1, random = ~Days.since.first.count2|Site,data=Data, na.action=na.exclude, correlation=corCAR1(form=~Days.since.first.count2|Site)))
summary(lme(Ndep.faeces~1, random = ~Days.since.first.count2|Site,data=Data, na.action=na.exclude, correlation=corCAR1(form=~Days.since.first.count2|Site)))


#Urine
summary(lme(Pdep.urine~1, random = ~Days.since.first.count2|Site,data=Data, na.action=na.exclude, correlation=corCAR1(form=~Days.since.first.count2|Site)))
summary(lme(Ndep.urine~1, random = ~Days.since.first.count2|Site,data=Data, na.action=na.exclude, correlation=corCAR1(form=~Days.since.first.count2|Site)))

#Total
summary(lme(totalP~1, random = ~Days.since.first.count2|Site,data=Data, na.action=na.exclude, correlation=corCAR1(form=~Days.since.first.count2|Site)))
summary(lme(totalN~1, random = ~Days.since.first.count2|Site,data=Data, na.action=na.exclude, correlation=corCAR1(form=~Days.since.first.count2|Site)))

#There is almost no temporal autocorrelation because the result (intercept) is nearly the same without the temporal autocorrelation structure
#Faeces
summary(lmer(Pdep.faeces~1+(1|Site),data=Data,na.action=na.exclude))
summary(lmer(Ndep.faeces~1+(1|Site),data=Data,na.action=na.exclude))

#Urine
summary(lmer(Pdep.urine~1+(1|Site),data=Data,na.action=na.exclude))
summary(lmer(Ndep.urine~1+(1|Site),data=Data,na.action=na.exclude))

#Total
summary(lmer(totalP~1+(1|Site),data=Data,na.action=na.exclude))
summary(lmer(totalN~1+(1|Site),data=Data,na.action=na.exclude))




###########################################################################################################
###############               simple calculations for main text                   ######################
###########################################################################################################
#How many dogs did we count in total?
sum(Data$Sum)

summary(Data$Off.ha.day)
summary(Data$On.ha.day)
summary(Data$Sum.ha.day)
tapply(Data$Sum.ha.day,Data$SiteID,mean.na) #Per site


#What percentage of dogs is off vs on leash?
mean.na(Data$Off.ha.day/Data$Sum.ha.day*100)
mean.na(Data$On.ha.day/Data$Sum.ha.day*100)

#Per site
tapply(Data$Off.ha.day/Data$Sum.ha.day*100,Data$SiteID,mean.na)

#How many days did we count on each day of the week?
tapply(Data$Sum.ha.day,Data$Day,length)

#How many dogs per year?
summary(Data$ryear)


#How many input of N and P in each site?
tapply(Data$totalP,Data$SiteID,mean)
tapply(Data$totalN,Data$SiteID,mean)






###########################################################################################################
###############                               figures                              ######################
###########################################################################################################
#Fig. 2 (dog densities)
#First rearrange the data a bit for calculation of mean and SD and for plotting
toplot<-matrix(0,3*nrow(Data),3)
colnames(toplot)<-c('off.on','counts','site')
toplot[1:nrow(Data),1]<-rep("off",nrow(Data))
toplot[(nrow(Data)+1):(2*nrow(Data)),1]<-rep("on",nrow(Data))
toplot[(nrow(Data)*2+1):(3*nrow(Data)),1]<-rep("sum",nrow(Data))
toplot[1:nrow(Data),2]<-Data$Off.ha.day
toplot[(nrow(Data)+1):(2*nrow(Data)),2]<-Data$On.ha.day
toplot[(nrow(Data)*2+1):(3*nrow(Data)),2]<-Data$Sum.ha.day
toplot[1:nrow(Data),3]<-Data$SiteID
toplot[(nrow(Data)+1):(2*nrow(Data)),3]<-Data$SiteID
toplot[(nrow(Data)*2+1):(3*nrow(Data)),3]<-Data$SiteID
toplot<-as.data.frame(toplot)
toplot$counts<-as.numeric(as.character(toplot$counts))
toplot$off.on<-factor(toplot$off.on)
toplot$site<-factor(toplot$site)

df2 <- data_summary(toplot, varname="counts", 
                    groupnames=c("off.on", "site"))


p2<- ggplot(df2, aes(x=off.on, y=counts, fill=site)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=counts-se, ymax=counts+se), width=.2,
                position=position_dodge(.9)) +
labs(x="Off vs on leash", y = "Number of dogs per hectare per day")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00','#FFDB6D','#00AFBB'))+
  geom_segment(aes(x=0.5,xend=1.5,y=mean(Data$Off.ha.day),yend=mean(Data$Off.ha.day)),color = "grey",lwd=1.5, lty="dashed")+
  geom_segment(aes(x=1.5,xend=2.5,y=mean(Data$On.ha.day),yend=mean(Data$On.ha.day)),color = "grey",lwd=1.5, lty="dashed")+
  geom_segment(aes(x=2.5,xend=3.5,y=mean(Data$Sum.ha.day),yend=mean(Data$Sum.ha.day)),color = "grey",lwd=1.5, lty="dashed")

print(p2)

#Fig. 2 in main text
pdf(file= "Fig2.pdf",width = 10, height = 5)  #save as pdf in wd
p2
dev.off()









##########################################################################################################
#Fig. 3 (fertilisation)
#Phosphorus
toplot<-matrix(0,3*nrow(Data),3)
colnames(toplot)<-c('source','deposition','site')
toplot[1:nrow(Data),1]<-rep("urine",nrow(Data))
toplot[(nrow(Data)+1):(2*nrow(Data)),1]<-rep("faeces",nrow(Data))
toplot[(nrow(Data)*2+1):(3*nrow(Data)),1]<-rep("sum",nrow(Data))
toplot[1:nrow(Data),2]<-Data$Pdep.urine
toplot[(nrow(Data)+1):(2*nrow(Data)),2]<-Data$Pdep.faeces
toplot[(nrow(Data)*2+1):(3*nrow(Data)),2]<-Data$totalP
toplot[1:nrow(Data),3]<-Data$SiteID
toplot[(nrow(Data)+1):(2*nrow(Data)),3]<-Data$SiteID
toplot[(nrow(Data)*2+1):(3*nrow(Data)),3]<-Data$SiteID
toplot<-as.data.frame(toplot)
toplot$deposition<-as.numeric(as.character(toplot$deposition))
toplot$source<-factor(toplot$source)
toplot$site<-factor(toplot$site)

df2 <- data_summary(toplot, varname="deposition", 
                    groupnames=c("source", "site"))


pP<- ggplot(df2, aes(x=ordered(source, levels=c('faeces','urine','sum')), y=deposition, fill=site)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=deposition-se, ymax=deposition+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Source of nutrients", y = "Deposition (kilogram per hectare per year)",title="Phosphorus")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00','#FFDB6D','#00AFBB'))+
  geom_segment(aes(x=0.5,xend=1.5,y=mean(Data$Pdep.faeces),yend=mean(Data$Pdep.faeces)),color = "grey",lwd=1.5, lty="dashed")+
  geom_segment(aes(x=1.5,xend=2.5,y=mean(Data$Pdep.urine),yend=mean(Data$Pdep.urine)),color = "grey",lwd=1.5, lty="dashed")+
  geom_segment(aes(x=2.5,xend=3.5,y=mean(Data$totalP),yend=mean(Data$totalP)),color = "grey",lwd=1.5, lty="dashed")+
      theme(legend.position = "none")
print(pP)


#Nitrogen
toplot<-matrix(0,3*nrow(Data),3)
colnames(toplot)<-c('source','deposition','site')
toplot[1:nrow(Data),1]<-rep("urine",nrow(Data))
toplot[(nrow(Data)+1):(2*nrow(Data)),1]<-rep("faeces",nrow(Data))
toplot[(nrow(Data)*2+1):(3*nrow(Data)),1]<-rep("sum",nrow(Data))
toplot[1:nrow(Data),2]<-Data$Ndep.urine
toplot[(nrow(Data)+1):(2*nrow(Data)),2]<-Data$Ndep.faeces
toplot[(nrow(Data)*2+1):(3*nrow(Data)),2]<-Data$totalN
toplot[1:nrow(Data),3]<-Data$SiteID
toplot[(nrow(Data)+1):(2*nrow(Data)),3]<-Data$SiteID
toplot[(nrow(Data)*2+1):(3*nrow(Data)),3]<-Data$SiteID
toplot<-as.data.frame(toplot)
toplot$deposition<-as.numeric(as.character(toplot$deposition))
toplot$source<-factor(toplot$source)
toplot$site<-factor(toplot$site)

df2 <- data_summary(toplot, varname="deposition", 
                    groupnames=c("source", "site"))

pN<- ggplot(df2, aes(x=ordered(source, levels=c('faeces','urine','sum')), y=deposition, fill=site)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=deposition-se, ymax=deposition+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Source of nutrients", y = "",title="Nitrogen")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00','#FFDB6D','#00AFBB'))+
  geom_segment(aes(x=0.5,xend=1.5,y=mean(Data$Ndep.faeces),yend=mean(Data$Ndep.faeces)),color = "grey",lwd=1.5, lty="dashed")+
  geom_segment(aes(x=1.5,xend=2.5,y=mean(Data$Ndep.urine),yend=mean(Data$Ndep.urine)),color = "grey",lwd=1.5, lty="dashed")+
  geom_segment(aes(x=2.5,xend=3.5,y=mean(Data$totalN),yend=mean(Data$totalN)),color = "grey",lwd=1.5, lty="dashed")
  
print(pN)


#Fig. 3 in main text
pdf(file= "Fig3.pdf",width = 10, height = 5)
grid.arrange(pP, pN, nrow = 1)
dev.off()






###########################################################################################################
###############           To obtain distributions with bootstrapping approach        ######################
###########################################################################################################
n = 999 #Number of bootstraps

results <- matrix(0,n,2)
colnames(results)<-c('P','N')

#make vectors to sample N and P concentrations from normal distribution 
Pconc.faeces.vector <- rnorm(1000, mean = Pconc.faeces, sd = 8.38)#(mg P / g dung)
Nconc.faeces.vector <- rnorm(1000, mean = Nconc.faeces, sd = 8.11) #(mg N / g dung)
Pconc.urine.vector <- rnorm(1000, mean = Pconc.urine, sd = 648.58)#(mg P / g dung)
#To avoid negative values in Pconc.urine.vector (because SD is > mean, only in this case)
Pconc.urine.vector <- Pconc.urine.vector[which(Pconc.urine.vector > 0)]
Nconc.urine.vector <- rnorm(1000, mean = Nconc.urine, sd = 3011.43) #(mg N / g dung)

#make vectors also for volume of urine and mass of faeces
uncertainty=0.2  #uncertainty, i.e. the standard deviation is equal to 20% of the mean only for the volume of the urine, mass of faeces and dog residence time (because no literature values available for the latter)
vol.urine.vector <- rnorm(1000, mean = vol.urine, sd = uncertainty*vol.urine) 
mass.faeces.vector <- rnorm(1000, mean = mass.faeces, sd = uncertainty*mass.faeces) 

for (i in 1:n){
  j<-sample(c(1:nrow(Data)),size=1,replace=T)
  r<-Data$Sum.ha.day[j]#Sample a dog count
  #Take into account uncertainty of duration of visit
  t.vector <- rnorm(1000, mean = Data$Hours[j], sd = uncertainty*Data$Hours[j]) #expressed in hours
  t<-sample(t.vector,size = 1, replace = T) #Sample a duration of the walk
  ryear = r*365.25 #(counts) average number of dogs per year per ha
  
  
  #Sample concentrations to take into account uncertainty in the nutrient concentrations in the faeces and urine
  Pconc.faeces <-sample(Pconc.faeces.vector,size = 1, replace = T)
  Nconc.faeces <- sample(Nconc.faeces.vector,size = 1, replace = T)
  Pconc.urine <- sample(Pconc.urine.vector,size = 1, replace = T)
  Nconc.urine <- sample(Nconc.urine.vector,size = 1, replace = T)
  
  vol.urine.sample <- sample(vol.urine.vector,size = 1, replace = T)
  mass.faeces.sample <- sample(mass.faeces.vector,size = 1, replace = T)
  
  
  #2) deposition from urine
  urine.dep = vol.urine.sample*ryear/4 #(L) volume of urine deposited in the area per year for all dogs, divide by four assuming dogs deposit one quarter of their urine during this walk
  Pdep.urine = Pconc.urine*urine.dep/1000/1000  #(kg P) total amount of P deposited in the area
  Ndep.urine = Nconc.urine*urine.dep/1000/1000  #(kg N) total amount of N deposited in the area
  
  
  #3) deposition from faeces
  total.mass.faeces = mass.faeces.sample*ryear #(kg) total mass of faeces deposited in the area
  Pdep.faeces = Pconc.faeces*total.mass.faeces/1000  #(kg P) total amount of P deposited in the area via faeces
  Ndep.faeces = Nconc.faeces*total.mass.faeces/1000  #(kg N) total amount of N deposited in the area from faeces
  
  
  #4) summed deposition
  results[i,1] <- Pdep.faeces + Pdep.urine ##(kg P / ha / yr)
  results[i,2] <- Ndep.faeces + Ndep.urine ##(kg N / ha / yr)
  
}

#View(results)

#11%iles and 89%iles
mean(results[,1]) #P
mean(results[,2]) #N

quantile(results[,1], probs = c(0.025, 0.05,0.11, 0.5,0.89,0.95,0.975)) #P
quantile(results[,2], probs = c(0.025, 0.05, 0.11,0.5,0.89,0.95, 0.975)) #N




########################################################################################################
#############################                FIGURES SI               ##################################
#########################################################################################################
#Supporting Information Fig. S1: uncertainty in model parameters
pdf(file= "FigS1.pdf")  #save as pdf in wd
par(mfrow=c(3,3))
hist(Pconc.faeces.vector, main="",xlab="P concentration in faeces (mg P / g)", col="#E69F00")
hist(Nconc.faeces.vector, main="",xlab="N concentration in faeces (mg N / g)", col="#E69F00")
hist(Pconc.urine.vector, main="",xlab="P concentration in urine (mg P / L)", col="#FFDB6D")
hist(Nconc.urine.vector, main="",xlab="N concentration in urine (mg N / L)", col="#FFDB6D")
hist(vol.urine.vector, main="",xlab="Volume of urine produced per day (L)", col="#56B4E9")
hist(mass.faeces.vector, main="",xlab="Mass of faeces produced per day (kg)", col="#56B4E9")
hist(results[,1], main="",xlab="Total phosphorus input (kg P/ha/yr)", col="#00AFBB")
abline(v=mean(results[,1]),col="grey",lty=1, lwd=2)
abline(v=quantile(results[,1], probs = c(0.95)),col="grey",lty=2, lwd=2)
abline(v=quantile(results[,1], probs = c(0.05)),col="grey",lty=2, lwd=2)
hist(results[,2], main="",xlab="Total nitrogen input (kg N/ha/yr)", col="#00AFBB")
abline(v=mean(results[,2]),col="grey",lty=1, lwd=2)
abline(v=quantile(results[,2], probs = c(0.95)),col="grey",lty=2, lwd=2)
abline(v=quantile(results[,2], probs = c(0.05)),col="grey",lty=2, lwd=2)
dev.off()
