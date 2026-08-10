##################################
#
# Overdose and Race
# By Drug
#
# Exploratory plots
#
##################################

library(ggplot2)
library(lemon)
library(plyr)
library(rgdal)
library(sf)
library(ggthemes)

# set data location
setwd("C:/Users/bnwhite/Desktop/bnw/projects/Kline/polyaddiction/data")

#load data
od.pres=read.csv('ODRacePrescription_07-19.csv')
od.hero=read.csv('ODRaceHeroin_07-19.csv')
od.fent=read.csv('ODRaceFentanyl_07-19.csv')
od.coca=read.csv('ODRaceCocaine_07-19.csv')
od.stim=read.csv('ODRacePsychostimulant_07-19.csv')
od.fent.hero=read.csv('ODRaceFentanylHeroin_07-19.csv')
od.fent.coca=read.csv('ODRaceFentanylCocaine_07-19.csv')
od.fent.stim=read.csv('ODRaceFentanylPsychostimulant_07-19.csv')

#rename
names(od.pres)=c('Race','County','Year','pDeaths','Population')
names(od.hero)=c('Race','County','Year','hDeaths','Population')
names(od.fent)=c('Race','County','Year','fDeaths','Population')
names(od.coca)=c('Race','County','Year','cDeaths','Population')
names(od.stim)=c('Race','County','Year','sDeaths','Population')
names(od.fent.hero)=c('Race','County','Year','fhDeaths','Population')
names(od.fent.coca)=c('Race','County','Year','fcDeaths','Population')
names(od.fent.stim)=c('Race','County','Year','fsDeaths','Population')

#join
od.drug=join(od.pres[,c('Race','Year','County','pDeaths','Population')],od.hero[,c('Race','Year','County','hDeaths','Population')],by=c('Race','Year','County','Population'),type='left')
od.drug=join(od.drug,od.fent[,c('Race','Year','County','fDeaths','Population')],by=c('Race','Year','County','Population'),type='left')
od.drug=join(od.drug,od.coca[,c('Race','Year','County','cDeaths','Population')],by=c('Race','Year','County','Population'),type='left')
od.drug=join(od.drug,od.stim[,c('Race','Year','County','sDeaths','Population')],by=c('Race','Year','County','Population'),type='left')
od.drug=join(od.drug,od.fent.hero[,c('Race','Year','County','fhDeaths','Population')],by=c('Race','Year','County','Population'),type='left')
od.drug=join(od.drug,od.fent.coca[,c('Race','Year','County','fcDeaths','Population')],by=c('Race','Year','County','Population'),type='left')
od.drug=join(od.drug,od.fent.stim[,c('Race','Year','County','fsDeaths','Population')],by=c('Race','Year','County','Population'),type='left')

#remove non-counties
od.drug=od.drug[!(od.drug$County %in% c('NonOH','Unknown','Total')),]

#limit to Black/White
od.drug=od.drug[od.drug$Race %in% c('Black','White'),]

#map merge ID
od.drug$NAME=od.drug$County

#adjust prescription deaths
od.drug$pDeathsAdj=od.drug$pDeaths-od.drug$fDeaths

#calculate rates
od.drug$pRate=od.drug$pDeathsAdj/od.drug$Population*100000
od.drug$hRate=od.drug$hDeaths/od.drug$Population*100000
od.drug$fRate=od.drug$fDeaths/od.drug$Population*100000
od.drug$cRate=od.drug$cDeaths/od.drug$Population*100000
od.drug$sRate=od.drug$sDeaths/od.drug$Population*100000
od.drug$fhRate=od.drug$fhDeaths/od.drug$Population*100000
od.drug$fcRate=od.drug$fcDeaths/od.drug$Population*100000
od.drug$fsRate=od.drug$fsDeaths/od.drug$Population*100000

#state counts
state.drug=aggregate(cbind(pDeathsAdj,hDeaths,fDeaths,cDeaths,sDeaths,fhDeaths,fcDeaths,fsDeaths,Population)~Year+Race,data=od.drug,FUN=sum)

#state rates
state.drug$pRate=state.drug$pDeathsAdj/state.drug$Population*100000
state.drug$hRate=state.drug$hDeaths/state.drug$Population*100000
state.drug$fRate=state.drug$fDeaths/state.drug$Population*100000
state.drug$cRate=state.drug$cDeaths/state.drug$Population*100000
state.drug$sRate=state.drug$sDeaths/state.drug$Population*100000
state.drug$fhRate=state.drug$fhDeaths/state.drug$Population*100000
state.drug$fcRate=state.drug$fcDeaths/state.drug$Population*100000
state.drug$fsRate=state.drug$fsDeaths/state.drug$Population*100000
state.drug$c.onlyRate=(state.drug$cDeaths-state.drug$fcDeaths)/state.drug$Population*100000
state.drug$s.onlyRate=(state.drug$sDeaths-state.drug$fsDeaths)/state.drug$Population*100000

#output directory
setwd("C:/Users/bnwhite/Desktop/bnw/projects/Kline/polyaddiction/output/exploratory")

#state time series
ts.pres=ggplot(data=state.drug,aes(x=Year,y=pRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Prescription Opioid Overdose',x='Year',y='Rate/100k')

ts.hero=ggplot(data=state.drug,aes(x=Year,y=hRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Heroin Overdose',x='Year',y='Rate/100k')

ts.fent=ggplot(data=state.drug,aes(x=Year,y=fRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Fentanyl Overdose',x='Year',y='Rate/100k')

ts.coca=ggplot(data=state.drug,aes(x=Year,y=cRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Cocaine Overdose',x='Year',y='Rate/100k')

ts.stim=ggplot(data=state.drug,aes(x=Year,y=sRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Psychostimulant Overdose',x='Year',y='Rate/100k')

ts.fent.hero=ggplot(data=state.drug,aes(x=Year,y=fhRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Fentanyl and Heroin Overdose',x='Year',y='Rate/100k')

ts.fent.coca=ggplot(data=state.drug,aes(x=Year,y=fcRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Fentanyl and Cocaine Overdose',x='Year',y='Rate/100k')

ts.fent.stim=ggplot(data=state.drug,aes(x=Year,y=fsRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Fentanyl and Psychostimulant Overdose',x='Year',y='Rate/100k')



ts.coca.only=ggplot(data=state.drug,aes(x=Year,y=c.onlyRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Cocaine Overdose (w/o Fentanyl)',x='Year',y='Rate/100k')

ts.stim.only=ggplot(data=state.drug,aes(x=Year,y=s.onlyRate,color=Race))+
  geom_path()+
  geom_vline(xintercept=2013,linetype='dashed')+
  labs(title='Psychostimulant Overdose (w/o Fentanyl)',x='Year',y='Rate/100k')

png('state_opioid_rates.png',width=6,height=9,res=300,units='in')

grid_arrange_shared_legend(ts.pres,ts.hero,ts.fent,nrow=3,ncol=1,position='right')

dev.off()


png('state_stimulant_rates.png',width=6,height=6,res=300,units='in')

grid_arrange_shared_legend(ts.coca,ts.stim,nrow=2,ncol=1,position='right')

dev.off()

png('state_poly_rates.png',width=6,height=9,res=300,units='in')

grid_arrange_shared_legend(ts.fent.hero,ts.fent.coca,ts.fent.stim,nrow=3,ncol=1,position='right')

dev.off()

png('state_stimulant_only_rates.png',width=6,height=6,res=300,units='in')

grid_arrange_shared_legend(ts.coca.only,ts.stim.only,nrow=2,ncol=1,position='right')

dev.off()


#maps
oh.map=st_read('C:/Users/dkline.MEDCTR/OneDrive - Wake Forest Baptist Health/Shapefiles/2020 County/cb_2020_us_county_500k.shp')
oh.map=oh.map[oh.map$STATEFP=='39',]

years=2013:2019

#stim cats
cats=c(0,0.1,1.1,2.3,4.3,7.9,14.2,20.8,180)
od.drug$sCats=cut(od.drug$sRate,cats,include.lowest=TRUE)

#cocaine cats
cats=c(0,0.1,0.4,3.1,6.9,12,174)
od.drug$cCats=cut(od.drug$cRate,cats,include.lowest=TRUE)


#stim/fent cats
cats=c(0,0.1,2.1,3.8,6.9,11.3,18.1,22,179)
od.drug$fsCats=cut(od.drug$sRate,cats,include.lowest=TRUE)

#cocaine/fent cats
cats=c(0,0.1,0.5,2.5,6.7,11,16.4,174)
od.drug$fcCats=cut(od.drug$cRate,cats,include.lowest=TRUE)

white.stim.map=list()
black.stim.map=list()

white.cocaine.map=list()
black.cocaine.map=list()

white.stim.fent.map=list()
black.stim.fent.map=list()

white.cocaine.fent.map=list()
black.cocaine.fent.map=list()


oh.map.data=merge(oh.map,od.drug,by='NAME')


for(t in 1:length(years)){
  
white.stim.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='White',],aes(fill=sCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' White'))


black.stim.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='Black',],aes(fill=sCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' Black'))

white.cocaine.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='White',],aes(fill=cCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' White'))

black.cocaine.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='Black',],aes(fill=cCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' Black'))

white.stim.fent.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='White',],aes(fill=fsCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' White'))


black.stim.fent.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='Black',],aes(fill=fsCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' Black'))

white.cocaine.fent.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='White',],aes(fill=fcCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' White'))

black.cocaine.fent.map[[t]]=ggplot(data=oh.map.data[oh.map.data$Year==years[t] & oh.map.data$Race=='Black',],aes(fill=fcCats),color='black')+
  geom_sf()+
  scale_fill_brewer(name='',palette = 'Reds',direction=1,drop=FALSE)+
  theme_map()+
  labs(title=paste0(years[t],' Black'))

}

png('stim_maps.png',height=5,width=12,units='in',res=300)

grid_arrange_shared_legend(white.stim.map[[1]],white.stim.map[[2]],white.stim.map[[3]],white.stim.map[[4]],white.stim.map[[5]],white.stim.map[[6]],white.stim.map[[7]],
                           black.stim.map[[1]],black.stim.map[[2]],black.stim.map[[3]],black.stim.map[[4]],black.stim.map[[5]],black.stim.map[[6]],black.stim.map[[7]],nrow=2,ncol=7,position='right')

dev.off()


png('cocaine_maps.png',height=5,width=12,units='in',res=300)

grid_arrange_shared_legend(white.cocaine.map[[1]],white.cocaine.map[[2]],white.cocaine.map[[3]],white.cocaine.map[[4]],white.cocaine.map[[5]],white.cocaine.map[[6]],white.cocaine.map[[7]],
                           black.cocaine.map[[1]],black.cocaine.map[[2]],black.cocaine.map[[3]],black.cocaine.map[[4]],black.cocaine.map[[5]],black.cocaine.map[[6]],black.cocaine.map[[7]],nrow=2,ncol=7,position='right')

dev.off()

png('stim_fent_maps.png',height=5,width=12,units='in',res=300)

grid_arrange_shared_legend(white.stim.fent.map[[1]],white.stim.fent.map[[2]],white.stim.fent.map[[3]],white.stim.fent.map[[4]],white.stim.fent.map[[5]],white.stim.fent.map[[6]],white.stim.fent.map[[7]],
                           black.stim.fent.map[[1]],black.stim.fent.map[[2]],black.stim.fent.map[[3]],black.stim.fent.map[[4]],black.stim.fent.map[[5]],black.stim.fent.map[[6]],black.stim.fent.map[[7]],nrow=2,ncol=7,position='right')

dev.off()


png('cocaine_fent_maps.png',height=5,width=12,units='in',res=300)

grid_arrange_shared_legend(white.cocaine.fent.map[[1]],white.cocaine.fent.map[[2]],white.cocaine.fent.map[[3]],white.cocaine.fent.map[[4]],white.cocaine.fent.map[[5]],white.cocaine.fent.map[[6]],white.cocaine.fent.map[[7]],
                           black.cocaine.fent.map[[1]],black.cocaine.fent.map[[2]],black.cocaine.fent.map[[3]],black.cocaine.fent.map[[4]],black.cocaine.fent.map[[5]],black.cocaine.fent.map[[6]],black.cocaine.fent.map[[7]],nrow=2,ncol=7,position='right')

dev.off()




