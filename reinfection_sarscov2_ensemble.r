library(pomp)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(plyr)
library(viridis)
library(modi)
library(bayestestR)
library(scales)



scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

scientific_10_0 <- function(x) {
  parse(text=gsub("\\+","",gsub("e", " %*% 10^", scales::scientific_format()(x))))
}

scientific_10_1 <- function(x) {
  parse(text=gsub("\\+","",gsub("1e", "10^", scientific_format()(x))))
}

scientific_10_2 <- function(x) {
  parse(text=c(gsub("\\+","",gsub("1e", "10^", scientific_format()(x[which(x<1)]))),
               x[which(x==1)],
               gsub("\\+","",gsub("1e", "10^", scientific_format()(x[which(x>1)])))))
}

scientific_10_3 <- function(x) {
  parse(text=c(gsub("\\+","",gsub("e", " %*% 10^", scales::scientific_format()(x[which(x<0)]))),
               x[which(x==0)],
               gsub("\\+","",gsub("e", " %*% 10^", scales::scientific_format()(x[which(x>0)])))))
}

#readin in data and model---------------------------
source("setup_sarscov2_ensemble.r")
#readin file for plotting
source("multiplot.r")

#Read in resulting fits------------

fit_all = read.csv("fit.csv")

fit_all$AIC = 2*fit_all$like+2*fit_all$enames

fit_all_trimmed = NULL

fit_all_new = NULL
for(y in 1:4){
  for(z in 1:160){
    sub = subset(fit_all,i==y&version==z&round==1)
    if(nrow(sub)>1){
      stop()
    }
    fit_all_new = rbind(fit_all_new,sub)
  }
}


for (k in 1:4){
  for(j in 1:160){
    sub = subset(fit_all,i==k&version==j)
    sub = arrange(sub,AIC)
    fit_all_trimmed = rbind(fit_all_trimmed,sub[1,])
  }
}

fit_all_trimmed = subset(fit_all_trimmed,yA==0)

fit_all_trimmed_ranked = NULL
for(k in 1:4){
  sub = subset(fit_all_trimmed,i==k)
  sub = arrange(sub,AIC)
  sub$rank = 1:nrow(sub)
  fit_all_trimmed_ranked = rbind(fit_all_trimmed_ranked,sub)
}
fit_all_trimmed = fit_all_trimmed_ranked

AIC_min = fit_all_trimmed%>%dplyr::group_by(version_dat)%>%dplyr::summarize(min = min(AIC))
AIC_min = arrange(AIC_min,min)
AIC_min$rank_min = as.character(1:nrow(AIC_min))

fit_all_trimmed = merge(fit_all_trimmed,AIC_min)


fit_all_trimmed$version_dat = factor(fit_all_trimmed$version_dat,levels = AIC_min$version_dat)
fit_all_trimmed$rank_min = factor(fit_all_trimmed$rank_min,levels =as.character(1:160))



#subset AIC to only include most successful of the 4 runs

fit_all_sub = NULL

for(k in 1:160){
  sub = subset(fit_all,version==k)
  sub = arrange(sub,AIC)
  fit_all_sub = rbind(fit_all_sub,sub[1,])
}

fit_all_sub = arrange(fit_all_sub,AIC)

fit_all_sub = subset(fit_all_sub,yA==0)

fit_all_sub$rank2 = 1:nrow(fit_all_sub)

fit_all_sub$diff = fit_all_sub$AIC-fit_all_sub$AIC[1]
fit_all_sub$weight = exp(-0.5*fit_all_sub$diff)/(sum(exp(-0.5*fit_all_sub$diff)))
fit_all_sub$summed_weight = fit_all_sub$weight[1]
for(i in 2:nrow(fit_all_sub)){
  fit_all_sub$summed_weight[i]=fit_all_sub$summed_weight[i-1]+fit_all_sub$weight[i]
}
fit_all_sub$evidenceratio = exp(-0.5*fit_all_sub$diff[1])/exp(-0.5*fit_all_sub$diff)

number_sims = which(fit_all_sub$summed_weight>=0.95)[1]

top_models = fit_all_sub[1:number_sims,]

#Simulate what unmodified top-ranked models look like
fit_all_og = NULL
for(z in 1:number_sims){
  sub = subset(fit_all_sub,rank2==z)
  print(z)
  params = as.vector(sub[1,1:42])
  tdat <- mcmv_NA %>%
    trajectory(params = params,format="data.frame")
  tdat = melt(tdat[,c(2:11)],id="time")
  colnames(tdat)=c("time","species","fit")
  fit_ode = merge(data.frame(tdat,type = "fit"),data.frame(data,type = "dat"),all=TRUE)
  fit_ode$species[fit_ode$id%in%c(paste0("N_V_",1:6))]="V"
  fit_ode$species[fit_ode$id%in%c(paste0("N_E_",1:6))]="E"
  fit_ode$species[fit_ode$id%in%c(paste0("N_T_",1:6))]="T"
  fit_ode$species[fit_ode$id%in%c(paste0("N_F1_",1:6))]="F1"
  fit_ode$species[fit_ode$id%in%c(paste0("N_F2_",1:6))]="F2"
  fit_ode$species[fit_ode$id%in%c(paste0("N_F3_",1:6))]="F3"
  fit_ode$species[fit_ode$id%in%c(paste0("N_A_",1:6))]="A"
  fit_ode$id[is.na(fit_ode$id)]=fit_ode$species[is.na(fit_ode$id)]
  species2 = as.vector(fit_ode$species)
  species2[which(species2=="F1")]="IFI27"
  species2[which(species2=="F2")]="IFI6"
  species2[which(species2=="F3")]="IFI16"
  fit_ode$names = species2
  fit_ode = data.frame(params,fit_ode,rank2=z,weight = sub$summed_weight)
  fit_all_og = rbind(fit_all_og,fit_ode)
  
}

fit_all_og = as.data.frame(fit_all_og%>%dplyr::group_by(time)%>%dplyr::filter(type =="fit",species=="V")%>%dplyr::summarize(val = weighted.quantile(fit,weight,prob=0.5)))


#Simulate top 28 runs that make up 95% likelihood----

for(z in 1:number_sims){
  sub = subset(fit_all_sub,rank2==z)
  print(z)
  params = as.vector(sub[1,1:42])
  for(four in seq(5,65,20)){
    for (eight in seq(5,65,20)){
      for(tau1 in c(0.1,0.2,0.3,0.4,0.5)){
        for(tau2 in c(0.1,0.2,0.3,0.4,0.5)){
        params$T_0 = four
        params$E_0 = eight
        params$tauE = tau1
        params$tauT = tau2
        tdat <- mcmv_NA %>%
          trajectory(params = params,format="data.frame")
        tdat = melt(tdat[,c(2:11)],id="time")
        colnames(tdat)=c("time","species","fit")
        
        fit_ode = merge(data.frame(tdat,type = "fit"),data.frame(data,type = "dat"),all=TRUE)
        fit_ode$species[fit_ode$id%in%c(paste0("N_V_",1:6))]="V"
        fit_ode$species[fit_ode$id%in%c(paste0("N_E_",1:6))]="E"
        fit_ode$species[fit_ode$id%in%c(paste0("N_T_",1:6))]="T"
        fit_ode$species[fit_ode$id%in%c(paste0("N_F1_",1:6))]="F1"
        fit_ode$species[fit_ode$id%in%c(paste0("N_F2_",1:6))]="F2"
        fit_ode$species[fit_ode$id%in%c(paste0("N_F3_",1:6))]="F3"
        fit_ode$species[fit_ode$id%in%c(paste0("N_A_",1:6))]="A"
        fit_ode$id[is.na(fit_ode$id)]=fit_ode$species[is.na(fit_ode$id)]
        species2 = as.vector(fit_ode$species)
        species2[which(species2=="F1")]="IFI27"
        species2[which(species2=="F2")]="IFI6"
        species2[which(species2=="F3")]="IFI16"
        fit_ode$names = species2
        fit_ode = subset(fit_ode,type=="fit"&species=="V")
        #fit_ode = fit_ode[which(fit_ode$fit==max(fit_ode$fit)),]
        fit_ode = data.frame(params,fit_ode,rank2=z,weight = sub$summed_weight)
        if(z==1&four==5&eight==5&tau1==0.1&tau2==0.1)
          write.table(fit_ode,"predictions.csv",sep = ",",row.names = FALSE,col.names = TRUE,append = FALSE)
        else
          write.table(fit_ode,"predictions.csv",sep = ",",row.names = FALSE,col.names = FALSE,append = TRUE)
        
      }
    }
  }
  }
}

fit_all_new = read.csv("predictions.csv")


fit_all_new$E_0_labels = paste0("j=",fit_all_new$E_0)
fit_all_new$T_0_labels = paste0("i=",fit_all_new$T_0)
fit_all_new$tauE_labels = paste0("i=",fit_all_new$tauE*10)
fit_all_new$tauT_labels = paste0("i=",fit_all_new$tauT*10)
fit_all_new$cd8_effect = "no"
fit_all_new$cd8_effect[which(fit_all_new$m!=0)]="yes"
fit_all_new$E_0_labels = factor(fit_all_new$E_0_labels,levels = c("j=5","j=25","j=45","j=65"))
fit_all_new$T_0_labels = factor(fit_all_new$T_0_labels,levels = c("i=5","i=25","i=45","i=65"))

fit_summary = as.data.frame(fit_all_new%>%dplyr::filter(type=="fit"&species=="V")%>%
                              dplyr::group_by(tauE,tauT,tauE_labels,tauT_labels,T_0,E_0,T_0_labels,E_0_labels,time)%>%
                              dplyr::summarize(value =weighted.quantile(fit,weight,prob = 0.5),
                                               twentyfive =weighted.quantile(fit,weight,prob = 0.25),
                                               seventyfive =weighted.quantile(fit,weight,prob = 0.75)
                                               ))

AUC = fit_summary%>%dplyr::group_by(tauE,tauT,T_0,E_0)%>%dplyr::summarize(AUC = area_under_curve(time, (value-3000), method = "trapezoid"),
                                                                          AUC_log = area_under_curve(time, (log10(value-3000+1)), method = "trapezoid")
                                                                          )
AUC$T_0_label = paste0("T*=",AUC$T_0)
AUC$T_0_label = factor(AUC$T_0_label,levels = paste0("T*=",seq(5,65,10)))
AUC$E_0_label = paste0("E*=",AUC$E_0)
AUC$E_0_label = factor(AUC$E_0_label,levels = paste0("E*=",seq(5,65,10)))

maxes = fit_summary%>%dplyr::group_by(tauE,tauT,T_0,E_0)%>%dplyr::slice(which.max(value))
maxes = maxes[,c(1,2,5,6,9,10)]
maxes$T_0_label = paste0("T*=",maxes$T_0)
maxes$T_0_label = factor(maxes$T_0_label,levels = paste0("T*=",seq(5,65,10)))
maxes$E_0_label = paste0("E*=",maxes$E_0)
maxes$E_0_label = factor(maxes$E_0_label,levels = paste0("E*=",seq(5,65,10)))

#Create linear model to determine relationship between the maximum viral load and starting characteristics of virus-specific C and CD8 T cells
hi2 = glm(log10(value)~E_0+T_0+tauE+tauT,dat = maxes)

log_change = 1/hi2$coefficients[2:5]
log_change[3:4] = log_change[3:4]*10
log_change

(max_plot = ggplot(subset(maxes,T_0%in%c(5,25,45,65)&E_0%in%c(5,25,45,65)),aes(x = tauT*10,y = tauE*10,fill = log10(value)))+facet_grid(E_0_label~T_0_label)+geom_tile(colour="black")+scale_fill_gradient(low="blue", high="red",breaks = c(4.2,5,6,7,8,9),limits =c(4.2,9))+
  scale_x_continuous(breaks = c(1,2,3,4,5))+scale_y_continuous(breaks = c(1,2,3,4,5))+
    theme_bw()+
  labs(x = "days post reinfection when proliferation\nof virus-specific CD4+ T cells begins",y = "days post reinfection when proliferation\nof virus-specific CD8+ T cells begins",
       title = "A\nT*=starting number of virus-specific CD4+ T cells/ml from BAL\nE*=starting number of virus-specific CD8+ T cells/ml from BAL",fill = "log10(maximum viral load)")+
  theme(text = element_text(size=10),
        title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        plot.title = element_text(size = 10),
        strip.text.y = element_text(size = 10)))


(AUC_plot = ggplot(subset(AUC,T_0%in%c(5,25,45,65)&E_0%in%c(5,25,45,65)&tauT%in%c(0.1,0.3,0.5)&tauE%in%c(0.1,0.3,0.5)),aes(x = tauT*10,y = tauE*10,fill = log10(AUC)))+facet_grid(E_0_label~T_0_label)+geom_tile(colour="black")+scale_fill_gradient(low="blue", high="red")+
    scale_x_continuous(breaks = c(1,3,5))+scale_y_continuous(breaks = c(1,3,5))+
    theme_bw()+
    labs(x = "days post reinfection when proliferation\nof virus-specific CD4+ T cells begins",y = "days post reinfection when proliferation\nof virus-specific CD8+ T cells begins",
         title = "A\nT*=starting number of virus-specific CD4+ T cells/ml from BAL\nE*=starting number of virus-specific CD8+ T cells/ml from BAL",fill = "log10(Area under\nviral load curve)")+
    theme(text = element_text(size=10),
          title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "right",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          plot.title = element_text(size = 10),
          strip.text.y = element_text(size = 10)))

(AUC_log_plot = ggplot(subset(AUC,T_0%in%c(5,25,45,65)&E_0%in%c(5,25,45,65)),aes(x = tauT*10,y = tauE*10,fill = AUC_log))+facet_grid(E_0_label~T_0_label)+geom_tile(colour="black")+scale_fill_gradient(low="blue", high="red")+
    scale_x_continuous(breaks = c(1,2,3,4,5))+scale_y_continuous(breaks = c(1,2,3,4,5))+
    theme_bw()+
    labs(x = "days post reinfection when proliferation\nof virus-specific CD4+ T cells begins",y = "days post reinfection when proliferation\nof virus-specific CD8+ T cells begins",
         title = "A\nT*=starting number of virus-specific CD4+ T cells/ml from BAL\nE*=starting number of virus-specific CD8+ T cells/ml from BAL",fill = "Area under\nlog10(viral load curve)")+
    theme(text = element_text(size=10),
          title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "right",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          plot.title = element_text(size = 10),
          strip.text.y = element_text(size = 10)))

#showing threshold of detection for qPCR
data2 = data.frame(time = seq(-1,12,0.1),value=3000)

(equal_ratio= ggplot()+
    geom_line(data = data2,aes(x = time,y = value),colour = "red",size = 1)+
    geom_line(data = fit_all_og,aes(x = time,y = val),colour = "#9933FF",size = 1)+
    coord_cartesian(xlim = c(0, 10))+
    scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    geom_line(data = subset(fit_all_new,species=="V"&type=="fit"&tauE==tauT&E_0==T_0&tauT%in%c(0.1,0.3,0.5)&tauE%in%c(0.1,0.3,0.5)),aes(y = fit,x = time,colour =cd8_effect,group = paste0(tauE,tauT,T_0_labels,E_0_labels,rank2)))+
    geom_line(data = subset(fit_summary,tauE==tauT&E_0==T_0&tauT%in%c(0.1,0.3,0.5)&tauE%in%c(0.1,0.3,0.5)),aes(x = time,y = value,group = paste0(tauE,tauT,T_0_labels,E_0_labels)),size = 1,colour = "#0000FF")+
    scale_y_log10(breaks = c(10^3,10^5,10^7,10^9),label=scientific_10_1,limits = c(10^3,10^9))+
    scale_colour_manual("virus-specific CD8+\nT cells clear infection",values =c("#FF6699","#99CC00"))+
    theme_bw()+
    labs(x = "days post infection",y = "SARS-CoV-2 copies/ml from BAL",title = "B\ni=starting value of virus-specific CD4+ and CD8+ T cells/ml from BAL\nj=days post reinfection when virus-specific CD4+ and CD8+ T cell proliferation begins")+
    facet_grid(paste0("j=",tauT*10)~T_0_labels)+
    theme(text = element_text(size=10),
          legend.text = element_text(size = 10),
          legend.position = "right",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          plot.title = element_text(size = 10),
          strip.text.y = element_text(size = 10)))

tiff("Fig7.tiff", height = 22.23, width =19.05, units = 'cm',
     compression = "lzw", res = 600)
multiplot(AUC_plot,equal_ratio)

dev.off()


