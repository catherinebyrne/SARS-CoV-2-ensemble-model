library(pomp)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(plyr)
library(viridis)
library(modi)
library(forcats)
library(cowplot)
library(spatstat.geom)
library(scales)
library(extrafont)
font_import()
fonts()


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

#readin in data and model
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

#select models that make 95% confidence set
number_sims = which(fit_all_sub$summed_weight>=0.95)[1]

top_models = fit_all_sub[1:number_sims,]
write.table(top_models,"top_models.csv",sep = ",",col.names = TRUE,row.names = FALSE)

#Begin making plots------------------------


#create heatmaps of runs
param_list = c("m","n","yV1","b1","yV2","b2","yV3","b3","yA","S")
bin = NULL
for (i in 1:length(param_list)){
  vec = rep(NA,nrow(fit_all_sub))
  chosen = grep(param_list[i],fit_all_sub$version_dat,fixed = TRUE)
  vec[chosen]=fit_all_sub$diff[chosen]
  bin = cbind(bin,vec)
}

colnames(bin) = param_list
fit_all_sub = cbind(fit_all_sub,bin)
heat = fit_all_sub[,c(51,56:65)]
heat$rank2 = as.character(heat$rank2)
heat$rank2 = factor(heat$rank2,levels = as.character(1:160))
heat = melt(heat, id.vars=c("rank2"))
heat$variable = factor(heat$variable,levels = parameter_importance$param)
colnames(heat) = revalue(colnames(heat),c("variable"="parameter"))
param_new = heat$parameter
param_new = revalue(param_new,c(
  "b1"="b [F1]","b2"="b [F2]","b3"="b [F3]",
  "yA"="y [A]",
  "yV1"="y [F1]","yV2"="y [F2]","yV3"="y [F3]"))
heat = data.frame(heat,param_new = param_new)

#Heatmap of only top-ranked models
heat_final = ggplot(subset(heat,rank2%in%c(1:number_sims)&parameter%in%c("S","m","n",paste0("yV",c(1,3)),paste0("b",c(1,3)))),aes(x = param_new,y = rank2,fill = value))+geom_tile(colour = "black")+theme_bw()+
  scale_fill_viridis( na.value="white")+ labs(fill=bquote(Delta~AIC)) +scale_x_discrete(labels = scales::parse_format())+labs(x = "parameter",y = "model rank",title = "A")+
  theme(text = element_text(size=10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10))

#heatmap of top 40 models
heat_final_all1 = ggplot(subset(heat,rank2%in%c(1:40)&param_new!="y [A]"),aes(x = param_new,y = rank2,fill = value))+geom_tile(colour = "black")+theme_bw()+
  scale_fill_viridis( na.value="white",begin =0,end = 0.6,breaks=seq(0,20,5),limits = c(0,22))+ labs(fill=bquote(Delta~AIC)) +scale_x_discrete(labels = scales::parse_format())+labs(x = "parameter",y = "model rank")+
  theme(text = element_text(size=10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10))

#heatmap of models ranked 41-70
heat_final_all2 = ggplot(subset(heat,rank2%in%c(41:70)&param_new!="y [A]"),aes(x = param_new,y = rank2,fill = value))+geom_tile(colour = "black")+theme_bw()+
  scale_fill_viridis( na.value="white",begin =0.6,end = 1,breaks = seq(20,100,20),limits = c(20,107))+ labs(fill=bquote(Delta~AIC)) +scale_x_discrete(labels = scales::parse_format())+labs(x = "parameter",y = "model rank")+
  theme(text = element_text(size=10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10))

#heatmap of models ranked 71-80
heat_final_all3 = ggplot(subset(heat,rank2%in%c(71:80)&param_new!="y [A]"),aes(x = param_new,y = rank2,fill = value))+geom_tile(colour = "black")+theme_bw()+
  scale_fill_gradient(low = "yellow", high = "red", na.value = "white",breaks = c(100,1000,2000,3000,4000),limits = c(100,4100))+
  labs(fill=bquote(Delta~AIC)) +scale_x_discrete(labels = scales::parse_format())+labs(x = "parameter",y = "model rank")+
  theme(text = element_text(size=10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10))

#Plot of range of AIC values for top-ranked models
one = ggplot(subset(fit_all_trimmed,rank_min%in%as.character(1:number_sims)),aes(x = rank_min,y = AIC-min(fit_all_trimmed$AIC)))+geom_point(aes(colour = as.character(i)),alpha = 0.5)+geom_boxplot(alpha = 0)+#coord_cartesian(ylim = c(1525, 1556))+
    theme_bw()+labs(colour = "Iteration",x = "Model Rank",y = bquote(Delta~AIC),title = "B")+
    #scale_y_continuous(breaks = seq(1525, 1555,5))+
    scale_y_continuous(breaks = c(0,5,10,15,20,25))+
    theme(text = element_text(size=10),
          legend.position = "right",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))

layout = matrix(c(1,3,
                  1,2,
                  1,2,
                  1,2,
                  4,4,
                  4,4),
                nrow = 6,byrow = TRUE)

tiff("Supp2.tiff", height = 21.23, width =19.5, units = 'cm',
     compression = "lzw", res = 600)
multiplot(heat_final_all1,heat_final_all2,heat_final_all3,one,layout = layout)

dev.off()

#Simulate top-ranked models that make up 95% likelihood----
fit_all_plot = NULL

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
  sub = data.frame(sub,fit_ode)
  fit_all_plot = rbind(fit_all_plot,sub)
}

#plot weighted median behaviour------

fit_all_plot2 = fit_all_plot

weighted_fit= as.data.frame(fit_all_plot2%>%dplyr::filter(type=="fit")%>%dplyr::group_by(time,names,species,type)%>%dplyr::summarize(fit_weighted = sum(fit*weight)/sum(weight),
                                                                                                                                     se = sum(weight/sum(weight)*abs(fit-fit_weighted)),
                                                                                                                                     upper = fit_weighted+se,
                                                                                                                                     lower = fit_weighted-se,
                                                                                                                                     median_weighted = weighted.quantile(fit,weight,prob = 0.5),
                                                                                                                                     twentyfive = weighted.quantile(fit,weight,prob = 0.25),
                                                                                                                                     seventyfive = weighted.quantile(fit,weight,prob = 0.75)))

#plot median
names = weighted_fit$species
names2 = revalue(names,c("E"="CD8+","T"="CD4+","V"="Genomic RNA","A"="IgG","F1"="IFI27","F2"="IFI6","F3"="IFI16"))

weighted_fit= data.frame(weighted_fit,names2 = names2)
weighted_fit$names2 = factor(weighted_fit$names2,levels = c("Genomic RNA","IgG","CD4+","CD8+","IFI27","IFI6","IFI16","S","I"))

names = fit_all_plot$species
names2 = revalue(names,c("E"="CD8+","T"="CD4+","V"="Genomic RNA","A"="IgG","F1"="IFI27","F2"="IFI6","F3"="IFI16"))

fit_all_plot2= data.frame(fit_all_plot2,names2 = names2)
fit_all_plot2$names2 = factor(fit_all_plot2$names2,levels = c("Genomic RNA","IgG","CD4+","CD8+","IFI27","IFI6","IFI16","S","I"))



colours3 = c("#CC79A7","#9933FF","#0000FF","#33CCFF","#009E73","#99CC00","#E69F00","#FF33CC", "red")

colour_dat = data.frame(names = levels(weighted_fit$names2),colours = colours3)


B2 = ggplot()+
    geom_hline(aes(yintercept = 3000),colour = "red")+
    geom_line(data = subset(fit_all_plot2,type == "fit"&species%in%c("V"))[,c("time","fit","rank2","names2")],aes(x = time,y = fit,group = rank2),alpha = 0.5,colour = "black")+
    geom_line(data = subset(weighted_fit,type == "fit"&species%in%c("V")), aes(x = time,y = median_weighted,colour = names2),size = 1,alpha =1)+
    
    geom_point(data =unique(subset(fit_all_plot2,type == "dat"&species%in%c("V"))[,c("time","fit","id","names2")]), aes(x = time,y = fit),alpha = 0.5)+
    scale_colour_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("V"))$names2))$colours))+
    scale_fill_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("V"))$names2))$colours))+
    theme_bw()+scale_y_log10(breaks = c(10^3,10^4,10^5,10^6,10^7,10^8,10^9),label=scientific_10_1,limits = c(10^3,10^9))+
    labs(x = "days post infection",y = "SARS-CoV-2 copies/ml from BAL",title = "B")+
    facet_wrap(~names2)+
    theme(legend.position = "none")+
    theme(text = element_text(size=10),
          legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))+scale_x_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))


C2 = ggplot()+
    geom_line(data = subset(fit_all_plot2,type == "fit"&species%in%c("A"))[,c("time","fit","rank2","names2")],aes(x = time,y = fit,group = rank2),alpha = 0.5,colour = "black")+
    geom_line(data = subset(weighted_fit,type == "fit"&species%in%c("A")), aes(x = time,y = median_weighted,colour = names2),size = 1,alpha =1)+
    geom_point(data =unique(subset(fit_all_plot2,type == "dat"&species%in%c("A"))[,c("time","fit","id","names2")]), aes(x = time,y = fit),alpha = 0.5)+
    scale_colour_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("A"))$names2))$colours))+
    scale_fill_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("A"))$names2))$colours))+
    theme_bw()+
    labs(x = "days post infection",y ="AUC for anti-spike titration curves from BAL",title ="C")+
    facet_wrap(~names2)+
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0,8*10^5),label=scientific_10_3)+
    theme(text = element_text(size=10),
          legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))+scale_x_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))


D2 = ggplot()+
    geom_line(data = subset(fit_all_plot2,type == "fit"&species%in%c("E","T"))[,c("time","fit","rank2","names2")],aes(x = time,y = fit,group = rank2),alpha = 0.5,colour = "black")+
    geom_line(data = subset(weighted_fit,type == "fit"&species%in%c("E","T")), aes(x = time,y = median_weighted,colour = names2),size = 1,alpha =1)+
    geom_point(data =unique(subset(fit_all_plot2,type == "dat"&species%in%c("E","T"))[,c("time","fit","id","names2")]), aes(x = time,y = fit),alpha = 0.5)+
    theme_bw()+
    scale_colour_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("E","T"))$names2))$colours))+
    scale_fill_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("E","T"))$names2))$colours))+
    labs(x = "days post infection",y =  "virus-specific\nCD4+ or CD8+ cells/ml from BAL",title = "D")+
    facet_wrap(~names2)+
    theme(legend.position = "none")+
    theme(text = element_text(size=10),
          legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))+scale_x_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))


E2 = ggplot()+
    geom_line(data = subset(fit_all_plot2,type == "fit"&species%in%c("F1","F2","F3"))[,c("time","fit","rank2","names2")],aes(x = time,y = fit,group = rank2),alpha = 0.5,colour = "black")+
    geom_line(data = subset(weighted_fit,type == "fit"&species%in%c("F1","F2","F3")), aes(x = time,y = median_weighted,colour = names2),size = 1,alpha =1)+
    geom_point(data =unique(subset(fit_all_plot2,type == "dat"&species%in%c("F1","F2","F3"))[,c("time","fit","id","names2")]), aes(x = time,y = fit),alpha = 0.5)+
    scale_colour_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("F1","F2","F3"))$names2))$colours))+
    scale_fill_manual(values = as.vector(subset(colour_dat,names%in%unique(subset(weighted_fit,species%in%c("F1","F2","F3"))$names2))$colours))+
    theme_bw()+
    labs(x = "days post infection",y = "RNA copies/cell from BAL",title = "E")+
    facet_wrap(~names2)+
    theme(legend.position = "none")+
    theme(text = element_text(size=10),
          legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))+scale_x_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))

layout = matrix(c(1,1,2,3,
                  1,1,4,4,
                  5,5,5,5),
                nrow = 3,byrow = TRUE)

tiff("Fig3.tiff", height = 22.23, width =19.05, units = 'cm',
     compression = "lzw", res = 600)
multiplot(heat_final,B2,C2,D2,E2,layout = layout)

dev.off()

#Calculate the parameter importance from all top runs----
parameter_importance = NULL
for(i in param_list){
  val = sum(fit_all_sub[c(1:number_sims),]$weight[grepl(paste0(i,"."), fit_all_sub[c(1:number_sims),]$version_dat, fixed=TRUE)])
  parameter_importance = rbind(parameter_importance,data.frame(param = i,importance = val))
}

(parameter_importance = parameter_importance%>% arrange(desc(importance)))

#plot fit parameter values-------------
param_vals = unique(subset(fit_all_sub,rank2<=number_sims)[,c(10:42,47,50,51,53)])
param_vals = melt(param_vals,id.vars = c("version_dat","weight","rank2","AIC"))
param_vals = subset(param_vals,variable!="o")
param_vals = subset(param_vals,value!=1.645000e+09)
param_vals = subset(param_vals,value!=2.900000e-02)

new = rep("name",nrow(param_vals))
new[param_vals$variable%in%c("tauT","tauE","tauA","tauF1","tauF2","tauF3")]="timing of proliferation"
new[param_vals$variable%in%c("alphaA")]="rate of anti-spike IgG\nproliferation"
new[param_vals$variable%in%c("yV1","yV2","yV3")]="damping of viral production\nby innate response"
new[param_vals$variable%in%c("b1","b2","b3")]="infected cell clearance\nby innate response"
new[param_vals$variable%in%c("m","n")]="infected cell clearance\nby T cells"
new[param_vals$variable%in%c("beta1","beta2","beta3")]="rate of innate response\nproliferation"
new[param_vals$variable%in%c("alphaE","alphaT")]="rate of virus-specific\nT cell proliferation"
new[param_vals$variable%in%c("delta","c","ds")]="clearance rates of virus,\ninfected cells, and susceptible cells"
new[param_vals$variable%in%c("z1","z2","z3")]="clearance rates of IFNs"
new[param_vals$variable%in%c("dE","dT")]="clearance rates of\nvirus-specific T cells"
new[param_vals$variable%in%c("p")] = "viral production\nrate"
new[param_vals$variable%in%c("eta")] = "infection rate"
new[param_vals$variable%in%c("yA")] = "damping\ninfection rate"
new[param_vals$variable%in%c("S_0")] = "number of\ntarget cells"

new2 = as.data.frame(param_vals%>%group_by(variable)%>%filter(value!=0)%>%tally())
param_vals = data.frame(param_vals,label = new)
param_vals = merge(param_vals,new2)

param_vals$n = as.character(param_vals$n)
param_vals$n = revalue(param_vals$n,c("1" = "1/12","6" = "6/12", "7"="7/12","8"="8/12","12"="12/12"))
param_vals$n = factor(param_vals$n,levels = c("1/12","6/12","7/12","8/12","12/12"))

param_vals[param_vals$label=="timing of proliferation",]$value = 10*(param_vals[param_vals$label=="timing of proliferation",]$value)

variable2 = param_vals$variable
variable2 = revalue(variable2,c("tauA"="tau [A]","tauT"="tau [T]","tauE"="tau [E]","tauF1"="tau [F1]","tauF2"="tau [F2]","tauF3"="tau [F3]",
                                "beta1"="alpha [F1]","beta2"="alpha [F2]","beta3"="alpha [F3]","alphaE"="alpha [E]","alphaT"="alpha [T]","alphaA"="alpha [A]",
                                "b1"="b [F1]","b2"="b [F2]","b3"="b [F3]",
                                "yA"="y [A]","S_0"="S [0]",
                                "yV1"="y [F1]","yV2"="y [F2]","yV3"="y [F3]",
                                "c"="d [V]","delta"="d [I]","z1"="d [F1]","z2"="d [F2]","z3"="d [F3]","dE"="d [E]","dT"="d [T]","ds"="d [S]"))
param_vals = data.frame(param_vals,variable2 = variable2)
param_vals$variable2 = factor(param_vals$variable2,levels = c("tau [F1]","tau [F2]","tau [F3]", "tau [T]","tau [E]","tau [A]","p" ,"y [F1]","y [F2]","y [F3]","d [V]" ,"eta","d [I]","alpha [F1]","d [F1]" ,    "alpha [F2]", "d [F2]"  ,   "alpha [F3]", "d [F3]",  "alpha [T]",   "alpha [E]" ,  
                                                              "d [E]" ,"d [T]","n","m","b [F1]","b [F2]","b [F3]","alpha [A]","y [A]",     
                                                              "S [0]","d [S]","o","inf" 
)) 
param_vals$grouping = NA
param_vals$grouping[param_vals$variable%in%c("b1","beta1","tauF1","yV1","z1")]="IFI27"
param_vals$grouping[param_vals$variable%in%c("b2","beta2","tauF2","yV2","z2")]="IFI6"
param_vals$grouping[param_vals$variable%in%c("b3","beta3","tauF3","yV3","z3")]="IFI16"
param_vals$grouping[param_vals$variable%in%c("alphaT","dT","n","tauT")]="virus-specific CD4+ T cells"
param_vals$grouping[param_vals$variable%in%c("alphaE","dE","m","tauE")]="virus-specific CD8+ T cells"
param_vals$grouping[param_vals$variable%in%c("alphaA","tauA")]="anti-spike IgG"
param_vals$grouping = factor(param_vals$grouping,levels =c("IFI27","IFI6","IFI16","virus-specific CD4+ T cells","virus-specific CD8+ T cells","anti-spike IgG"))
param_vals = arrange(param_vals,grouping)


weighted_param_vals = as.data.frame(param_vals%>%dplyr::group_by(variable2)%>%dplyr::summarize(value =weighted.quantile(value,weight,prob = 0.5)))
median_param_vals = as.data.frame(param_vals%>%dplyr::group_by(variable2)%>%dplyr::filter(value!=0)%>%dplyr::summarize(median =median(value)))

colours5 = c("grey20","grey40","grey60","grey75","grey90")
colours6 = c("red","#E69F00","#99CC00","#33CCFF","#9933FF","#FF33CC")

#Use these plots just to get legends
start = ggplot(subset(param_vals,!is.na(grouping)),aes(x = variable2,y = value,colour = grouping))+geom_boxplot(outlier.shape = NA,alpha = 0)+geom_jitter(alpha = 0.2,width = 0.25,height = 0)+
  scale_y_log10()+facet_wrap(~label,scales = "free")+theme_bw()+theme(legend.position  = "right")+labs(y = "Value",x = "Parameter")+scale_fill_manual("Appearance\nFrequency",values = colours5)+scale_colour_manual("Compartment",values = colours6)+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size =10))

start2 = ggplot(subset(param_vals,!is.na(grouping)),aes(x = variable2,y = max(value),fill = n),colour = "black")+geom_bar(stat = "identity",alpha = 0.75)+
  facet_wrap(~label,scales = "free")+theme_bw()+theme(legend.position  = "right")+labs(y = "Value",x = "Parameter")+scale_fill_manual("Appearance Frequency",values = colours5)+scale_colour_manual("Compartment",values = colours6)+
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size =10))

legend <- cowplot::get_legend(start)
legend2 <- cowplot::get_legend(start2)

tiff("Fig4_legend1.tiff", height = 22.23, width =19.05, units = 'cm',
     compression = "lzw", res = 600)
grid.newpage()
grid.draw(legend)
dev.off()

tiff("Fig4_legend2.tiff", height = 22.23, width =19.05, units = 'cm',
     compression = "lzw", res = 600)
grid.newpage()
grid.draw(legend2)
dev.off()

#Create function to create each parameter grouping plot
plot_fun  = function(label_chosen,ylab,plot_name,log,chosen_breaks,colour_chosen){
  dat_plot = subset(param_vals,label==label_chosen&value!=0)
  chosen = colour_chosen
  
  chosen2 = data.frame(names = c("IFI27","IFI6","IFI16","virus-specific CD4+ T cells","virus-specific CD8+ T cells","anti-spike IgG"),colours = colours6)
  chosen2 = subset(chosen2,names%in%unique(dat_plot$grouping))
  chosen2 = as.vector(chosen2$colours)
  
  weighted_dat = subset(weighted_param_vals,variable2%in%unique(dat_plot$variable2))
  
  if(length(unique(dat_plot$variable2))==1){
    shading <- data.frame(min =0.5,
                          max = 1.5,
                          col = seq(1,length(unique(dat_plot$n)),1))
  }else{
    min = 0.5
    max = 1.5
    min_col = min
    max_col = max
    for(i in 2:length(unique(dat_plot$variable2))){
      if (unique(subset(dat_plot,variable2==unique(dat_plot$variable2)[i])$n)==unique(subset(dat_plot,variable2==unique(dat_plot$variable2)[i-1])$n)){
        max_col[length(max_col)] = max_col[length(max_col)]+1
      }else{
        min_col = c(min_col,max_col[length(max_col)])
        max_col = c(max_col,max_col[length(max_col)]+1)
      }
    }
    shading = data.frame(min = min_col,
                         max = max_col,
                         col = seq(1,length(unique(dat_plot$n)),1))
  }
  
  X = ggplot()+
    geom_rect(data = shading,
              aes(xmin = min, xmax = max, ymin = min(min(weighted_dat$value),min(dat_plot$value))*0.75, ymax = Inf,
                  fill = forcats::fct_inorder(factor(col))),alpha = 0.75, inherit.aes = F)+
    scale_fill_manual(values = chosen)+
    geom_boxplot(data = dat_plot,aes(x = variable2,y = value,colour = forcats::fct_inorder(grouping)),outlier.shape = NA,alpha = 0)+
    geom_jitter(data = dat_plot,aes(x = variable2,y = value,colour = forcats::fct_inorder(grouping)),alpha = 0.2,width = 0.25,height = 0)+
    facet_wrap(~label,scales = "free")+
    geom_point(data = weighted_dat,aes(x = variable2,y = value),colour = "black")+
    theme_bw()+
    scale_colour_manual(values = chosen2)+labs(y = ylab)+
    theme(text = element_text(size=10),
          legend.position = "none",
          axis.title.x = element_blank(),axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))+
    scale_x_discrete(labels = scales::parse_format()) +labs(title = plot_name)
  if(log=="yes"&chosen_breaks=="NA"){
    X = X+scale_y_log10(ylab,
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x)))
  }else if(log=="mod"&chosen_breaks== "NA"){
    X = X+scale_y_continuous(label=scientific_10_0)
  }else if(log=="yes"&chosen_breaks!="NA"){
    X = X+scale_y_log10(label = scientific_10,breaks = c(0.01,0.015,0.02,0.025),limits = c(0.01,0.026))
  }else if(log=="mod"&chosen_breaks!= "NA"){
    X = X+scale_y_continuous(label=scientific_10_0,breaks = c(0.01,0.015,0.02,0.025),limits = c(0.01,0.026))
  }
  
  assign(plot_name,X,envir = .GlobalEnv)
  return(plot_name)
}

plot_fun("timing of proliferation","days post infection","A","no","NA","grey90")
plot_fun("rate of innate response\nproliferation","per infected cell\nper day","B","yes","NA","grey90")
plot_fun("infected cell clearance\nby innate response","per IFN-stimulating gene copy/cell\nper infected cell\nper day","C","no","NA",c("grey60","grey20"))
plot_fun("damping of viral production\nby innate response","per IFN-stimulating gene copy/cell","D","no","NA",c("grey60","grey20"))
plot_fun("rate of virus-specific\nT cell proliferation","per infected cell\nper day","E","yes","NA","grey90")
plot_fun("infected cell clearance\nby T cells","per virus-specific CD4+ or CD8+ cell/ml\nper infected cell\nper day","F","no","NA",c("grey90","grey40"))
plot_fun("rate of anti-spike IgG\nproliferation","per infected cell\nper day","G","yes","NA","grey90")


layout = matrix(c(1,1,NA,
                  2,3,4,
                  5,6,7
), nrow=3, byrow=TRUE)

tiff("Fig4.tiff", height = 22.23, width =19.05, units = 'cm',
     compression = "lzw", res = 600)
multiplot(A,B,C,D,E,F,G,layout =layout)
dev.off()

plot_fun2  = function(label_chosen,ylab,plot_name,log,chosen_breaks,colour_chosen){
  dat_plot = subset(param_vals,label==label_chosen&value!=0)
  chosen = colour_chosen
  
  weighted_dat = subset(weighted_param_vals,variable2%in%unique(dat_plot$variable2))
  
  if(length(unique(dat_plot$variable2))==1){
    shading <- data.frame(min =0.5,
                          max = 1.5,
                          col = seq(1,length(unique(dat_plot$n)),1))
  }else{
    min = 0.5
    max = 1.5
    min_col = min
    max_col = max
    for(i in 2:length(unique(dat_plot$variable2))){
      if (unique(subset(dat_plot,variable2==unique(dat_plot$variable2)[i])$n)==unique(subset(dat_plot,variable2==unique(dat_plot$variable2)[i-1])$n)){
        max_col[length(max_col)] = max_col[length(max_col)]+1
      }else{
        min_col = c(min_col,max_col[length(max_col)])
        max_col = c(max_col,max_col[length(max_col)]+1)
      }
    }
    shading = data.frame(min = min_col,
                         max = max_col,
                         col = seq(1,length(unique(dat_plot$n)),1))
  }

  X = ggplot()+
    geom_rect(data = shading,
              aes(xmin = min, xmax = max, ymin = min(min(weighted_dat$value),min(dat_plot$value))*0.8, ymax = Inf,
                  fill = forcats::fct_inorder(factor(col))),alpha = 0.75, inherit.aes = F)+
    scale_fill_manual(values = chosen)+
    geom_boxplot(data = dat_plot,aes(x = variable2,y = value),outlier.shape = NA,alpha = 0,colour = "blue")+
    geom_jitter(data = dat_plot,aes(x = variable2,y = value),alpha = 0.2,width = 0.25,height = 0,colour = "blue")+
    facet_wrap(~label,scales = "free")+
    geom_point(data = weighted_dat,aes(x = variable2,y = value),colour = "black")+
    theme_bw()+labs(y = ylab)+
    theme(text = element_text(size=10),
          legend.position = "none",
          axis.title.x = element_blank(),axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))+
    scale_x_discrete(labels = scales::parse_format()) +labs(title = plot_name)
  if(log=="yes"&chosen_breaks=="NA"){
    X = X+scale_y_log10(ylab,
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x)))
  }else if(log=="mod"&chosen_breaks== "NA"){
    X = X+scale_y_continuous(label=scientific_10_0)
  }else if(log=="yes"&chosen_breaks!="NA"){
    X = X+scale_y_log10(label = scientific_10,breaks = c(0.01,0.015,0.02,0.025),limits = c(0.01,0.026))
  }else if(log=="mod"&chosen_breaks!= "NA"){
    X = X+scale_y_continuous(label=scientific_10_0,breaks = c(0.01,0.015,0.02,0.025),limits = c(0.01,0.026))
  }
  
  assign(plot_name,X,envir = .GlobalEnv)
  return(plot_name)
}

plot_fun2("clearance rates of virus,\ninfected cells, and susceptible cells","per virus or cell\nper day","A","yes","NA",c("grey75","grey90"))
plot_fun2("clearance rates of IFNs","per gene copy/cell\nper day","B","no","NA","grey90")
plot_fun2("clearance rates of\nvirus-specific T cells","per virus-specific CD4+ or CD8+ cell/ml\nper day","C","yes","NA","grey90")
plot_fun2("viral production\nrate","per-infected-cell\nper-day","D","mod","NA","grey90")
plot_fun2("infection rate","per-gRNA copy/ml\nper-day","E","mod","YES","grey90")

layout = matrix(c(1,1,1,1,4,4,
                  2,2,2,2,5,5,
                  3,3,3,NA,NA,NA
), nrow=3, byrow=TRUE)

tiff("Supp3.tiff", height = 22.23, width =19.05, units = 'cm',
     compression = "lzw", res = 600)
multiplot(A,B,C,D,E,layout =layout)
dev.off()


#look at correlation between parameter values of top-ranked models ----
corr_dat = subset(fit_all_sub,rank2<=number_sims)[,c(10:42)]
corr_dat = corr_dat[,which(!colnames(corr_dat)%in%c("S_0","b2","b3","yV2","yV3"))]
corr_dat$ds[corr_dat$o==0]=NA
corr_dat[corr_dat==0]=NA


colnames(corr_dat) = revalue(colnames(corr_dat),c("tauA"="tau [A]","tauT"="tau [T]","tauE"="tau [E]","tauF1"="tau [F1]","tauF2"="tau [F2]","tauF3"="tau [F3]",
                                                  "beta1"="alpha [F1]","beta2"="alpha [F2]","beta3"="alpha [F3]","alphaE"="alpha [E]","alphaT"="alpha [T]","alphaA"="alpha [A]",
                                                  "b1"="b [F1]","b2"="b [F2]","b3"="b [F3]",
                                                  "yA"="y [A]","S_0"="S [0]",
                                                  "yV1"="y [F1]","yV2"="y [F2]","yV3"="y [F3]",
                                                  "c"="d [V]","delta"="d [I]","z1"="d [F1]","z2"="d [F2]","z3"="d [F3]","dE"="d [E]","dT"="d [T]","ds"="d [S]"))

levels = c("tau [F1]","tau [F2]","tau [F3]","tau [T]","tau [E]","tau [A]",
           "alpha [F1]","alpha [F2]","alpha [F3]","alpha [T]","alpha [E]","alpha [A]",
           "b [F1]","n","m","y [F1]",
           "d [F1]","d [F2]","d [F3]","d [T]","d [E]","d [V]","d [I]","d [S]","eta","p")

options(show.error.messages = FALSE)
corr_vals = NULL
for (i in 1:(length(levels)-1)){
  p1 = levels[i]
  for(j in (i+1):length(levels)){
    p2 = levels[j]
    test = try(cor.test(corr_dat[,which(colnames(corr_dat)==levels[i])],corr_dat[,which(colnames(corr_dat)==levels[j])]))
    if("try-error" %in% class(test)){
      val = NA
      pval = NA
    }else{
      val = as.numeric(test$estimate)
      pval =as.numeric(test$p.val)
    }
    corr_vals = rbind(corr_vals,data.frame(p1 = p1,p2 = p2,val = val,pval = pval))
  }
}

corr_vals$na = "FALSE"
corr_vals$na[corr_vals$pval<0.05]=TRUE


corr_plot = ggplot(corr_vals,aes(x  = p2,y = p1,fill = val))+geom_tile(colour="white", size=0.25)+scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                                                                                                                       high="red",name = "Correlation",breaks = c(-1,-0.5,0,0.5,1),limits = c(-1,1))+
  scale_x_discrete(labels = scales::parse_format(),guide = guide_axis(angle = 0))+
  scale_y_discrete(labels = scales::parse_format()) +labs(x = "",y = "")+
  theme(text = element_text(size=10),
        legend.position = "right",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        title = element_text(size = 10))+
  geom_point(aes(size=ifelse(na, "dot", "no_dot")),shape = 4) +
  scale_size_manual(values=c(dot=3, no_dot=NA), guide="none")

tiff("Supp4.tiff", height = 17, width =19.5, units = 'cm',
     compression = "lzw", res = 600)
corr_plot

dev.off()


#Plot different immune components impact on clearing infection-----
mE = subset(fit_all_plot,species=="E"&type=="fit")
mE$new = mE$fit*mE$m

nT = subset(fit_all_plot,species=="T"&type=="fit")
nT$new= nT$fit*nT$n

bF1 = subset(fit_all_plot,species%in%c("F1")&type=="fit")
bF1$bF1 = bF1$fit*bF1$b1

bF2 = subset(fit_all_plot,species%in%c("F2")&type=="fit")
bF2$bF2 = bF2$fit*bF2$b2

bF3 = subset(fit_all_plot,species%in%c("F3")&type=="fit")
bF3$bF3 = bF3$fit*bF3$b3

bF = merge(bF3[,c("time","weight","rank2","bF3")],bF2[,c("time","weight","rank2","bF2")])
bF = merge(bF,bF1[,c("time","weight","rank2","bF1")])
bF$new = bF$bF1+bF$bF2+bF$bF3

I = subset(fit_all_plot,species=="I"&type=="fit")[,c("time","fit","rank2")]
I = arrange(I,time,rank2)

immune = rbind(data.frame(bF[,c(1:3,7)],name = "bF"),
               data.frame(nT[,c("time","weight","rank2","new")],name = "nT"),
               data.frame(mE[,c("time","weight","rank2","new")],name = "mE"))

immune_mult = immune
immune_mult = arrange(immune_mult,time,rank2)
immune_mult = data.frame(immune_mult,I = rep(I$fit,each = 3))
immune_mult$new_mult = immune_mult$new*immune_mult$I


weighted_immune =  as.data.frame(immune%>%dplyr::group_by(time,name)%>%dplyr::summarize(fit_weighted = sum(new*weight)/sum(weight),
                                                                                        se = sum(weight/sum(weight)*abs(new-fit_weighted)),
                                                                                        upper = fit_weighted+se,
                                                                                        lower = fit_weighted-se,
                                                                                        median_weighted = weighted.quantile(new,weight,probs = 0.5),
                                                                                        twentyfive = weighted.quantile(new,weight,probs = 0.25),
                                                                                        seventyfive = weighted.quantile(new,weight,probs = 0.75),
                                                                                        zero = weighted.quantile(new,weight,probs = 0),
                                                                                        hundred = weighted.quantile(new,weight,probs = 1)
))


weighted_immune_mult =  as.data.frame(immune_mult%>%dplyr::group_by(time,name)%>%dplyr::summarize(fit_weighted = sum(new_mult*weight)/sum(weight),
                                                                                                  se = sum(weight/sum(weight)*abs(new_mult-fit_weighted)),
                                                                                                  upper = fit_weighted+se,
                                                                                                  lower = fit_weighted-se,
                                                                                                  median_weighted = weighted.quantile(new_mult,weight,probs = 0.5),
                                                                                                  twentyfive = weighted.quantile(new_mult,weight,probs = 0.25),
                                                                                                  seventyfive = weighted.quantile(new_mult,weight,probs = 0.75)
))


#plotting weighted median

A3 = ggplot()+
    geom_line(data = weighted_immune,aes(x = time,y = median_weighted,colour = name),size = 1,alpha =1)+
    geom_ribbon(data = weighted_immune,aes(x = time,ymin = twentyfive,ymax = seventyfive,fill = name),alpha = 0.5)+
    theme_bw()+
    #scale_y_log10(label = scientific_10_1)+
    labs(x = "day",y = "measure",colour = "Term",fill = "Term",title = "A")+
    scale_y_continuous(limits = c(0,25))+scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    theme(text = element_text(size=12),
          legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          title = element_text(size = 12))+
    labs(x ="days since infection",y = "per-cell-per-day clearance rates")


B3 = ggplot()+
    geom_line(data = weighted_immune_mult,aes(x = time,y = median_weighted+1,colour = name),size = 1,alpha =1)+
    geom_ribbon(data = weighted_immune_mult,aes(x = time,ymin = twentyfive+1,ymax = seventyfive+1,fill = name),alpha = 0.5)+
    theme_bw()+
    labs(x = "day",y = "measure",fill = "Mediated by:",colour = "Mediated by:",title ="B")+
    theme(text = element_text(size=12),
          legend.position = "right",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          title = element_text(size = 12))+
    labs(x ="days since infection",y = "per-day clearance rates+1")+scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    scale_y_log10(label=scientific_10_2,breaks= c(10^0,10^2,10^4,10^6,10^8),limits = c(1,10^8))+
    
    scale_color_discrete(labels = c("Innate immune response",
                                    "virus-specific CD4+ T cells",
                                    "virus-specific CD8+ T cells"))+
    scale_fill_discrete(labels = c("Innate immune response",
                                   "virus-specific CD4+ T cells",
                                   "virus-specific CD8+ T cells"))

legend <- cowplot::get_legend(B3)
grid.newpage()
grid.draw(legend)

B3 = B3+theme(legend.position = "none")

#plot impact of interferon on viral production

yF1 = subset(fit_all_plot,species%in%c("F1")&type=="fit")
yF1$yF1 = exp(-yF1$fit*yF1$yV1)

yF2 = subset(fit_all_plot,species%in%c("F2")&type=="fit")
yF2$yF2 = exp(-yF2$fit*yF2$yV2)

yF3 = subset(fit_all_plot,species%in%c("F3")&type=="fit")
yF3$yF3 = exp(-yF3$fit*yF3$yV3)
yF = merge(yF3[,c("time","p","weight","rank2","yF3")],yF2[,c("time","p","weight","rank2","yF2")])
yF = merge(yF,yF1[,c("time","p","weight","rank2","yF1")])
new =c(yF$yF1*yF$yF2*yF$yF3*yF$p,yF$yF1*yF$yF2*yF$yF3)

v_damp = data.frame(yF[,c(1,3,4)],new = new,name = rep(c("yFp","yF"),each=(length(new)/2)))

weighted_vdamp =  as.data.frame(v_damp%>%dplyr::group_by(time,name)%>%dplyr::summarize(fit_weighted = sum(new*weight)/sum(weight),
                                                                                       se = sum(weight/sum(weight)*abs(new-fit_weighted)),
                                                                                       upper = fit_weighted+se,
                                                                                       lower = fit_weighted-se,
                                                                                       median_weighted = weighted.quantile(new,weight,probs = 0.5),
                                                                                       twentyfive = weighted.quantile(new,weight,probs = 0.25),
                                                                                       seventyfive = weighted.quantile(new,weight,probs = 0.75),
                                                                                       five = weighted.quantile(new,weight,probs = 0.05),
                                                                                       ninetyfive = weighted.quantile(new,weight,probs = 0.95)
))

B4 = ggplot(subset(weighted_vdamp,name=="yFp"))+
    geom_line(aes(x = time,y = median_weighted,colour = name),size = 1,alpha =1)+
    geom_ribbon(aes(x = time,ymin = twentyfive,ymax = seventyfive,fill = name),alpha = 0.5)+
    #geom_ribbon(aes(x = time,ymin = five,ymax = ninetyfive,fill = name),alpha = 0.25)+
    theme_bw()+
    labs(x = "day",y = "measure",colour = "Term",fill = "Term",title = "C")+
    scale_x_continuous(breaks = c(0,2,4,6,8,10))+scale_y_continuous(label=scientific_10_3,breaks = c(0,10^4,2*10^4,3*10^4,4*10^4),limits = c(0,4*10^4))+
    theme(text = element_text(size=12),
          legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          title = element_text(size = 12))+
    labs(x ="days since infection",y = "per-cell-per-day virus production")

layout = matrix(c(1,2,3,NA),
                nrow = 2,byrow = TRUE)
multiplot(A3,B3,B4,layout = layout)


tiff("Fig5.tiff", height = 19, width =13.2, units = 'cm',
     compression = "lzw", res = 600)
multiplot(A3,B3,B4,layout = layout)

dev.off()

tiff("Fig5_legend.tiff", height = 19, width =13.2, units = 'cm',
     compression = "lzw", res = 600)
grid.newpage()
grid.draw(legend)

dev.off()


#Try changing parameter values and seeing what happens-----------

all_fit = top_models

param_list = c("m","n","yV1","b1","yV2","b2","yV3","b3")

all_dat = NULL
for(i in param_list){
    for(k in c(1,0.5,0.1,0.01,0)){
      chosen = all_fit#[grepl(paste0(i,"."), all_fit$version_dat, fixed=TRUE),]
      chosen[,which(colnames(chosen)==i)]=chosen[,which(colnames(chosen)==i)]*k
      for (j in 1:nrow(chosen)){
        rank2 = chosen[j,"rank2"]
        version = chosen[j,c("version","version_dat","weight")]
        like_og = chosen[j,"like"]
        params_chosen = as.vector(chosen[j,c(1:42)])
        
        tdat <- mcmv %>%
          trajectory(params = params_chosen,format="array")
        
        ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                        params=params_chosen)
        like <- apply(ell,1,sum)
        
        tdat <- mcmv_NA %>%
          trajectory(params = params_chosen,format="data.frame")
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
        
        all_dat = rbind(all_dat,data.frame(fit_ode,test_param =i,version,rank = rank2,change = k,like = -like,like_og = like_og,like_diff = abs(-like-like_og)))
      }
    }
  }


  for(k in c(1,0)){
    chosen = all_fit#[grepl("S.", all_fit$version_dat, fixed=TRUE),]
    chosen[,which(colnames(chosen)=="o")]=chosen[,which(colnames(chosen)=="o")]*k
    for (j in 1:nrow(chosen)){
      rank2 = chosen[j,"rank2"]
      version = chosen[j,c("version","version_dat","weight")]
      like_og = chosen[j,"like"]
      params_chosen = as.vector(chosen[j,c(1:42)])
      
      tdat <- mcmv %>%
        trajectory(params = params_chosen,format="array")
      
      ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                      params=params_chosen)
      like <- apply(ell,1,sum)
      
      tdat <- mcmv_NA %>%
        trajectory(params = params_chosen,format="data.frame")
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
      
      all_dat = rbind(all_dat,data.frame(fit_ode,test_param ="S",version,rank = rank2,change = k,like = -like,like_og = like_og,like_diff = abs(-like-like_og)))
    }
  }

param_list2 = c("yV1","b1","yV2","b2","yV3","b3")

for(k in c(1,0.5,0.1,0.01,0)){
  chosen = all_fit
  chosen[,which(colnames(chosen)%in%param_list2)]=chosen[,which(colnames(chosen)%in%param_list2)]*k
  for (j in 1:nrow(chosen)){
    rank2 = chosen[j,"rank2"]
    version = chosen[j,c("version","version_dat","weight")]
    like_og = chosen[j,"like"]
    params_chosen = as.vector(chosen[j,c(1:42)])
    
    tdat <- mcmv %>%
      trajectory(params = params_chosen,format="array")
    
    ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                    params=params_chosen)
    like <- apply(ell,1,sum)
    
    tdat <- mcmv_NA %>%
      trajectory(params = params_chosen,format="data.frame")
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
    
    all_dat = rbind(all_dat,data.frame(fit_ode,test_param ="innate",version,rank = rank2,change = k,like = -like,like_og = like_og,like_diff = abs(-like-like_og)))
  }
}

param_list3 = c("b1","b2")

for(k in c(1,0.5,0.1,0.01,0)){
  chosen = all_fit
  chosen[,which(colnames(chosen)%in%param_list3)]=chosen[,which(colnames(chosen)%in%param_list3)]*k
  for (j in 1:nrow(chosen)){
    rank2 = chosen[j,"rank2"]
    version = chosen[j,c("version","version_dat","weight")]
    like_og = chosen[j,"like"]
    params_chosen = as.vector(chosen[j,c(1:42)])
    
    tdat <- mcmv %>%
      trajectory(params = params_chosen,format="array")
    
    ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                    params=params_chosen)
    like <- apply(ell,1,sum)
    
    tdat <- mcmv_NA %>%
      trajectory(params = params_chosen,format="data.frame")
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
    
    all_dat = rbind(all_dat,data.frame(fit_ode,test_param ="bi",version,rank = rank2,change = k,like = -like,like_og = like_og,like_diff = abs(-like-like_og)))
  }
}

param_list4 = c("yV1","yV3")

for(k in c(1,0.5,0.1,0.01,0)){
  chosen = all_fit
  chosen[,which(colnames(chosen)%in%param_list4)]=chosen[,which(colnames(chosen)%in%param_list4)]*k
  for (j in 1:nrow(chosen)){
    rank2 = chosen[j,"rank2"]
    version = chosen[j,c("version","version_dat","weight")]
    like_og = chosen[j,"like"]
    params_chosen = as.vector(chosen[j,c(1:42)])
    
    tdat <- mcmv %>%
      trajectory(params = params_chosen,format="array")
    
    ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                    params=params_chosen)
    like <- apply(ell,1,sum)
    
    tdat <- mcmv_NA %>%
      trajectory(params = params_chosen,format="data.frame")
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
    
    all_dat = rbind(all_dat,data.frame(fit_ode,test_param ="yi",version,rank = rank2,change = k,like = -like,like_og = like_og,like_diff = abs(-like-like_og)))
  }
}

change = as.data.frame(all_dat%>%
                         dplyr::filter(type=="fit")%>%
                         dplyr::group_by(test_param,change,time,species)%>%
                         dplyr::summarize(fit_weighted = sum(fit*weight)/sum(weight),
                                          se = sum(weight/sum(weight)*abs(fit-fit_weighted)),
                                          upper = fit_weighted+se,
                                          lower = fit_weighted-se,
                                          fifty = weighted.quantile(fit,weight,probs = 0.5),
                                          twentyfive = weighted.quantile(fit,weight,probs = 0.25),
                                          seventyfive = weighted.quantile(fit,weight,probs = 0.75))
)



colours4 = c("#0000FF","#CC79A7","#99CC00","#33CCFF","#E69F00")

change$test_param =  revalue(change$test_param,c(
  "b1"="infected cell clearance\nby IFI27","b2"="infected cell clearance\nby IFI6","b3"="infected cell clearance\nby IFI16","bi" = "infected cell clearance\nby IFNs",
  "yA"="hinderance of infection\nby anti-spike IgG",
  "yV1"="hinderance of viral\nproduction by IFI27","yV2"="hinderance of viral\nproduction by IFI6","yV3"="hinderance of viral\nproduction by IFI16","yi" = "hinderance of viral\nproduction by IFNs",
  "n"="infected cell clearance\nby CD4s",
  "m"="infected cell clearance\nby CD8s",
  "S"="target cell limitation",
  "innate"="all aspects of\ninnate response"))
change = subset(change,test_param%in%c("all aspects of\ninnate response","infected cell clearance\nby CD4s","infected cell clearance\nby IFI27","infected cell clearance\nby CD8s","hinderance of viral\nproduction by IFI27","target cell limitation","hinderance of viral\nproduction by IFI16","infected cell clearance\nby IFI16"))

change$test_param =factor(change$test_param,levels = c("all aspects of\ninnate response","infected cell clearance\nby CD4s","infected cell clearance\nby IFI27","infected cell clearance\nby CD8s","hinderance of viral\nproduction by IFI27","target cell limitation","hinderance of viral\nproduction by IFI16","infected cell clearance\nby IFI16"))

data2 = data.frame(time = seq(-1,12,0.1),value=3000)

one = ggplot()+
    geom_line(data = data2,aes(x = time,y = value),colour = "red",size = 1)+
    coord_cartesian(xlim = c(0, 10))+
    geom_line(data = subset(change,change%in%c(0,0.01,0.1,0.5,1)&species=="V"&!test_param%in%c("infected cell clearance\nby IFI27","infected cell clearance\nby CD4s","all aspects of\ninnate response")),aes(x = time,y = fifty,colour = as.character(change),group = paste0(test_param,change)),size = 1,alpha = 1)+
    facet_wrap(~test_param)+scale_y_log10(breaks = c(10^3,10^5,10^7,10^9,10^11),limits = c(10^3,10^11),label=scientific_10_1)+
    theme_bw()+scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    #scale_colour_manual("Scaling Factor",values = colours4)+
    #scale_fill_manual("Scaling Factor",values = colours4)+
    labs(x = "days post infection",y = "Genomic RNA copies/ml from BAL",colour = "Scaling Factor")+
    theme(text = element_text(size=10),
          legend.position = "right",
          legend.text = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))

one_data = subset(change,change%in%c(1)&species=="V"&test_param%in%c("infected cell clearance\nby IFI27","infected cell clearance\nby CD4s","all aspects of\ninnate response"))
one_data$test_param =factor(one_data$test_param,levels = c("all aspects of\ninnate response","infected cell clearance\nby IFI27","infected cell clearance\nby CD4s"))
one_data = one_data[,c(1,3:11)]

change2  = subset(change,change%in%c(0.01,0.1,0.5)&species=="V"&test_param%in%c("infected cell clearance\nby IFI27","infected cell clearance\nby CD4s","all aspects of\ninnate response"))
change2$test_param =factor(change2$test_param,levels = c("all aspects of\ninnate response","infected cell clearance\nby IFI27","infected cell clearance\nby CD4s"))
change2$change2 = as.character(change2$change)
change2$change2 = factor(change2$change2,levels = c("0.5","0.1","0.01"))

two = ggplot()+
    geom_line(data = data2,aes(x = time,y = value),colour = "red",size = 1)+
    coord_cartesian(xlim = c(0, 10))+
    geom_line(data = one_data,
              aes(x = time,y = fifty),size = 1,alpha = 1,colour = "purple")+
    geom_line(data = change2,
              aes(x = time,y = fifty,colour = as.character(test_param),group = paste0(test_param,change)),size = 1,alpha = 1)+
    
    geom_ribbon(data = change2,
                aes(x = time,ymax = seventyfive,ymin = twentyfive,fill = as.character(test_param),group = paste0(test_param,change)),size = 1,alpha = 0.5)+
    facet_grid(change2~test_param)+scale_y_log10(breaks = c(10^3,10^5,10^7,10^9,10^11),limits = c(10^3,10^11),label=scientific_10_1)+
    theme_bw()+scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    #scale_colour_manual("Scaling Factor",values = colours4)+
    #scale_fill_manual("Scaling Factor",values = colours4)+
    labs(x = "days post infection",y = "Genomic RNA copies/ml from BAL",colour = "Scaling Factor",fill = "Scaling Factor")+
    theme(text = element_text(size=10),
          legend.position = "none",
          legend.text = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10))



tiff("Fig6.tiff", height = 15, width =18, units = 'cm',
     compression = "lzw", res = 600)
two
dev.off()

tiff("Supp5.tiff", height = 15, width =18, units = 'cm',
     compression = "lzw", res = 600)
one
dev.off()



