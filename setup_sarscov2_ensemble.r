library(pomp)
library(tidyr)
library(plyr)
library(dplyr)

#read in all data

data_all = read.csv("all_dat_extrapolated.csv")
data_all = subset(data_all,day%in%c(0,1,2,4,7,10))
data_all[data_all$type=="gRNA"&data_all$day==0,"day"]=-1
colnames(data_all)[4]="species"
colnames(data_all)[3]="fit"
colnames(data_all)[1]="time"

#turn v specific percent to a number
data_all$fit[data_all$species=="CD8_unadjusted"]=round((data_all$fit[data_all$species=="CD8_unadjusted"]/100*1100)/(1-data_all$fit[data_all$species=="CD8_unadjusted"]/100),0)
data_all$fit[data_all$species=="CD4_unadjusted"]=round((data_all$fit[data_all$species=="CD4_unadjusted"]/100*814)/(1-data_all$fit[data_all$species=="CD4_unadjusted"]/100),0)

monkey_list = c("DG3V", "DGCX", "DHGF", "DHKM", "DGRX", "DG4i")
data_all$monkey=revalue(data_all$monkey,c("DG3V"="1","DGCX"="2","DHGF"="3","DHKM"="4","DGRX"="5","DG4i"="6"))

data_1 = subset(data_all,species%in%c("CD8_unadjusted","CD4_unadjusted","gRNA","IgG"))
data_1$species = revalue(data_1$species,c("CD8_unadjusted"="E","CD4_unadjusted"="T","gRNA"="V","IgG"="A"))

#calculate starting value for IgG
azero = median(subset(data_1,species=="A"&time==0)$fit)

#choose IFNs included in model and define mapping for each iteration-----------------------------
ifn="IFI27"
ifn2 = "IFI6"
ifn3 = 'IFI16'

ifn_dat = subset(data_all,species==ifn)
f1zero = round(median(subset(ifn_dat,time==0)$fit),0)

ifn2_dat = subset(data_all,species==ifn2)
f2zero = round(median(subset(ifn2_dat,time==0)$fit),0)

ifn3_dat = subset(data_all,species==ifn3)
f3zero = round(median(subset(ifn3_dat,time==0)$fit),0)

#get all possible variable assignments-----------------------------------------------------------
variables = NULL
for (m_chosen in c(0.1,0)){
  for (n_chosen in c(0.1,0)){
    for (yV1_chosen in c(0.1,0)){
      for(b1_chosen in c(0.1,0)){ 
        for (yV2_chosen in c(0)){
          for(b2_chosen in c(0)){ 
            for (yV3_chosen in c(0)){
              for(b3_chosen in c(0)){ 
                for(yA_chosen in c(10^-6,0)){ #0.1
                  for(o_chosen in c(1,0)){
                    variables = rbind(variables,data.frame(m = m_chosen,
                                                           n = n_chosen,
                                                           yV1 = yV1_chosen,
                                                           yV2 = yV2_chosen,
                                                           yV3 = yV3_chosen,
                                                           b1 = b1_chosen,
                                                           b2 = b2_chosen,
                                                           b3 = b3_chosen,
                                                           yA = yA_chosen,
                                                           o = o_chosen))
                  }}}}}}}}}}

for (m_chosen in c(0.1,0)){
  for (n_chosen in c(0.1,0)){
    for (yV1_chosen in c(0)){
      for(b1_chosen in c(0)){ 
        for (yV2_chosen in c(0.1,0)){
          for(b2_chosen in c(0.1,0)){ 
            for (yV3_chosen in c(0)){
              for(b3_chosen in c(0)){ 
                for(yA_chosen in c(10^-6,0)){ #0.1
                  for(o_chosen in c(1,0)){
                    variables = rbind(variables,data.frame(m = m_chosen,
                                                           n = n_chosen,
                                                           yV1 = yV1_chosen,
                                                           yV2 = yV2_chosen,
                                                           yV3 = yV3_chosen,
                                                           b1 = b1_chosen,
                                                           b2 = b2_chosen,
                                                           b3 = b3_chosen,
                                                           yA = yA_chosen,
                                                           o = o_chosen))
                  }}}}}}}}}}

for (m_chosen in c(0.1,0)){
  for (n_chosen in c(0.1,0)){
    for (yV1_chosen in c(0)){
      for(b1_chosen in c(0)){ 
        for (yV2_chosen in c(0)){
          for(b2_chosen in c(0)){ 
            for (yV3_chosen in c(0.1,0)){
              for(b3_chosen in c(0.1,0)){ 
                for(yA_chosen in c(10^-6,0)){ #0.1
                  for(o_chosen in c(1,0)){
                    variables = rbind(variables,data.frame(m = m_chosen,
                                                           n = n_chosen,
                                                           yV1 = yV1_chosen,
                                                           yV2 = yV2_chosen,
                                                           yV3 = yV3_chosen,
                                                           b1 = b1_chosen,
                                                           b2 = b2_chosen,
                                                           b3 = b3_chosen,
                                                           yA = yA_chosen,
                                                           o = o_chosen))
                  }}}}}}}}}}

variables = distinct(variables)

variables$version = 1:nrow(variables)

#select correct ordering of IFNs based on ifn_combos
data = rbind(data_1,
             data.frame(ifn_dat[,c(1:3)],species="F1"),
             data.frame(ifn2_dat[,c(1:3)],species="F2"),
             data.frame(ifn3_dat[,c(1:3)],species="F3"))
data_plot = data
#select correct initial conditions for IFNs
inits = data.frame(f1zero,f2zero,f3zero)
colnames(inits) = c("F1","F2","F3")

#format data appropriately
data$fit = round(data$fit,0)
data = subset(data,species%in%c("V","F1","F2","F3","F4","E","T","A"))
data$id = paste0("N_",data$species,"_",data$monkey)
data = data[,c(1,3,5)]
data_chosen = spread(data,id,fit)
no_dat = data.frame(time = seq(0,10,1))
data_chosen = merge(no_dat,data_chosen)

init <- Csnippet("
                  V_r = V_0;
                   V = V_0+3000;
                   I = 0;
                   F1 = F1_0;
                   F2 = F2_0;
                   F3 = F3_0;
                   E = E_0;
                   T = T_0;
                   A = A_0;
                   S = S_0;
                   ")

dmeas <- Csnippet("
      lik=0;
      if (!ISNA(N_F1_1)){
      lik+=dpois(N_F1_1,F1*rho2+1e-6,1);
      }
      if (!ISNA(N_F1_2)){
      lik+=dpois(N_F1_2,F1*rho2+1e-6,1);
      }
      if (!ISNA(N_F1_3)){
      lik+=dpois(N_F1_3,F1*rho2+1e-6,1);
      }
      if (!ISNA(N_F1_4)){
      lik+=dpois(N_F1_4,F1*rho2+1e-6,1);
      }
      
      if (!ISNA(N_F2_1)){
      lik+=dpois(N_F2_1,F2*rho2+1e-6,1);
      }
      if (!ISNA(N_F2_2)){
      lik+=dpois(N_F2_2,F2*rho2+1e-6,1);
      }
      if (!ISNA(N_F2_3)){
      lik+=dpois(N_F2_3,F2*rho2+1e-6,1);
      }
      if (!ISNA(N_F2_4)){
      lik+=dpois(N_F2_4,F2*rho2+1e-6,1);
      }
      
      if (!ISNA(N_F3_1)){
      lik+=dpois(N_F3_1,F3*rho2+1e-6,1);
      }
      if (!ISNA(N_F3_2)){
      lik+=dpois(N_F3_2,F3*rho2+1e-6,1);
      }
      if (!ISNA(N_F3_3)){
      lik+=dpois(N_F3_3,F3*rho2+1e-6,1);
      }
      if (!ISNA(N_F3_4)){
      lik+=dpois(N_F3_4,F3*rho2+1e-6,1);
      }
      
      if (!ISNA(N_E_1)){
      lik+=dpois(N_E_1,E*rho2+1e-6,1);
      }
      if (!ISNA(N_E_2)){
      lik+=dpois(N_E_2,E*rho2+1e-6,1);
      }
      if (!ISNA(N_E_3)){
      lik+=dpois(N_E_3,E*rho2+1e-6,1);
      }
      if (!ISNA(N_E_4)){
      lik+=dpois(N_E_4,E*rho2+1e-6,1);
      }
      if (!ISNA(N_E_5)){
      lik+=dpois(N_E_5,E*rho2+1e-6,1);
      }
      if (!ISNA(N_E_6)){
      lik+=dpois(N_E_6,E*rho2+1e-6,1);
      }
      
      if (!ISNA(N_T_1)){
      lik+=dpois(N_T_1,T*rho2+1e-6,1);
      }
      if (!ISNA(N_T_2)){
      lik+=dpois(N_T_2,T*rho2+1e-6,1);
      }
      if (!ISNA(N_T_3)){
      lik+=dpois(N_T_3,T*rho2+1e-6,1);
      }
      if (!ISNA(N_T_4)){
      lik+=dpois(N_T_4,T*rho2+1e-6,1);
      }
      if (!ISNA(N_T_5)){
      lik+=dpois(N_T_5,T*rho2+1e-6,1);
      }
      if (!ISNA(N_T_6)){
      lik+=dpois(N_T_6,T*rho2+1e-6,1);
      }
      
      if(!ISNA(N_V_1)){
      lik+=dnorm(log(N_V_1),log(V),rho1,1);
      }
      if(!ISNA(N_V_2)){
      lik+=dnorm(log(N_V_2),log(V),rho1,1);
      }
      if(!ISNA(N_V_3)){
      lik+=dnorm(log(N_V_3),log(V),rho1,1);
      }
      if(!ISNA(N_V_4)){
      lik+=dnorm(log(N_V_4),log(V),rho1,1);
      }
      if(!ISNA(N_V_5)){
      lik+=dnorm(log(N_V_5),log(V),rho1,1);
      }
      if(!ISNA(N_V_6)){
      lik+=dnorm(log(N_V_6),log(V),rho1,1);
      }
      
      if(!ISNA(N_A_1)){
      lik+=dnorm(log(N_A_1),log(A),rho1,1);
      }
      if(!ISNA(N_A_2)){
      lik+=dnorm(log(N_A_2),log(A),rho1,1);
      }
      if(!ISNA(N_A_3)){
      lik+=dnorm(log(N_A_3),log(A),rho1,1);
      }
      if(!ISNA(N_A_4)){
      lik+=dnorm(log(N_A_4),log(A),rho1,1);
      }
      if(!ISNA(N_A_5)){
      lik+=dnorm(log(N_A_5),log(A),rho1,1);
      }
      if(!ISNA(N_A_6)){
      lik+=dnorm(log(N_A_6),log(A),rho1,1);
      }
      
      lik = (give_log) ? lik: exp(lik);
      ")


model <- vectorfield(
  Csnippet("   double wT;
              if(t<tauT*10)
              wT = 0;
              else
              wT = alphaT;
              
              double wE;
              if(t<tauE*10)
              wE = 0;
              else
              wE = alphaE;
              
              double wA;
              if(t<tauA*10)
              wA = 0;
              else
              wA = alphaA;
              
              double wF1;
              if(t<tauF1*10)
              wF1 = 0;
              else
              wF1 = beta1;
              
              double wF2;
              if(t<tauF2*10)
              wF2 = 0;
              else
              wF2 = beta2;
              
              double wF3;
              if(t<tauF3*10)
              wF3 = 0;
              else
              wF3 = beta3;
              
              
              DV_r = p*exp(-yV1*F1)*exp(-yV2*F2)*exp(-yV3*F3)*I-c*V_r;
              DV = DV_r;
              
              DI =eta/S_0*exp(-yA*(A))*S*V_r-delta*I-m*I*E-n*I*T-b1*I*F1-b2*I*F2-b3*I*F3;
              
              DF1 = F1_0*z1+wF1*I-z1*F1;
              
              DF2 = F2_0*z2+wF2*I-z2*F2;
              
              DF3 = F3_0*z3+wF3*I-z3*F3;
              
              DE = -dE*E+wE*I;
              
              DT =-dT*T+wT*I;

              DA = wA*I;
              
              DS = S_0*ds-S*ds-o*eta/S_0*exp(-yA*(A))*S*V_r;
           
           "
  ))


param_names = c("V_0","A_0","F1_0","F2_0","F3_0","S_0","T_0","E_0",
                "rho1","rho2",
                "p","c",
                "yV1","yV2","yV3","eta","delta",
                "z1","z2","z3","beta1","beta2","beta3",
                "m","n",
                "alphaE","dE","tauE",
                "alphaT","dT","tauT","tauF1","tauF2","tauF3",
                "b1","b2","b3","alphaA","tauA","yA","ds","o")

state_names = c("V_r","V","I","F1","F2","F3","E","T","A","S")

par_trans = parameter_trans(log = c("rho1","p","yV1","yV2","yV3","c","eta","delta","alphaE","alphaT","dE","dT","m","n",
                                    "b1","b2","b3","alphaA","S_0","beta1","beta2","beta3","ds"),
                            logit = c("rho2","tauA","tauT","tauE","tauF1","tauF2","tauF3","z1","z2","z3","yA"))

mcmv <- pomp(
  data=data_chosen,
  times="time",t0=0,
  skeleton = model,
  dmeasure=dmeas, 
  rinit = init,
  partrans=par_trans,
  paramnames=param_names,
  statenames=state_names)

data_NA = merge(data.frame(time = seq(0,10,0.1)),data_chosen,all = TRUE)
data_NA[1:nrow(data_NA),2:ncol(data_NA)]=NA

mcmv_NA <- pomp(
  data=data_NA,
  times="time",t0=0,
  skeleton = model,
  dmeasure=dmeas, 
  rinit = init,
  partrans=par_trans,
  paramnames=param_names,
  statenames=state_names)


