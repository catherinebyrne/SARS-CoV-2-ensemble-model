source("setup_sarscov2_ensemble.r")

#run model---------------------------------------------------------------------------------------
for (i in 1:4){ #repeat fitting 4 times
  for(j in 1:nrow(variables)){
    print(j)
    input = as.numeric(variables[j,])
    m_chosen = input[1]
    n_chosen = input[2]
    yV1_chosen = input[3]
    yV2_chosen = input[4]
    yV3_chosen = input[5]
    b1_chosen = input[6]
    b2_chosen = input[7]
    b3_chosen = input[8]
    yA_chosen = input[9]
    o_chosen = input[10]
    version = input[11]

    #parameters we want to fit
    enames = c("p","c",
               "eta","delta",
               "beta1","z1","tauF1","beta2","z2","tauF2","beta3","z3","tauF3",
               "alphaE","dE",
               "tauE",
               "alphaT","dT",
               "tauT",
               "alphaA",
               "tauA")
    version_dat = NULL
    if(n_chosen!=0){
      enames = c(enames,"n")
      version_dat = paste0(version_dat,"n.")
    }
    if(yA_chosen!=0){
      enames = c(enames,"yA")
      version_dat = paste0(version_dat,"yA.")
    }
    if(m_chosen!=0){
      enames = c(enames,"m")
      version_dat = paste0(version_dat,"m.")
    }
    if(yV1_chosen!=0){
      enames = c(enames,"yV1")
      version_dat = paste0(version_dat,"yV1.")
    }
    if(yV2_chosen!=0){
      enames = c(enames,"yV2")
      version_dat = paste0(version_dat,"yV2.")
    }
    if(yV3_chosen!=0){
      enames = c(enames,"yV3")
      version_dat = paste0(version_dat,"yV3.")
    }
    if(b1_chosen!=0){
      enames = c(enames,"b1")
      version_dat = paste0(version_dat,"b1.")
    }
    if(b2_chosen!=0){
      enames = c(enames,"b2")
      version_dat = paste0(version_dat,"b2.")
    }
    if(b3_chosen!=0){
      enames = c(enames,"b3")
      version_dat = paste0(version_dat,"b3.")
    }
    if(o_chosen!=0){
      enames = c(enames,"ds")
      version_dat = paste0(version_dat,"S.")
    }
    if(is.null(version_dat)){
      version_dat = "none"
    }


    if (i!=1){
      set.seed(i)
      rand = runif(31,0.9,1.1)
    }else{
      rand = runif(31,1,1)
    }
    params_init = c(V_0 = 14710,A_0 = 366592,T_0 = 0,E_0 = 0,F1_0=f1zero,F2_0=f2zero,F3_0=f3zero,rho1=0.5,rho2=0.95,
                    p=29312.99*rand[1],
                    yV1=yV1_chosen*rand[2],
                    yV2=yV2_chosen*rand[3],
                    yV3=yV3_chosen*rand[4],
                    c=36.66984*rand[5],
                    eta=0.01816938*rand[6],
                    delta=2.373594*rand[7],
                    beta1=0.6173115*10^-4*rand[8],
                    z1=0.1649753*rand[9],
                    beta2=0.99*10^-4*rand[10],
                    z2=0.5973211*rand[11],
                    beta3=0.99*10^-4*rand[12],
                    z3=0.3697688*rand[13],
                    alphaE = 0.5443405*rand[14],
                    alphaT = 0.007150051*rand[15],
                    dE = 0.01107405*rand[16],
                    dT =0.0300499*rand[17],
                    m = m_chosen*rand[18],
                    n = n_chosen*rand[19],
                    b1 = b1_chosen*rand[20],
                    b2 = b2_chosen*rand[21],
                    b3 = b3_chosen*rand[22],
                    alphaA = 10^(-8)*rand[23],
                    tauA = 0.9*rand[24],
                    tauT = 0.5411309*rand[25],
                    tauE = 0.6934965*rand[26],
                    tauF1 = 0.003836362*rand[27],
                    tauF2 = 0.003836362*rand[28],
                    tauF3 = 0.003836362*rand[29],
                    yA = yA_chosen*rand[30],
                    S_0 = 1.645*10^9,#taken from Leander et al 2021 - not fitting
                    ds = 0.029*rand[31],
                    o = o_chosen
    )

    pars_fit = NULL
    error = NULL
    k=0
    likelihood = 10^8
    diff = 1

    while(is.null(error)&diff>0.5&length(getLoadedDLLs())<614){
      k=k+1
      #create objective function -  quantifies the mismatch between model predictions and data
      #ofun saves information each time it is evaluated
      tryCatch({
        ofun <- mcmv %>%
          traj_objfun(
            est=enames,
            dmeasure=dmeas,
            paramnames=param_names,
            statenames=state_names,
            params =params_init
          )
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
      #fit - optimizer searches parameter space to find parameters under which the likelihood of the data, given a trajectory of the deterministic skeleton, is maximized.
    
      tryCatch({
        fit <- optim(
          fn=ofun,
          par=coef(ofun,enames,transform=TRUE),
          method="Nelder-Mead",
          control=list(trace=0)
        )
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
      #parameters of fit
      tryCatch({
        pars_fit = rbind(pars_fit,coef(ofun))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
      tryCatch({
        likelihood_new = fit$value
        diff = likelihood-likelihood_new
        likelihood = likelihood_new
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
      tryCatch({
        params_init = pars_fit[k,]
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    
    if(length(getLoadedDLLs())>=614) #R Can get over-loaded - need to restart R if it does
      stop()
    
    #now get likelihood from fits
    options(show.error.messages = FALSE)
    tdat <- try(mcmv %>%
                  trajectory(params = params_init,format="array"))
    if("try-error" %in% class(tdat)){
      params_init = pars_fit[k-1,]
      tdat <- try(mcmv %>%
                    trajectory(params =params_init,format="array"))
    }
    if("try-error" %in% class(tdat)){
      params_init = pars_fit[k-2,]
      tdat <- try(mcmv %>%
                    trajectory(params =params_init,format="array"))
    }
    options(show.error.messages = TRUE)
    
    ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                    params=params_init)
    like <- -apply(ell,1,sum)
    
    while(is.nan(like)){
      k=k-1
      params_init = pars_fit[k,]
      tdat <- try(mcmv %>%
                    trajectory(params =params_init,format="array"))
      ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                      params=params_init)
      like <- -apply(ell,1,sum)
    }
    
    fit_all = data.frame(t(params_init),like =like,n_data = nrow(data),enames = length(enames),version = version,version_dat = version_dat,i=i,round = 1)
        if(j==1&i==1){
          write.table(fit_all,"fit.csv",col.names = TRUE,row.names = FALSE,append = FALSE,sep = ",")
        }else{
          write.table(fit_all,"fit.csv",col.names = FALSE,row.names = FALSE,append = TRUE,sep = ",")
        }
  }
}

#read in fits-------------------------------------------------------------------------
fit_all = read.csv("fit.csv")

#do additional fitting----------------------------------------------------------------
#want to fit models where yA=0, starting from best-fitting place where yA!=0
yAs = subset(fit_all,yA!=0)[,c(1:42,47,48)]
yAs$version_dat = gsub("yA.","",yAs$version_dat)
yAs$version_dat[which(yAs$version_dat == "")]="none"
yAs$eta = yAs$eta*exp(-yAs$A_0*yAs$yA)
yAs$yA = 0
names = unique(subset(fit_all,version_dat%in%yAs$version_dat)[,c(46,47)])
yAs = merge(names,yAs)
yAs$version = as.numeric(yAs$version)
yAs = arrange(yAs,version,i)

write.table(yAs,"yAs.csv",col.names = TRUE,row.names = FALSE,sep= ",")

#want to fit models where S0=0, starting from best-fitting place where S0!=0
S0s = subset(fit_all,o==1)[,c(1:42,47,48)]
S0s$version_dat = gsub("S.","",S0s$version_dat)
S0s$version_dat[which(S0s$version_dat == "")]="none"
S0s$ds= 0.029 # set to initial parameter value 
S0s$o=0 
names = unique(subset(fit_all,version_dat%in%S0s$version_dat)[,c(46,47)])
S0s = merge(names,S0s)
S0s$version = as.numeric(S0s$version)
S0s = arrange(S0s,version)

write.table(S0s,"S0s.csv",col.names= TRUE,row.names=FALSE,sep = ",")

#Rerun model for yA and S0 scenarios 
for(file in c("yAs.csv","S0s.csv")){
  file_chosen = read.csv(file)
  if(file=="yAs.csv"){
    round_chosen=2
  }else{
    round_chosen=3
  }
  for (j in 1:nrow(file_chosen)){
    i = file_chosen$i[j]
    version = file_chosen$version[j]
    version_dat = file_chosen$version_dat[j]
    input = as.numeric(variables[version,])
    m_chosen = input[1]
    n_chosen = input[2]
    yV1_chosen = input[3]
    yV2_chosen = input[4]
    yV3_chosen = input[5]
    b1_chosen = input[6]
    b2_chosen = input[7]
    b3_chosen = input[8]
    yA_chosen = input[9]
    o_chosen = input[10]
      
      #parameters we want to fit
      enames = c("p","c",
                 "eta","delta",
                 "beta1","z1","tauF1","beta2","z2","tauF2","beta3","z3","tauF3",
                 "alphaE","dE",
                 "tauE",
                 "alphaT","dT",
                 "tauT",
                 "alphaA",
                 "tauA")
      version_dat = NULL
      if(n_chosen!=0){
        enames = c(enames,"n")
        version_dat = paste0(version_dat,"n.")
      }
      if(yA_chosen!=0){
        enames = c(enames,"yA")
        version_dat = paste0(version_dat,"yA.")
      }
      if(m_chosen!=0){
        enames = c(enames,"m")
        version_dat = paste0(version_dat,"m.")
      }
      if(yV1_chosen!=0){
        enames = c(enames,"yV1")
        version_dat = paste0(version_dat,"yV1.")
      }
      if(yV2_chosen!=0){
        enames = c(enames,"yV2")
        version_dat = paste0(version_dat,"yV2.")
      }
      if(yV3_chosen!=0){
        enames = c(enames,"yV3")
        version_dat = paste0(version_dat,"yV3.")
      }
      if(b1_chosen!=0){
        enames = c(enames,"b1")
        version_dat = paste0(version_dat,"b1.")
      }
      if(b2_chosen!=0){
        enames = c(enames,"b2")
        version_dat = paste0(version_dat,"b2.")
      }
      if(b3_chosen!=0){
        enames = c(enames,"b3")
        version_dat = paste0(version_dat,"b3.")
      }
      if(o_chosen!=0){
        enames = c(enames,"ds")
        version_dat = paste0(version_dat,"S.")
      }
      if(is.null(version_dat)){
        version_dat = "none"
      }
      
    
      params_init = as.vector(file_chosen[j,3:(ncol(file_chosen)-1)][1,])
      
      pars_fit = NULL
      error = NULL
      k=0
      likelihood = 10^8
      diff = 1
      
      while(is.null(error)&k<=0&length(getLoadedDLLs())<614){
        k=k+1
        #create objective function -  quantifies the mismatch between model predictions and data
        #ofun saves information each time it is evaluated
        tryCatch({
          ofun <- mcmv %>%
            traj_objfun(
              est=enames,
              dmeasure=dmeas,
              paramnames=param_names,
              statenames=state_names,
              params =params_init
            )
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        #fit - optimizer searches parameter space to find parameters under which the likelihood of the data, given a trajectory of the deterministic skeleton, is maximized.
        
        tryCatch({
          fit <- optim(
            fn=ofun,
            par=coef(ofun,enames,transform=TRUE),
            method="Nelder-Mead",
            control=list(trace=0)
          )
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        #parameters of fit
        tryCatch({
          pars_fit = rbind(pars_fit,coef(ofun))
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        tryCatch({
          likelihood_new = fit$value
          diff = likelihood-likelihood_new
          likelihood = likelihood_new
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        tryCatch({
          params_init = pars_fit[k,]
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
      
      if(length(getLoadedDLLs())>=614)
        stop()
      
      #now get likelihood from fits
      options(show.error.messages = FALSE)
      tdat <- try(mcmv %>%
                    trajectory(params = params_init,format="array"))
      if("try-error" %in% class(tdat)){
        params_init = pars_fit[k-1,]
        tdat <- try(mcmv %>%
                      trajectory(params =params_init,format="array"))
      }
      if("try-error" %in% class(tdat)){
        params_init = pars_fit[k-2,]
        tdat <- try(mcmv %>%
                      trajectory(params =params_init,format="array"))
      }
      options(show.error.messages = TRUE)
      
      ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                      params=params_init)
      like <- -apply(ell,1,sum)
      
      while(is.nan(like)){
        k=k-1
        params_init = pars_fit[k,]
        tdat <- try(mcmv %>%
                      trajectory(params =params_init,format="array"))
        ell <- dmeasure(mcmv,y=obs(mcmv),x=tdat,times=time(mcmv),log=TRUE,
                        params=params_init)
        like <- -apply(ell,1,sum)
      }
      
      fit_all = data.frame(t(params_init),like =like,n_data = nrow(data),enames = length(enames),version = version,version_dat = version_dat,i=i,round = round_chosen)
      write.table(fit_all,"fit.csv",col.names = FALSE,row.names = FALSE,append = TRUE,sep = ",")
      
  }
}


