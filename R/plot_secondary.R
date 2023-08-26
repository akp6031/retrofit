### Plots for paper
# 20220929
# Error in match.names(clabs, names(xi)) : 
# names do not match previous names
# **sequences problem
# colnames(temp)=c("rmse","model")
# temp$method="N=10,M=3"
#
# 20220929
# Insufficient values in manual scale. 7 needed but only 4 provided.
# match models
# temp=data.frame(model=c(1,7,8,9),rmse=0.9,
## Simulation Plots
RetrofitPlotSecondary <- function(dir, iter, lambdas, file, save=FALSE){
  library(ggplot2)
  library(cowplot)
  library(matrixStats)
  library(DescTools)
  path_grid = paste(dir, file, ".pdf", sep="")
  path_grid_col3 = paste(dir, file,"_col3", ".pdf", sep="")
  Method3 = "Extra 5 cell types"
  ### Read synthetic data
  # iter = 4000
  folder1 = paste(dir, "/iter_",iter,"_lamb", lambdas[[1]],"_seed1-12-11",sep="")
  folder2 = paste(dir, "/iter_",iter,"_lamb", lambdas[[2]],"_seed1-12-11",sep="")
  folder3 = paste(dir, "/iter_",iter,"_lamb", lambdas[[3]],"_seed1-12-11",sep="")
  folder4 = paste(dir, "/iter_",iter,"_lamb", lambdas[[4]],"_seed1-12-11",sep="")
  folder5 = paste(dir, "/iter_",iter,"_lamb", lambdas[[5]],"_seed1-12-11",sep="")
  folder6 = paste(dir, "/iter_",iter,"_lamb", lambdas[[6]],"_seed1-12-11",sep="")
  folder7 = paste(dir, "/iter_",iter,"_lamb", lambdas[[7]],"_seed1-12-11",sep="")
  folder_source = paste(dir, "/Source/", sep="")
  setwd(folder_source)
  # hard code => 0.01
  folder_reference = paste(dir, "/iter_",iter,"_lamb", 0.01,"_seed1-12-11",sep="")
  ### plot for cor(W)
  W=as.matrix(read.csv("Cerebellum_W_K=10.csv")[,-1])
  W_mod1=read.csv(paste(folder_reference, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod2=read.csv(paste(folder_reference, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod3=read.csv(paste(folder_reference, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod4=W_mod3
  W_mod5=read.csv(paste(folder1, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod6=read.csv(paste(folder1, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod7=read.csv(paste(folder1, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod8=W_mod7
  W_mod9=read.csv(paste(folder2, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod10=read.csv(paste(folder2, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod11=read.csv(paste(folder2, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod12=W_mod11
  W_mod13=read.csv(paste(folder3, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod14=read.csv(paste(folder3, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod15=read.csv(paste(folder3, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod16=W_mod15
  W_mod17=read.csv(paste(folder4, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod18=read.csv(paste(folder4, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod19=read.csv(paste(folder4, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod20=W_mod19
  W_mod21=read.csv(paste(folder5, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod22=read.csv(paste(folder5, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod23=read.csv(paste(folder5, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod24=W_mod23
  W_mod25=read.csv(paste(folder6, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod26=read.csv(paste(folder6, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod27=read.csv(paste(folder6, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod28=W_mod27
  W_mod29=read.csv(paste(folder7, "/mapped/N=20,M=5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod30=read.csv(paste(folder7, "/mapped/N=10,M=3_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod31=read.csv(paste(folder7, "/mapped/extra5_loc_W_mapped_cor.csv", sep=""))[,-1]
  W_mod32=W_mod31
  cor_W1=sort(diag(cor(W,W_mod1)),decreasing=T,na.last=T)
  cor_W2=sort(diag(cor(W,W_mod2)),decreasing=T,na.last=T)
  cor_W3=sort(diag(cor(W,W_mod3)),decreasing=T,na.last=T)[1:5]
  cor_W4=sort(diag(cor(W,W_mod4)),decreasing=T,na.last=T)[1:5]
  cor_W5=sort(diag(cor(W,W_mod5)),decreasing=T,na.last=T)
  cor_W6=sort(diag(cor(W,W_mod6)),decreasing=T,na.last=T)
  cor_W7=sort(diag(cor(W,W_mod7)),decreasing=T,na.last=T)[1:5]
  cor_W8=sort(diag(cor(W,W_mod8)),decreasing=T,na.last=T)[1:5]
  cor_W9=sort(diag(cor(W,W_mod9)),decreasing=T,na.last=T)
  cor_W10=sort(diag(cor(W,W_mod10)),decreasing=T,na.last=T)
  cor_W11=sort(diag(cor(W,W_mod11)),decreasing=T,na.last=T)[1:5]
  cor_W12=sort(diag(cor(W,W_mod12)),decreasing=T,na.last=T)[1:5]
  cor_W13=sort(diag(cor(W,W_mod13)),decreasing=T,na.last=T)
  cor_W14=sort(diag(cor(W,W_mod14)),decreasing=T,na.last=T)
  cor_W15=sort(diag(cor(W,W_mod15)),decreasing=T,na.last=T)[1:5]
  cor_W16=sort(diag(cor(W,W_mod16)),decreasing=T,na.last=T)[1:5]
  cor_W17=sort(diag(cor(W,W_mod17)),decreasing=T,na.last=T)
  cor_W18=sort(diag(cor(W,W_mod18)),decreasing=T,na.last=T)
  cor_W19=sort(diag(cor(W,W_mod19)),decreasing=T,na.last=T)[1:5]
  cor_W20=sort(diag(cor(W,W_mod20)),decreasing=T,na.last=T)[1:5]
  cor_W21=sort(diag(cor(W,W_mod21)),decreasing=T,na.last=T)
  cor_W22=sort(diag(cor(W,W_mod22)),decreasing=T,na.last=T)
  cor_W23=sort(diag(cor(W,W_mod23)),decreasing=T,na.last=T)[1:5]
  cor_W24=sort(diag(cor(W,W_mod24)),decreasing=T,na.last=T)[1:5]
  cor_W25=sort(diag(cor(W,W_mod25)),decreasing=T,na.last=T)
  cor_W26=sort(diag(cor(W,W_mod26)),decreasing=T,na.last=T)
  cor_W27=sort(diag(cor(W,W_mod27)),decreasing=T,na.last=T)[1:5]
  cor_W28=sort(diag(cor(W,W_mod28)),decreasing=T,na.last=T)[1:5]
  cor_W29=sort(diag(cor(W,W_mod29)),decreasing=T,na.last=T)
  cor_W30=sort(diag(cor(W,W_mod30)),decreasing=T,na.last=T)
  cor_W31=sort(diag(cor(W,W_mod31)),decreasing=T,na.last=T)[1:5]
  cor_W32=sort(diag(cor(W,W_mod32)),decreasing=T,na.last=T)[1:5]
  
  
  df=as.data.frame(rbind(# cbind(1:length(cor_W1), cor_W1),
                         cbind(1:length(cor_W5), cor_W5),
                         cbind(1:length(cor_W9), cor_W9),
                         cbind(1:length(cor_W13), cor_W13),
                         cbind(1:length(cor_W17), cor_W17),
                         cbind(1:length(cor_W21), cor_W21),
                         cbind(1:length(cor_W25), cor_W25),
                         cbind(1:length(cor_W29), cor_W29)))
  df$model=c(# rep("MAIN", length(cor_W1)),
             rep(paste("LAMBDA=", lambdas[[1]], sep=""), length(cor_W5)),
             rep(paste("LAMBDA=", lambdas[[2]], sep=""),length(cor_W9)),
             rep(paste("LAMBDA=", lambdas[[3]], sep=""),length(cor_W13)),
             rep(paste("LAMBDA=", lambdas[[4]], sep=""),length(cor_W17)),
             rep(paste("LAMBDA=", lambdas[[5]], sep=""),length(cor_W21)),
             rep(paste("LAMBDA=", lambdas[[6]], sep=""),length(cor_W25)),
             rep(paste("LAMBDA=", lambdas[[7]], sep=""),length(cor_W29))
             )
  df$method="N=20,M=5 (Exact 10 cell types)"
  colnames(df)=c("x","cor","model","method")
  
  temp=as.data.frame(rbind(# cbind(1:length(cor_W2), cor_W2),
                           cbind(1:length(cor_W6), cor_W6),
                           cbind(1:length(cor_W10), cor_W10),
                           cbind(1:length(cor_W14), cor_W14),
                           cbind(1:length(cor_W18), cor_W18),
                           cbind(1:length(cor_W22), cor_W22),
                           cbind(1:length(cor_W26), cor_W26),
                           cbind(1:length(cor_W30), cor_W30)))
  temp$model=c(# rep("MAIN", length(cor_W2)),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""), length(cor_W6)),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),length(cor_W10)),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),length(cor_W14)),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),length(cor_W18)),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),length(cor_W22)),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),length(cor_W26)),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),length(cor_W30)))
  temp$method="N=10,M=3 (Exact 10 cell types)"
  colnames(temp)=c("x","cor","model","method")
  df=as.data.frame(rbind(df,temp))
  
  temp=as.data.frame(rbind(# cbind(1:length(cor_W3), cor_W3),
                           cbind(1:length(cor_W7), cor_W7),
                           cbind(1:length(cor_W11), cor_W11),
                           cbind(1:length(cor_W15), cor_W15),
                           cbind(1:length(cor_W19), cor_W19),
                           cbind(1:length(cor_W23), cor_W23),
                           cbind(1:length(cor_W27), cor_W27),
                           cbind(1:length(cor_W31), cor_W31)))
  temp$model=c(# rep("MAIN", length(cor_W3)),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""), length(cor_W7)),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),length(cor_W11)),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),length(cor_W15)),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),length(cor_W19)),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),length(cor_W23)),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),length(cor_W27)),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),length(cor_W31))
               )
  temp$method=Method3
  colnames(temp)=c("x","cor","model","method")
  df=as.data.frame(rbind(df,temp))
  
  temp=as.data.frame(rbind(# cbind(1:length(cor_W4), cor_W4),
                           cbind(1:length(cor_W8), cor_W8),
                           cbind(1:length(cor_W12), cor_W12),
                           cbind(1:length(cor_W16), cor_W16),
                           cbind(1:length(cor_W20), cor_W20),
                           cbind(1:length(cor_W24), cor_W24),
                           cbind(1:length(cor_W28), cor_W28),
                           cbind(1:length(cor_W32), cor_W32)))
  temp$model=c(# rep("MAIN", length(cor_W4)),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""),length(cor_W8)),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),length(cor_W12)),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),length(cor_W16)),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),length(cor_W20)),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),length(cor_W24)),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),length(cor_W28)),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),length(cor_W32)))
  temp$method="N=10,M=3 (Missing 2 cell types)"
  colnames(temp)=c("x","cor","model","method")
  df=as.data.frame(rbind(df,temp))
  
  df$model=factor(df$model,levels=c(# "MAIN",
                                    paste("LAMBDA=", lambdas[[1]], sep=""), 
                                    paste("LAMBDA=", lambdas[[2]], sep=""), 
                                    paste("LAMBDA=", lambdas[[3]], sep=""),
                                    paste("LAMBDA=", lambdas[[4]], sep=""),
                                    paste("LAMBDA=", lambdas[[5]], sep=""),
                                    paste("LAMBDA=", lambdas[[6]], sep=""),
                                    paste("LAMBDA=", lambdas[[7]], sep="")))
  df$method=factor(df$method,levels=c("N=20,M=5 (Exact 10 cell types)","N=10,M=3 (Exact 10 cell types)", Method3,"N=10,M=3 (Missing 2 cell types)"))
  
  
  
  ## Plot for RMSE(H)
  H=read.csv("N=20,M=5_H.csv")[,-1]
  # H_mod1=read.csv("N=20,M=5_loc_H_retrofit.csv")[,-1]
  H_mod1=read.csv(paste(folder_reference, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod2=read.csv("N=20,M=5_loc_RCTD_H.csv")[,-1]
  H_mod3=read.csv("N=20,M=5_loc_RCTD-D_H.csv")[,-1]
  H_mod4=read.csv("N=20,M=5_loc_Stereoscope.csv")[,-1]
  H_mod5=t(read.csv("N=20,M=5_loc_NMFreg.csv")[,-11])
  H_mod6=read.csv("N=20,M=5_loc_SPOTlight.csv")[,-1]
  H_mod7=read.csv(paste(folder1, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod8=read.csv(paste(folder2, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod9=read.csv(paste(folder3, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod10=read.csv(paste(folder4, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod11=read.csv(paste(folder5, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod12=read.csv(paste(folder6, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod13=read.csv(paste(folder7, "/mapped/N=20,M=5_loc_H_mapped_cor.csv", sep=""))[,-1]
  
  X=read.csv("N=20,M=5_loc_X.csv")[,-1]
  W=as.matrix(read.csv("Cerebellum_W_K=10.csv")[,-1])
  
  K=nrow(H)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H,H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H,H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H,H_mod7)),decreasing=T,na.last=T)
  cor_H8=sort(diag(cor(H,H_mod8)),decreasing=T,na.last=T)
  cor_H9=sort(diag(cor(H,H_mod9)),decreasing=T,na.last=T)
  cor_H10=sort(diag(cor(H,H_mod10)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H,H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H,H_mod12)),decreasing=T,na.last=T)
  cor_H13=sort(diag(cor(H,H_mod13)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_Model9=matrix(NA, ncol=S, nrow=nrow(H_mod9))
  H_Model10=matrix(NA, ncol=S, nrow=nrow(H_mod10))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  H_Model13=matrix(NA, ncol=S, nrow=nrow(H_mod13))
  H_True=matrix(NA, ncol=S, nrow=nrow(H))
  
  rmseh1=NULL
  rmseh2=NULL
  rmseh3=NULL
  rmseh4=NULL
  rmseh5=NULL
  rmseh6=NULL
  rmseh7=NULL
  rmseh8=NULL
  rmseh9=NULL
  rmseh10=NULL
  rmseh11=NULL
  rmseh12=NULL
  rmseh13=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model8[,i]=H_mod8[,i]/sum(H_mod8[,i])
    H_Model9[,i]=H_mod9[,i]/sum(H_mod9[,i])
    H_Model10[,i]=H_mod10[,i]/sum(H_mod10[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    H_Model13[,i]=H_mod13[,i]/sum(H_mod13[,i])
    
    rmseh1=c(rmseh1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmseh2=c(rmseh2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmseh3=c(rmseh3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmseh4=c(rmseh4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmseh5=c(rmseh5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmseh6=c(rmseh6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmseh7=c(rmseh7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmseh8=c(rmseh8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
    rmseh9=c(rmseh9,sqrt(mean((H_True[,i]-H_Model9[,i])^2)))
    rmseh10=c(rmseh10,sqrt(mean((H_True[,i]-H_Model10[,i])^2)))
    rmseh11=c(rmseh11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmseh12=c(rmseh12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
    rmseh13=c(rmseh13,sqrt(mean((H_True[,i]-H_Model13[,i])^2)))
  }
  
  ## prop variance explained
  
  # W_hat, H_hat and Theta_hat from retrofit (with all components before mapping)
  # W_hat=read.csv("N=20,M=5_loc_W_hat_L=20.csv")[,-1]
  # H_hat=read.csv("N=20,M=5_loc_H_hat_L=20.csv")[,-1]
  # Thet_hat=read.csv("N=20,M=5_loc_Theta_hat_L=20.csv")[,-1]
  W_hat=read.csv(paste(folder_reference, "/decomposed/N=20,M=5_loc_W_decomposed.csv", sep=""))[,-1]
  H_hat=read.csv(paste(folder_reference, "/decomposed/N=20,M=5_loc_H_decomposed.csv", sep=""))[,-1]
  Thet_hat=read.csv(paste(folder_reference, "/decomposed/N=20,M=5_loc_TH_decomposed.csv", sep=""))[,-1]
  X_hat=(as.matrix(W_hat) %*% diag(Thet_hat) + 0.01) %*% as.matrix(H_hat)
  temp=rep(NA,length(Thet_hat))
  for(k in 1:length(temp)){
    temp[k]=sum(as.matrix(W_hat[,k] * Thet_hat[k] + 0.01) %*% as.matrix(H_hat[k,]))
  }
  prop1=sort(temp/sum(temp),decreasing=T)
  ind1=sort(temp/sum(temp),decreasing=T,index.return=T)$ix
  
  X_hat1=W %*% H_Model1
  X_hat2=W %*% H_Model2
  X_hat3=W %*% H_Model3
  X_hat4=W %*% H_Model4
  X_hat5=W %*% H_Model5
  X_hat6=W %*% H_Model6
  X_hat7=W %*% H_Model7
  X_hat8=W %*% H_Model8
  X_hat9=W %*% H_Model9
  X_hat10=W %*% H_Model10
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  X_hat13=W %*% H_Model13
  X_hat99=(as.matrix(W_hat[,ind1[1:12]]) %*% diag(Thet_hat[ind1[1:12]]) + 0.01) %*% as.matrix(H_hat[ind1[1:12],]) #no. of components that explain about 80%
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  cor_X10=sort(diag(cor(X,X_hat10)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X13=sort(diag(cor(X,X_hat13)),decreasing=T,na.last=T)
  cor_X99=sort(diag(cor(X,X_hat99)),decreasing=T,na.last=T)
  
  nrmse1=NULL
  nrmse2=NULL
  nrmse3=NULL
  nrmse4=NULL
  nrmse5=NULL
  nrmse6=NULL
  nrmse7=NULL
  nrmse8=NULL
  nrmse9=NULL
  nrmse10=NULL
  nrmse11=NULL
  nrmse12=NULL
  nrmse13=NULL
  nrmse99=NULL
  for(i in 1:S){
    nrmse1=c(nrmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    nrmse2=c(nrmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    nrmse3=c(nrmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    nrmse4=c(nrmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    nrmse5=c(nrmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    nrmse6=c(nrmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    nrmse7=c(nrmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    nrmse8=c(nrmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    nrmse9=c(nrmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
    nrmse10=c(nrmse10,sqrt(mean(((X[,i]-X_hat10[,i])/sd(X[,i]))^2)))
    nrmse11=c(nrmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    nrmse12=c(nrmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    nrmse13=c(nrmse13,sqrt(mean(((X[,i]-X_hat13[,i])/sd(X[,i]))^2)))
    nrmse99=c(nrmse99,sqrt(mean(((X[,i]-X_hat99[,i])/sd(X[,i]))^2)))
  }

  temp=as.data.frame(rbind(cbind(rmseh7,7),cbind(rmseh8,8), cbind(rmseh9,9),cbind(rmseh10,10), cbind(rmseh11,11), cbind(rmseh12,12), cbind(rmseh13,13)))
  colnames(temp)=c("rmse","model")
  temp$method="N=20,M=5"
  df1=as.data.frame(rbind(temp))
  df2=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
                          cbind(seq(0,1,length.out = 1000),cor_H7),
                          cbind(seq(0,1,length.out = 1000),cor_H8),
                          cbind(seq(0,1,length.out = 1000),cor_H9),
                          cbind(seq(0,1,length.out = 1000),cor_H10),
                          cbind(seq(0,1,length.out = 1000),cor_H11),
                          cbind(seq(0,1,length.out = 1000),cor_H12),
                          cbind(seq(0,1,length.out = 1000),cor_H13)))
  df2$model=c(# rep("MAIN",1000), 
              rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000))
  df2$method="N=20,M=5"
  colnames(df2)=c("x","cor","model","method")
  
  df3=as.data.frame(rbind(# cbind(nrmse1,1),
                          # cbind(nrmse99,99),
                          cbind(nrmse7,7), cbind(nrmse8,8),cbind(nrmse9,9),
                          cbind(nrmse10,10),cbind(nrmse11,11),
                          cbind(nrmse12,12),cbind(nrmse13,13)
  ))
  colnames(df3)=c("rmse","model")
  df3$method="N=20,M=5"
  
  df4=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
                          # cbind(seq(0,1,length.out = 1000),cor_X99),
                          cbind(seq(0,1,length.out = 1000),cor_X7),
                          cbind(seq(0,1,length.out = 1000),cor_X8),
                          cbind(seq(0,1,length.out = 1000),cor_X9),
                          cbind(seq(0,1,length.out = 1000),cor_X10),
                          cbind(seq(0,1,length.out = 1000),cor_X11),
                          cbind(seq(0,1,length.out = 1000),cor_X12),
                          cbind(seq(0,1,length.out = 1000),cor_X13)))
  df4$model=c(# rep("MAIN",1000),
              # rep("PVE80%",1000),
              rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
              rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000)
              )
  df4$method="N=20,M=5"
  colnames(df4)=c("x","cor","model","method")
  
  label = c(
    format(ks.test(rmseh1,rmseh1, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh6, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh7, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh8, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh9, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh10, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh11, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh12, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh13, alternative="greater")$p.value, digits=3)
  )
  
  ## add p-values from above KS tests
  ann.text1=data.frame(model=c(7,8,9,10,11,12,13),rmse=0.9,
                       method=factor("N=20,M=5",
                                     levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")),
                       label=label[c(7,8,9,10,11,12,13)])
  # label=c("1.7e-207","8.4e-113","3.1e-10","1.8e-281","2.4e-121"))
  
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H13),digits=3)
  )
  
  ann.text2=data.frame(label=c(# paste("MAIN:", label[[1]]),
                               paste(paste("LAMBDA=", lambdas[[1]], sep=""), label[[7]]),
                               paste(paste("LAMBDA=", lambdas[[2]], sep=""), label[[8]]),
                               paste(paste("LAMBDA=", lambdas[[3]], sep=""), label[[9]]),
                               paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                               paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                               paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                               paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                               ),
                       model=c(# "MAIN",
                               paste("LAMBDA=", lambdas[[1]], sep=""),
                               paste("LAMBDA=", lambdas[[2]], sep=""), 
                               paste("LAMBDA=", lambdas[[3]], sep=""),
                               paste("LAMBDA=", lambdas[[4]], sep=""),
                               paste("LAMBDA=", lambdas[[5]], sep=""),
                               paste("LAMBDA=", lambdas[[6]], sep=""),
                               paste("LAMBDA=", lambdas[[7]], sep="")
                               ),
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")))
  
  label = c(
    format(ks.test(nrmse1,nrmse1, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse2, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse3, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse4, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse5, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse6, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse7, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse8, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse9, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse10, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse11, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse12, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse13, alternative="greater")$p.value, digits=3)
    # format(ks.test(nrmse1,nrmse7, alternative="greater")$p.value, digits=3)
  )
  label99 = format(ks.test(nrmse1,nrmse99, alternative="greater")$p.value, digits=3)
  
  ann.text3=data.frame(model=c(# 99,
                               7,8,9,10,11,12,13),
                       rmse=0.9,
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3", Method3,"Missing 2 cell types")),
                       label=c(# label99,
                               label[[7]],label[[8]],label[[9]],label[[10]],label[[11]],label[[12]],label[[13]]))
  # label=c("5.1e-48","5.5e-10","0.186","5.9e-86","1.6e-28"))#,"2.1e-28"))
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X13),digits=3)
  )
  label99 = round(AUC(x=seq(0,1,length.out = 1000),y=cor_X99),digits=3)
  
  ann.text4=data.frame(label = c(# paste("MAIN:", label[[1]]),
                                 # paste("PVE80%:", label99),
                                 paste(paste("LAMBDA=", lambdas[[1]], sep=""), label[[7]]),
                                 paste(paste("LAMBDA=", lambdas[[2]], sep=""), label[[8]]),
                                 paste(paste("LAMBDA=", lambdas[[3]], sep=""), label[[9]]),
                                 paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                                 paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                                 paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                                 paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                                 ),
                      model=c(# "MAIN",
                              # "PVE80%",
                              paste("LAMBDA=", lambdas[[1]], sep=""),
                              paste("LAMBDA=", lambdas[[2]], sep=""), 
                              paste("LAMBDA=", lambdas[[3]], sep=""),
                              paste("LAMBDA=", lambdas[[4]], sep=""),
                              paste("LAMBDA=", lambdas[[5]], sep=""),
                              paste("LAMBDA=", lambdas[[6]], sep=""),
                              paste("LAMBDA=", lambdas[[7]], sep="")
                              ),
                      method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types"))
  )
  
  
  # N=10,M=3
  H=read.csv("N=10,M=3_H.csv")[,-1]
  # H_mod1=read.csv("N=10,M=3_loc_H_retrofit.csv")[,-1]
  H_mod1=read.csv(paste(folder_reference, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod2=read.csv("N=10,M=3_loc_RCTD_H.csv")[,-1]
  H_mod3=read.csv("N=10,M=3_loc_RCTD-D_H.csv")[,-1]
  H_mod4=read.csv("N=10,M=3_loc_Stereoscope.csv")[,-1]
  H_mod5=t(read.csv("N=10,M=3_loc_NMFreg.csv")[,-11])
  H_mod6=read.csv("N=10,M=3_loc_SPOTlight.csv")[,-1]
  H_mod7=read.csv(paste(folder1, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod8=read.csv(paste(folder2, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod9=read.csv(paste(folder3, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod10=read.csv(paste(folder4, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod11=read.csv(paste(folder5, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod12=read.csv(paste(folder6, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod13=read.csv(paste(folder7, "/mapped/N=10,M=3_loc_H_mapped_cor.csv", sep=""))[,-1]
  
  X=read.csv("N=10,M=3_loc_X.csv")[,-1]
  
  K=nrow(H)
  S=ncol(H)
  
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H,H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H,H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H,H_mod7)),decreasing=T,na.last=T)
  cor_H8=sort(diag(cor(H,H_mod8)),decreasing=T,na.last=T)
  cor_H9=sort(diag(cor(H,H_mod9)),decreasing=T,na.last=T)
  cor_H10=sort(diag(cor(H,H_mod10)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H,H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H,H_mod12)),decreasing=T,na.last=T)
  cor_H13=sort(diag(cor(H,H_mod13)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_Model9=matrix(NA, ncol=S, nrow=nrow(H_mod9))
  H_Model10=matrix(NA, ncol=S, nrow=nrow(H_mod10))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  H_Model13=matrix(NA, ncol=S, nrow=nrow(H_mod13))
  H_True=matrix(NA, ncol=S, nrow=nrow(H))
  
  rmseh1=NULL
  rmseh2=NULL
  rmseh3=NULL
  rmseh4=NULL
  rmseh5=NULL
  rmseh6=NULL
  rmseh7=NULL
  rmseh8=NULL
  rmseh9=NULL
  rmseh10=NULL
  rmseh11=NULL
  rmseh12=NULL
  rmseh13=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model8[,i]=H_mod8[,i]/sum(H_mod8[,i])
    H_Model9[,i]=H_mod9[,i]/sum(H_mod9[,i])
    H_Model10[,i]=H_mod10[,i]/sum(H_mod10[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    H_Model13[,i]=H_mod13[,i]/sum(H_mod13[,i])
    
    rmseh1=c(rmseh1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmseh2=c(rmseh2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmseh3=c(rmseh3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmseh4=c(rmseh4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmseh5=c(rmseh5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmseh6=c(rmseh6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmseh7=c(rmseh7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmseh8=c(rmseh8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
    rmseh9=c(rmseh9,sqrt(mean((H_True[,i]-H_Model9[,i])^2)))
    rmseh10=c(rmseh10,sqrt(mean((H_True[,i]-H_Model10[,i])^2)))
    rmseh11=c(rmseh11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmseh12=c(rmseh12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
    rmseh13=c(rmseh13,sqrt(mean((H_True[,i]-H_Model13[,i])^2)))
  }
  
  # prop variance explained
  # W_hat=read.csv("N=10,M=3_loc_W_hat_L=20.csv")[,-1]
  # H_hat=read.csv("N=10,M=3_loc_H_hat_L=20.csv")[,-1]
  # Thet_hat=read.csv("N=10,M=3_loc_Theta_hat_L=20.csv")[,-1]
  W_hat=read.csv(paste(folder_reference, "/decomposed/N=10,M=3_loc_W_decomposed.csv", sep=""))[,-1]
  H_hat=read.csv(paste(folder_reference, "/decomposed/N=10,M=3_loc_H_decomposed.csv", sep=""))[,-1]
  Thet_hat=read.csv(paste(folder_reference, "/decomposed/N=10,M=3_loc_TH_decomposed.csv", sep=""))[,-1]
  X_hat=(as.matrix(W_hat) %*% diag(Thet_hat) + 0.01) %*% as.matrix(H_hat)
  temp=rep(NA,length(Thet_hat))
  for(k in 1:length(temp)){
    temp[k]=sum(as.matrix(W_hat[,k] * Thet_hat[k] + 0.01) %*% as.matrix(H_hat[k,]))
  }
  prop2=sort(temp/sum(temp),decreasing=T)
  ind1=sort(temp/sum(temp),decreasing=T,index.return=T)$ix
  
  X_hat1=W %*% H_Model1
  X_hat2=W %*% H_Model2
  X_hat3=W %*% H_Model3
  X_hat4=W %*% H_Model4
  X_hat5=W %*% H_Model5
  X_hat6=W %*% H_Model6
  X_hat7=W %*% H_Model7
  X_hat8=W %*% H_Model8
  X_hat9=W %*% H_Model9
  X_hat10=W %*% H_Model10
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  X_hat13=W %*% H_Model13
  X_hat99=(as.matrix(W_hat[,ind1[1:13]]) %*% diag(Thet_hat[ind1[1:13]]) + 0.01) %*% as.matrix(H_hat[ind1[1:13],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  cor_X10=sort(diag(cor(X,X_hat10)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X13=sort(diag(cor(X,X_hat13)),decreasing=T,na.last=T)
  cor_X99=sort(diag(cor(X,X_hat99)),decreasing=T,na.last=T)
  
  nrmse1=NULL
  nrmse2=NULL
  nrmse3=NULL
  nrmse4=NULL
  nrmse5=NULL
  nrmse6=NULL
  nrmse7=NULL
  nrmse8=NULL
  nrmse9=NULL
  nrmse10=NULL
  nrmse11=NULL
  nrmse12=NULL
  nrmse13=NULL
  nrmse99=NULL
  for(i in 1:S){
    nrmse1=c(nrmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    nrmse2=c(nrmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    nrmse3=c(nrmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    nrmse4=c(nrmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    nrmse5=c(nrmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    nrmse6=c(nrmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    nrmse7=c(nrmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    nrmse8=c(nrmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    nrmse9=c(nrmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
    nrmse10=c(nrmse10,sqrt(mean(((X[,i]-X_hat10[,i])/sd(X[,i]))^2)))
    nrmse11=c(nrmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    nrmse12=c(nrmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    nrmse13=c(nrmse13,sqrt(mean(((X[,i]-X_hat13[,i])/sd(X[,i]))^2)))
    nrmse99=c(nrmse99,sqrt(mean(((X[,i]-X_hat99[,i])/sd(X[,i]))^2)))
  }
  
  # temp=as.data.frame(rbind(cbind(rmseh1,1),cbind(rmseh7,7),cbind(rmseh8,8), cbind(rmseh9,9), cbind(rmseh10,10), cbind(rmseh11,11)))
  temp=as.data.frame(rbind(cbind(rmseh7,7),cbind(rmseh8,8), cbind(rmseh9,9), cbind(rmseh10,10), cbind(rmseh11,11), cbind(rmseh12,12), cbind(rmseh13,13)))
  colnames(temp)=c("rmse","model")
  temp$method="N=10,M=3"
  df1=as.data.frame(rbind(df1,temp))
  
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
                           cbind(seq(0,1,length.out = 1000),cor_H7),
                           cbind(seq(0,1,length.out = 1000),cor_H8),
                           cbind(seq(0,1,length.out = 1000),cor_H9),
                           cbind(seq(0,1,length.out = 1000),cor_H10),
                           cbind(seq(0,1,length.out = 1000),cor_H11),
                           cbind(seq(0,1,length.out = 1000),cor_H12),
                           cbind(seq(0,1,length.out = 1000),cor_H13)))
  temp$model=c(# rep("MAIN",1000), 
               rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000)
  )
  colnames(temp)=c("x","cor","model")
  temp$method="N=10,M=3"
  df2=as.data.frame(rbind(df2,temp))
  

  temp=as.data.frame(rbind(# cbind(nrmse1,1),
                           # cbind(nrmse99,99),
                           cbind(nrmse7,7), cbind(nrmse8,8),cbind(nrmse9,9),
                           cbind(nrmse10,10),cbind(nrmse11,11),
                           cbind(nrmse12,12),cbind(nrmse13,13)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="N=10,M=3"
  df3=as.data.frame(rbind(df3,temp))
  
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
                            # cbind(seq(0,1,length.out = 1000),cor_X99),
                            cbind(seq(0,1,length.out = 1000),cor_X7),
                            cbind(seq(0,1,length.out = 1000),cor_X8),
                            cbind(seq(0,1,length.out = 1000),cor_X9),
                            cbind(seq(0,1,length.out = 1000),cor_X10),
                            cbind(seq(0,1,length.out = 1000),cor_X11),
                            cbind(seq(0,1,length.out = 1000),cor_X12),
                            cbind(seq(0,1,length.out = 1000),cor_X13)))
  
  temp$model=c(# rep("MAIN",1000),
               # rep("PVE80%",1000),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000)
  )
  temp$method="N=10,M=3"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  label = c(
    format(ks.test(rmseh1,rmseh1, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh6, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh7, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh8, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh9, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh10, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh11, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh12, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh13, alternative="greater")$p.value, digits=3)
  )
  
  temp=data.frame(model=c(7,8,9,10,11,12,13),rmse=0.9,
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3", Method3,"Missing 2 cell types")),
                  label=label[c(7,8,9,10,11,12,13)])
  # label=c("1.0e-4","0.003","0.975","2.0e-242","6.6e-83"))
  ann.text1=as.data.frame(rbind(ann.text1,temp))
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H13),digits=3)
  )
  
  temp=data.frame(label=c(# paste("MAIN:", label[[1]]),
                          paste(paste("LAMBDA=", lambdas[[1]], sep=""), label[[7]]),
                          paste(paste("LAMBDA=", lambdas[[2]], sep=""), label[[8]]),
                          paste(paste("LAMBDA=", lambdas[[3]], sep=""), label[[9]]),
                          paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                          paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                          paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                          paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                          ),
                  model=c(# "MAIN",
                          paste("LAMBDA=", lambdas[[1]], sep=""),
                          paste("LAMBDA=", lambdas[[2]], sep=""), 
                          paste("LAMBDA=", lambdas[[3]], sep=""),
                          paste("LAMBDA=", lambdas[[4]], sep=""),
                          paste("LAMBDA=", lambdas[[5]], sep=""),
                          paste("LAMBDA=", lambdas[[6]], sep=""),
                          paste("LAMBDA=", lambdas[[7]], sep="")
                          ),
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  
  label=c(
    format(ks.test(nrmse1,nrmse1, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse2, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse3, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse4, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse5, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse6, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse7, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse8, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse9, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse10, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse11, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse12, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse13, alternative="greater")$p.value, digits=3)
  )
  label99 = format(ks.test(nrmse1,nrmse99, alternative="greater")$p.value, digits=3)
  
  
  temp=data.frame(model=c(# 99,
                          7,8,9,10,11,12,13),rmse=0.9,
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",
                                                    Method3,"Missing 2 cell types")),
                  label=c(# label99,
                          label[[7]],label[[8]],label[[9]],label[[10]],label[[11]], label[[12]], label[[13]]))
  # label=c("0.171","4.0e-4","0.294","0.844","0.254"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X13),digits=3)
  )
  label99 = round(AUC(x=seq(0,1,length.out = 1000),y=cor_X99),digits=3)
  
  temp=data.frame(label = c(# paste("MAIN:", label[[1]]),
                            # paste("PVE80%:", label99),
                            paste(paste("LAMBDA=", lambdas[[1]], sep=""), label[[7]]),
                            paste(paste("LAMBDA=", lambdas[[2]], sep=""), label[[8]]),
                            paste(paste("LAMBDA=", lambdas[[3]], sep=""), label[[9]]),
                            paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                            paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                            paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                            paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                            ),
                  model=c(# "MAIN",
                          # "PVE80%:",
                          paste("LAMBDA=", lambdas[[1]], sep=""),
                          paste("LAMBDA=", lambdas[[2]], sep=""), 
                          paste("LAMBDA=", lambdas[[3]], sep=""),
                          paste("LAMBDA=", lambdas[[4]], sep=""),
                          paste("LAMBDA=", lambdas[[5]], sep=""),
                          paste("LAMBDA=", lambdas[[6]], sep=""),
                          paste("LAMBDA=", lambdas[[7]], sep="")
                          ),
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")))
  
  
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  
  ## extra 5
  H=read.csv("extra5_H.csv")[,-1]
  # H_mod1=read.csv("extra5_loc_H_retrofit.csv")[,-1]
  H_mod1=read.csv(paste(folder_reference, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod2=read.csv("extra5_loc_RCTD_H.csv")[,-1]
  H_mod3=read.csv("extra5_loc_RCTD-D_H.csv")[,-1]
  H_mod4=read.csv("extra5_loc_Stereoscope.csv")[,-1]
  H_mod5=t(read.csv("extra5_loc_NMFreg.csv")[,-11])
  H_mod6=read.csv("extra5_loc_SPOTlight.csv")[,-1]
  H_mod7=read.csv(paste(folder1, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod8=read.csv(paste(folder2, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod9=read.csv(paste(folder3, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod10=read.csv(paste(folder4, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod11=read.csv(paste(folder5, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod12=read.csv(paste(folder6, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod13=read.csv(paste(folder7, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  
  
  X=read.csv("extra5_loc_X.csv")[,-1]
  
  H_mod1=H_mod1[1:5,]
  H_mod2=H_mod2[1:5,]
  H_mod3=H_mod3[1:5,]
  H_mod4=H_mod4[1:5,]
  H_mod5=H_mod5[1:5,]
  H_mod6=H_mod6[1:5,]
  H_mod7=H_mod7[1:5,]
  H_mod8=H_mod8[1:5,]
  H_mod9=H_mod9[1:5,]
  H_mod10=H_mod10[1:5,]
  H_mod11=H_mod11[1:5,]
  H_mod12=H_mod12[1:5,]
  H_mod13=H_mod13[1:5,]
  
  K=nrow(H)
  G=nrow(W)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H,H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H,H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H,H_mod7)),decreasing=T,na.last=T)
  cor_H8=sort(diag(cor(H,H_mod8)),decreasing=T,na.last=T)
  cor_H9=sort(diag(cor(H,H_mod9)),decreasing=T,na.last=T)
  cor_H10=sort(diag(cor(H,H_mod10)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H,H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H,H_mod12)),decreasing=T,na.last=T)
  cor_H13=sort(diag(cor(H,H_mod13)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_Model9=matrix(NA, ncol=S, nrow=nrow(H_mod9))
  H_Model10=matrix(NA, ncol=S, nrow=nrow(H_mod10))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  H_Model13=matrix(NA, ncol=S, nrow=nrow(H_mod13))
  H_True=matrix(NA, ncol=S, nrow=nrow(H))
  
  rmseh1=NULL
  rmseh2=NULL
  rmseh3=NULL
  rmseh4=NULL
  rmseh5=NULL
  rmseh6=NULL
  rmseh7=NULL
  rmseh8=NULL
  rmseh9=NULL
  rmseh10=NULL
  rmseh11=NULL
  rmseh12=NULL
  rmseh13=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model8[,i]=H_mod8[,i]/sum(H_mod8[,i])
    H_Model9[,i]=H_mod9[,i]/sum(H_mod9[,i])
    H_Model10[,i]=H_mod10[,i]/sum(H_mod10[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    H_Model13[,i]=H_mod13[,i]/sum(H_mod13[,i])
    
    rmseh1=c(rmseh1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmseh2=c(rmseh2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmseh3=c(rmseh3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmseh4=c(rmseh4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmseh5=c(rmseh5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmseh6=c(rmseh6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmseh7=c(rmseh7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmseh8=c(rmseh8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
    rmseh9=c(rmseh9,sqrt(mean((H_True[,i]-H_Model9[,i])^2)))
    rmseh10=c(rmseh10,sqrt(mean((H_True[,i]-H_Model10[,i])^2)))
    rmseh11=c(rmseh11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmseh12=c(rmseh12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
    rmseh13=c(rmseh13,sqrt(mean((H_True[,i]-H_Model13[,i])^2)))
  }
  
  # prop variance explained
  # W_hat=read.csv("extra5_loc_W_hat_L=10.csv")[,-1]
  # H_hat=read.csv("extra5_loc_H_hat_L=10.csv")[,-1]
  # Thet_hat=read.csv("extra5_loc_Theta_hat_L=10.csv")[,-1]
  W_hat=read.csv(paste(folder_reference, "/decomposed/extra5_loc_W_decomposed.csv", sep=""))[,-1]
  H_hat=read.csv(paste(folder_reference, "/decomposed/extra5_loc_H_decomposed.csv", sep=""))[,-1]
  Thet_hat=read.csv(paste(folder_reference, "/decomposed/extra5_loc_TH_decomposed.csv", sep=""))[,-1]
  X_hat=(as.matrix(W_hat) %*% diag(Thet_hat) + 0.01) %*% as.matrix(H_hat)
  temp=rep(NA,length(Thet_hat))
  for(k in 1:length(temp)){
    temp[k]=sum(as.matrix(W_hat[,k] * Thet_hat[k] + 0.01) %*% as.matrix(H_hat[k,]))
  }
  prop3=sort(temp/sum(temp),decreasing=T)
  ind1=sort(temp/sum(temp),decreasing=T,index.return=T)$ix
  
  W=W[,1:5]
  X_hat1=W %*% H_Model1
  X_hat2=W %*% H_Model2
  X_hat3=W %*% H_Model3
  X_hat4=W %*% H_Model4
  X_hat5=W %*% H_Model5
  X_hat6=W %*% H_Model6
  X_hat7=W %*% H_Model7
  X_hat8=W %*% H_Model8
  X_hat9=W %*% H_Model9
  X_hat10=W %*% H_Model10
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  X_hat13=W %*% H_Model13
  X_hat99=(as.matrix(W_hat[,ind1[1:6]]) %*% diag(Thet_hat[ind1[1:6]]) + 0.01) %*% as.matrix(H_hat[ind1[1:6],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  cor_X10=sort(diag(cor(X,X_hat10)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X13=sort(diag(cor(X,X_hat13)),decreasing=T,na.last=T)
  cor_X99=sort(diag(cor(X,X_hat99)),decreasing=T,na.last=T)
  
  nrmse1=NULL
  nrmse2=NULL
  nrmse3=NULL
  nrmse4=NULL
  nrmse5=NULL
  nrmse6=NULL
  nrmse7=NULL
  nrmse8=NULL
  nrmse9=NULL
  nrmse10=NULL
  nrmse11=NULL
  nrmse12=NULL
  nrmse13=NULL
  nrmse99=NULL
  for(i in 1:S){
    nrmse1=c(nrmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    nrmse2=c(nrmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    nrmse3=c(nrmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    nrmse4=c(nrmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    nrmse5=c(nrmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    nrmse6=c(nrmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    nrmse7=c(nrmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    nrmse8=c(nrmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    nrmse9=c(nrmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
    nrmse10=c(nrmse10,sqrt(mean(((X[,i]-X_hat10[,i])/sd(X[,i]))^2)))
    nrmse11=c(nrmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    nrmse12=c(nrmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    nrmse13=c(nrmse13,sqrt(mean(((X[,i]-X_hat13[,i])/sd(X[,i]))^2)))
    nrmse99=c(nrmse99,sqrt(mean(((X[,i]-X_hat99[,i])/sd(X[,i]))^2)))
  }
  
  temp=as.data.frame(rbind(cbind(rmseh7,7),cbind(rmseh8,8),cbind(rmseh9,9), cbind(rmseh10,10), cbind(rmseh11,11), cbind(rmseh12,12), cbind(rmseh13,13)))
  colnames(temp)=c("rmse","model")
  temp$method=Method3
  df1=as.data.frame(rbind(df1,temp))
  
  
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out=1000),cor_H1),
                           cbind(seq(0,1,length.out=1000),cor_H7),
                           cbind(seq(0,1,length.out=1000),cor_H8),
                           cbind(seq(0,1,length.out=1000),cor_H9),
                           cbind(seq(0,1,length.out = 1000),cor_H10),
                           cbind(seq(0,1,length.out = 1000),cor_H11),
                           cbind(seq(0,1,length.out = 1000),cor_H12),
                           cbind(seq(0,1,length.out = 1000),cor_H13)
                           ))
  temp$model=c(# rep("MAIN",1000),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000)
               )
  # temp$model=c(rep("RETROFIT",1000),
  #              rep("RCTD",1000),
  #              rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
  #              rep("NMFreg",1000),rep("SPOTlight",1000))
  colnames(temp)=c("x","cor","model")
  temp$method=Method3
  df2=as.data.frame(rbind(df2,temp))
  
  temp=as.data.frame(rbind(# cbind(nrmse1,1),
                           # cbind(nrmse99,99),
                           cbind(nrmse7,7),cbind(nrmse8,8),
                           cbind(nrmse9,9),
                           cbind(nrmse10,10),cbind(nrmse11,11),
                           cbind(nrmse12,12),cbind(nrmse13,13)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=Method3
  df3=as.data.frame(rbind(df3,temp))
  
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out=1000),cor_X1),
                           # cbind(seq(0,1,length.out=1000),cor_X99),
                           cbind(seq(0,1,length.out=1000),cor_X7),
                           cbind(seq(0,1,length.out=1000),cor_X8),
                           cbind(seq(0,1,length.out=1000),cor_X9),
                           cbind(seq(0,1,length.out = 1000),cor_X10),
                           cbind(seq(0,1,length.out = 1000),cor_X11),
                           cbind(seq(0,1,length.out = 1000),cor_X12),
                           cbind(seq(0,1,length.out = 1000),cor_X13)
                           ))
  temp$model=c(# rep("MAIN",1000),
               # rep("PVE80%",1000),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000)
               )
  temp$method=Method3
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  label=c(
    format(ks.test(rmseh1,rmseh1, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh6, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh7, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh8, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh9, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh10, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh11, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh12, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh13, alternative="greater")$p.value, digits=3)
  )
  
  temp=data.frame(model=c(7,8,9,10,11,12,13),rmse=0.9,
                  method=factor(Method3,levels=c("N=20,M=5","N=10,M=3", Method3,"Missing 2 cell types")),
                  label=label[c(7,8,9,10,11,12,13)])
  # label=c("2.9e-35","3.5e-29","0.002","2.3e-150","1e-81"))
  ann.text1=as.data.frame(rbind(ann.text1,temp))
  
  cor_H6[is.na(cor_H6)]=0
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H13),digits=3)
  )
  
  temp=data.frame(label=c(# paste("MAIN:", label[[1]]),
                          paste(paste("LAMBDA=", lambdas[[1]], sep=""), label[[7]]),
                          paste(paste("LAMBDA=", lambdas[[2]], sep=""), label[[8]]),
                          paste(paste("LAMBDA=", lambdas[[3]], sep=""), label[[9]]),
                          paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                          paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                          paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                          paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                          ),
                  model=c(# "MAIN",
                          paste("LAMBDA=", lambdas[[1]], sep=""),
                          paste("LAMBDA=", lambdas[[2]], sep=""), 
                          paste("LAMBDA=", lambdas[[3]], sep=""),
                          paste("LAMBDA=", lambdas[[4]], sep=""),
                          paste("LAMBDA=", lambdas[[5]], sep=""),
                          paste("LAMBDA=", lambdas[[6]], sep=""),
                          paste("LAMBDA=", lambdas[[7]], sep="")
                  ),
                  method=factor(Method3,levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  
  label=c(
    format(ks.test(nrmse1,nrmse1, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse2, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse3, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse4, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse5, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse6, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse7, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse8, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse9, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse10, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse11, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse12, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse13, alternative="greater")$p.value, digits=3)
  )
  label99 = format(ks.test(nrmse1,nrmse99, alternative="greater")$p.value, digits=3)
  temp=data.frame(model=c(# 99,
                          7,8,9,10,11,12,13),rmse=0.9,
                  method=factor(Method3,levels=c("N=20,M=5","N=10,M=3",
                                                              Method3,"Missing 2 cell types")),
                  label=c(# label99,
                          label[[7]],label[[8]],label[[9]],label[[10]],label[[11]],label[[12]],label[[13]]))
  # label=c("9.9e-10","1.6e-08","0.294","5.3e-37","1.4e-31"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  cor_X6[is.na(cor_X6)]=0
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X13),digits=3)
  )
  label99 = round(AUC(x=seq(0,1,length.out = 1000),y=cor_X99),digits=3)
  
  temp=data.frame(label = c(# paste("MAIN:", label[[1]]),
                            # paste("PVE80%:", label99),
                            paste(paste("LAMBDA=", lambdas[[1]], sep=""), label[[7]]),
                            paste(paste("LAMBDA=", lambdas[[2]], sep=""), label[[8]]),
                            paste(paste("LAMBDA=", lambdas[[3]], sep=""), label[[9]]),
                            paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                            paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                            paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                            paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                            ),
                  model=c(# "MAIN",
                          # "PVE80%",
                          paste("LAMBDA=", lambdas[[1]], sep=""),
                          paste("LAMBDA=", lambdas[[2]], sep=""), 
                          paste("LAMBDA=", lambdas[[3]], sep=""),
                          paste("LAMBDA=", lambdas[[4]], sep=""),
                          paste("LAMBDA=", lambdas[[5]], sep=""),
                          paste("LAMBDA=", lambdas[[6]], sep=""),
                          paste("LAMBDA=", lambdas[[7]], sep="")
                          ),
                  method=factor(Method3,levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")))
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  
  
  #### missing 2
  H=read.csv("extra5_H.csv")[,-1]
  # H_mod1=read.csv("extra5_loc_H_retrofit.csv")[,-1]
  H_mod1=read.csv(paste(folder_reference, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod2=read.csv("miss2_loc_RCTD_H.csv")[,-1]
  H_mod3=read.csv("miss2_loc_RCTD-D_H.csv")[,-1]
  H_mod4=read.csv("miss2_loc_Stereoscope.csv")[,-1]
  H_mod5=t(read.csv("miss2_loc_NMFreg.csv")[,-11])
  H_mod6=read.csv("miss2_loc_SPOTlight.csv")[,-1]
  H_mod7=read.csv(paste(folder1, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod8=read.csv(paste(folder2, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod9=read.csv(paste(folder3, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod10=read.csv(paste(folder4, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod11=read.csv(paste(folder5, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod12=read.csv(paste(folder6, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  H_mod13=read.csv(paste(folder7, "/mapped/extra5_loc_H_mapped_cor.csv", sep=""))[,-1]
  
  
  X=read.csv("extra5_loc_X.csv")[,-1]
  
  H_mod1=H_mod1[1:5,]
  H_mod2=H_mod2[1:3,]
  H_mod3=H_mod3[1:3,]
  H_mod4=H_mod4[1:3,]
  H_mod5=H_mod5[1:3,]
  H_mod6=H_mod6[1:3,]
  H_mod7=H_mod7[1:5,]
  H_mod8=H_mod8[1:5,]
  H_mod9=H_mod9[1:5,]
  H_mod10=H_mod10[1:5,]
  H_mod11=H_mod11[1:5,]
  H_mod12=H_mod12[1:5,]
  H_mod13=H_mod13[1:5,]
  
  K=nrow(H)
  G=nrow(W)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H[1:3,],H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H[1:3,],H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H[1:3,],H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H[1:3,],H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H[1:3,],H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H,H_mod7)),decreasing=T,na.last=T)
  cor_H8=sort(diag(cor(H,H_mod8)),decreasing=T,na.last=T)
  cor_H9=sort(diag(cor(H,H_mod9)),decreasing=T,na.last=T)
  cor_H10=sort(diag(cor(H,H_mod10)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H,H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H,H_mod12)),decreasing=T,na.last=T)
  cor_H13=sort(diag(cor(H,H_mod13)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_Model9=matrix(NA, ncol=S, nrow=nrow(H_mod9))
  H_Model10=matrix(NA, ncol=S, nrow=nrow(H_mod10))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  H_Model13=matrix(NA, ncol=S, nrow=nrow(H_mod13))
  H_True=matrix(NA, ncol=S, nrow=nrow(H))
  
  rmseh1=NULL
  rmseh2=NULL
  rmseh3=NULL
  rmseh4=NULL
  rmseh5=NULL
  rmseh6=NULL
  rmseh7=NULL
  rmseh8=NULL
  rmseh9=NULL
  rmseh10=NULL
  rmseh11=NULL
  rmseh12=NULL
  rmseh13=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model8[,i]=H_mod8[,i]/sum(H_mod8[,i])
    H_Model9[,i]=H_mod9[,i]/sum(H_mod9[,i])
    H_Model10[,i]=H_mod10[,i]/sum(H_mod10[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    H_Model13[,i]=H_mod13[,i]/sum(H_mod13[,i])
    
    rmseh1=c(rmseh1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmseh2=c(rmseh2,sqrt(mean((H_True[1:3,i]-H_Model2[,i])^2)))
    rmseh3=c(rmseh3,sqrt(mean((H_True[1:3,i]-H_Model3[,i])^2)))
    rmseh4=c(rmseh4,sqrt(mean((H_True[1:3,i]-H_Model4[,i])^2)))
    rmseh5=c(rmseh5,sqrt(mean((H_True[1:3,i]-H_Model5[,i])^2)))
    rmseh6=c(rmseh6,sqrt(mean((H_True[1:3,i]-H_Model6[,i])^2)))
    rmseh7=c(rmseh7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmseh8=c(rmseh8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
    rmseh9=c(rmseh9,sqrt(mean((H_True[,i]-H_Model9[,i])^2)))
    rmseh10=c(rmseh10,sqrt(mean((H_True[,i]-H_Model10[,i])^2)))
    rmseh11=c(rmseh11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmseh12=c(rmseh12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
    rmseh13=c(rmseh13,sqrt(mean((H_True[,i]-H_Model13[,i])^2)))
  }
  
  X_hat=(as.matrix(W_hat) %*% diag(Thet_hat) + 0.01) %*% as.matrix(H_hat)
  temp=rep(NA,length(Thet_hat))
  for(k in 1:length(temp)){
    temp[k]=sum(as.matrix(W_hat[,k] * Thet_hat[k] + 0.01) %*% as.matrix(H_hat[k,]))
  }
  prop=sort(temp/sum(temp),decreasing=T)
  ind1=sort(temp/sum(temp),decreasing=T,index.return=T)$ix
  
  
  X_hat1=W %*% H_Model1
  X_hat2=W[,1:3] %*% H_Model2
  X_hat3=W[,1:3] %*% H_Model3
  X_hat4=W[,1:3] %*% H_Model4
  X_hat5=W[,1:3] %*% H_Model5
  X_hat6=W[,1:3] %*% H_Model6
  X_hat7=W %*% H_Model7
  X_hat8=W %*% H_Model8
  X_hat9=W %*% H_Model9
  X_hat10=W %*% H_Model10
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  X_hat13=W %*% H_Model13
  X_hat99=(as.matrix(W_hat[,ind1[1:6]]) %*% diag(Thet_hat[ind1[1:6]]) + 0.01) %*% as.matrix(H_hat[ind1[1:6],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  cor_X10=sort(diag(cor(X,X_hat10)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X13=sort(diag(cor(X,X_hat13)),decreasing=T,na.last=T)
  cor_X99=sort(diag(cor(X,X_hat99)),decreasing=T,na.last=T)
  
  nrmse1=NULL
  nrmse2=NULL
  nrmse3=NULL
  nrmse4=NULL
  nrmse5=NULL
  nrmse6=NULL
  nrmse7=NULL
  nrmse8=NULL
  nrmse9=NULL
  nrmse10=NULL
  nrmse11=NULL
  nrmse12=NULL
  nrmse13=NULL
  nrmse99=NULL
  for(i in 1:S){
    nrmse1=c(nrmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    nrmse2=c(nrmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    nrmse3=c(nrmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    nrmse4=c(nrmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    nrmse5=c(nrmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    nrmse6=c(nrmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    nrmse7=c(nrmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    nrmse8=c(nrmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    nrmse9=c(nrmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
    nrmse10=c(nrmse10,sqrt(mean(((X[,i]-X_hat10[,i])/sd(X[,i]))^2)))
    nrmse11=c(nrmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    nrmse12=c(nrmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    nrmse13=c(nrmse13,sqrt(mean(((X[,i]-X_hat13[,i])/sd(X[,i]))^2)))
    nrmse99=c(nrmse99,sqrt(mean(((X[,i]-X_hat99[,i])/sd(X[,i]))^2)))
  }
  
  temp=as.data.frame(rbind(cbind(rmseh7,7),cbind(rmseh8,8),cbind(rmseh9,9),cbind(rmseh10,10),cbind(rmseh11,11),cbind(rmseh12,12),cbind(rmseh13,13)))
  colnames(temp)=c("rmse","model")
  temp$method="Missing 2 cell types"
  df1=as.data.frame(rbind(df1,temp))
  
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
                           cbind(seq(0,1,length.out = 1000),cor_H7),
                           cbind(seq(0,1,length.out = 1000),cor_H8),
                           cbind(seq(0,1,length.out = 1000),cor_H9),
                           cbind(seq(0,1,length.out = 1000),cor_H10),
                           cbind(seq(0,1,length.out = 1000),cor_H11),
                           cbind(seq(0,1,length.out = 1000),cor_H12),
                           cbind(seq(0,1,length.out = 1000),cor_H13)
                           ))
  
  temp$model=c(# rep("MAIN",1000),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000)
               )
  colnames(temp)=c("x","cor","model")
  temp$method="Missing 2 cell types"
  df2=as.data.frame(rbind(df2,temp))
  
  temp=as.data.frame(rbind(# cbind(nrmse1,1),
                           # cbind(nrmse99,99),
                           cbind(nrmse7,7),cbind(nrmse8,8),
                           cbind(nrmse9,9),
                           cbind(nrmse10,10),cbind(nrmse11,11),
                           cbind(nrmse12,12),cbind(nrmse13,13)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="Missing 2 cell types"
  df3=as.data.frame(rbind(df3,temp))
  
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
                           # cbind(seq(0,1,length.out = 1000),cor_X99),
                           cbind(seq(0,1,length.out = 1000),cor_X7),
                           cbind(seq(0,1,length.out = 1000),cor_X8),
                           cbind(seq(0,1,length.out = 1000),cor_X9),
                           cbind(seq(0,1,length.out = 1000),cor_X10),
                           cbind(seq(0,1,length.out = 1000),cor_X11),
                           cbind(seq(0,1,length.out = 1000),cor_X12),
                           cbind(seq(0,1,length.out = 1000),cor_X13)))
  temp$model=c(# rep("MAIN",1000),
               # rep("PVE80%",1000),
               rep(paste("LAMBDA=", lambdas[[1]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[2]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[3]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[4]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[5]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[6]], sep=""),1000),
               rep(paste("LAMBDA=", lambdas[[7]], sep=""),1000)
               )
  temp$method="Missing 2 cell types"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  label=c(
    format(ks.test(rmseh1,rmseh1, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh6, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh7, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh8, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh9, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh10, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh11, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh12, alternative="greater")$p.value, digits=3),
    format(ks.test(rmseh1,rmseh13, alternative="greater")$p.value, digits=3)
  )
  
  temp=data.frame(model=c(7,8,9,10,11,12,13),rmse=0.9,
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                                Method3,"Missing 2 cell types")),
                  label=label[c(7,8,9,10,11,12,13)])
  # label=c("4e-4","0.265","5.8e-8","3.9e-7","1.3e-14"))
  ann.text1=as.data.frame(rbind(ann.text1,temp))
  
  
  cor_H1[is.na(cor_H1)]=0
  cor_H2[is.na(cor_H2)]=0
  cor_H3[is.na(cor_H3)]=0
  cor_H4[is.na(cor_H4)]=0
  cor_H5[is.na(cor_H5)]=0
  cor_H6[is.na(cor_H6)]=0
  cor_H7[is.na(cor_H7)]=0
  cor_H8[is.na(cor_H8)]=0
  cor_H9[is.na(cor_H9)]=0
  cor_H10[is.na(cor_H10)]=0
  cor_H11[is.na(cor_H11)]=0
  cor_H12[is.na(cor_H12)]=0
  cor_H13[is.na(cor_H13)]=0
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H13),digits=3)
  )
  
  temp=data.frame(label=c(# paste("MAIN:", label[[1]]),
                          paste(paste("LAMBDA=", lambdas[[1]], ":", sep=""), label[[7]]),
                          paste(paste("LAMBDA=", lambdas[[2]], ":", sep=""), label[[8]]),
                          paste(paste("LAMBDA=", lambdas[[3]], ":", sep=""), label[[9]]),
                          paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                          paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                          paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                          paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                          ),
                  model=c(# "MAIN",
                          paste("LAMBDA=", lambdas[[1]], sep=""),
                          paste("LAMBDA=", lambdas[[2]], sep=""),
                          paste("LAMBDA=", lambdas[[3]], sep=""),
                          paste("LAMBDA=", lambdas[[4]], sep=""),
                          paste("LAMBDA=", lambdas[[5]], sep=""),
                          paste("LAMBDA=", lambdas[[6]], sep=""),
                          paste("LAMBDA=", lambdas[[7]], sep="")
                          ),
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  
  label=c(
    format(ks.test(nrmse1,nrmse1, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse2, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse3, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse4, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse5, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse6, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse7, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse8, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse9, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse10, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse11, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse12, alternative="greater")$p.value, digits=3),
    format(ks.test(nrmse1,nrmse13, alternative="greater")$p.value, digits=3)
  )
  label99 = format(ks.test(nrmse1,nrmse99, alternative="greater")$p.value, digits=3)
  
  temp=data.frame(model=c(# 99,
                          7,8,9,10,11,12,13),rmse=0.9,
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                                Method3,"Missing 2 cell types")),
                  label=c(# label99,
                          label[[7]],label[[8]],label[[9]],label[[10]],label[[11]],label[[12]],label[[13]]))
  # label=c("0.171","0.047","5.4e-9","0.254","2e-4"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  cor_X1[is.na(cor_X1)]=0
  cor_X2[is.na(cor_X2)]=0
  cor_X3[is.na(cor_X3)]=0
  cor_X4[is.na(cor_X4)]=0
  cor_X5[is.na(cor_X5)]=0
  cor_X6[is.na(cor_X6)]=0
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X10),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X12),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X13),digits=3)
  )
  label99 = round(AUC(x=seq(0,1,length.out = 1000),y=cor_X99),digits=3)
  
  temp=data.frame(label = c(# paste("MAIN:", label[[1]]),
                            # paste("PVE80%:", label99),
                            paste(paste("LAMBDA=", lambdas[[1]], ":", sep=""), label[[7]]),
                            paste(paste("LAMBDA=", lambdas[[2]], ":", sep=""), label[[8]]),
                            paste(paste("LAMBDA=", lambdas[[3]], ":", sep=""), label[[9]]),
                            paste(paste("LAMBDA=", lambdas[[4]], sep=""), label[[10]]),
                            paste(paste("LAMBDA=", lambdas[[5]], sep=""), label[[11]]),
                            paste(paste("LAMBDA=", lambdas[[6]], sep=""), label[[12]]),
                            paste(paste("LAMBDA=", lambdas[[7]], sep=""), label[[13]])
                            ),
                  model=c(# "MAIN",
                          # "PVE80%",
                          paste("LAMBDA=", lambdas[[1]], sep=""),
                          paste("LAMBDA=", lambdas[[2]], sep=""), 
                          paste("LAMBDA=", lambdas[[3]], sep=""),
                          paste("LAMBDA=", lambdas[[4]], sep=""),
                          paste("LAMBDA=", lambdas[[5]], sep=""),
                          paste("LAMBDA=", lambdas[[6]], sep=""),
                          paste("LAMBDA=", lambdas[[7]], sep="")
                          ),
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types")))
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  
  
  df1$model=as.factor(df1$model)
  # df2$model=factor(df2$model,levels=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight"))
  df2$model=factor(df2$model,levels=c(# "MAIN",
                                      paste("LAMBDA=", lambdas[[1]], sep=""), 
                                      paste("LAMBDA=", lambdas[[2]], sep=""), 
                                      paste("LAMBDA=", lambdas[[3]], sep=""),
                                      paste("LAMBDA=", lambdas[[4]], sep=""),
                                      paste("LAMBDA=", lambdas[[5]], sep=""),
                                      paste("LAMBDA=", lambdas[[6]], sep=""),
                                      paste("LAMBDA=", lambdas[[7]], sep="")
                                      ))
  df1$method=factor(df1$method,levels=c("N=20,M=5","N=10,M=3", Method3,"Missing 2 cell types"))
  df2$method=factor(df2$method,levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types"))
  
  # df3$model=factor(df3$model,levels=c(1,7,2,3,4,5,6))
  df3$model=factor(df3$model,levels=c(# 99,
                                      7,8,9,10,11,12,13))
  # df4$model=factor(df4$model,levels=c("RETROFIT","RETROFIT (PVE=80%)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight"))
  df4$model=factor(df4$model,levels=c(# "MAIN",
                                      # "PVE80%",
                                      paste("LAMBDA=", lambdas[[1]], sep=""),
                                      paste("LAMBDA=", lambdas[[2]], sep=""), 
                                      paste("LAMBDA=", lambdas[[3]], sep=""),
                                      paste("LAMBDA=", lambdas[[4]], sep=""),
                                      paste("LAMBDA=", lambdas[[5]], sep=""),
                                      paste("LAMBDA=", lambdas[[6]], sep=""),
                                      paste("LAMBDA=", lambdas[[7]], sep="")
                                      ))
  df3$method=factor(df3$method,levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types"))
  df4$method=factor(df4$method,levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types"))
  
  ann.text1$model=as.factor(ann.text1$model)
  ann.text2$model=as.factor(ann.text2$model)
  ann.text3$model=as.factor(ann.text3$model)
  ann.text4$model=factor(ann.text4$model,levels=c("N=20,M=5","N=10,M=3",Method3,"Missing 2 cell types"))
  
  ann.text1$rmse=c(rep(1,7),rep(1,7),rep(1,7),rep(1,7))
  # ann.text2$cor=rep(c(0.55,0.45,0.35,0.25,0.15,0.05),4)
  ann.text2$cor=rep(c(0.70,0.60,0.50,0.40, 0.30, 0.20,0.10),4)
  ann.text2$x=0.05
  ann.text3$rmse=c(rep(3.5,7),rep(3.5,7),rep(3.5,7),rep(3.5,7))
  ann.text4$cor=rep(c(0.70,0.60,0.50,0.40, 0.30, 0.20,0.10),4)
  ann.text4$x=0.05
  
  p1=ggplot(df,aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + 
    geom_point(aes(color=model))+
    # xlab("Cell Types") +
    ylab(expression(paste("Correlation (",W^0,",",widetilde(W),")"))) +
    theme_bw()+
    scale_color_manual(values = c('grey30', "#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),
                       labels=c(# "MAIN",
                                paste("LAMBDA=", lambdas[[1]], sep=""), 
                                paste("LAMBDA=", lambdas[[2]], sep=""), 
                                paste("LAMBDA=", lambdas[[3]], sep=""),
                                paste("LAMBDA=", lambdas[[4]], sep=""),
                                paste("LAMBDA=", lambdas[[5]], sep=""),
                                paste("LAMBDA=", lambdas[[6]], sep=""),
                                paste("LAMBDA=", lambdas[[7]], sep="")
                                )) +
    ylim(0,1)+
    theme(legend.position="none", axis.title=element_text(size=9), axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm")) +
    facet_grid(cols=vars(method),scales = "free_x")
  
  ## plots for RMSE(H), Cor(H), RMSE(X), Cor(X)
  p2=ggplot(df1, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("RMSE (",H,",",widetilde(H),")"))) +
    xlab("") +
    geom_violin() + theme_bw()+ ylim(0,1)+
    # geom_boxplot(width=0.1)+
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(0, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c('gray30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"))+
    scale_x_discrete(labels=c(# "MAIN",
                              paste("LAMBDA=", lambdas[[1]], sep=""),
                              paste("LAMBDA=", lambdas[[2]], sep=""), 
                              paste("LAMBDA=", lambdas[[3]], sep=""),
                              paste("LAMBDA=", lambdas[[4]], sep=""),
                              paste("LAMBDA=", lambdas[[5]], sep=""),
                              paste("LAMBDA=", lambdas[[6]], sep=""),
                              paste("LAMBDA=", lambdas[[7]], sep="")
    )) +facet_grid(cols=vars(method)) + 
    geom_text(data = ann.text1,aes(label = label),size=2.3, color=rep("black", 7*4))
  
  p3=ggplot(df2, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + xlab("Percentile") +
    ylab(expression(paste("Correlation (",H,",",widetilde(H),")"))) +
    theme_bw()+
    scale_color_manual(values = c('gray30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),
                       labels=c(# "MAIN",
                                paste("LAMBDA=", lambdas[[1]], sep=""), 
                                paste("LAMBDA=", lambdas[[2]], sep=""), 
                                paste("LAMBDA=", lambdas[[3]], sep=""),
                                paste("LAMBDA=", lambdas[[4]], sep=""),
                                paste("LAMBDA=", lambdas[[5]], sep=""),
                                paste("LAMBDA=", lambdas[[6]], sep=""),
                                paste("LAMBDA=", lambdas[[7]], sep="")
    )) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.title=element_text(size=9),axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm")) +
    geom_text(data = ann.text2,aes(label = label),size=2.2,
              hjust=0, color=rep(c('gray30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),4))
  
  
  p4=ggplot(df3, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("NRMSE (",X,",",widetilde(X),")"))) +
    xlab("") +
    geom_violin() + theme_bw() + # ylim(0,3.5) +
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          axis.text.y = element_text(size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c(# 'gray35',
                               'grey20',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"))+
    scale_x_discrete(labels=c(# "MAIN",
                              # "PVE80%",
                              paste("LAMBDA=", lambdas[[1]], sep=""), 
                              paste("LAMBDA=", lambdas[[2]], sep=""), 
                              paste("LAMBDA=", lambdas[[3]], sep=""),
                              paste("LAMBDA=", lambdas[[4]], sep=""),
                              paste("LAMBDA=", lambdas[[5]], sep=""),
                              paste("LAMBDA=", lambdas[[6]], sep=""),
                              paste("LAMBDA=", lambdas[[7]], sep="")
    )) +facet_wrap(~method,scales="free",nrow=1,ncol=4)+
    geom_text(data = ann.text3,aes(label = label),size=2.3, color=rep("black", 7*4))
  # facet_wrap(~method,scales="free",nrow=1,ncol=4) +
  
  p5=ggplot(df4, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model,linetype=model)) + xlab("Percentile") +
    ylab(expression(paste("Correlation (",X,",",widetilde(X),")"))) +
    theme_bw()+
    scale_color_manual(values = c(# 'gray30',
                                  'grey30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),
                       labels=c(# "MAIN",
                                # 'PVE80%',
                                paste("LAMBDA=", lambdas[[1]], sep=""), 
                                paste("LAMBDA=", lambdas[[2]], sep=""), 
                                paste("LAMBDA=", lambdas[[3]], sep=""),
                                paste("LAMBDA=", lambdas[[4]], sep=""),
                                paste("LAMBDA=", lambdas[[5]], sep=""),
                                paste("LAMBDA=", lambdas[[6]], sep=""),
                                paste("LAMBDA=", lambdas[[7]], sep="")
    )) +
    scale_linetype_manual(values = c('solid','dashed','solid','solid', 'solid','solid','solid','solid'),
                          labels=c(# "MAIN",
                            # 'PVE80%',
                            paste("LAMBDA=", lambdas[[1]], sep=""), 
                            paste("LAMBDA=", lambdas[[2]], sep=""), 
                            paste("LAMBDA=", lambdas[[3]], sep=""),
                            paste("LAMBDA=", lambdas[[4]], sep=""),
                            paste("LAMBDA=", lambdas[[5]], sep=""),
                            paste("LAMBDA=", lambdas[[6]], sep=""),
                            paste("LAMBDA=", lambdas[[7]], sep="")
                          )) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.text.x = element_text(size=7),axis.title.x = element_text(vjust=10),
          #axis.title.x = element_text(vjust = 12),
          axis.title=element_text(size=9),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm")) +
    geom_text(data = ann.text4,aes(label = label),size=2.2,hjust=0,
              color=rep(c(# 'gray30',
                          'grey30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),4))
  
  plot_grid(
    p1,NULL,p2,NULL,p3,NULL, p4,NULL,p5,
    rel_heights = c(1.5, -0.25, 1, -0.05, 1.2, -0.1, 1, -0.05, 1.2),
    align='hv',
    labels=c("A","","B","","C","","D","","E"),# vjust=12,
    #label_colour="red",
    nrow=9, ncol=1)
  
  if(save){
    print('ggsave')
    ggsave(path_grid, width=35, height=45, units="cm")
  }
  
  ################################################################################
  ################################################################################
  # method3 only
  ################################################################################
  ################################################################################
  
  df = df[df$method == Method3, ]
  df1 = df1[df1$method == Method3, ]
  df2 = df2[df2$method == Method3, ]
  df3 = df3[df3$method == Method3, ]
  df4 = df4[df4$method == Method3, ]
  
  ann.text1 = ann.text1[ann.text1$method == Method3, ]
  ann.text2 = ann.text2[ann.text2$method == Method3, ]
  ann.text3 = ann.text3[ann.text3$method == Method3, ]
  ann.text4 = ann.text4[ann.text4$method == Method3, ]
  
  
  p1=ggplot(df,aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + 
    geom_point(aes(color=model))+
    xlab("Cell Types") +
    ylab(expression(paste("Correlation (",W^0,",",widetilde(W),")"))) +
    theme_bw()+
    scale_color_manual(values = c('grey30', "#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),
                       labels=c(# "MAIN",
                         paste("LAMBDA=", lambdas[[1]], sep=""), 
                         paste("LAMBDA=", lambdas[[2]], sep=""), 
                         paste("LAMBDA=", lambdas[[3]], sep=""),
                         paste("LAMBDA=", lambdas[[4]], sep=""),
                         paste("LAMBDA=", lambdas[[5]], sep=""),
                         paste("LAMBDA=", lambdas[[6]], sep=""),
                         paste("LAMBDA=", lambdas[[7]], sep="")
                       )) +
    ylim(0,1)+
    theme(legend.position="none", axis.title=element_text(size=9), axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm"))
  
  ## plots for RMSE(H), Cor(H), RMSE(X), Cor(X)
  p2=ggplot(df1, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("RMSE (",H,",",widetilde(H),")"))) +
    xlab("") +
    geom_violin() + theme_bw()+ ylim(0,1)+
    # geom_boxplot(width=0.1)+
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(0, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c('gray30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"))+
    scale_x_discrete(labels=c(# "MAIN",
      paste("LAMBDA=", lambdas[[1]], sep=""),
      paste("LAMBDA=", lambdas[[2]], sep=""), 
      paste("LAMBDA=", lambdas[[3]], sep=""),
      paste("LAMBDA=", lambdas[[4]], sep=""),
      paste("LAMBDA=", lambdas[[5]], sep=""),
      paste("LAMBDA=", lambdas[[6]], sep=""),
      paste("LAMBDA=", lambdas[[7]], sep="")
    )) +
    geom_text(data = ann.text1,aes(label = label),size=2.3, color=rep("black", 7*1))
  
  p3=ggplot(df2, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + xlab("Normalized Rank") +
    ylab(expression(paste("Correlation (",H,",",widetilde(H),")"))) +
    theme_bw()+
    scale_color_manual(values = c('gray30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),
                       labels=c(# "MAIN",
                         paste("LAMBDA=", lambdas[[1]], sep=""), 
                         paste("LAMBDA=", lambdas[[2]], sep=""), 
                         paste("LAMBDA=", lambdas[[3]], sep=""),
                         paste("LAMBDA=", lambdas[[4]], sep=""),
                         paste("LAMBDA=", lambdas[[5]], sep=""),
                         paste("LAMBDA=", lambdas[[6]], sep=""),
                         paste("LAMBDA=", lambdas[[7]], sep="")
                       )) +
    ylim(0,1)+
    theme(legend.position = "none",strip.text = element_blank(),
          axis.title=element_text(size=9),axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm")) +
    geom_text(data = ann.text2,aes(label = label),size=2.2,
              hjust=0, color=rep(c('gray30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),1))
  
  
  p4=ggplot(df3, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("NRMSE (",X,",",widetilde(X),")"))) +
    xlab("") +
    geom_violin() + theme_bw() + ylim(0,3.5) +
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          axis.text.y = element_text(size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c(# 'gray35',
      'grey20',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"))+
    scale_x_discrete(labels=c(# "MAIN",
      # "PVE80%",
      paste("LAMBDA=", lambdas[[1]], sep=""), 
      paste("LAMBDA=", lambdas[[2]], sep=""), 
      paste("LAMBDA=", lambdas[[3]], sep=""),
      paste("LAMBDA=", lambdas[[4]], sep=""),
      paste("LAMBDA=", lambdas[[5]], sep=""),
      paste("LAMBDA=", lambdas[[6]], sep=""),
      paste("LAMBDA=", lambdas[[7]], sep="")
    )) +
    geom_text(data = ann.text3,aes(label = label),size=2.3, color=rep("black", 7*1))
  # facet_wrap(~method,scales="free",nrow=1,ncol=4) +
  
  p5=ggplot(df4, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model,linetype=model)) + xlab("Normalized Rank") +
    ylab(expression(paste("Correlation (",X,",",widetilde(X),")"))) +
    theme_bw()+
    scale_color_manual(values = c(# 'gray30',
      'grey30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),
      labels=c(# "MAIN",
        # 'PVE80%',
        paste("LAMBDA=", lambdas[[1]], sep=""), 
        paste("LAMBDA=", lambdas[[2]], sep=""), 
        paste("LAMBDA=", lambdas[[3]], sep=""),
        paste("LAMBDA=", lambdas[[4]], sep=""),
        paste("LAMBDA=", lambdas[[5]], sep=""),
        paste("LAMBDA=", lambdas[[6]], sep=""),
        paste("LAMBDA=", lambdas[[7]], sep="")
      )) +
    scale_linetype_manual(values = c('solid','dashed','solid','solid', 'solid','solid','solid','solid'),
                          labels=c(# "MAIN",
                            # 'PVE80%',
                            paste("LAMBDA=", lambdas[[1]], sep=""), 
                            paste("LAMBDA=", lambdas[[2]], sep=""), 
                            paste("LAMBDA=", lambdas[[3]], sep=""),
                            paste("LAMBDA=", lambdas[[4]], sep=""),
                            paste("LAMBDA=", lambdas[[5]], sep=""),
                            paste("LAMBDA=", lambdas[[6]], sep=""),
                            paste("LAMBDA=", lambdas[[7]], sep="")
                          )) +
    ylim(0,1)+
    theme(legend.position = "none",strip.text = element_blank(),
          axis.text.x = element_text(size=7),axis.title.x = element_text(vjust=10),
          #axis.title.x = element_text(vjust = 12),
          axis.title=element_text(size=9),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm")) +
    geom_text(data = ann.text4,aes(label = label),size=2.2,hjust=0,
              color=rep(c(# 'gray30',
                'grey30',"#D55E00",'#56B4E9','#0072B2', '#009E73', "#CC79A7", "#4B0092"),1))
  
  plot_grid(
    p2,NULL,p3,NULL, p4,NULL,p5,NULL,p1,
    rel_widths = c(1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1),
    # rel_heights = c(1.5, -0.25, 1, -0.05, 1.2, -0.1, 1, -0.05, 1.2),
    align='hv',
    # labels=c("A","","B","","C","","D","","E"),# vjust=12,
    #label_colour="red",
    nrow=1, ncol=9)
  
  if(save){
    print('ggsave')
    ggsave(path_grid_col3, width=55, height=10, units="cm")
  }
}