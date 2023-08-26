### Plots for paper

## Simulation Plots
RetrofitPlotComparison <- function(dir1, dir2, dir3, dir4, file, save=FALSE){
  path_grid = paste("../", file, ".pdf", sep="")
  original_dir=getwd()
  setwd(dir1)
  library(ggplot2)
  library(cowplot)
  library(matrixStats)
  library(DescTools)
  ### Read synthetic data
  
  ### plot for cor(W)
  W=as.matrix(read.csv(paste(dir1, "datasets/HippoI_W.csv", sep="/"),   row.names = 1, check.names = FALSE))
  # Mapped W_hat from RETROFIT (no. of genes x no. of cell types)
  W_mod1=read.csv(paste(dir1, "mapped/HippoI__map_cor_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  W_mod2=read.csv(paste(dir2, "mapped/HippoI__map_cor_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  W_mod3=read.csv(paste(dir3, "mapped/HippoI__map_cor_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  W_mod4=read.csv(paste(dir4, "mapped/HippoI__map_cor_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  
  cor_W1=sort(diag(cor(W,W_mod1)),decreasing=T,na.last=T)
  cor_W2=sort(diag(cor(W,W_mod2)),decreasing=T,na.last=T)
  cor_W3=sort(diag(cor(W,W_mod3)),decreasing=T,na.last=T)
  cor_W4=sort(diag(cor(W,W_mod4)),decreasing=T,na.last=T)
  
  
  df=as.data.frame(c(cor_W1,cor_W1,cor_W1,cor_W1))
  df$x=as.factor(c(1:length(cor_W1),1:length(cor_W1),
                   1:length(cor_W1),1:length(cor_W1)))
  df$method=factor(c(rep("N=20,M=5 (Exact 10 cell types)",length(cor_W1)),
                     rep("N=10,M=3 (Exact 10 cell types)",length(cor_W1)),
                     rep("N=10,M=3 (Extra 5 cell types)",length(cor_W1)),
                     rep("N=10,M=3 (Missing 2 cell types)",length(cor_W1))),
                   c("N=20,M=5 (Exact 10 cell types)","N=10,M=3 (Exact 10 cell types)",
                     "N=10,M=3 (Extra 5 cell types)",
                     "N=10,M=3 (Missing 2 cell types)"))
  colnames(df)=c("cor","x","method")
  
  
  ## Plot for RMSE(H)
  X=read.csv(paste(dir1, "datasets/HippoI_X.csv", sep="/"), row.names = 1, check.names = FALSE)
  H=read.csv(paste(dir1, "datasets/HippoI_H.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod1=read.csv(paste(dir1, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod2=read.csv(paste(dir2, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod3=read.csv(paste(dir3, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod4=read.csv(paste(dir4, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  
  W=as.matrix(read.csv(paste(dir1, "datasets/HippoI_W.csv", sep="/"),   row.names = 1, check.names = FALSE))
  
  K=nrow(H)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_True=matrix(NA, ncol=S, nrow=nrow(H))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
  }
  
  # W_hat, H_hat and Theta_hat from retrofit (with all components before mapping)
  W_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  Thet_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_th.csv", sep="/"))[,-1]
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
  X_hat80=(as.matrix(W_hat[,ind1[1:12]]) %*% diag(Thet_hat[ind1[1:12]]) +
             0.01) %*% as.matrix(H_hat[ind1[1:12],]) #no. of components that explain about 80%
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  
  #B1
  df1=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
                          cbind(rmse_Model3,3),cbind(rmse_Model4,4)
  ))
  df1$method="N=20,M=5"
  colnames(df1)=c("rmse","model","method")
  
  #C1
  df2=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H2),
    cbind(seq(0,1,length.out = 1000),cor_H3),
    cbind(seq(0,1,length.out = 1000),cor_H4)))
  df2$model=c(rep("m1",1000),rep("m2",1000),
              rep("m3",1000),rep("m4",1000))
  df2$method="N=20,M=5"
  colnames(df2)=c("x","cor","model","method")
  
  #D1
  df3=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse1,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),
    cbind(rmse4,4)
  ))
  colnames(df3)=c("rmse","model")
  df3$method="N=20,M=5"
  
  #E1
  df4=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X80),
    cbind(seq(0,1,length.out = 1000),cor_X2),
    cbind(seq(0,1,length.out = 1000),cor_X3),
    cbind(seq(0,1,length.out = 1000),cor_X4)))
  df4$model=c(rep("m1",1000),
              rep("m1 (PVE=80%)",1000),rep("m2",1000),
              rep("m3",1000),rep("m4",1000))
  df4$method="N=20,M=5"
  colnames(df4)=c("x","cor","model","method")
  
  # B1
  digit_length=2
  label = c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=digit_length),
    # format(ks.test(rmse_Model1,rmse_Model11, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=digit_length)
  )
  
  ## add p-values from above KS tests
  ann.text1=data.frame(model=2:4,rmse=0.9,
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",
                                                         "Extra 5 cell types","Missing 2 cell types")),
                       label=label)
  # label=c("1.7e-207","8.4e-113","3.1e-10","1.8e-281","2.4e-121"))
  
  # C1
  digit_length=3
  label = c(
    # round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length)
  )
  
  H_mod1=read.csv(paste(dir1, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod2=read.csv(paste(dir2, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod3=read.csv(paste(dir3, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod4=read.csv(paste(dir4, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  
  ## add AUC from above results
  ann.text2=data.frame(label=c(paste("m1:", label[[1]]),
                               paste("m2:", label[[2]]),
                               paste("m3:", label[[3]]),
                               paste("m4:", label[[4]])),
                       # label=c("RETROFIT: 0.952","RCTD: 0.651",
                       #         "RCTD-D: 0.810","Stereoscope: 0.939","NMFreg: 0.628",
                       #         "SPOTlight: 0.762"),
                       model=c("m1","m2","m3","m4"),
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",
                                                         "Extra 5 cell types","Missing 2 cell types")))
  
  # D1
  digit_length=2
  label = c(
    format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=digit_length)
  )
  
  ann.text3=data.frame(model=2:4,rmse=0.9,
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",
                                                         "Extra 5 cell types","Missing 2 cell types")),
                       label=label)
  # label=c("5.1e-48","5.5e-10","0.186","5.9e-86","1.6e-28"))#,"2.1e-28"))
  
  # E1
  digit_length=3
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length)
  )
  ann.text4=data.frame(label = c(paste("m1:", label[[1]]),
                                 paste("m1 (PVE=80%):", label[[2]]),
                                 paste("m2:", label[[3]]),
                                 paste("m3:", label[[4]]),
                                 paste("m4:", label[[5]])),
                       # label=c("RETROFIT: 0.966","RETROFIT (PVE=80%): 0.912",
                       #         "RCTD: 0.821",
                       #         "RCTD-D: 0.850","Stereoscope: 0.973",
                       #         "NMFreg: 741","SPOTlight: 0.873"),
                       model=c("m1","m1 (PVE=80%)","m2","m3","m4"),
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",
                                                         "Extra 5 cell types","Missing 2 cell types")))
  
  
  # N=10,M=3
  X=read.csv(paste(dir1, "datasets/HippoI_X.csv", sep="/"), row.names = 1, check.names = FALSE)
  H=read.csv(paste(dir1, "datasets/HippoI_H.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod1=read.csv(paste(dir1, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod2=read.csv(paste(dir2, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod3=read.csv(paste(dir3, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod4=read.csv(paste(dir4, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  
  K=nrow(H)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
  }
  
  # prop variance explained
  W_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  Thet_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_th.csv", sep="/"))[,-1]
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
  X_hat80=(as.matrix(W_hat[,ind1[1:13]]) %*% diag(Thet_hat[ind1[1:13]]) + 0.01) %*% as.matrix(H_hat[ind1[1:13],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  
  # B2
  temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
                           cbind(rmse_Model3,3),cbind(rmse_Model4,4)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="N=10,M=3"
  df1=as.data.frame(rbind(df1,temp))
  
  # C2
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H2),
    cbind(seq(0,1,length.out = 1000),cor_H3),
    cbind(seq(0,1,length.out = 1000),cor_H4)))
  temp$model=c(rep("m1",1000),rep("m2",1000),
               rep("m3",1000),rep("m4",1000))
  colnames(temp)=c("x","cor","model")
  temp$method="N=10,M=3"
  df2=as.data.frame(rbind(df2,temp))
  
  # D2
  temp=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse1,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),cbind(rmse4,4)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="N=10,M=3"
  df3=as.data.frame(rbind(df3,temp))
  
  # E2
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X80),
    cbind(seq(0,1,length.out = 1000),cor_X2),
    cbind(seq(0,1,length.out = 1000),cor_X3),
    cbind(seq(0,1,length.out = 1000),cor_X4)))
  temp$model=c(rep("m1",1000),
               rep("m1 (PVE=80%)",1000),rep("m2",1000),
               rep("m3",1000),rep("m4",1000))
  temp$method="N=10,M=3"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  # B2
  digit_length=2
  label = c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=digit_length)
  )
  
  temp=data.frame(model=2:4,rmse=0.9,
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("1.0e-4","0.003","0.975","2.0e-242","6.6e-83"))
  ann.text1=as.data.frame(rbind(ann.text1,temp))
  
  # C2
  digit_length=3
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length)
  )
  temp=data.frame(label = c(paste("m1:", label[[1]]),
                            paste("m2:", label[[2]]),
                            paste("m3:", label[[3]]),
                            paste("m4:", label[[4]])),
                  # label=c("RETROFIT: 0.954","RCTD: 0.926",
                  #         "RCTD-D: 0.937","Stereoscope: 0.958",
                  #         "NMFreg: 0.764",
                  #         "SPOTlight: 0.769"),
                  model=c("m1","m2","m3","m4"),
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  
  # D2
  digit_length=2
  label=c(
    format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=digit_length)
  )
  
  
  temp=data.frame(model=2:4,rmse=0.9,
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("0.171","4.0e-4","0.294","0.844","0.254"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  # E2
  digit_length=3
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length)
  )
  temp=data.frame(label = c(paste("m1:", label[[1]]),
                            paste("m1 (PVE=80%):", label[[2]]),
                            paste("m2:", label[[3]]),
                            paste("m3:", label[[4]]),
                            paste("m4:", label[[5]])),
                  # label=c("RETROFIT: 0.968",
                  #         "RETROFIT (PVE=80%): 0.933","RCTD: 0.930",
                  #         "RCTD-D: 0.907","Stereoscope: 0.946",
                  #         "NMFreg: 0.802","SPOTlight: 0.832"),
                  model=c("m1","m1 (PVE=80%)","m2","m3","m4"),
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")))
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  
  ## extra 5
  X=read.csv(paste(dir1, "datasets/HippoI_X.csv", sep="/"), row.names = 1, check.names = FALSE)
  H=read.csv(paste(dir1, "datasets/HippoI_H.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod1=read.csv(paste(dir1, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod2=read.csv(paste(dir2, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod3=read.csv(paste(dir3, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod4=read.csv(paste(dir4, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  
  K=nrow(H)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
  }
  
  # prop variance explained
  W_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  Thet_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_th.csv", sep="/"))[,-1]
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
  X_hat80=(as.matrix(W_hat[,ind1[1:13]]) %*% diag(Thet_hat[ind1[1:13]]) + 0.01) %*% as.matrix(H_hat[ind1[1:13],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  
  # B2
  temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
                           cbind(rmse_Model3,3),cbind(rmse_Model4,4)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="Extra 5 cell types"
  df1=as.data.frame(rbind(df1,temp))
  
  # C2
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H2),
    cbind(seq(0,1,length.out = 1000),cor_H3),
    cbind(seq(0,1,length.out = 1000),cor_H4)))
  temp$model=c(rep("m1",1000),rep("m2",1000),
               rep("m3",1000),rep("m4",1000))
  colnames(temp)=c("x","cor","model")
  temp$method="Extra 5 cell types"
  df2=as.data.frame(rbind(df2,temp))
  
  # D2
  temp=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse1,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),cbind(rmse4,4)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="Extra 5 cell types"
  df3=as.data.frame(rbind(df3,temp))
  
  # E2
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X80),
    cbind(seq(0,1,length.out = 1000),cor_X2),
    cbind(seq(0,1,length.out = 1000),cor_X3),
    cbind(seq(0,1,length.out = 1000),cor_X4)))
  temp$model=c(rep("m1",1000),
               rep("m1 (PVE=80%)",1000),rep("m2",1000),
               rep("m3",1000),rep("m4",1000))
  temp$method="Extra 5 cell types"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  # B2
  digit_length=2
  label = c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=digit_length)
  )
  
  temp=data.frame(model=2:4,rmse=0.9,
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("1.0e-4","0.003","0.975","2.0e-242","6.6e-83"))
  ann.text1=as.data.frame(rbind(ann.text1,temp))
  
  # C2
  digit_length=3
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length)
  )
  temp=data.frame(label = c(paste("m1:", label[[1]]),
                            paste("m2:", label[[2]]),
                            paste("m3:", label[[3]]),
                            paste("m4:", label[[4]])),
                  # label=c("RETROFIT: 0.954","RCTD: 0.926",
                  #         "RCTD-D: 0.937","Stereoscope: 0.958",
                  #         "NMFreg: 0.764",
                  #         "SPOTlight: 0.769"),
                  model=c("m1","m2","m3","m4"),
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  
  # D2
  digit_length=2
  label=c(
    format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=digit_length)
  )
  
  
  temp=data.frame(model=2:4,rmse=0.9,
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("0.171","4.0e-4","0.294","0.844","0.254"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  # E2
  digit_length=3
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length)
  )
  temp=data.frame(label = c(paste("m1:", label[[1]]),
                            paste("m1 (PVE=80%):", label[[2]]),
                            paste("m2:", label[[3]]),
                            paste("m3:", label[[4]]),
                            paste("m4:", label[[5]])),
                  # label=c("RETROFIT: 0.968",
                  #         "RETROFIT (PVE=80%): 0.933","RCTD: 0.930",
                  #         "RCTD-D: 0.907","Stereoscope: 0.946",
                  #         "NMFreg: 0.802","SPOTlight: 0.832"),
                  model=c("m1","m1 (PVE=80%)","m2","m3","m4"),
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")))
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  
  
  #### missing 2
  X=read.csv(paste(dir1, "datasets/HippoI_X.csv", sep="/"), row.names = 1, check.names = FALSE)
  H=read.csv(paste(dir1, "datasets/HippoI_H.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod1=read.csv(paste(dir1, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod2=read.csv(paste(dir2, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod3=read.csv(paste(dir3, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_mod4=read.csv(paste(dir4, "mapped/HippoI__map_cor_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  
  K=nrow(H)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
  }
  
  # prop variance explained
  W_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_w.csv", sep="/"),  row.names = 1, check.names = FALSE)
  H_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_h.csv", sep="/"),  row.names = 1, check.names = FALSE)
  Thet_hat=read.csv(paste(dir1, "decomposed/HippoI__decomp_th.csv", sep="/"))[,-1]
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
  X_hat80=(as.matrix(W_hat[,ind1[1:13]]) %*% diag(Thet_hat[ind1[1:13]]) + 0.01) %*% as.matrix(H_hat[ind1[1:13],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  
  # B2
  temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
                           cbind(rmse_Model3,3),cbind(rmse_Model4,4)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="Missing 2 cell types"
  df1=as.data.frame(rbind(df1,temp))
  
  # C2
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H2),
    cbind(seq(0,1,length.out = 1000),cor_H3),
    cbind(seq(0,1,length.out = 1000),cor_H4)))
  temp$model=c(rep("m1",1000),rep("m2",1000),
               rep("m3",1000),rep("m4",1000))
  colnames(temp)=c("x","cor","model")
  temp$method="Missing 2 cell types"
  df2=as.data.frame(rbind(df2,temp))
  
  # D2
  temp=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse1,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),cbind(rmse4,4)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="Missing 2 cell types"
  df3=as.data.frame(rbind(df3,temp))
  
  # E2
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X80),
    cbind(seq(0,1,length.out = 1000),cor_X2),
    cbind(seq(0,1,length.out = 1000),cor_X3),
    cbind(seq(0,1,length.out = 1000),cor_X4)))
  temp$model=c(rep("m1",1000),
               rep("m1 (PVE=80%)",1000),rep("m2",1000),
               rep("m3",1000),rep("m4",1000))
  temp$method="Missing 2 cell types"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  # B2
  digit_length=2
  label = c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=digit_length)
  )
  
  temp=data.frame(model=2:4,rmse=0.9,
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("1.0e-4","0.003","0.975","2.0e-242","6.6e-83"))
  ann.text1=as.data.frame(rbind(ann.text1,temp))
  
  # C2
  digit_length=3
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length)
  )
  temp=data.frame(label = c(paste("m1:", label[[1]]),
                            paste("m2:", label[[2]]),
                            paste("m3:", label[[3]]),
                            paste("m4:", label[[4]])),
                  # label=c("RETROFIT: 0.954","RCTD: 0.926",
                  #         "RCTD-D: 0.937","Stereoscope: 0.958",
                  #         "NMFreg: 0.764",
                  #         "SPOTlight: 0.769"),
                  model=c("m1","m2","m3","m4"),
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  
  # D2
  digit_length=2
  label=c(
    format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=digit_length),
    format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=digit_length)
  )
  
  
  temp=data.frame(model=2:4,rmse=0.9,
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("0.171","4.0e-4","0.294","0.844","0.254"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  # E2
  digit_length=3
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length)
  )
  temp=data.frame(label = c(paste("m1:", label[[1]]),
                            paste("m1 (PVE=80%):", label[[2]]),
                            paste("m2:", label[[3]]),
                            paste("m3:", label[[4]]),
                            paste("m4:", label[[5]])),
                  # label=c("RETROFIT: 0.968",
                  #         "RETROFIT (PVE=80%): 0.933","RCTD: 0.930",
                  #         "RCTD-D: 0.907","Stereoscope: 0.946",
                  #         "NMFreg: 0.802","SPOTlight: 0.832"),
                  model=c("m1","m1 (PVE=80%)","m2","m3","m4"),
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")))
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  
  # B-
  df1$model=as.factor(df1$model)
  # C-
  df2$model=factor(df2$model,levels=c("m1","m2","m3","m4"))
  # B-
  df1$method=factor(df1$method,levels=c("N=20,M=5","N=10,M=3", "Extra 5 cell types","Missing 2 cell types"))
  # C-
  df2$method=factor(df2$method,levels=c("N=20,M=5","N=10,M=3", "Extra 5 cell types","Missing 2 cell types"))
  # D-
  df3$model=factor(df3$model,levels=c(1,8,2,3,4))
  # E-
  df4$model=factor(df4$model,levels=c("m1","m1 (PVE=80%)","m2","m3","m4"))
  # D-
  df3$method=factor(df3$method,levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types"))
  # E-
  df4$method=factor(df4$method,levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types"))
  # B-
  ann.text1$model=as.factor(ann.text1$model)
  # C-
  ann.text2$model=as.factor(ann.text2$model)
  # D-
  ann.text3$model=as.factor(ann.text3$model)
  # E-
  ann.text4$model=factor(ann.text4$model,levels=c("m1","m1 (PVE=80%)","m2","m3","m4"))
  # C-
  ann.text2$cor=rep(c(0.65,0.55,0.45,0.35),4)
  ann.text2$x=0.05
  
  # D-
  cnt = 3
  ann.text3$rmse=c(rep(3,cnt),rep(58,cnt),rep(40,cnt),rep(15,cnt))
  
  # E-
  ann.text4$cor=rep(c(0.75,0.65,0.55,0.45,0.35),4)
  ann.text4$x=0.05
  
  # plot for cor(W)
  # A-
  p1=ggplot(df,aes(x=x,y=cor)) +
    geom_line(aes(x=x,y=cor,group=method)) + geom_point(aes(x=x,y=cor,group=method))+
    xlab("Cell Types") +
    ylab(expression(paste("Correlation (",W^0,",",widetilde(W),")"))) +
    theme_bw()+
    ylim(0,1)+
    theme(axis.title=element_text(size=9), axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm")) +
    facet_grid(cols=vars(method),scales = "free_x")
  
  ## plots for RMSE(H), Cor(H), RMSE(X), Cor(X)
  # B-
  p2=ggplot(df1, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("RMSE (",H,",",widetilde(H),")"))) +
    xlab("") +
    geom_violin() + theme_bw()+ ylim(0,1)+
    # geom_boxplot(width=0.1)+
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(0, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c('gray30','#56B4E9','#0072B2','#009E73'))+
    scale_x_discrete(labels=c("m1","m2","m3","m4"))+
    facet_grid(cols=vars(method))+
    geom_text(data = ann.text1,aes(label = label),size=1.8, color=rep(c("black"),3*4))
  
  # C-
  p3=ggplot(df2, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + xlab("Percentile") +
    ylab(expression(paste("Correlation (",H,",",widetilde(H),")"))) +
    theme_bw()+
    scale_color_manual(values = c('gray30','#56B4E9','#0072B2','#009E73'),
                       labels=c("m1","m2","m3","m4")) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.title=element_text(size=9),axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    geom_text(data = ann.text2,aes(label = label),size=2.2,
              hjust=0,color=rep(c('gray30','#56B4E9','#0072B2','#009E73'),4))
  
  # D-
  # p4=ggplot(df3, aes(x=model, y=rmse, fill=model)) +
  #   ylab(expression(paste("NRMSE (",X,",",widetilde(X),")"))) +
  #   xlab("") +
  #   geom_violin() + theme_bw()+ #ylim(0,60)+
  #   # geom_boxplot(width=0.1)+
  #   theme(legend.position = "none",axis.title=element_text(size=9),
  #         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
  #         strip.text = element_blank(),
  #         plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
  #   scale_fill_manual(values=c('gray35','grey20','#56B4E9','#0072B2',
  #                              '#009E73',"#D55E00", "#CC79A7", "#E69F00"))+
  #   scale_x_discrete(labels=c("RETROFIT","RETROFIT\n(PVE=80%)","RCTD","RCTD-D","Stereoscope","NMFreg","SPOTlight","STdeconvolve"
  #   )) +facet_wrap(~method,scales="free",nrow=1,ncol=4)+
  #   geom_text(data = ann.text3,aes(label = label),size=2.2, color=rep(c("black"),6*4))
  
  # E-
  p5=ggplot(df4, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model,linetype=model)) + xlab("Percentile") +
    ylab(expression(paste("Correlation (",X,",",widetilde(X),")"))) +
    theme_bw()+
    scale_color_manual(values = c('gray30','grey30','#56B4E9','#0072B2','#009E73'),
                       labels=c("m1","m1\n(PVE=80%)","m2","m3","m4")) +
    scale_linetype_manual(values = c('solid','dashed','solid','solid', 'solid'),
                          labels=c("m1","m1\n(PVE=80%)","m2","m3","m4")) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.text.x = element_text(size=7),axis.title.x = element_text(vjust=10),
          #axis.title.x = element_text(vjust = 12),
          axis.title=element_text(size=9),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    geom_text(data = ann.text4,aes(label = label),size=2.2,hjust=0,
              color=rep(c('gray35','grey20','#56B4E9','#0072B2', '#009E73'),4))
  
  ## plot for supplementary
  # plot_grid(
  #   p1,NULL,p2,NULL,p3,NULL, p4,NULL,p5,
  #   rel_heights = c(1, -0.25, 1, -0.05, 1, -0.1, 1, -0.05, 1),
  #   align='hv',
  #   labels=c("A","","B","","C","","D","","E"),# vjust=12,
  #   #label_colour="red",
  #   nrow=9, ncol=1)
  
  
  ## plot for main paper
  quantile(df3$rmse[df3$method=="N=10,M=3"],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df3$rmse[df3$method=="Extra 5 cell types"],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df3$rmse[df3$method=="Missing 2 cell types"],c(0.75,0.8,0.9), na.rm = TRUE)
  
  # D-
  ann.text3$rmse=rep(3.2,3*4)
  df5=df3[complete.cases(df3),]
  quantile(df5$rmse[df5$method=="N=10,M=3"],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df5$rmse[df5$method=="Extra 5 cell types"],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df5$rmse[df5$method=="Missing 2 cell types"],c(0.75,0.8,0.9), na.rm = TRUE)
  
  df5=df5[df5$rmse<3,]
  p4=ggplot(df5, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("NRMSE (",X,",",hat(X),")"))) +
    xlab("") +
    geom_violin() + theme_bw()+ ylim(0,3.5)+
    # geom_boxplot(width=0.1)+
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                     size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c('gray35','grey20','#56B4E9','#0072B2', '#009E73'))+
    scale_x_discrete(labels=c("m1","m1\n(PVE=80%)","m2","m3","m4"
    )) + facet_grid(cols=vars(method)) +
    geom_text(data = ann.text3,aes(label = label),size=1.8,color=rep(c("black"),3*4))
  
  
  plot_grid(
    p1,NULL,p2,NULL,p3,NULL, p4,NULL,p5,
    rel_heights = c(1, -0.25, 1, -0.05, 1, -0.1, 1, -0.05, 1),
    align='hv',
    labels=c("A","","B","","C","","D","","E"),# vjust=12,
    #label_colour="red",
    nrow=9, ncol=1)
  
  if(save){
    ggsave(path_grid, width=25, height=25, units="cm")  
  }
  setwd(original_dir)
}
