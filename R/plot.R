### Plots for paper

## Simulation Plots
RetrofitPlot <- function(dir, file, save=FALSE){
  
  path_A = paste("../", file, "_A", ".pdf", sep="")
  path_BCDE = paste("../", file, "_BCDE", ".pdf", sep="")
  original_dir=getwd()
  setwd(dir)
  library(ggplot2)
  library(cowplot)
  library(matrixStats)
  library(DescTools)
  
  color_func <- function(label){
    if(label > 0.05){
      return ("gray30")
    } else {
      return ("black")
    }
  }
  
  MethodTitle1="10 cells per spot \n up to 3 cell types per spot \n exact reference (10 true)"
  MethodTitle2="20 cells per spot \n up to 5 cell types per spot \n exact reference (10 true)"
  MethodTitle3="10 cells per spot \n up to 3 cell types per spot \n imperfect reference (5 true, 5 extra)"
  MethodTitle4="10 cells per spot \n up to 3 cell types per spot \n imperfect reference (3 true, 5 extra)"
  # MethodTitle1="N20M5"
  # MethodTitle2="N10M3"
  # MethodTitle3="Extra5"
  # MethodTitle4="Missing2"
  ### Read synthetic data
  ### plot for cor(W)
  W=as.matrix(read.csv("Cerebellum_W_K=10.csv",   row.names = 1, check.names = FALSE))
  # Mapped W_hat from RETROFIT (no. of genes x no. of cell types)
  W_mod1=read.csv("N=20,M=5_loc_W_retrofit.csv",  row.names = 1, check.names = FALSE)
  W_mod2=read.csv("N=10,M=3_loc_W_retrofit.csv",  row.names = 1, check.names = FALSE)
  W_mod3=read.csv("extra5_loc_W_retrofit.csv",    row.names = 1, check.names = FALSE)
  W_mod4=W_mod3
  W_mod5=read.csv("../mapped/N=20,M=5_loc_W_mapped_cor.csv",    row.names = 1, check.names = FALSE)
  W_mod6=read.csv("../mapped/N=20,M=5_loc_W_mapped_marker.csv", row.names = 1, check.names = FALSE)
  W_mod7=read.csv("../mapped/N=10,M=3_loc_W_mapped_cor.csv",    row.names = 1, check.names = FALSE)
  W_mod8=read.csv("../mapped/N=10,M=3_loc_W_mapped_marker.csv", row.names = 1, check.names = FALSE)
  W_mod9=read.csv("../mapped/extra5_loc_W_mapped_cor.csv",      row.names = 1, check.names = FALSE)
  W_mod10=read.csv("../mapped/extra5_loc_W_mapped_marker.csv",  row.names = 1, check.names = FALSE)
  W_mod11=W_mod9
  W_mod12=W_mod10
  
  cell_types = colnames(W)[1:5]
  W_mod3=W_mod3[,cell_types]
  W_mod4=W_mod4[,cell_types]
  W_mod9=W_mod9[,cell_types]
  W_mod10=W_mod10[,cell_types]
  W_mod11=W_mod11[,cell_types]
  W_mod12=W_mod12[,cell_types]
  cor_W1=sort(diag(cor(W,W_mod1)),decreasing=T,na.last=T)
  cor_W2=sort(diag(cor(W,W_mod2)),decreasing=T,na.last=T)
  cor_W3=sort(diag(cor(W[,cell_types],W_mod3)),decreasing=T,na.last=T)
  cor_W4=sort(diag(cor(W[,cell_types],W_mod4)),decreasing=T,na.last=T)
  cor_W5=sort(diag(cor(W,W_mod5)),decreasing=T,na.last=T)
  cor_W6=sort(diag(cor(W,W_mod6)),decreasing=T,na.last=T)
  cor_W7=sort(diag(cor(W,W_mod7)),decreasing=T,na.last=T)
  cor_W8=sort(diag(cor(W,W_mod8)),decreasing=T,na.last=T)
  cor_W9=sort(diag(cor(W[,cell_types],W_mod3)),decreasing=T,na.last=T)
  cor_W10=sort(diag(cor(W[,cell_types],W_mod4)),decreasing=T,na.last=T)
  cor_W11=sort(diag(cor(W[,cell_types],W_mod9)),decreasing=T,na.last=T)
  cor_W12=sort(diag(cor(W[,cell_types],W_mod11)),decreasing=T,na.last=T)
  
  
  df_corrW=as.data.frame(c(cor_W7,cor_W5,cor_W9,cor_W11))
  df_corrW$x=as.factor(c(1:length(cor_W7),1:length(cor_W5),1:length(cor_W9),1:length(cor_W11)))
  df_corrW$method=factor(c(rep(MethodTitle1,length(cor_W7)),
                           rep(MethodTitle2,length(cor_W5)),
                           rep(MethodTitle3,length(cor_W9)),
                           rep(MethodTitle4,length(cor_W11))),
                         c(MethodTitle1,
                           MethodTitle2,
                           MethodTitle3,
                           MethodTitle4))
  colnames(df_corrW)=c("cor","x","method")
  ############################################################################################################
  ################################# N=10,M=3  ################################################################
  ############################################################################################################
  H=read.csv("N=10,M=3_H.csv", row.names = 1, check.names = FALSE)
  H_mod1=read.csv("N=10,M=3_loc_H_retrofit.csv", row.names = 1, check.names = FALSE)
  H_mod2=read.csv("N=10,M=3_loc_RCTD_H.csv", row.names = 1, check.names = FALSE)
  H_mod3=read.csv("N=10,M=3_loc_RCTD-D_H.csv", row.names = 1, check.names = FALSE)
  H_mod4=read.csv("N=10,M=3_loc_Stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("N=10,M=3_loc_NMFreg.csv", check.names = FALSE)[,-11])
  H_mod6=read.csv("N=10,M=3_loc_SPOTlight.csv", row.names = 1, check.names = FALSE)
  H_mod7=read.csv("N=10,M=3_loc_SPOTlight.csv", row.names = 1, check.names = FALSE)
  H_txt7=read.table("N=10,M=3_loc_STdeconvolve.txt", header = TRUE)
  H_mod11=read.csv("../mapped/N=10,M=3_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod12=read.csv("../mapped/N=10,M=3_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  X=read.csv("N=10,M=3_loc_X.csv", row.names = 1, check.names = FALSE)
  
  K=nrow(H)
  S=ncol(H)
  
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H,H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H,H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H,H_mod7)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H,H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H,H_mod12)),decreasing=T,na.last=T)
  # hijacking STdeconvolve
  cor_H7=H_txt7[,"cor_H"]
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  # H_True=matrix(NA, ncol=S, nrow=nrow(H))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model11=NULL
  rmse_Model12=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model11=c(rmse_Model11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmse_Model12=c(rmse_Model12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
  }
  # hijacking STdeconvolve
  rmse_Model7=H_txt7[,"rmse_H"]
  
  # prop variance explained
  W_hat=read.csv("N=10,M=3_loc_W_hat_L=20.csv", row.names = 1, check.names = FALSE)
  H_hat=read.csv("N=10,M=3_loc_H_hat_L=20.csv", row.names = 1, check.names = FALSE)
  Thet_hat=read.csv("N=10,M=3_loc_Theta_hat_L=20.csv")[,-1]
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
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  X_hat80=(as.matrix(W_hat[,ind1[1:13]]) %*% diag(Thet_hat[ind1[1:13]]) + 0.01) %*% as.matrix(H_hat[ind1[1:13],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  # Hijacking STdeconvolve
  cor_X7=H_txt7[,"cor_X"]
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse11=NULL
  rmse12=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse11=c(rmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    rmse12=c(rmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  # Hijacking STdeconvolve
  rmse7=H_txt7[,"rmse_X"]
  
  # temp=as.data.frame(rbind(# cbind(rmse_Model1,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model11,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model3,3),cbind(rmse_Model4,4),
  #                          cbind(rmse_Model5,5),cbind(rmse_Model6,6)
  # ))
  
  # rmseH2
  temp=as.data.frame(rbind(cbind(rmse_Model11,1),cbind(rmse_Model2,2),
                           cbind(rmse_Model3,3),cbind(rmse_Model4,4),
                           cbind(rmse_Model5,5),cbind(rmse_Model6,6),
                           cbind(rmse_Model7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle1
  temp$title="10 cells/spot \n 3 cell types/spot \n Extact Ref (10 true)"
  df_rmseH=as.data.frame(rbind(temp))
  
  # corrH
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H11),
    cbind(seq(0,1,length.out = 1000),cor_H2),
    cbind(seq(0,1,length.out = 1000),cor_H3),
    cbind(seq(0,1,length.out = 1000),cor_H4),
    cbind(seq(0,1,length.out = 1000),cor_H5),
    cbind(seq(0,1,length.out = 1000),cor_H6),
    cbind(seq(0,1,length.out = 1000),cor_H7)))
  temp$model=c(rep("RETROFIT",1000),rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  colnames(temp)=c("x","cor","model")
  temp$method=MethodTitle1
  df_corrH=as.data.frame(rbind(temp))
  
  # nrmseX2
  temp=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse11,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),cbind(rmse4,4),
    cbind(rmse5,5),cbind(rmse6,6),
    cbind(rmse7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle1
  df_nrmseX=as.data.frame(rbind(temp))
  
  # corrX2
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X11),
    cbind(seq(0,1,length.out = 1000),cor_X80),
    cbind(seq(0,1,length.out = 1000),cor_X2),
    cbind(seq(0,1,length.out = 1000),cor_X3),
    cbind(seq(0,1,length.out = 1000),cor_X4),
    cbind(seq(0,1,length.out = 1000),cor_X5),
    cbind(seq(0,1,length.out = 1000),cor_X6),
    cbind(seq(0,1,length.out = 1000),cor_X7)))
  temp$model=c(rep("RETROFIT",1000),
               rep("RETROFIT (PVE)",1000),rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  temp$method=MethodTitle1
  colnames(temp)=c("x","cor","model","method")
  df_corrX=as.data.frame(rbind(temp))
  
  # rmseH2
  digit_length=2
  label = c(
    ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value,
    # ks.test(rmse_Model1,rmse_Model11, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model7, alternative="greater")$p.value
  )
  
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle1,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func)
  )
  # label=c("1.0e-4","0.003","0.975","2.0e-242","6.6e-83"))
  ann.text_rmseH=as.data.frame(rbind(temp))
  
  # corrH
  digit_length=3
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=digit_length), nsmall = 3)
  )
  temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
                            paste("RCTD:", label[[2]]),
                            paste("RCTD-D:", label[[3]]),
                            paste("Stereoscope:", label[[4]]),
                            paste("NMFreg:", label[[5]]),
                            paste("SPOTlight:", label[[6]]),
                            paste("STdeconvolve:", label[[7]])),
                  # label=c("RETROFIT: 0.954","RCTD: 0.926",
                  #         "RCTD-D: 0.937","Stereoscope: 0.958",
                  #         "NMFreg: 0.764",
                  #         "SPOTlight: 0.769"),
                  model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight", "STdeconvolve"),
                  method=factor(MethodTitle1,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrH=as.data.frame(rbind(temp))
  
  # nrmseX2
  digit_length=2
  label=c(
    ks.test(rmse1,rmse2, alternative="greater")$p.value,
    # ks.test(rmse1,rmse80, alternative="greater")$p.value,
    ks.test(rmse1,rmse3, alternative="greater")$p.value,
    ks.test(rmse1,rmse4, alternative="greater")$p.value,
    ks.test(rmse1,rmse5, alternative="greater")$p.value,
    ks.test(rmse1,rmse6, alternative="greater")$p.value,
    ks.test(rmse1,rmse7, alternative="greater")$p.value
  )
  
  
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle1,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func))
  # label=c("0.171","4.0e-4","0.294","0.844","0.254"))
  ann.text_nrmseX=as.data.frame(rbind(temp))
  
  # corrX2
  digit_length=3
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=digit_length), nsmall = 3)
  )
  temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
                            paste("RETROFIT (PVE):", label[[2]]),
                            paste("RCTD:", label[[3]]),
                            paste("RCTD-D:", label[[4]]),
                            paste("Stereoscope:", label[[5]]),
                            paste("NMFreg:", label[[6]]),
                            paste("SPOTlight:", label[[7]]),
                            paste("STdeconvolve:", label[[7]])),
                  # label=c("RETROFIT: 0.968",
                  #         "RETROFIT (PVE): 0.933","RCTD: 0.930",
                  #         "RCTD-D: 0.907","Stereoscope: 0.946",
                  #         "NMFreg: 0.802","SPOTlight: 0.832"),
                  model=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight", "STdeconvolve"),
                  method=factor(MethodTitle1,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrX=as.data.frame(rbind(temp))
  
  
  ############################################################################################################
  ################################# N=20,M=5 ################################################################
  ############################################################################################################
  ## Plot for RMSE(H)
  H=read.csv("N=20,M=5_H.csv",                    row.names = 1, check.names = FALSE)
  H_mod1=read.csv("N=20,M=5_loc_H_retrofit.csv",  row.names = 1, check.names = FALSE)
  H_mod2=read.csv("N=20,M=5_loc_RCTD_H.csv",      row.names = 1, check.names = FALSE)
  H_mod3=read.csv("N=20,M=5_loc_RCTD-D_H.csv",    row.names = 1, check.names = FALSE)
  H_mod4=read.csv("N=20,M=5_loc_Stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("N=20,M=5_loc_NMFreg.csv",check.names = FALSE)[,-11])
  H_mod6=read.csv("N=20,M=5_loc_SPOTlight.csv",   row.names = 1, check.names = FALSE)
  H_mod7=read.csv("N=20,M=5_loc_SPOTlight.csv", row.names = 1, check.names = FALSE)
  H_txt7=read.table("N=20,M=5_loc_STdeconvolve.txt", header = TRUE)
  X=read.csv("N=20,M=5_loc_X.csv",                row.names = 1, check.names = FALSE)
  W=as.matrix(read.csv("Cerebellum_W_K=10.csv",   row.names = 1, check.names = FALSE))
  H_mod11=read.csv("../mapped/N=20,M=5_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod12=read.csv("../mapped/N=20,M=5_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  
  K=nrow(H)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H,H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H,H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H,H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H,H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H,H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H,H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H,H_mod7)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H,H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H,H_mod12)),decreasing=T,na.last=T)
  # hijacking STdeconvolve
  cor_H7=H_txt7[,"cor_H"]
  
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  H_True=matrix(NA, ncol=S, nrow=nrow(H))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model11=NULL
  rmse_Model12=NULL
  
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model11=c(rmse_Model11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmse_Model12=c(rmse_Model12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
  }
  # hijacking STdeconvolve
  rmse_Model7=H_txt7[,"rmse_H"]
  ## prop variance explained
  
  # W_hat, H_hat and Theta_hat from retrofit (with all components before mapping)
  W_hat=read.csv("N=20,M=5_loc_W_hat_L=20.csv", row.names = 1, check.names = FALSE)
  H_hat=read.csv("N=20,M=5_loc_H_hat_L=20.csv", row.names = 1, check.names = FALSE)
  Thet_hat=read.csv("N=20,M=5_loc_Theta_hat_L=20.csv")[,-1]
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
  X_hat80=(as.matrix(W_hat[,ind1[1:12]]) %*% diag(Thet_hat[ind1[1:12]]) +
             0.01) %*% as.matrix(H_hat[ind1[1:12],]) #no. of components that explain about 80%
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  # Hijacking STdeconvolve
  cor_X7=H_txt7[,"cor_X"]
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse11=NULL
  rmse12=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse11=c(rmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    rmse12=c(rmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  # Hijacking STdeconvolve
  rmse7=H_txt7[,"rmse_X"]
  
  # rmseH
  temp=as.data.frame(rbind(cbind(rmse_Model11,1),cbind(rmse_Model2,2),
                           cbind(rmse_Model3,3),cbind(rmse_Model4,4),
                           cbind(rmse_Model5,5),cbind(rmse_Model6,6),
                           cbind(rmse_Model7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle2
  temp$title = "20 cells/spot \n 5 cell types/spot \n Extact Ref (10 true)"
  df_rmseH = as.data.frame(rbind(df_rmseH, temp))
  
  # corrH1
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H11),
    cbind(seq(0,1,length.out = 1000),cor_H2),
    cbind(seq(0,1,length.out = 1000),cor_H3),
    cbind(seq(0,1,length.out = 1000),cor_H4),
    cbind(seq(0,1,length.out = 1000),cor_H5),
    cbind(seq(0,1,length.out = 1000),cor_H6),
    cbind(seq(0,1,length.out = 1000),cor_H7)))
  temp$model=c(rep("RETROFIT",1000),rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  temp$method=MethodTitle2
  colnames(temp)=c("x","cor","model","method")
  df_corrH = as.data.frame(rbind(df_corrH, temp))
  
  # nrmseX1
  temp=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse11,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),cbind(rmse4,4),
    cbind(rmse5,5),cbind(rmse6,6),
    cbind(rmse7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle2
  df_nrmseX = as.data.frame(rbind(df_nrmseX, temp))
  
  # corrX1
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X11),
    cbind(seq(0,1,length.out = 1000),cor_X80),
    cbind(seq(0,1,length.out = 1000),cor_X2),
    cbind(seq(0,1,length.out = 1000),cor_X3),
    cbind(seq(0,1,length.out = 1000),cor_X4),
    cbind(seq(0,1,length.out = 1000),cor_X5),
    cbind(seq(0,1,length.out = 1000),cor_X6),
    cbind(seq(0,1,length.out = 1000),cor_X7)))
  temp$model=c(rep("RETROFIT",1000),
               rep("RETROFIT (PVE)",1000),rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  temp$method=MethodTitle2
  colnames(temp)=c("x","cor","model","method")
  df_corrX = as.data.frame(rbind(df_corrX, temp))
  
  # rmseH
  digit_length=2
  label = c(
    ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value,
    # ks.test(rmse_Model1,rmse_Model11, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model7, alternative="greater")$p.value
  )
  
  ## add p-values from above KS tests
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle2,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func))
  ann.text_rmseH=as.data.frame(rbind(ann.text_rmseH,temp))
  # label=c("1.7e-207","8.4e-113","3.1e-10","1.8e-281","2.4e-121"))
  
  # corrH1
  digit_length=3
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=digit_length), nsmall = 3)
  )
  
  H_mod1=read.csv("N=20,M=5_loc_H_retrofit.csv", row.names = 1, check.names = FALSE)
  H_mod2=read.csv("N=20,M=5_loc_RCTD_H.csv", row.names = 1, check.names = FALSE)
  H_mod3=read.csv("N=20,M=5_loc_RCTD-D_H.csv", row.names = 1, check.names = FALSE)
  H_mod4=read.csv("N=20,M=5_loc_Stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("N=20,M=5_loc_NMFreg.csv", check.names = FALSE)[,-11])
  H_mod6=read.csv("N=20,M=5_loc_SPOTlight.csv", row.names = 1, check.names = FALSE)
  H_mod7=read.csv("N=20,M=5_loc_SPOTlight.csv", row.names = 1, check.names = FALSE)
  H_txt7=read.table("N=20,M=5_loc_STdeconvolve.txt", header = TRUE)
  
  ## add AUC from above results
  temp=data.frame(label=c(paste("RETROFIT:", label[[1]]),
                          paste("RCTD:", label[[2]]),
                          paste("RCTD-D:", label[[3]]),
                          paste("Stereoscope:", label[[4]]),
                          paste("NMFreg:", label[[5]]),
                          paste("SPOTlight:", label[[6]]),
                          paste("STdeconvolve:", label[[7]])),
                  # label=c("RETROFIT: 0.952","RCTD: 0.651",
                  #         "RCTD-D: 0.810","Stereoscope: 0.939","NMFreg: 0.628",
                  #         "SPOTlight: 0.762"),
                  model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight", "STdeconvolve"),
                  method=factor(MethodTitle2,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrH=as.data.frame(rbind(ann.text_corrH,temp))
  
  # nrmseX1
  digit_length=2
  label = c(
    ks.test(rmse1,rmse2, alternative="greater")$p.value,
    # ks.test(rmse1,rmse11, alternative="greater")$p.value,
    ks.test(rmse1,rmse3, alternative="greater")$p.value,
    ks.test(rmse1,rmse4, alternative="greater")$p.value,
    ks.test(rmse1,rmse5, alternative="greater")$p.value,
    ks.test(rmse1,rmse6, alternative="greater")$p.value,
    ks.test(rmse1,rmse7, alternative="greater")$p.value
    # ks.test(rmse1,rmse80, alternative="greater")$p.value
  )
  
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle2,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func))
  ann.text_nrmseX=as.data.frame(rbind(ann.text_nrmseX,temp))
  # label=c("5.1e-48","5.5e-10","0.186","5.9e-86","1.6e-28"))#,"2.1e-28"))
  
  # corrX1
  digit_length=3
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=digit_length), nsmall = 3)
  )
  temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
                            paste("RETROFIT (PVE):", label[[2]]),
                            paste("RCTD:", label[[3]]),
                            paste("RCTD-D:", label[[4]]),
                            paste("Stereoscope:", label[[5]]),
                            paste("NMFreg:", label[[6]]),
                            paste("SPOTlight:", label[[7]]),
                            paste("STdeconvolve:", label[[8]])),
                  # label=c("RETROFIT: 0.966","RETROFIT (PVE): 0.912",
                  #         "RCTD: 0.821",
                  #         "RCTD-D: 0.850","Stereoscope: 0.973",
                  #         "NMFreg: 741","SPOTlight: 0.873"),
                  model=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight", "STdeconvolve"),
                  method=factor(MethodTitle2,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrX=as.data.frame(rbind(ann.text_corrX,temp))
  
  ############################################################################################################
  ################################# extra5  ################################################################
  ############################################################################################################
  H=read.csv("extra5_H.csv",                    row.names = 1, check.names = FALSE)
  H_mod1=read.csv("extra5_loc_H_retrofit.csv",  row.names = 1, check.names = FALSE)
  H_mod2=read.csv("extra5_loc_RCTD_H.csv",      row.names = 1, check.names = FALSE)
  H_mod3=read.csv("extra5_loc_RCTD-D_H.csv",    row.names = 1, check.names = FALSE)
  H_mod4=read.csv("extra5_loc_Stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("extra5_loc_NMFreg.csv",    check.names = FALSE)[,-11])
  H_mod6=read.csv("extra5_loc_SPOTlight.csv",   row.names = 1, check.names = FALSE)
  H_mod7=read.csv("extra5_loc_SPOTlight.csv",   row.names = 1, check.names = FALSE)
  H_txt7=read.table("extra5_loc_STdeconvolve.txt", header = TRUE)
  H_mod11=read.csv("../mapped/extra5_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod12=read.csv("../mapped/extra5_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  X=read.csv("extra5_loc_X.csv", row.names = 1, check.names = FALSE)
  
  cell_types = rownames(H)[1:5]
  H_mod1=H_mod1[cell_types,]
  H_mod2=H_mod2[cell_types,]
  H_mod3=H_mod3[cell_types,]
  H_mod4=H_mod4[cell_types,]
  H_mod5=H_mod5[cell_types,]
  H_mod6=H_mod6[cell_types,]
  H_mod7=H_mod7[cell_types,]
  H_mod11=H_mod11[cell_types,]
  H_mod12=H_mod12[cell_types,]
  
  K=nrow(H)
  G=nrow(W)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H[cell_types,],H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H[cell_types,],H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H[cell_types,],H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H[cell_types,],H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H[cell_types,],H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H[cell_types,],H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H[cell_types,],H_mod7)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H[cell_types,],H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H[cell_types,],H_mod12)),decreasing=T,na.last=T)
  # hijacking STdeconvolve
  cor_H7=H_txt7[,"cor_H"]
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  # H_True=matrix(NA, ncol=S, nrow=nrow(H))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model11=NULL
  rmse_Model12=NULL
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model11=c(rmse_Model11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmse_Model12=c(rmse_Model12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
  }
  # hijacking STdeconvolve
  rmse_Model7=H_txt7[,"rmse_H"]
  
  # prop variance explained
  W_hat=read.csv("extra5_loc_W_hat_L=10.csv", row.names = 1, check.names = FALSE)
  H_hat=read.csv("extra5_loc_H_hat_L=10.csv", row.names = 1, check.names = FALSE)
  Thet_hat=read.csv("extra5_loc_Theta_hat_L=10.csv")[,-1]
  X_hat=(as.matrix(W_hat) %*% diag(Thet_hat) + 0.01) %*% as.matrix(H_hat)
  temp=rep(NA,length(Thet_hat))
  for(k in 1:length(temp)){
    temp[k]=sum(as.matrix(W_hat[,k] * Thet_hat[k] + 0.01) %*% as.matrix(H_hat[k,]))
  }
  prop3=sort(temp/sum(temp),decreasing=T)
  ind1=sort(temp/sum(temp),decreasing=T,index.return=T)$ix
  
  W=W[,cell_types]
  X_hat1=W %*% H_Model1
  X_hat2=W %*% H_Model2
  X_hat3=W %*% H_Model3
  X_hat4=W %*% H_Model4
  X_hat5=W %*% H_Model5
  X_hat6=W %*% H_Model6
  X_hat7=W %*% H_Model7
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  X_hat80=(as.matrix(W_hat[,ind1[1:6]]) %*% diag(Thet_hat[ind1[1:6]]) + 0.01) %*% as.matrix(H_hat[ind1[1:6],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  # Hijacking STdeconvolve
  cor_X7=H_txt7[,"cor_X"]
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse11=NULL
  rmse12=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse11=c(rmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    rmse12=c(rmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  # Hijacking STdeconvolve
  rmse7=H_txt7[,"rmse_X"]
  # temp=as.data.frame(rbind(# cbind(rmse_Model1,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model11,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model3,3),cbind(rmse_Model4,4),
  #                          cbind(rmse_Model5,5),cbind(rmse_Model6,6)
  # ))
  
  # rmseH3
  temp=as.data.frame(rbind(cbind(rmse_Model11,1),cbind(rmse_Model2,2),
                           cbind(rmse_Model3,3),cbind(rmse_Model4,4),
                           cbind(rmse_Model5,5),cbind(rmse_Model6,6),
                           cbind(rmse_Model7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle3
  temp$title = "10 cells/spot \n 3 cell types/spot \n Imperfect Ref (5 true, 5 extra)"
  df_rmseH=as.data.frame(rbind(df_rmseH,temp))
  
  # corrH3
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out=1000),cor_H1),
    cbind(seq(0,1,length.out=1000),cor_H11),
    cbind(seq(0,1,length.out=1000),cor_H2),
    cbind(seq(0,1,length.out=1000),cor_H3),
    cbind(seq(0,1,length.out=1000),cor_H4),
    cbind(seq(0,1,length.out=1000),cor_H5),
    cbind(seq(0,1,length.out=1000),cor_H6),
    cbind(seq(0,1,length.out=1000),cor_H7)))
  temp$model=c(rep("RETROFIT",1000),
               rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  colnames(temp)=c("x","cor","model")
  temp$method=MethodTitle3
  df_corrH=as.data.frame(rbind(df_corrH,temp))
  
  # nrmseX3
  temp=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse11,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),cbind(rmse4,4),
    cbind(rmse5,5),cbind(rmse6,6),
    cbind(rmse7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle3
  df_nrmseX=as.data.frame(rbind(df_nrmseX,temp))
  
  # corrX3
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out=1000),cor_X1),
    cbind(seq(0,1,length.out=1000),cor_X11),
    cbind(seq(0,1,length.out=1000),cor_X80),
    cbind(seq(0,1,length.out=1000),cor_X2),
    cbind(seq(0,1,length.out=1000),cor_X3),
    cbind(seq(0,1,length.out=1000),cor_X4),
    cbind(seq(0,1,length.out=1000),cor_X5),
    cbind(seq(0,1,length.out=1000),cor_X6),
    cbind(seq(0,1,length.out=1000),cor_X7)))
  temp$model=c(rep("RETROFIT",1000),rep("RETROFIT (PVE)",1000),
               rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  temp$method=MethodTitle3
  colnames(temp)=c("x","cor","model","method")
  df_corrX=as.data.frame(rbind(df_corrX,temp))
  
  # rmseH3
  digit_length=2
  label=c(
    ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value,
    # ks.test(rmse_Model1,rmse_Model11, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model7, alternative="greater")$p.value
  )
  
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle3,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func))
  # label=c("2.9e-35","3.5e-29","0.002","2.3e-150","1e-81"))
  ann.text_rmseH=as.data.frame(rbind(ann.text_rmseH,temp))
  
  # corrH3
  digit_length=3
  cor_H6[is.na(cor_H6)]=0
  cor_H7[is.na(cor_H7)]=0
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=digit_length), nsmall = 3)
  )
  temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
                            paste("RCTD:", label[[2]]),
                            paste("RCTD-D:", label[[3]]),
                            paste("Stereoscope:", label[[4]]),
                            paste("NMFreg:", label[[5]]),
                            paste("SPOTlight:", label[[6]]),
                            paste("STdeconvolve:", label[[7]])),
                  # label=c("RETROFIT: 0.977","RCTD: 0.891",
                  #         "RCTD-D: 0.915","Stereoscope: 0.967","NMFreg: 0.904",
                  #         "SPOTlight: 0.750"),
                  model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight","STdeconvolve"),
                  method=factor(MethodTitle3,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrH=as.data.frame(rbind(ann.text_corrH,temp))
  
  # nrmseX3
  digit_length=2
  label=c(
    ks.test(rmse1,rmse2, alternative="greater")$p.value,
    # ks.test(rmse1,rmse80, alternative="greater")$p.value,
    ks.test(rmse1,rmse3, alternative="greater")$p.value,
    ks.test(rmse1,rmse4, alternative="greater")$p.value,
    ks.test(rmse1,rmse5, alternative="greater")$p.value,
    ks.test(rmse1,rmse6, alternative="greater")$p.value,
    ks.test(rmse1,rmse7, alternative="greater")$p.value
  )
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle3,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func))
  # label=c("9.9e-10","1.6e-08","0.294","5.3e-37","1.4e-31"))
  ann.text_nrmseX=as.data.frame(rbind(ann.text_nrmseX,temp))
  
  # corrX3
  digit_length=3
  cor_X6[is.na(cor_X6)]=0
  cor_X7[is.na(cor_X7)]=0
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=digit_length), nsmall = 3)
  )
  temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
                            paste("RETROFIT (PVE):", label[[2]]),
                            paste("RCTD:", label[[3]]),
                            paste("RCTD-D:", label[[4]]),
                            paste("Stereoscope:", label[[5]]),
                            paste("NMFreg:", label[[6]]),
                            paste("SPOTlight:", label[[7]]),
                            paste("STdeconvolve:", label[[8]])),
                  # label=c("RETROFIT: 0.986","RETROFIT (PVE): 0.963",
                  #         "RCTD: 0.941",
                  #         "RCTD-D: 0.920","Stereoscope: 0.983",
                  #         "NMFreg: 0.807","SPOTlight: 0.805"),
                  model=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight", "STdeconvolve"),
                  method=factor(MethodTitle3,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrX=as.data.frame(rbind(ann.text_corrX,temp))
  
  
  ############################################################################################################
  ################################# missing5  ################################################################
  ############################################################################################################
  H=read.csv("extra5_H.csv",                  row.names = 1, check.names = FALSE)
  H_mod1=read.csv("extra5_loc_H_retrofit.csv",row.names = 1, check.names = FALSE)
  H_mod2=read.csv("miss2_loc_RCTD_H.csv",     row.names = 1, check.names = FALSE)
  H_mod3=read.csv("miss2_loc_RCTD-D_H.csv",   row.names = 1, check.names = FALSE)
  H_mod4=read.csv("miss2_loc_Stereoscope.csv",row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("miss2_loc_NMFreg.csv",   check.names = FALSE)[,-11])
  H_mod6=read.csv("miss2_loc_SPOTlight.csv",  row.names = 1, check.names = FALSE)
  H_mod7=read.csv("miss2_loc_SPOTlight.csv",  row.names = 1, check.names = FALSE)
  H_txt7=read.table("extra5_loc_STdeconvolve.txt", header = TRUE)
  H_mod11=read.csv("../mapped/extra5_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod12=read.csv("../mapped/extra5_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  X=read.csv("extra5_loc_X.csv", row.names = 1, check.names = FALSE)
  
  cell_types = rownames(H)[1:5]
  cell_types_missing = rownames(H)[c(3,5)]
  H_mod2[c(4,5),]=0
  H_mod3[c(4,5),]=0
  H_mod4[c(4,5),]=0
  H_mod5[c(4,5),]=0
  H_mod6[c(4,5),]=0
  H_mod7[c(4,5),]=0
  rownames(H_mod2)[c(4,5)] = cell_types_missing
  rownames(H_mod3)[c(4,5)] = cell_types_missing
  rownames(H_mod4)[c(4,5)] = cell_types_missing
  rownames(H_mod5)[c(4,5)] = cell_types_missing
  rownames(H_mod6)[c(4,5)] = cell_types_missing
  rownames(H_mod7)[c(4,5)] = cell_types_missing
  # rownames(H_mod7)[c(4,5)] = cell_types_missing # stdconvolve => same structure as extra5_loc_H_retrofit.csv
  H_mod1=H_mod1[cell_types,]
  H_mod2=H_mod2[cell_types,]
  H_mod3=H_mod3[cell_types,]
  H_mod4=H_mod4[cell_types,]
  H_mod5=H_mod5[cell_types,]
  H_mod6=H_mod6[cell_types,]
  H_mod7=H_mod7[cell_types,]
  H_mod11=H_mod11[cell_types,]
  H_mod12=H_mod12[cell_types,]
  
  K=nrow(H)
  G=nrow(W)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H[cell_types,],H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H[cell_types,],H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H[cell_types,],H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H[cell_types,],H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H[cell_types,],H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H[cell_types,],H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H[cell_types,],H_mod7)),decreasing=T,na.last=T)
  cor_H11=sort(diag(cor(H[cell_types,],H_mod11)),decreasing=T,na.last=T)
  cor_H12=sort(diag(cor(H[cell_types,],H_mod12)),decreasing=T,na.last=T)
  # hijacking STdeconvolve
  cor_H7=H_txt7[,"cor_H"]
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model11=matrix(NA, ncol=S, nrow=nrow(H_mod11))
  H_Model12=matrix(NA, ncol=S, nrow=nrow(H_mod12))
  # H_True=matrix(NA, ncol=S, nrow=nrow(H))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model11=NULL
  rmse_Model12=NULL
  for(i in 1:S){
    H_True[,i]=H[,i]/sum(H[,i])
    H_Model1[,i]=H_mod1[,i]/sum(H_mod1[,i])
    H_Model2[,i]=H_mod2[,i]/sum(H_mod2[,i])
    H_Model3[,i]=H_mod3[,i]/sum(H_mod3[,i])
    H_Model4[,i]=H_mod4[,i]/sum(H_mod4[,i])
    H_Model5[,i]=H_mod5[,i]/sum(H_mod5[,i])
    H_Model6[,i]=H_mod6[,i]/sum(H_mod6[,i])
    H_Model7[,i]=H_mod7[,i]/sum(H_mod7[,i])
    H_Model11[,i]=H_mod11[,i]/sum(H_mod11[,i])
    H_Model12[,i]=H_mod12[,i]/sum(H_mod12[,i])
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model11=c(rmse_Model11,sqrt(mean((H_True[,i]-H_Model11[,i])^2)))
    rmse_Model12=c(rmse_Model12,sqrt(mean((H_True[,i]-H_Model12[,i])^2)))
  }
  # hijacking STdeconvolve
  rmse_Model7=H_txt7[,"rmse_H"]
  
  X_hat=(as.matrix(W_hat) %*% diag(Thet_hat) + 0.01) %*% as.matrix(H_hat)
  temp=rep(NA,length(Thet_hat))
  for(k in 1:length(temp)){
    temp[k]=sum(as.matrix(W_hat[,k] * Thet_hat[k] + 0.01) %*% as.matrix(H_hat[k,]))
  }
  prop=sort(temp/sum(temp),decreasing=T)
  ind1=sort(temp/sum(temp),decreasing=T,index.return=T)$ix
  
  
  X_hat1=W %*% H_Model1
  X_hat2=W %*% H_Model2
  X_hat3=W %*% H_Model3
  X_hat4=W %*% H_Model4
  X_hat5=W %*% H_Model5
  X_hat6=W %*% H_Model6
  X_hat7=W %*% H_Model7
  X_hat11=W %*% H_Model11
  X_hat12=W %*% H_Model12
  X_hat80=(as.matrix(W_hat[,ind1[1:6]]) %*% diag(Thet_hat[ind1[1:6]]) + 0.01) %*% as.matrix(H_hat[ind1[1:6],])
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X11=sort(diag(cor(X,X_hat11)),decreasing=T,na.last=T)
  cor_X12=sort(diag(cor(X,X_hat12)),decreasing=T,na.last=T)
  cor_X80=sort(diag(cor(X,X_hat80)),decreasing=T,na.last=T)
  # Hijacking STdeconvolve
  cor_X7=H_txt7[,"cor_X"]
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse11=NULL
  rmse12=NULL
  rmse80=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse11=c(rmse11,sqrt(mean(((X[,i]-X_hat11[,i])/sd(X[,i]))^2)))
    rmse12=c(rmse12,sqrt(mean(((X[,i]-X_hat12[,i])/sd(X[,i]))^2)))
    rmse80=c(rmse80,sqrt(mean(((X[,i]-X_hat80[,i])/sd(X[,i]))^2)))
  }
  # Hijacking STdeconvolve
  rmse7=H_txt7[,"rmse_X"]
  
  # temp=as.data.frame(rbind(# cbind(rmse_Model1,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model11,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model3,3),cbind(rmse_Model4,4),
  #                          cbind(rmse_Model5,5),cbind(rmse_Model6,6)
  # ))
  
  # rmseH4
  temp=as.data.frame(rbind(cbind(rmse_Model11,1),cbind(rmse_Model2,2),
                           cbind(rmse_Model3,3),cbind(rmse_Model4,4),
                           cbind(rmse_Model5,5),cbind(rmse_Model6,6),
                           cbind(rmse_Model7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle4
  temp$title="10 cells/spot \n 3 cell types/spot \n Imperfect Ref (3 true, 5 extra)"
  df_rmseH=as.data.frame(rbind(df_rmseH,temp))
  
  # corrH4
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_H1),
    cbind(seq(0,1,length.out = 1000),cor_H11),
    cbind(seq(0,1,length.out = 1000),cor_H2),
    cbind(seq(0,1,length.out = 1000),cor_H3),
    cbind(seq(0,1,length.out = 1000),cor_H4),
    cbind(seq(0,1,length.out = 1000),cor_H5),
    cbind(seq(0,1,length.out = 1000),cor_H6),
    cbind(seq(0,1,length.out = 1000),cor_H7)))
  temp$model=c(rep("RETROFIT",1000),rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  colnames(temp)=c("x","cor","model")
  temp$method=MethodTitle4
  df_corrH=as.data.frame(rbind(df_corrH,temp))
  
  # nrmseX4
  temp=as.data.frame(rbind(# cbind(rmse1,1),cbind(rmse80,7),
    cbind(rmse11,1),cbind(rmse80,8),
    cbind(rmse2,2),
    cbind(rmse3,3),cbind(rmse4,4),
    cbind(rmse5,5),cbind(rmse6,6),
    cbind(rmse7,7)
  ))
  colnames(temp)=c("rmse","model")
  temp$method=MethodTitle4
  df_nrmseX=as.data.frame(rbind(df_nrmseX,temp))
  
  # corrX4
  temp=as.data.frame(rbind(# cbind(seq(0,1,length.out = 1000),cor_X1),
    cbind(seq(0,1,length.out = 1000),cor_X11),
    cbind(seq(0,1,length.out = 1000),cor_X80),
    cbind(seq(0,1,length.out = 1000),cor_X2),
    cbind(seq(0,1,length.out = 1000),cor_X3),
    cbind(seq(0,1,length.out = 1000),cor_X4),
    cbind(seq(0,1,length.out = 1000),cor_X5),
    cbind(seq(0,1,length.out = 1000),cor_X6),
    cbind(seq(0,1,length.out = 1000),cor_X7)))
  temp$model=c(rep("RETROFIT",1000),rep("RETROFIT (PVE)",1000),rep("RCTD",1000),
               rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
               rep("NMFreg",1000),rep("SPOTlight",1000),
               rep("STdeconvolve",1000))
  temp$method=MethodTitle4
  colnames(temp)=c("x","cor","model","method")
  df_corrX=as.data.frame(rbind(df_corrX,temp))
  
  # rmseH4
  digit_length=2
  label=c(
    ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value,
    # ks.test(rmse_Model1,rmse_Model11, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value,
    ks.test(rmse_Model1,rmse_Model7, alternative="greater")$p.value
  )
  
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle4,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func))
  # label=c("4e-4","0.265","5.8e-8","3.9e-7","1.3e-14"))
  ann.text_rmseH=as.data.frame(rbind(ann.text_rmseH,temp))
  
  # corrH4
  digit_length=3
  cor_H1[is.na(cor_H1)]=0
  cor_H2[is.na(cor_H2)]=0
  cor_H3[is.na(cor_H3)]=0
  cor_H4[is.na(cor_H4)]=0
  cor_H5[is.na(cor_H5)]=0
  cor_H6[is.na(cor_H6)]=0
  cor_H7[is.na(cor_H7)]=0
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=digit_length), nsmall = 3)
  )
  temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
                            paste("RCTD:", label[[2]]),
                            paste("RCTD-D:", label[[3]]),
                            paste("Stereoscope:", label[[4]]),
                            paste("NMFreg:", label[[5]]),
                            paste("SPOTlight:", label[[6]]),
                            paste("STdeconvolve:", label[[7]])),
                  
                  # label=c("RETROFIT: 0.613",
                  #         "RCTD: 0.505",
                  #         "RCTD-D: 0.509","Stereoscope: 0.440",
                  #         "NMFreg: 0.576",
                  #         "SPOTlight: 0.411"),
                  model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight", "STdeconvolve"),
                  method=factor(MethodTitle4,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrH=as.data.frame(rbind(ann.text_corrH,temp))
  
  # nrmseX4
  digit_length=2
  label=c(
    ks.test(rmse1,rmse2, alternative="greater")$p.value,
    # ks.test(rmse1,rmse80, alternative="greater")$p.value,
    ks.test(rmse1,rmse3, alternative="greater")$p.value,
    ks.test(rmse1,rmse4, alternative="greater")$p.value,
    ks.test(rmse1,rmse5, alternative="greater")$p.value,
    ks.test(rmse1,rmse6, alternative="greater")$p.value,
    ks.test(rmse1,rmse7, alternative="greater")$p.value
  )
  
  temp=data.frame(model=2:7,rmse=0.9,
                  method=factor(MethodTitle4,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)),
                  label=format(label, digits=digit_length),
                  color=sapply(label, color_func))
  # label=c("0.171","0.047","5.4e-9","0.254","2e-4"))
  ann.text_nrmseX=as.data.frame(rbind(ann.text_nrmseX,temp))
  
  # corrX4
  digit_length=3
  cor_X1[is.na(cor_X1)]=0
  cor_X2[is.na(cor_X2)]=0
  cor_X3[is.na(cor_X3)]=0
  cor_X4[is.na(cor_X4)]=0
  cor_X5[is.na(cor_X5)]=0
  cor_X6[is.na(cor_X6)]=0
  cor_X7[is.na(cor_X7)]=0
  label = c(
    # format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X11),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X80),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=digit_length), nsmall = 3),
    format(round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=digit_length), nsmall = 3)
  )
  temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
                            paste("RETROFIT (PVE):", label[[2]]),
                            paste("RCTD:", label[[3]]),
                            paste("RCTD-D:", label[[4]]),
                            paste("Stereoscope:", label[[5]]),
                            paste("NMFreg:", label[[6]]),
                            paste("SPOTlight:", label[[7]]),
                            paste("STdeconvolve:", label[[8]])),
                  # label=c("RETROFIT: 0.552","RETROFIT (PVE): 0.963",
                  #         "RCTD: 0.541",
                  #         "RCTD-D: 0.541","Stereoscope: 0.537",
                  #         "NMFreg: 0.546","SPOTlight: 0.421"),
                  model=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope",
                          "NMFreg","SPOTlight","STdeconvolve"),
                  method=factor(MethodTitle4,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4)))
  ann.text_corrX=as.data.frame(rbind(ann.text_corrX,temp))
  
  # df_corrW
  # df_corrW$method=factor(df_corrW$method,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4))
  # rmseH-
  df_rmseH$model=as.factor(df_rmseH$model)
  # corrH-
  df_corrH$model=factor(df_corrH$model,levels=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight", "STdeconvolve"))
  # rmseH-
  df_rmseH$method=factor(df_rmseH$method,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4))
  # corrH-
  df_corrH$method=factor(df_corrH$method,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4))
  # nrmseX-
  df_nrmseX$model=factor(df_nrmseX$model,levels=c(1,8,2,3,4,5,6,7))
  # corrX-
  df_corrX$model=factor(df_corrX$model,levels=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight", "STdeconvolve"))
  # nrmseX-
  df_nrmseX$method=factor(df_nrmseX$method,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4))
  # corrX-
  df_corrX$method=factor(df_corrX$method,levels=c(MethodTitle1,MethodTitle2,MethodTitle3,MethodTitle4))
  # rmseH-
  ann.text_rmseH$model=as.factor(ann.text_rmseH$model)
  # corrH-
  ann.text_corrH$model=as.factor(ann.text_corrH$model)
  # nrmseX-
  ann.text_nrmseX$model=as.factor(ann.text_nrmseX$model)
  # corrX-
  ann.text_corrX$model=factor(ann.text_corrX$model,levels=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight", "STdeconvolve"))
  # corrH-
  ann.text_corrH$cor=rep(c(0.65,0.55,0.45,0.35,0.25,0.15,0.05),4)
  ann.text_corrH$x=0.05
  ann.text_corrH$color=rep(c('gray30','#56B4E9','#0072B2',
                             '#009E73',"#D55E00", "#CC79A7", "#4B0092"),4)
  
  # corrX-
  ann.text_corrX$cor=rep(c(0.75,0.65,0.55,0.45,0.35,0.25,0.15,0.05),4)
  ann.text_corrX$x=0.05
  
  # nrmseX- label location
  cnt = 6
  ann.text_rmseH$rmse=c(rep(0.75,cnt),rep(0.75,cnt),rep(0.75,cnt),rep(0.75,cnt))
  ann.text_nrmseX$rmse=c(rep(3.5,cnt),rep(3.5,cnt),rep(3.5,cnt),rep(3.5,cnt))
  
  
  # rmseH-
  plot_rmseH=ggplot(df_rmseH, aes(x=model, y=rmse, fill=model)) +
    geom_violin(lwd=0) + theme_bw()+ ylim(0,0.75)+
    xlab("") +
    ylab(expression(paste("RMSE (",H,",",widetilde(H),")"))) +
    theme(legend.position = "None",
          axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          axis.text.y = element_text(size=7),
          # strip.text = element_blank(),
          plot.margin = unit(c(0.2, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c('gray30','#56B4E9','#0072B2','#009E73',"#D55E00", "#CC79A7", "#4B0092"))+
    scale_x_discrete(labels=c("RETROFIT","RCTD","RCTD-D","Stereoscope","NMFreg","SPOTlight", "STdeconvolve"))+
    geom_text(data = ann.text_rmseH,aes(label = label),size=1.8, 
              color=ann.text_rmseH$color) +
    facet_grid(cols=vars(method))
  
  # corrH-
  plot_corrH=ggplot(df_corrH, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + 
    xlab("Normalized Rank") +
    ylab(expression(paste("Correlation (",H,",",widetilde(H),")"))) +
    theme_bw()+
    scale_color_manual(values = c('gray30','#56B4E9','#0072B2',
                                  '#009E73',"#D55E00", "#CC79A7", "#4B0092"),
                       labels=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight", "STdeconvolve")) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.title=element_text(size=9),axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          axis.text.y = element_text(size=7),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    geom_text(data=ann.text_corrH,aes(label = label),size=2.2,hjust=0,
              color=ann.text_corrH$color) 
  
  # corrX-
  plot_corrX=ggplot(df_corrX, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model,linetype=model)) + 
    xlab("Normalized Rank") +
    ylab(expression(paste("Correlation (",X,",",widetilde(X),")"))) +
    theme_bw()+
    scale_color_manual(values = c('gray30','grey30','#56B4E9','#0072B2',
                                  '#009E73',"#D55E00", "#CC79A7", "#4B0092"),
                       labels=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight", "STdeconvolve")) +
    scale_linetype_manual(values = c('solid','dashed','solid','solid', 'solid',"solid", "solid", "solid"),
                          labels=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight", "STdeconvolve")) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.text.x = element_text(size=7),
          axis.title.x = element_text(vjust=10),
          axis.text.y = element_text(size=7),
          #axis.title.x = element_text(vjust = 12),
          axis.title=element_text(size=9),
          plot.margin = unit(c(-0.5, 0.1, 0, 0.1), "cm"))+
    geom_text(data = ann.text_corrX,aes(label = label),size=2.2,hjust=0,
              color=rep(c('gray35','grey20','#56B4E9','#0072B2',
                          '#009E73',"#D55E00", "#CC79A7", "#4B0092"),4))
  
  ## plot for main paper
  quantile(df_nrmseX$rmse[df_nrmseX$method==MethodTitle1],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df_nrmseX$rmse[df_nrmseX$method==MethodTitle3],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df_nrmseX$rmse[df_nrmseX$method==MethodTitle4],c(0.75,0.8,0.9), na.rm = TRUE)
  
  # filter
  df_nrmseX$rmse[df_nrmseX$rmse>3]<-3
  
  # df5=df_nrmseX[complete.cases(df_nrmseX),]
  # quantile(df5$rmse[df5$method==MethodTitle1],c(0.75,0.8,0.9), na.rm = TRUE)
  # quantile(df5$rmse[df5$method==MethodTitle3],c(0.75,0.8,0.9), na.rm = TRUE)
  # quantile(df5$rmse[df5$method==MethodTitle4],c(0.75,0.8,0.9), na.rm = TRUE)
  
  # df5=df5[df5$rmse<3,]
  plot_nrmseX=ggplot(df_nrmseX, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("NRMSE (",X,",",widetilde(X),")"))) +
    xlab("") +
    geom_violin(lwd=0) + theme_bw()+ ylim(0,3.5)+
    # geom_boxplot(width=0.1)+
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          axis.text.y = element_text(size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(-0.5, 0.1, 0, 0.1), "cm"))+
    scale_fill_manual(values=c('gray35','grey20','#56B4E9','#0072B2', '#009E73',"#D55E00", "#CC79A7", "#4B0092"))+
    scale_x_discrete(labels=c("RETROFIT","RETROFIT (PVE)","RCTD","RCTD-D","Stereoscope","NMFreg","SPOTlight", "STdeconvolve"
    )) + facet_grid(cols=vars(method)) +
    geom_text(data=ann.text_nrmseX,
              aes(label = label),size=1.8,
              color=ann.text_nrmseX$color)
  # color=c("black","black","grey60","black","black","black",
  #         "grey60","black","grey60","grey60","grey60","black",
  #         "black","black","grey60","black","black","black",
  #         "grey60","black","black","grey60","black","black"))
  
  # plot for cor(W)
  plot_corrW=ggplot(df_corrW,aes(x=x,y=cor)) +
    geom_line(aes(x=x,y=cor,group=method)) + geom_point(aes(x=x,y=cor,group=method))+
    xlab("Cell Types") +
    ylab(expression(paste("Correlation (",W^0,",",widetilde(W),")"))) +
    theme_bw()+
    ylim(0,1.05)+
    theme(axis.title=element_text(size=9), axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm")) +
    facet_grid(cols=vars(method),scales = "free_x")
  
  
  plot_grid(
    NULL, plot_rmseH,NULL,NULL,NULL,NULL,NULL,NULL, NULL,NULL,
    rel_heights = c(0.5, 1.2, 0, 1, 0, 1, 0, 1, 0, 1),
    align='hv',
    # label_x = 0,
    label_y = 1,
    # hjust = -0.5,
    # vjust = 0.8,
    labels=c("","A","","B","","C","","D","","E", ""),
    nrow=10, ncol=1)
  
  if(save){
    ggsave(path_A, width=25, height=35, units="cm")  
  }
  plot_grid(
    NULL, NULL,NULL,plot_corrH,NULL, plot_nrmseX,NULL,plot_corrX, NULL, plot_corrW,
    rel_heights = c(0.5, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    align='hv',
    # label_x = 0,
    label_y = 1.0,
    # hjust = -0.5,
    # vjust = 0.8,
    labels=c("","A","","B","","C","","D","","E", ""),
    nrow=10, ncol=1)
  if(save){
    ggsave(path_BCDE, width=25, height=35, units="cm")  
  }
  setwd(original_dir)
}
