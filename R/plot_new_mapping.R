### Plots for paper
Normalize <- function(mat){
  mat_norm <- mat
  rowsums = rowSums(mat)
  for (i in 1:length(rowsums)){
    if(rowsums[i] != 0){
      mat_norm[i,] = mat[i,]/rowsums[i]
    }
  }
  return(mat_norm)
}


## Simulation Plots
RetrofitPlotNewMapping <- function(dir, file, save=FALSE){
  path_grid = paste("../", file, ".pdf", sep="")
  setwd(dir)
  library(ggplot2)
  library(cowplot)
  library(matrixStats)
  library(DescTools)
  ### Read synthetic data
  labels = list()
  
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
  cor_W10=sort(diag(cor(W[,cell_types],W_mod3)),decreasing=T,na.last=T)
  cor_W11=sort(diag(cor(W[,cell_types],W_mod3)),decreasing=T,na.last=T)
  cor_W12=sort(diag(cor(W[,cell_types],W_mod3)),decreasing=T,na.last=T)
  
  # df=as.data.frame(c(pmax(0, cor_W5),pmax(0, cor_W7), pmax(0, cor_W9), pmax(0, cor_W11)))
  # 
  # df$x=as.factor(c(1:length(cor_W5),1:length(cor_W7),1:length(cor_W9[]),1:length(cor_W11)))
  # df$method=factor(c(rep("N=20,M=5 (Exact 10 cell types)",length(cor_W5)),
  #                    rep("N=10,M=3 (Exact 10 cell types)",length(cor_W7)),
  #                    rep("N=10,M=3 (Extra 5 cell types)",length(cor_W9)),
  #                    rep("N=10,M=3 (Missing 2 cell types)",length(cor_W11))),
  #                  c("N=20,M=5 (Exact 10 cell types)",
  #                    "N=10,M=3 (Exact 10 cell types)",
  #                    "N=10,M=3 (Extra 5 cell types)",
  #                    "N=10,M=3 (Missing 2 cell types)"))
  # colnames(df)=c("cor","x","method")
  df=as.data.frame(rbind(cbind(1:length(cor_W1), cor_W1),
                         cbind(1:length(cor_W5), cor_W5),
                         cbind(1:length(cor_W6), cor_W6)))
  df$model=c(rep("CORRELATION_OLD", length(cor_W1)),
             rep("CORRELATION_NEW", length(cor_W5)),
             rep("MARKER",length(cor_W6)))
  df$method="N=20,M=5 (Exact 10 cell types)"
  colnames(df)=c("x","cor","model","method")
  
  temp=as.data.frame(rbind(cbind(1:length(cor_W2), cor_W2),
                           cbind(1:length(cor_W7), cor_W7),
                           cbind(1:length(cor_W8), cor_W8)))
  temp$model=c(rep("CORRELATION_OLD", length(cor_W2)),
               rep("CORRELATION_NEW", length(cor_W7)),
               rep("MARKER",length(cor_W8)))
  temp$method="N=10,M=3 (Exact 10 cell types)"
  colnames(temp)=c("x","cor","model","method")
  df=as.data.frame(rbind(df,temp))
  
  temp=as.data.frame(rbind(cbind(1:length(cor_W3), cor_W3),
                           cbind(1:length(cor_W9), cor_W9),
                           cbind(1:length(cor_W10), cor_W10)))
  temp$model=c(rep("CORRELATION_OLD", length(cor_W3)),
               rep("CORRELATION_NEW", length(cor_W9)),
               rep("MARKER",length(cor_W10)))
  temp$method="N=10,M=3 (Extra 5 cell types)"
  colnames(temp)=c("x","cor","model","method")
  df=as.data.frame(rbind(df,temp))
  
  temp=as.data.frame(rbind(cbind(1:length(cor_W4), cor_W4),
                           cbind(1:length(cor_W11), cor_W11),
                           cbind(1:length(cor_W12), cor_W12)))
  temp$model=c(rep("CORRELATION_OLD", length(cor_W4)),
               rep("CORRELATION_NEW", length(cor_W11)),
               rep("MARKER",length(cor_W12)))
  temp$method="N=10,M=3 (Missing 2 cell types)"
  colnames(temp)=c("x","cor","model","method")
  df=as.data.frame(rbind(df,temp))
  
  df$model=factor(df$model,levels=c("CORRELATION_OLD","CORRELATION_NEW", "MARKER"))
  df$method=factor(df$method,levels=c("N=20,M=5 (Exact 10 cell types)","N=10,M=3 (Exact 10 cell types)", "N=10,M=3 (Extra 5 cell types)","N=10,M=3 (Missing 2 cell types)"))
  

  
  ## Plot for RMSE(H)
  H=read.csv("N=20,M=5_H.csv",                    row.names = 1, check.names = FALSE)
  H_mod1=read.csv("N=20,M=5_loc_H_retrofit.csv",  row.names = 1, check.names = FALSE)
  H_mod2=read.csv("N=20,M=5_loc_RCTD_H.csv",      row.names = 1, check.names = FALSE)
  H_mod3=read.csv("N=20,M=5_loc_RCTD-D_H.csv",    row.names = 1, check.names = FALSE)
  H_mod4=read.csv("N=20,M=5_loc_stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("N=20,M=5_loc_NMFreg.csv",check.names = FALSE)[,-11])
  H_mod6=read.csv("N=20,M=5_loc_SPOTlight.csv",   row.names = 1, check.names = FALSE)
  H_mod7=read.csv("../mapped/N=20,M=5_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod8=read.csv("../mapped/N=20,M=5_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  
  X=read.csv("N=20,M=5_loc_X.csv",                row.names = 1, check.names = FALSE)
  W=as.matrix(read.csv("Cerebellum_W_K=10.csv",   row.names = 1, check.names = FALSE))
  
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
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model8=NULL
  
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
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model8=c(rmse_Model8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
    
  }
  
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
  X_hat7=(as.matrix(W_hat[,ind1[1:12]]) %*% diag(Thet_hat[ind1[1:12]]) +
            0.01) %*% as.matrix(H_hat[ind1[1:12],]) #no. of components that explain about 80%
  X_hat8=W %*% H_Model7
  X_hat9=W %*% H_Model8
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse8=NULL
  rmse9=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse8=c(rmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    rmse9=c(rmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
  }
  
  # df1=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
  #                         cbind(rmse_Model3,3),cbind(rmse_Model4,4),
  #                         cbind(rmse_Model5,5),cbind(rmse_Model6,6),
  #                         cbind(rmse_Model7,7),cbind(rmse_Model8,8)
  # ))
  df1=as.data.frame(rbind(cbind(rmse_Model1,1), cbind(rmse_Model7,7),cbind(rmse_Model8,8)))
  df1$method="N=20,M=5"
  colnames(df1)=c("rmse","model","method")
  
  # df2=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_H1),
  #                         cbind(seq(0,1,length.out = 1000),cor_H2),
  #                         cbind(seq(0,1,length.out = 1000),cor_H3),
  #                         cbind(seq(0,1,length.out = 1000),cor_H4),
  #                         cbind(seq(0,1,length.out = 1000),cor_H5),
  #                         cbind(seq(0,1,length.out = 1000),cor_H6)))
  df2=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_H1),
                          cbind(seq(0,1,length.out = 1000),cor_H7),
                          cbind(seq(0,1,length.out = 1000),cor_H8)))
  # df2$model=c(rep("RETROFIT",1000),rep("RCTD",1000),
  #             rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
  #             rep("NMFreg",1000),rep("SPOTlight",1000))
  df2$model=c(rep("CORRELATION_OLD",1000), 
              rep("CORRELATION_NEW",1000),
              rep("MARKER",1000))
  df2$method="N=20,M=5"
  colnames(df2)=c("x","cor","model","method")
  
  # df3=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
  #                         cbind(rmse2,2),
  #                         cbind(rmse3,3),cbind(rmse4,4),
  #                         cbind(rmse5,5),cbind(rmse6,6)
  # ))
  df3=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
                          cbind(rmse8,8),cbind(rmse9,9)
  ))
  colnames(df3)=c("rmse","model")
  df3$method="N=20,M=5"
  
  # df4=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_X1),
  #                         cbind(seq(0,1,length.out = 1000),cor_X7),
  #                         cbind(seq(0,1,length.out = 1000),cor_X2),
  #                         cbind(seq(0,1,length.out = 1000),cor_X3),
  #                         cbind(seq(0,1,length.out = 1000),cor_X4),
  #                         cbind(seq(0,1,length.out = 1000),cor_X5),
  #                         cbind(seq(0,1,length.out = 1000),cor_X6)))
  df4=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_X1),
                          cbind(seq(0,1,length.out = 1000),cor_X7),
                          cbind(seq(0,1,length.out = 1000),cor_X8),
                          cbind(seq(0,1,length.out = 1000),cor_X9)))
  df4$model=c(rep("CORRELATION_OLD",1000),
              rep("PVE=80%",1000),
              rep("CORRELATION_NEW",1000),
              rep("MARKER",1000))
  df4$method="N=20,M=5"
  colnames(df4)=c("x","cor","model","method")
  
  label = c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value, digits=3)
  )
  
  ## add p-values from above KS tests
  ann.text1=data.frame(model=2:6,rmse=0.9,
                       method=factor("N=20,M=5",
                                     levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")),
                       label=label)
  # label=c("1.7e-207","8.4e-113","3.1e-10","1.8e-281","2.4e-121"))
  
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3)
  )
  
  H_mod1=read.csv("N=20,M=5_loc_H_retrofit.csv", row.names = 1, check.names = FALSE)
  H_mod2=read.csv("N=20,M=5_loc_RCTD_H.csv", row.names = 1, check.names = FALSE)
  H_mod3=read.csv("N=20,M=5_loc_RCTD-D_H.csv", row.names = 1, check.names = FALSE)
  H_mod4=read.csv("N=20,M=5_loc_stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("N=20,M=5_loc_NMFreg.csv", check.names = FALSE)[,-11])
  H_mod6=read.csv("N=20,M=5_loc_SPOTlight.csv", row.names = 1, check.names = FALSE)
  H_mod7=read.csv("../mapped/N=20,M=5_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod8=read.csv("../mapped/N=20,M=5_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  
  ## add AUC from above results
  # ann.text2=data.frame(label=c(paste("RETROFIT:", label[[1]]),
  #                              paste("RCTD:", label[[2]]),
  #                              paste("RCTD-D:", label[[3]]),
  #                              paste("Stereoscope:", label[[4]]),
  #                              paste("NMFreg:", label[[5]]),
  #                              paste("SPOTlight:", label[[6]])),
  #                      model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
  #                              "NMFreg","SPOTlight"),
  #                      method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",
  #                                                        "Extra 5 cell types","Missing 2 cell types")))
  ann.text2=data.frame(label=c(paste("CORRELATION_OLD:", label[[1]]),
                               paste("CORRELATION_NEW:", label[[7]]),
                               paste("MARKER:", label[[8]])),
                       model=c("CORRELATION_OLD","CORRELATION_NEW","MARKER"),
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  labels[[length(labels)+1]] = c(label[[1]], label[[7]], label[[8]])
  
  label = c(
    format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=3),
    # format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=3),
    # format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=3),
    # format(ks.test(rmse1,rmse5, alternative="greater")$p.value, digits=3),
    # format(ks.test(rmse1,rmse6, alternative="greater")$p.value, digits=3)
    format(ks.test(rmse1,rmse7, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse1,rmse8, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse1,rmse9, alternative="greater")$p.value, digits=3)
  )
  
  ann.text3=data.frame(rmse=0.9,
                       model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")),
                       label=label)
  # label=c("5.1e-48","5.5e-10","0.186","5.9e-86","1.6e-28"))#,"2.1e-28"))
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3)
  )
  # ann.text4=data.frame(label = c(paste("RETROFIT:", label[[1]]),
  #                                paste("RETROFIT (PVE=80%):", label[[2]]),
  #                                paste("RCTD:", label[[3]]),
  #                                paste("RCTD-D:", label[[4]]),
  #                                paste("Stereoscope:", label[[5]]),
  #                                paste("NMFreg:", label[[6]]),
  #                                paste("SPOTlight:", label[[7]])),
  #                      model=c("RETROFIT","RETROFIT (PVE=80%)","RCTD","RCTD-Doublet","Stereoscope",
  #                              "NMFreg","SPOTlight"),
  #                      method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3",
  #                                                        "Extra 5 cell types","Missing 2 cell types")))
  ann.text4=data.frame(label = c(paste("CORRELATION_OLD:", label[[1]]),
                                 paste("PVE=80%:", label[[2]]),
                                 paste("CORRELATION_NEW:", label[[8]]),
                                 paste("MARKER:", label[[9]])
                                 ),
                       model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                       method=factor("N=20,M=5",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  labels[[length(labels)+1]] = c(label[[1]], label[[2]], label[[8]], label[[9]])
  
  # N=10,M=3
  H=read.csv("N=10,M=3_H.csv", row.names = 1, check.names = FALSE)
  H_mod1=read.csv("N=10,M=3_loc_H_retrofit.csv", row.names = 1, check.names = FALSE)
  H_mod2=read.csv("N=10,M=3_loc_RCTD_H.csv", row.names = 1, check.names = FALSE)
  H_mod3=read.csv("N=10,M=3_loc_RCTD-D_H.csv", row.names = 1, check.names = FALSE)
  H_mod4=read.csv("N=10,M=3_loc_stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("N=10,M=3_loc_NMFreg.csv", check.names = FALSE)[,-11])
  H_mod6=read.csv("N=10,M=3_loc_SPOTlight.csv", row.names = 1, check.names = FALSE)
  H_mod7=read.csv("../mapped/N=10,M=3_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod8=read.csv("../mapped/N=10,M=3_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  
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
  cor_H8=sort(diag(cor(H,H_mod8)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model8=NULL
  
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
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model8=c(rmse_Model8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
    
  }
  
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
  X_hat7=(as.matrix(W_hat[,ind1[1:13]]) %*% diag(Thet_hat[ind1[1:13]]) +
            0.01) %*% as.matrix(H_hat[ind1[1:13],])
  X_hat8=W %*% H_Model7
  X_hat9=W %*% H_Model8
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse8=NULL
  rmse9=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse8=c(rmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    rmse9=c(rmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
  }
  
  
  # temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model3,3),cbind(rmse_Model4,4),
  #                          cbind(rmse_Model5,5),cbind(rmse_Model6,6),
  #                          cbind(rmse_Model7,7),cbind(rmse_Model8,8)
  # ))
  temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model7,7),cbind(rmse_Model8,8)))
  colnames(temp)=c("rmse","model")
  temp$method="N=10,M=3"
  df1=as.data.frame(rbind(df1,temp))
  
  # temp=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_H1),
  #                          cbind(seq(0,1,length.out = 1000),cor_H2),
  #                          cbind(seq(0,1,length.out = 1000),cor_H3),
  #                          cbind(seq(0,1,length.out = 1000),cor_H4),
  #                          cbind(seq(0,1,length.out = 1000),cor_H5),
  #                          cbind(seq(0,1,length.out = 1000),cor_H6)))
  temp=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_H1),
                           cbind(seq(0,1,length.out = 1000),cor_H7),
                           cbind(seq(0,1,length.out = 1000),cor_H8)))
  # temp$model=c(rep("RETROFIT",1000),rep("RCTD",1000),
  #              rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
  #              rep("NMFreg",1000),rep("SPOTlight",1000))
  temp$model=c(rep("CORRELATION_OLD",1000),
               rep("CORRELATION_NEW",1000),
               rep("MARKER",1000))
  colnames(temp)=c("x","cor","model")
  temp$method="N=10,M=3"
  df2=as.data.frame(rbind(df2,temp))
  
  # temp=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
  #                          cbind(rmse2,2),
  #                          cbind(rmse3,3),cbind(rmse4,4),
  #                          cbind(rmse5,5),cbind(rmse6,6)
  # ))
  temp=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
                           cbind(rmse8,8),cbind(rmse9,9)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="N=10,M=3"
  df3=as.data.frame(rbind(df3,temp))
  
  # temp=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_X1),
  #                          cbind(seq(0,1,length.out = 1000),cor_X7),
  #                          cbind(seq(0,1,length.out = 1000),cor_X2),
  #                          cbind(seq(0,1,length.out = 1000),cor_X3),
  #                          cbind(seq(0,1,length.out = 1000),cor_X4),
  #                          cbind(seq(0,1,length.out = 1000),cor_X5),
  #                          cbind(seq(0,1,length.out = 1000),cor_X6)))
  temp=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_X1),
                           cbind(seq(0,1,length.out = 1000),cor_X7),
                           cbind(seq(0,1,length.out = 1000),cor_X8),
                           cbind(seq(0,1,length.out = 1000),cor_X9)))
  # temp$model=c(rep("RETROFIT",1000),
  #              rep("RETROFIT (PVE=80%)",1000),rep("RCTD",1000),
  #              rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
  #              rep("NMFreg",1000),rep("SPOTlight",1000))
  temp$model=c(rep("CORRELATION_OLD",1000),
               rep("PVE=80%",1000),
               rep("CORRELATION_NEW",1000),
               rep("MARKER",1000))
  temp$method="N=10,M=3"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  label = c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value, digits=3)
  )
  
  temp=data.frame(model=2:6,rmse=0.9,
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",
                                                    "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
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
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3)
  )
  # temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
  #                           paste("RCTD:", label[[2]]),
  #                           paste("RCTD-D:", label[[3]]),
  #                           paste("Stereoscope:", label[[4]]),
  #                           paste("NMFreg:", label[[5]]),
  #                           paste("SPOTlight:", label[[6]])),
  #                 model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
  #                         "NMFreg","SPOTlight"),
  #                 method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3",
  #                                                   "Extra 5 cell types","Missing 2 cell types")))
  temp=data.frame(label=c(paste("CORRELATION_OLD:", label[[1]]),
                               paste("CORRELATION_NEW:", label[[7]]),
                               paste("MARKER:", label[[8]])),
                  model=c("CORRELATION_OLD","CORRELATION_NEW","MARKER"),
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  labels[[length(labels)+1]] = c(label[[1]], label[[7]], label[[8]])
  # label=c(
  #   format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse5, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse6, alternative="greater")$p.value, digits=3)
  # )
  label=c(0,0,0,0)
  
  
  temp=data.frame(rmse=0.9,
                  model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3", "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("0.171","4.0e-4","0.294","0.844","0.254"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3)
  )
  # temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
  #                           paste("RETROFIT (PVE=80%):", label[[2]]),
  #                           paste("RCTD:", label[[3]]),
  #                           paste("RCTD-D:", label[[4]]),
  #                           paste("Stereoscope:", label[[5]]),
  #                           paste("NMFreg:", label[[6]]),
  #                           paste("SPOTlight:", label[[7]])),
  #                 model=c("RETROFIT","RETROFIT (PVE=80%)","RCTD","RCTD-Doublet","Stereoscope",
  #                         "NMFreg","SPOTlight"),
  #                 method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  temp=data.frame(label = c(paste("CORRELATION_OLD:", label[[1]]),
                            paste("PVE=80%:", label[[2]]),
                            paste("CORRELATION_NEW:", label[[8]]),
                            paste("MARKER:", label[[9]])),
                  model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                  method=factor("N=10,M=3",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  
  
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  labels[[length(labels)+1]] = c(label[[1]], label[[2]], label[[8]], label[[9]])
  
  ## extra 5
  H=read.csv("extra5_H.csv",                    row.names = 1, check.names = FALSE)
  H_mod1=read.csv("extra5_loc_H_retrofit.csv",  row.names = 1, check.names = FALSE)
  H_mod2=read.csv("extra5_loc_RCTD_H.csv",      row.names = 1, check.names = FALSE)
  H_mod3=read.csv("extra5_loc_RCTD-D_H.csv",    row.names = 1, check.names = FALSE)
  H_mod4=read.csv("extra5_loc_stereoscope.csv", row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("extra5_loc_NMFreg.csv",    check.names = FALSE)[,-11])
  H_mod6=read.csv("extra5_loc_SPOTlight.csv",   row.names = 1, check.names = FALSE)
  H_mod7=read.csv("../mapped/extra5_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod8=read.csv("../mapped/extra5_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  
  X=read.csv("extra5_loc_X.csv", row.names = 1, check.names = FALSE)
  
  cell_types = rownames(H)[1:5]
  H_mod1=H_mod1[cell_types,]
  H_mod2=H_mod2[cell_types,]
  H_mod3=H_mod3[cell_types,]
  H_mod4=H_mod4[cell_types,]
  H_mod5=H_mod5[cell_types,]
  H_mod6=H_mod6[cell_types,]
  H_mod7=H_mod7[cell_types,]
  H_mod8=H_mod8[cell_types,]
  
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
  cor_H8=sort(diag(cor(H[cell_types,],H_mod8)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model8=NULL
  
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
    
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model8=c(rmse_Model8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
  }
  
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
  
  W=W[,1:5]
  X_hat1=W %*% H_Model1
  X_hat2=W %*% H_Model2
  X_hat3=W %*% H_Model3
  X_hat4=W %*% H_Model4
  X_hat5=W %*% H_Model5
  X_hat6=W %*% H_Model6
  X_hat7=(as.matrix(W_hat[,ind1[1:6]]) %*% diag(Thet_hat[ind1[1:6]]) +0.01) %*% as.matrix(H_hat[ind1[1:6],])
  X_hat8=W %*% H_Model7
  X_hat9=W %*% H_Model8
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse8=NULL
  rmse9=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse8=c(rmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    rmse9=c(rmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
  }
  
  # temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model3,3),cbind(rmse_Model4,4),
  #                          cbind(rmse_Model5,5),cbind(rmse_Model6,6),
  #                          cbind(rmse_Model7,7),cbind(rmse_Model8,8)
  # ))
  temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model7,7),cbind(rmse_Model8,8)))
  colnames(temp)=c("rmse","model")
  temp$method="Extra 5 cell types"
  df1=as.data.frame(rbind(df1,temp))
  
  # temp=as.data.frame(rbind(cbind(seq(0,1,length.out=1000),cor_H1),
  #                          cbind(seq(0,1,length.out=1000),cor_H2),
  #                          cbind(seq(0,1,length.out=1000),cor_H3),
  #                          cbind(seq(0,1,length.out=1000),cor_H4),
  #                          cbind(seq(0,1,length.out=1000),cor_H5),
  #                          cbind(seq(0,1,length.out=1000),cor_H6)))
  temp=as.data.frame(rbind(cbind(seq(0,1,length.out=1000),cor_H1),
                           cbind(seq(0,1,length.out=1000),cor_H7),
                           cbind(seq(0,1,length.out=1000),cor_H8)))
  temp$model=c(rep("CORRELATION_OLD",1000),
               rep("CORRELATION_NEW",1000),
               rep("MARKER",1000))
  # temp$model=c(rep("RETROFIT",1000),
  #              rep("RCTD",1000),
  #              rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
  #              rep("NMFreg",1000),rep("SPOTlight",1000))
  colnames(temp)=c("x","cor","model")
  temp$method="Extra 5 cell types"
  df2=as.data.frame(rbind(df2,temp))
  
  # temp=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
  #                          cbind(rmse2,2),
  #                          cbind(rmse3,3),cbind(rmse4,4),
  #                          cbind(rmse5,5),cbind(rmse6,6)
  # ))
  temp=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
                           cbind(rmse8,8),cbind(rmse9,9)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="Extra 5 cell types"
  df3=as.data.frame(rbind(df3,temp))
  
  # temp=as.data.frame(rbind(cbind(seq(0,1,length.out=1000),cor_X1),
  #                          cbind(seq(0,1,length.out=1000),cor_X7),
  #                          cbind(seq(0,1,length.out=1000),cor_X2),
  #                          cbind(seq(0,1,length.out=1000),cor_X3),
  #                          cbind(seq(0,1,length.out=1000),cor_X4),
  #                          cbind(seq(0,1,length.out=1000),cor_X5),
  #                          cbind(seq(0,1,length.out=1000),cor_X6)))
  temp=as.data.frame(rbind(cbind(seq(0,1,length.out=1000),cor_X1),
                           cbind(seq(0,1,length.out=1000),cor_X7),
                           cbind(seq(0,1,length.out=1000),cor_X8),
                           cbind(seq(0,1,length.out=1000),cor_X9)))
  temp$model=c(rep("CORRELATION_OLD",1000),rep("PVE=80%",1000),
               rep("CORRELATION_NEW",1000),rep("MARKER",1000))
  temp$method="Extra 5 cell types"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  label=c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value, digits=3)
  )
  
  temp=data.frame(model=2:6,rmse=0.9,
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                              "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
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
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3)
  )
  # temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
  #                           paste("RCTD:", label[[2]]),
  #                           paste("RCTD-D:", label[[3]]),
  #                           paste("Stereoscope:", label[[4]]),
  #                           paste("NMFreg:", label[[5]]),
  #                           paste("SPOTlight:", label[[6]])),
  #                 model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
  #                         "NMFreg","SPOTlight"),
  #                 method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3",
  #                                                             "Extra 5 cell types","Missing 2 cell types")))
  temp=data.frame(label=c(paste("CORRELATION_OLD:", label[[1]]),
                          paste("CORRELATION_NEW:", label[[7]]),
                          paste("MARKER:", label[[8]])),
                  model=c("CORRELATION_OLD","CORRELATION_NEW","MARKER"),
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  labels[[length(labels)+1]] = c(label[[1]], label[[7]], label[[8]])
  # label=c(
  #   format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse5, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse6, alternative="greater")$p.value, digits=3)
  # )
  label=c(0,0,0,0)
  temp=data.frame(rmse=0.9,
                  model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")),
                  label=label)
  # label=c("9.9e-10","1.6e-08","0.294","5.3e-37","1.4e-31"))
  ann.text3=as.data.frame(rbind(ann.text3,temp))
  
  cor_X6[is.na(cor_X6)]=0
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3)
  )
  # temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
  #                           paste("RETROFIT (PVE=80%):", label[[2]]),
  #                           paste("RCTD:", label[[3]]),
  #                           paste("RCTD-D:", label[[4]]),
  #                           paste("Stereoscope:", label[[5]]),
  #                           paste("NMFreg:", label[[6]]),
  #                           paste("SPOTlight:", label[[7]])),
  #                 model=c("RETROFIT","RETROFIT (PVE=80%)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight"),
  #                 method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  temp=data.frame(label = c(paste("CORRELATION_OLD:", label[[1]]),
                            paste("PVE=80%:", label[[2]]),
                            paste("CORRELATION_NEW:", label[[8]]),
                            paste("MARKER:", label[[9]])),
                  model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                  method=factor("Extra 5 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  labels[[length(labels)+1]] = c(label[[1]], label[[2]], label[[8]], label[[9]])
  
  #### missing 2
  H=read.csv("extra5_H.csv",                  row.names = 1, check.names = FALSE)
  H_mod1=read.csv("extra5_loc_H_retrofit.csv",row.names = 1, check.names = FALSE)
  H_mod2=read.csv("miss2_loc_RCTD_H.csv",     row.names = 1, check.names = FALSE)
  H_mod3=read.csv("miss2_loc_RCTD-D_H.csv",   row.names = 1, check.names = FALSE)
  H_mod4=read.csv("miss2_loc_stereoscope.csv",row.names = 1, check.names = FALSE)
  H_mod5=t(read.csv("miss2_loc_NMFreg.csv",   check.names = FALSE)[,-11])
  H_mod6=read.csv("miss2_loc_SPOTlight.csv",  row.names = 1, check.names = FALSE)
  H_mod7=read.csv("../mapped/extra5_loc_H_mapped_cor.csv", row.names = 1, check.names = FALSE)
  H_mod8=read.csv("../mapped/extra5_loc_H_mapped_marker.csv", row.names = 1, check.names = FALSE)
  
  X=read.csv("extra5_loc_X.csv", row.names = 1, check.names = FALSE)
  
  cell_types = rownames(H)[1:5]
  cell_types_missing = rownames(H)[c(1,2,4)]
  H_mod1=H_mod1[cell_types,]
  H_mod2=H_mod2[cell_types_missing,]
  H_mod3=H_mod3[cell_types_missing,]
  H_mod4=H_mod4[cell_types_missing,]
  H_mod5=H_mod5[cell_types_missing,]
  H_mod6=H_mod6[cell_types_missing,]
  H_mod7=H_mod7[cell_types,]
  H_mod8=H_mod8[cell_types,]
  
  K=nrow(H)
  G=nrow(W)
  S=ncol(H)
  
  cor_H1=sort(diag(cor(H[cell_types,],H_mod1)),decreasing=T,na.last=T)
  cor_H2=sort(diag(cor(H[cell_types_missing,],H_mod2)),decreasing=T,na.last=T)
  cor_H3=sort(diag(cor(H[cell_types_missing,],H_mod3)),decreasing=T,na.last=T)
  cor_H4=sort(diag(cor(H[cell_types_missing,],H_mod4)),decreasing=T,na.last=T)
  cor_H5=sort(diag(cor(H[cell_types_missing,],H_mod5)),decreasing=T,na.last=T)
  cor_H6=sort(diag(cor(H[cell_types_missing,],H_mod6)),decreasing=T,na.last=T)
  cor_H7=sort(diag(cor(H[cell_types,],H_mod7)),decreasing=T,na.last=T)
  cor_H8=sort(diag(cor(H[cell_types,],H_mod8)),decreasing=T,na.last=T)
  
  H_Model1=matrix(NA, ncol=S, nrow=nrow(H_mod1))
  H_Model2=matrix(NA, ncol=S, nrow=nrow(H_mod2))
  H_Model3=matrix(NA, ncol=S, nrow=nrow(H_mod3))
  H_Model4=matrix(NA, ncol=S, nrow=nrow(H_mod4))
  H_Model5=matrix(NA, ncol=S, nrow=nrow(H_mod5))
  H_Model6=matrix(NA, ncol=S, nrow=nrow(H_mod6))
  H_Model7=matrix(NA, ncol=S, nrow=nrow(H_mod7))
  H_Model8=matrix(NA, ncol=S, nrow=nrow(H_mod8))
  H_True=H
  
  rmse_Model1=NULL
  rmse_Model2=NULL
  rmse_Model3=NULL
  rmse_Model4=NULL
  rmse_Model5=NULL
  rmse_Model6=NULL
  rmse_Model7=NULL
  rmse_Model8=NULL
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
    rmse_Model1=c(rmse_Model1,sqrt(mean((H_True[,i]-H_Model1[,i])^2)))
    rmse_Model2=c(rmse_Model2,sqrt(mean((H_True[cell_types_missing,i]-H_Model2[,i])^2)))
    rmse_Model3=c(rmse_Model3,sqrt(mean((H_True[cell_types_missing,i]-H_Model3[,i])^2)))
    rmse_Model4=c(rmse_Model4,sqrt(mean((H_True[cell_types_missing,i]-H_Model4[,i])^2)))
    rmse_Model5=c(rmse_Model5,sqrt(mean((H_True[cell_types_missing,i]-H_Model5[,i])^2)))
    rmse_Model6=c(rmse_Model6,sqrt(mean((H_True[cell_types_missing,i]-H_Model6[,i])^2)))
    rmse_Model7=c(rmse_Model7,sqrt(mean((H_True[,i]-H_Model7[,i])^2)))
    rmse_Model8=c(rmse_Model8,sqrt(mean((H_True[,i]-H_Model8[,i])^2)))
  }
  
  X_hat=(as.matrix(W_hat) %*% diag(Thet_hat) + 0.01) %*% as.matrix(H_hat)
  temp=rep(NA,length(Thet_hat))
  for(k in 1:length(temp)){
    temp[k]=sum(as.matrix(W_hat[,k] * Thet_hat[k] + 0.01) %*% as.matrix(H_hat[k,]))
  }
  prop=sort(temp/sum(temp),decreasing=T)
  ind1=sort(temp/sum(temp),decreasing=T,index.return=T)$ix
  
  
  X_hat1=W %*% H_Model1
  X_hat2=W[,cell_types_missing] %*% H_Model2
  X_hat3=W[,cell_types_missing] %*% H_Model3
  X_hat4=W[,cell_types_missing] %*% H_Model4
  X_hat5=W[,cell_types_missing] %*% H_Model5
  X_hat6=W[,cell_types_missing] %*% H_Model6
  X_hat7=(as.matrix(W_hat[,ind1[1:6]]) %*% diag(Thet_hat[ind1[1:6]]) +
            0.01) %*% as.matrix(H_hat[ind1[1:6],])
  X_hat8=W %*% H_Model7
  X_hat9=W %*% H_Model8
  
  cor_X1=sort(diag(cor(X,X_hat1)),decreasing=T,na.last=T)
  cor_X2=sort(diag(cor(X,X_hat2)),decreasing=T,na.last=T)
  cor_X3=sort(diag(cor(X,X_hat3)),decreasing=T,na.last=T)
  cor_X4=sort(diag(cor(X,X_hat4)),decreasing=T,na.last=T)
  cor_X5=sort(diag(cor(X,X_hat5)),decreasing=T,na.last=T)
  cor_X6=sort(diag(cor(X,X_hat6)),decreasing=T,na.last=T)
  cor_X7=sort(diag(cor(X,X_hat7)),decreasing=T,na.last=T)
  cor_X8=sort(diag(cor(X,X_hat8)),decreasing=T,na.last=T)
  cor_X9=sort(diag(cor(X,X_hat9)),decreasing=T,na.last=T)
  
  rmse1=NULL
  rmse2=NULL
  rmse3=NULL
  rmse4=NULL
  rmse5=NULL
  rmse6=NULL
  rmse7=NULL
  rmse8=NULL
  rmse9=NULL
  for(i in 1:S){
    rmse1=c(rmse1,sqrt(mean(((X[,i]-X_hat1[,i])/sd(X[,i]))^2)))
    rmse2=c(rmse2,sqrt(mean(((X[,i]-X_hat2[,i])/sd(X[,i]))^2)))
    rmse3=c(rmse3,sqrt(mean(((X[,i]-X_hat3[,i])/sd(X[,i]))^2)))
    rmse4=c(rmse4,sqrt(mean(((X[,i]-X_hat4[,i])/sd(X[,i]))^2)))
    rmse5=c(rmse5,sqrt(mean(((X[,i]-X_hat5[,i])/sd(X[,i]))^2)))
    rmse6=c(rmse6,sqrt(mean(((X[,i]-X_hat6[,i])/sd(X[,i]))^2)))
    rmse7=c(rmse7,sqrt(mean(((X[,i]-X_hat7[,i])/sd(X[,i]))^2)))
    rmse8=c(rmse8,sqrt(mean(((X[,i]-X_hat8[,i])/sd(X[,i]))^2)))
    rmse9=c(rmse9,sqrt(mean(((X[,i]-X_hat9[,i])/sd(X[,i]))^2)))
  }
  
  # temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model2,2),
  #                          cbind(rmse_Model3,3),cbind(rmse_Model4,4),
  #                          cbind(rmse_Model5,5),cbind(rmse_Model6,6),
  #                          cbind(rmse_Model7,7),cbind(rmse_Model8,8)
  # ))
  temp=as.data.frame(rbind(cbind(rmse_Model1,1),cbind(rmse_Model7,7),cbind(rmse_Model8,8)))
  colnames(temp)=c("rmse","model")
  temp$method="Missing 2 cell types"
  df1=as.data.frame(rbind(df1,temp))
  
  # temp=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_H1),
  #                          cbind(seq(0,1,length.out = 1000),cor_H2),
  #                          cbind(seq(0,1,length.out = 1000),cor_H3),
  #                          cbind(seq(0,1,length.out = 1000),cor_H4),
  #                          cbind(seq(0,1,length.out = 1000),cor_H5),
  #                          cbind(seq(0,1,length.out = 1000),cor_H6)))
  temp=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_H1),
                           cbind(seq(0,1,length.out = 1000),cor_H7),
                           cbind(seq(0,1,length.out = 1000),cor_H8)))
  # temp$model=c(rep("RETROFIT",1000),rep("RCTD",1000),
  #              rep("RCTD-Doublet",1000),rep("Stereoscope",1000),
  #              rep("NMFreg",1000),rep("SPOTlight",1000))
  temp$model=c(rep("CORRELATION_OLD",1000),
               rep("CORRELATION_NEW",1000),
               rep("MARKER",1000))
  colnames(temp)=c("x","cor","model")
  temp$method="Missing 2 cell types"
  df2=as.data.frame(rbind(df2,temp))
  
  # temp=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
  #                          cbind(rmse2,2),
  #                          cbind(rmse3,3),cbind(rmse4,4),
  #                          cbind(rmse5,5),cbind(rmse6,6)
  # ))
  temp=as.data.frame(rbind(cbind(rmse1,1),cbind(rmse7,7),
                           cbind(rmse8,8),cbind(rmse9,9)
  ))
  colnames(temp)=c("rmse","model")
  temp$method="Missing 2 cell types"
  df3=as.data.frame(rbind(df3,temp))
  
  temp=as.data.frame(rbind(cbind(seq(0,1,length.out = 1000),cor_X1),
                           cbind(seq(0,1,length.out = 1000),cor_X7),
                           cbind(seq(0,1,length.out = 1000),cor_X8),
                           cbind(seq(0,1,length.out = 1000),cor_X9)))
  temp$model=c(rep("CORRELATION_OLD",1000),rep("PVE=80%",1000),
               rep("CORRELATION_NEW",1000),rep("MARKER",1000))
  temp$method="Missing 2 cell types"
  colnames(temp)=c("x","cor","model","method")
  df4=as.data.frame(rbind(df4,temp))
  
  label=c(
    format(ks.test(rmse_Model1,rmse_Model2, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model3, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model4, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model5, alternative="greater")$p.value, digits=3),
    format(ks.test(rmse_Model1,rmse_Model6, alternative="greater")$p.value, digits=3)
  )
  
  temp=data.frame(model=2:6,rmse=0.9,
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
                                                                "Extra 5 cell types","Missing 2 cell types")),
                  label=label)
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
  label = c(
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H1),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3)
  )
  # temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
  #                           paste("RCTD:", label[[2]]),
  #                           paste("RCTD-D:", label[[3]]),
  #                           paste("Stereoscope:", label[[4]]),
  #                           paste("NMFreg:", label[[5]]),
  #                           paste("SPOTlight:", label[[6]])),
  #                 model=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope",
  #                         "NMFreg","SPOTlight"),
  #                 method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3",
  #                                                               "Extra 5 cell types","Missing 2 cell types")))
  temp=data.frame(label=c(paste("CORRELATION_OLD:", label[[1]]),
                          paste("CORRELATION_NEW:", label[[7]]),
                          paste("MARKER:", label[[8]])),
                  model=c("CORRELATION_OLD","CORRELATION_NEW","MARKER"),
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  ann.text2=as.data.frame(rbind(ann.text2,temp))
  labels[[length(labels)+1]] = c(label[[1]], label[[7]], label[[8]])
  # label=c(
  #   format(ks.test(rmse1,rmse2, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse3, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse4, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse5, alternative="greater")$p.value, digits=3),
  #   format(ks.test(rmse1,rmse6, alternative="greater")$p.value, digits=3)
  # )
  label=c(0,0,0,0)
  temp=data.frame(rmse=0.9,
                  model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")),
                  label=label)
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
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X7),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X2),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X3),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X4),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X5),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X6),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X8),digits=3),
    round(AUC(x=seq(0,1,length.out = 1000),y=cor_X9),digits=3)
  )
  # temp=data.frame(label = c(paste("RETROFIT:", label[[1]]),
  #                           paste("RETROFIT (PVE=80%):", label[[2]]),
  #                           paste("RCTD:", label[[3]]),
  #                           paste("RCTD-D:", label[[4]]),
  #                           paste("Stereoscope:", label[[5]]),
  #                           paste("NMFreg:", label[[6]]),
  #                           paste("SPOTlight:", label[[7]])),
  #                 model=c("RETROFIT","RETROFIT (PVE=80%)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight"),
  #                 method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  temp=data.frame(label = c(paste("CORRELATION_OLD:", label[[1]]),
                            paste("PVE=80%:", label[[2]]),
                            paste("CORRELATION_NEW:", label[[8]]),
                            paste("MARKER:", label[[9]])),
                  model=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"),
                  method=factor("Missing 2 cell types",levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types")))
  ann.text4=as.data.frame(rbind(ann.text4,temp))
  labels[[length(labels)+1]] = c(label[[1]], label[[2]], label[[8]], label[[9]])
  
  df1$model=as.factor(df1$model)
  # df2$model=factor(df2$model,levels=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight"))
  df2$model=factor(df2$model,levels=c("CORRELATION_OLD","CORRELATION_NEW", "MARKER"))
  df1$method=factor(df1$method,levels=c("N=20,M=5","N=10,M=3", "Extra 5 cell types","Missing 2 cell types"))
  df2$method=factor(df2$method,levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types"))
  
  # df3$model=factor(df3$model,levels=c(1,7,2,3,4,5,6))
  df3$model=factor(df3$model,levels=c(1,7,8,9))
  # df4$model=factor(df4$model,levels=c("RETROFIT","RETROFIT (PVE=80%)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight"))
  df4$model=factor(df4$model,levels=c("CORRELATION_OLD","PVE=80%","CORRELATION_NEW","MARKER"))
  df3$method=factor(df3$method,levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types"))
  df4$method=factor(df4$method,levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types"))
  
  ann.text1$model=as.factor(ann.text1$model)
  ann.text2$model=as.factor(ann.text2$model)
  ann.text3$model=as.factor(ann.text3$model)
  ann.text4$model=factor(ann.text4$model,levels=c("N=20,M=5","N=10,M=3","Extra 5 cell types","Missing 2 cell types"))
  
  # ann.text2$cor=rep(c(0.55,0.45,0.35,0.25,0.15,0.05),4)
  ann.text2$cor=rep(c(0.55,0.35,0.15),4)
  ann.text2$x=0.05
  
  # ann.text3$rmse=c(rep(3,5),rep(58,5),rep(40,5),rep(15,5))
  ann.text3$rmse=c(rep(3,4),rep(58,4),rep(40,4),rep(15,4))
  
  ann.text4$cor=rep(c(0.65,0.45,0.25,0.05),4)
  ann.text4$x=0.05
  
  # plot for cor(W)
  # p1=ggplot(df,aes(x=x,y=cor)) +
  #   geom_line(aes(x=x,y=cor,group=method)) + geom_point(aes(x=x,y=cor,group=method))+
  #   xlab("Cell Types") +
  #   ylab(expression(paste("Correlation (",W^0,",",widetilde(W),")"))) +
  #   theme_bw()+
  #   ylim(0,1)+
  #   theme(axis.title=element_text(size=9), axis.title.x = element_text(vjust=10),
  #         axis.text.x = element_text(size=7),
  #         plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm")) +
  #   facet_grid(cols=vars(method),scales = "free_x")
  p1=ggplot(df,aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + 
    geom_point(aes(color=model))+
    # xlab("Cell Types") +
    ylab(expression(paste("Correlation (",W^0,",",widetilde(W),")"))) +
    theme_bw()+
    scale_color_manual(values = c('grey30', '#56B4E9', '#D55E00'),
                       labels=c("CORRELATION_OLD","CORRELATION_NEW", "MARKER")) +
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
    scale_fill_manual(values=c('grey30','#56B4E9','#D55E00'))+
    scale_x_discrete(labels=c("CORRELATION_OLD","CORRELATION_NEW","MARKER"
    )) +facet_grid(cols=vars(method))
  # +
  #   geom_text(data = ann.text1,aes(label = label),size=2.3,
  #             color=c("black","black","black","black","black",
  #                     "black","black","grey60","black","black",
  #                     "black","black","black","black","black",
  #                     "black","grey60","black","black","black"))

  # p2=ggplot(df1, aes(x=model, y=rmse, fill=model)) +
  #   ylab(expression(paste("RMSE (",H,",",widetilde(H),")"))) +
  #   xlab("") +
  #   geom_violin() + theme_bw()+ ylim(0,1)+
  #   # geom_boxplot(width=0.1)+
  #   theme(legend.position = "none",axis.title=element_text(size=9),
  #         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
  #                                    size=7),
  #         strip.text = element_blank(),
  #         plot.margin = unit(c(0, 0.1, -0.1, 0.1), "cm"))+
  #   scale_fill_manual(values=c('grey30','#56B4E9','#0072B2',
  #                              '#009E73',"#D55E00", "#CC79A7", "#CC79A7", "#CC79A7"))+
  #   scale_x_discrete(labels=c("RETROFIT","RCTD","RCTD-D","Stereoscope","NMFreg","SPOTlight", "temp1", "temp2"
  #   )) +facet_grid(cols=vars(method))+
  #   geom_text(data = ann.text1,aes(label = label),size=2.3,
  #             color=c("black","black","black","black","black",
  #                     "black","black","grey60","black","black",
  #                     "black","black","black","black","black",
  #                     "black","grey60","black","black","black"))
  
  # p3=ggplot(df2, aes(x=x,y=cor, group=model)) +
  #   geom_line(aes(color=model)) + xlab("Percentile") +
  #   ylab(expression(paste("Correlation (",H,",",widetilde(H),")"))) +
  #   theme_bw()+
  #   scale_color_manual(values = c('grey30','#56B4E9','#0072B2',
  #                                 '#009E73',"#D55E00", "#CC79A7"),
  #                      labels=c("RETROFIT","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight")) +
  #   ylim(0,1)+facet_grid(cols=vars(method)) +
  #   theme(legend.position = "none",strip.text = element_blank(),
  #         axis.title=element_text(size=9),axis.title.x = element_text(vjust=10),
  #         axis.text.x = element_text(size=7),
  #         plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
  #   geom_text(data = ann.text2,aes(label = label),size=2.2,
  #             hjust=0,color=rep(c('grey30','#56B4E9','#0072B2',
  #                                 '#009E73',"#D55E00", "#CC79A7"),4))
  p3=ggplot(df2, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model)) + xlab("Percentile") +
    ylab(expression(paste("Correlation (",H,",",widetilde(H),")"))) +
    theme_bw()+
    scale_color_manual(values = c('grey30','#56B4E9','#D55E00'),
                       labels=c("CORRELATION_OLD","CORRELATION_NEW","MARKER")) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.title=element_text(size=9),axis.title.x = element_text(vjust=10),
          axis.text.x = element_text(size=7),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm")) +
    geom_text(data = ann.text2,aes(label = label),size=2.2, hjust=0, 
              color=rep(c('grey30','#56B4E9','#D55E00'),4))
  
  # p4=ggplot(df3, aes(x=model, y=rmse, fill=model)) +
  #   ylab(expression(paste("NRMSE (",X,",",widetilde(X),")"))) +
  #   xlab("") +
  #   geom_violin() + theme_bw()+ #ylim(0,60)+
  #   # geom_boxplot(width=0.1)+
  #   theme(legend.position = "none",axis.title=element_text(size=9),
  #         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
  #                                    size=7),
  #         strip.text = element_blank(),
  #         plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
  #   scale_fill_manual(values=c('grey35','grey20','#56B4E9','#0072B2',
  #                              '#009E73',"#D55E00", "#CC79A7"))+
  #   scale_x_discrete(labels=c("RETROFIT","RETROFIT\n(PVE=80%)","RCTD","RCTD-D","Stereoscope","NMFreg","SPOTlight"
  #   )) +facet_wrap(~method,scales="free",nrow=1,ncol=4)+
  #   geom_text(data = ann.text3,aes(label = label),size=2.3,
  #             color=c("black","black","grey60","black","black",
  #                     "grey60","black","grey60","grey60","grey60",
  #                     "black","black","grey60","black","black",
  #                     "grey60","black","black","grey60","black"))
  
  p4=ggplot(df3, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("NRMSE (",X,",",widetilde(X),")"))) +
    xlab("") +
    geom_violin() + theme_bw() +
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c('grey35','grey20','#56B4E9','#0072B2'))+
    scale_x_discrete(labels=c("CORRELATION_OLD","PVE80%","CORRELATION_NEW","MARKER"
    )) +facet_wrap(~method,scales="free",nrow=1,ncol=4) # +
    geom_text(data = ann.text3,aes(label = label),size=2.3,
              color=rep(c('grey30','grey30','#56B4E9','#D55E00'),4))
              # color=c("black","black","grey60","black","black",
              #         "grey60","black","grey60","grey60","grey60",
              #         "black","black","grey60","black","black",
              #         "grey60","black","black","grey60","black"))
  
  # p5=ggplot(df4, aes(x=x,y=cor, group=model)) +
  #   geom_line(aes(color=model,linetype=model)) + xlab("Percentile") +
  #   ylab(expression(paste("Correlation (",X,",",widetilde(X),")"))) +
  #   theme_bw()+
  #   scale_color_manual(values = c('grey30','grey30','#56B4E9','#0072B2',
  #                                 '#009E73',"#D55E00", "#CC79A7"),
  #                      labels=c("RETROFIT","RETROFIT\n(PVE=80%)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight")) +
  #   scale_linetype_manual(values = c('solid','dashed','solid','solid',
  #                                    'solid',"solid", "solid"),
  #                         labels=c("RETROFIT","RETROFIT\n(PVE=80%)","RCTD","RCTD-Doublet","Stereoscope","NMFreg","SPOTlight")) +
  #   ylim(0,1)+facet_grid(cols=vars(method)) +
  #   theme(legend.position = "none",strip.text = element_blank(),
  #         axis.text.x = element_text(size=7),axis.title.x = element_text(vjust=10),
  #         #axis.title.x = element_text(vjust = 12),
  #         axis.title=element_text(size=9),
  #         plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
  #   geom_text(data = ann.text4,aes(label = label),size=2.2,hjust=0,
  #             color=rep(c('grey35','grey20','#56B4E9','#0072B2',
  #                         '#009E73',"#D55E00", "#CC79A7"),4))
  p5=ggplot(df4, aes(x=x,y=cor, group=model)) +
    geom_line(aes(color=model,linetype=model)) + xlab("Percentile") +
    ylab(expression(paste("Correlation (",X,",",widetilde(X),")"))) +
    theme_bw()+
    scale_color_manual(values = c('grey30','grey30','#56B4E9','#D55E00'),
                       labels=c("CORRELATION_OLD","PVE80%","CORRELATION_NEW","MARKER")) +
    scale_linetype_manual(values = c('solid','dashed','solid','solid'),
                          labels=c("CORRELATION_OLD","PVE80%","CORRELATION_NEW","MARKER")) +
    ylim(0,1)+facet_grid(cols=vars(method)) +
    theme(legend.position = "none",strip.text = element_blank(),
          axis.text.x = element_text(size=7),axis.title.x = element_text(vjust=10),
          #axis.title.x = element_text(vjust = 12),
          axis.title=element_text(size=9),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm")) +
    geom_text(data = ann.text4,aes(label = label),size=2.2,hjust=0,
              color=rep(c('grey30','grey30','#56B4E9','#D55E00'),4))
  
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
  
  ann.text3$rmse=rep(3.2,16)
  df5=df3[complete.cases(df3),]
  quantile(df5$rmse[df5$method=="N=10,M=3"],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df5$rmse[df5$method=="Extra 5 cell types"],c(0.75,0.8,0.9), na.rm = TRUE)
  quantile(df5$rmse[df5$method=="Missing 2 cell types"],c(0.75,0.8,0.9), na.rm = TRUE)
  
  df5=df5[df5$rmse<3,]
  # p4=ggplot(df5, aes(x=model, y=rmse, fill=model)) +
  #   ylab(expression(paste("NRMSE (",X,",",hat(X),")"))) +
  #   xlab("") +
  #   geom_violin() + theme_bw()+ ylim(0,3.5)+
  #   # geom_boxplot(width=0.1)+
  #   theme(legend.position = "none",axis.title=element_text(size=9),
  #         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=7),
  #         strip.text = element_blank(),
  #         plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
  #   scale_fill_manual(values=c('grey35','grey20','#56B4E9','#0072B2',
  #                              '#009E73',"#D55E00", "#CC79A7"))+
  #   scale_x_discrete(labels=c("RETROFIT","RETROFIT\n(PVE=80%)","RCTD","RCTD-D","Stereoscope","NMFreg","SPOTlight"
  #   )) + facet_grid(cols=vars(method)) +
  #   geom_text(data = ann.text3,aes(label = label),size=2.3,
  #             color=c("black","black","grey60","black","black",
  #                     "grey60","black","grey60","grey60","grey60",
  #                     "black","black","grey60","black","black",
  #                     "grey60","black","black","grey60","black"))
  p4=ggplot(df5, aes(x=model, y=rmse, fill=model)) +
    ylab(expression(paste("NRMSE (",X,",",hat(X),")"))) +
    xlab("") +
    geom_violin() + theme_bw()+ ylim(0,3.5)+
    theme(legend.position = "none",axis.title=element_text(size=9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=7),
          strip.text = element_blank(),
          plot.margin = unit(c(-0.5, 0.1, -0.1, 0.1), "cm"))+
    scale_fill_manual(values=c('grey30','grey30','#56B4E9','#D55E00'))+
    scale_x_discrete(labels=c("CORRELATION_OLD","PVE80%","CORRELATION_NEW","MARKER"
    )) + facet_grid(cols=vars(method)) # +
    geom_text(data = ann.text3,aes(label = label),size=2.3,
              color=rep(c('grey30','grey30','#56B4E9','#D55E00'),4))
  
  
  plot_grid(
    p1,NULL,p2,NULL,p3,NULL, p4,NULL,p5,
    rel_heights = c(1, -0.25, 1, -0.05, 1, -0.1, 1, -0.05, 1),
    align='hv',
    labels=c("A","","B","","C","","D","","E"),# vjust=12,
    #label_colour="red",
    nrow=9, ncol=1)
  
  if(save){
    print('ggsave')
    ggsave(path_grid, width=25, height=25, units="cm")
  }
  
  return (labels)
}
