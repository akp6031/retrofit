MouseBrainSimulation <- function(folder) {
  setwd("~/Research/retrofit/retrofit/R")
  # setwd("/storage/home/akp6031/retrofit/R")
  
  # file_source = "N=20,M=5_loc_X"
  # file_source = "N=10,M=3_loc_X"
  # file_source = "extra5_loc_X"
  file_source = "Cerebellum_X"
  ref_file_source = "Cerebellum_W_X_filtered"
  ref_marker_file_source = "Cerebellum_pseudomarkers"
  L = 20
  K = 11
  in_dir = "../results"
  out_dir = "../../Mouse/results"
  iterations=5000
  # alpha_w_0 = 0.03
  # beta_w_0 = 0.0002
  # alpha_th_0 = 1
  lambda=0
  seed=1
  
  if(!is.null(folder)){
    out_dir = paste(out_dir, folder, sep="/")
    dir.create(out_dir)
  }
  ref_w_path = paste(in_dir, paste(ref_file_source, ".csv", sep=""), sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  
  ref_marker_path = paste(in_dir, paste(ref_marker_file_source, ".csv", sep=""), sep="/")
  print(ref_marker_path)
  ref_marker_d=read.csv(ref_marker_path, row.names = 1)
  ref_marker = list()
  for(r in 1:nrow(ref_marker_d)){
    gene = ref_marker_d[[1]][r]
    cell_type = ref_marker_d[[2]][r]
    if(is.null(ref_marker[[cell_type]])){
      ref_marker[[cell_type]] = c()
    }
    ref_marker[[cell_type]] = c(ref_marker[[cell_type]], gene)
  }
  
  in_path = paste(in_dir, paste(file_source, ".csv", sep=""), sep="/")
  x=read.csv(in_path, row.names = 1)
  
  result = RetrofitMain(x, ref_cor=ref_w, ref_marker=ref_marker, iterations=iterations, L=L, K=K,
                        # alpha_w_0=alpha_w_0, beta_w_0=beta_w_0, alpha_th_0=alpha_th_0,
                        lambda=lambda, seed=seed)
  
  decomp_h_path = paste(out_dir, paste(file_source,"__decomp_h.csv", sep=""), sep="/")
  decomp_w_path = paste(out_dir, paste(file_source,"__decomp_w.csv", sep=""), sep="/")
  decomp_th_path = paste(out_dir, paste(file_source,"__decomp_th.csv", sep=""), sep="/")
  match_w_path = paste(out_dir, paste(file_source,"__map_cor_w.csv", sep=""), sep="/")
  match_h_path = paste(out_dir, paste(file_source,"__map_cor_h.csv", sep=""), sep="/")
  match_marker_w_path = paste(out_dir, paste(file_source,"__map_marker_w.csv", sep=""), sep="/")
  match_marker_h_path = paste(out_dir, paste(file_source,"__map_marker_h.csv", sep=""), sep="/")
  cor_w_path = paste(out_dir, paste(file_source,"__cor_w.csv", sep=""), sep="/")
  
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_cor_map, match_h_path)
  write.csv(result$w_cor_map, match_w_path)
  write.csv(result$h_marker_map, match_marker_h_path)
  write.csv(result$w_marker_map, match_marker_w_path)
  
  paths <- list(decomp_w_path, decomp_h_path, decomp_th_path, match_w_path, match_h_path)
  print(paths)
  return(paths)
}


RetrofitMainPlotSimulation <- function(iterations, seed) {
  setwd("~/Research/retrofit/retrofit/results/local")
  folder = paste('iter_', iterations, "_seed", seed[[1]], "-",seed[[2]], "-", seed[[3]], sep="")
  dir.create(folder)
  dir.create("supplementary")
  dir.create("datasets")
  dir.create("decomposed")
  dir.create("mapped")
  file.copy("plot", folder, recursive=TRUE)
  
  # file_source = "N=10,M=3_loc_X"
  # file_source = "extra5_loc_X"
  
  # common
  output_dir = paste("~/Research/retrofit/retrofit/results/local", folder, "plot", sep="/")
  iterations=iterations
  ref_file_source = "Cerebellum_W_K=10.csv"
  ref_w_path = paste(output_dir, ref_file_source, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  
  # scenario 
  file_source = "N=20,M=5_loc_"
  L = 20
  K = 10
  
  # run main
  in_path = paste(output_dir, paste(file_source,"X.csv", sep=""), sep="/")
  x=read.csv(in_path, row.names = 1)
  setwd("~/Research/retrofit/retrofit/R")
  result = RetrofitMain(x, ref_w, iterations=iterations, L=L, K=K, seed=seed[[1]])
  # w=w, 
  # h=h, 
  # th=th, 
  # w_cor_map=w_cor_map, 
  # h_cor_map=h_cor_map, 
  # w_cor_map_correlation=w_cor_map_correlation,
  # w_marker_map=w_marker_map,
  # h_marker_map=h_marker_map
  # save files
  decomp_h_path   = paste(output_dir, "N=20,M=5_loc_H_hat_L=20.csv", sep="/")
  decomp_w_path   = paste(output_dir, "N=20,M=5_loc_W_hat_L=20.csv", sep="/")
  decomp_th_path  = paste(output_dir, "N=20,M=5_loc_Theta_hat_L=20.csv", sep="/")
  match_w_path    = paste(output_dir, "N=20,M=5_loc_W_retrofit.csv", sep="/")
  match_h_path    = paste(output_dir, "N=20,M=5_loc_H_retrofit.csv", sep="/")
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_match, match_h_path)
  write.csv(result$w_match, match_w_path)
  out_data__x_path = paste("datasets/",file_source,"X.csv", sep="")
  out_data__ref_cor_path = "datasets/ref_cor__Cerebellum_W_K=10.csv"
  out_data__ref_marker_path = "datasets/ref_marker__Cerebellum_pseudomarkers.csv"
  
  
  
  #   out_sup__cor_w_path = paste("supplementary/", file_source,"W_mapped_cor__not_normed.csv", sep="")
  #   out_sup__cor_h_path = paste("supplementary/",file_source,"H_mapped_cor__not_normed.csv", sep="")
  #   out_sup__cor_w_check_old_path = paste("supplementary/",file_source,"W_mapped_cor__check_not_normed.csv", sep="")
  #   out_sup__cor_w_check_path = paste("supplementary/",file_source,"W_mapped_cor__check.csv", sep="")
  
  
  
  
  
  #   out_decomposed__w_path = paste("decomposed/",file_source,"W_decomposed.csv", sep="")
  #   out_decomposed__h_path = paste("decomposed/",file_source,"H_decomposed.csv", sep="")
  
  #   out_mapped__cor_w_path = paste("mapped/",file_source,"W_mapped_cor.csv", sep="")
  #   out_mapped__cor_h_path = paste("mapped/",file_source,"H_mapped_cor.csv", sep="")
  #   out_mapped__marker_w_path = paste("mapped/",file_source,"W_mapped_marker.csv", sep="")
  #   out_mapped__marker_h_path = paste("mapped/",file_source,"H_mapped_marker.csv", sep="")
  # scenario 
  file_source = "N=10,M=3_loc_X.csv"
  L = 20
  K = 10
  
  # run main
  in_path = paste(output_dir, file_source, sep="/")
  x=read.csv(in_path, row.names = 1)
  setwd("~/Research/retrofit/retrofit/R")
  result = RetrofitMain(x, ref_w, iterations=iterations, L=L, K=K, seed=seed[[2]])
  # save files
  decomp_h_path   = paste(output_dir, "N=10,M=3_loc_H_hat_L=20.csv", sep="/")
  decomp_w_path   = paste(output_dir, "N=10,M=3_loc_W_hat_L=20.csv", sep="/")
  decomp_th_path  = paste(output_dir, "N=10,M=3_loc_Theta_hat_L=20.csv", sep="/")
  match_w_path    = paste(output_dir, "N=10,M=3_loc_W_retrofit.csv", sep="/")
  match_h_path    = paste(output_dir, "N=10,M=3_loc_H_retrofit.csv", sep="/")
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_match, match_h_path)
  write.csv(result$w_match, match_w_path)
  
  # scenario 
  file_source = "extra5_loc_X.csv"
  L = 10
  K = 5
  
  # run main
  in_path = paste(output_dir, file_source, sep="/")
  x=read.csv(in_path, row.names = 1)
  setwd("~/Research/retrofit/retrofit/R")
  result = RetrofitMain(x, ref_w, iterations=iterations, L=L, K=K, seed=seed[[3]])
  # save files
  decomp_h_path   = paste(output_dir, "extra5_loc_H_hat_L=10.csv", sep="/")
  decomp_w_path   = paste(output_dir, "extra5_loc_W_hat_L=10.csv", sep="/")
  decomp_th_path  = paste(output_dir, "extra5_loc_Theta_hat_L=10.csv", sep="/")
  match_w_path    = paste(output_dir, "extra5_loc_W_retrofit.csv", sep="/")
  match_h_path    = paste(output_dir, "extra5_loc_H_retrofit.csv", sep="/")
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_match, match_h_path)
  write.csv(result$w_match, match_w_path)
  
  # plot
  RetrofitPlot(dir=output_dir, file=folder)
}

RunScenarios <- function(){
  # for(i in 16:30){
  #   RetrofitMainPlotSimulation(4000, i)  
  # }
  # RetrofitMainPlotSimulation(4000, c(28, 12, 4))
  SecondaryPlotSimulation(20,0)
  SecondaryPlotSimulation(20,0.02)
  SecondaryPlotSimulation(16,0.01)
  SecondaryPlotSimulation(16,0.02)
  SecondaryPlotSimulation(16,0.03)
  SecondaryPlotSimulation(13,0.01)
  SecondaryPlotSimulation(13,0.02)
  SecondaryPlotSimulation(13,0.03)
}

ColonSimulation <- function(name, iterations) {
  config = list(A1=c("A1_X", "A1_X_ref"),
                A2=c("A2_X", "A2_X_ref"),
                A3=c("A3_X", "A3_X_ref"),
                A4=c("A4_X", "A4_X_ref"),
                A8=c("A8_X", "A8_X_ref"),
                A9=c("A9_X", "A9_X_ref"))
  
  w = c()
  iterations = iterations
  
  for (seed in 1:2){
    file_source = config[[name]][1]
    ref_file_source = config[[name]][2]
    in_path = paste(file_source, ".csv", sep="")
    in_ref_path = paste(ref_file_source, ".csv", sep="")
    
    # load
    setwd("~/Research/retrofit/retrofit/results")
    ref_w=read.csv(in_ref_path, row.names = 1)
    x=read.csv(in_path, row.names = 1)
    
    # run main
    L = 16
    K = 8
    result = RetrofitMain(x, ref_w, iterations=iterations, L=L, K=K, seed=seed)  
    
    # save files
    setwd("~/Research/retrofit/retrofit/results/local")
    output_dir = paste("colon_", file_source, "_iter_", iterations, "_seed", seed, sep="")
    dir.create(output_dir)
    decomp_h_path   = paste(output_dir, paste(file_source, "_decomposed_H", ".csv", sep=""), sep="/")
    decomp_w_path   = paste(output_dir, paste(file_source, "_decomposed_W", ".csv", sep=""), sep="/")
    decomp_th_path  = paste(output_dir, paste(file_source, "_decomposed_TH", ".csv", sep=""), sep="/")
    match_w_path    = paste(output_dir, paste(file_source, "_matched_W", ".csv", sep=""), sep="/")
    match_h_path    = paste(output_dir, paste(file_source, "_matched_H", ".csv", sep=""), sep="/")
    cor_w_path      = paste(output_dir, paste(file_source, "_cor_W", ".csv", sep=""), sep="/")
    write.csv(result$h, decomp_h_path)
    write.csv(result$w, decomp_w_path)
    write.csv(result$th, decomp_th_path)
    write.csv(result$h_match, match_h_path)
    write.csv(result$w_match, match_w_path)
    write.csv(result$w_cor, cor_w_path)
    
    w[[seed]] = result$w_match 
  }
  
  cor_seed1_v_2 = cor(w[[1]], w[[2]])
  cor_seed1_v_2_path = paste(output_dir, "cor_seed1_v_2.csv", sep="/")
  write.csv(cor_seed1_v_2, cor_seed1_v_2_path)
}


SecondaryPlotSimulation <- function(L, lambda) {
  setwd("~/Research/retrofit/retrofit/results/local")
  iterations = 4000
  seed = c(1,12,11)
  folder = paste('secondary_L_', L, '_lambda_', lambda, '_iter_', iterations, "_seed", seed[[1]], "-",seed[[2]], "-", seed[[3]], sep="")
  dir.create(folder)
  file.copy("plot", folder, recursive=TRUE)
  
  # common
  input_dir = paste(folder, "plot", sep="/")
  output_dir = input_dir
  ref_file_source = "Cerebellum_W_K=10.csv"
  ref_w_path = paste(input_dir, ref_file_source, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1, check.names = FALSE)
  
  # scenario 
  file_source = "N=20,M=5_loc_X.csv"
  L = L
  K = 10
  
  # run main
  in_path = paste(input_dir, file_source, sep="/")
  x=read.csv(in_path, row.names = 1)
  result = RetrofitMain(x, ref_cor = ref_w, iterations=iterations, L=L, K=K, seed=seed[[1]], lambda = lambda)
  # save files
  decomp_h_path   = paste(output_dir, "N=20,M=5_loc_H_hat_L=20.csv", sep="/")
  decomp_w_path   = paste(output_dir, "N=20,M=5_loc_W_hat_L=20.csv", sep="/")
  decomp_th_path  = paste(output_dir, "N=20,M=5_loc_Theta_hat_L=20.csv", sep="/")
  match_w_path    = paste(output_dir, "N=20,M=5_loc_W_retrofit.csv", sep="/")
  match_h_path    = paste(output_dir, "N=20,M=5_loc_H_retrofit.csv", sep="/")
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_cor_map, match_h_path)
  write.csv(result$w_cor_map, match_w_path)
  
  # scenario 
  file_source = "N=10,M=3_loc_X.csv"
  L = L
  K = 10
  
  # run main
  in_path = paste(input_dir, file_source, sep="/")
  x=read.csv(in_path, row.names = 1)
  result = RetrofitMain(x, ref_cor = ref_w, iterations=iterations, L=L, K=K, seed=seed[[2]], lambda=lambda)
  # save files
  decomp_h_path   = paste(output_dir, "N=10,M=3_loc_H_hat_L=20.csv", sep="/")
  decomp_w_path   = paste(output_dir, "N=10,M=3_loc_W_hat_L=20.csv", sep="/")
  decomp_th_path  = paste(output_dir, "N=10,M=3_loc_Theta_hat_L=20.csv", sep="/")
  match_w_path    = paste(output_dir, "N=10,M=3_loc_W_retrofit.csv", sep="/")
  match_h_path    = paste(output_dir, "N=10,M=3_loc_H_retrofit.csv", sep="/")
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_cor_map, match_h_path)
  write.csv(result$w_cor_map, match_w_path)
  
  # scenario 
  file_source = "extra5_loc_X.csv"
  L = max(L-5, 10)
  K = 10
  
  # run main
  in_path = paste(input_dir, file_source, sep="/")
  x=read.csv(in_path, row.names = 1)
  result = RetrofitMain(x, ref_cor = ref_w, iterations=iterations, L=L, K=K, seed=seed[[3]], lambda=lambda)
  # save files
  decomp_h_path   = paste(output_dir, "extra5_loc_H_hat_L=10.csv", sep="/")
  decomp_w_path   = paste(output_dir, "extra5_loc_W_hat_L=10.csv", sep="/")
  decomp_th_path  = paste(output_dir, "extra5_loc_Theta_hat_L=10.csv", sep="/")
  match_w_path    = paste(output_dir, "extra5_loc_W_retrofit.csv", sep="/")
  match_h_path    = paste(output_dir, "extra5_loc_H_retrofit.csv", sep="/")
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_cor_map, match_h_path)
  write.csv(result$w_cor_map, match_w_path)
  
  # plot
  # RetrofitPlot(dir=output_dir, file = folder, save=TRUE)
}


MatchSimulation <- function(dir, file) {
  in_file = paste(file, ".csv", sep="")
  in_path = paste(dir, in_file, sep="/")
  out_h_file = paste(file,"_out_H.csv", sep="")
  out_h_path = paste(dir, out_h_file, sep="/")
  out_w_file = paste(file,"_out_W.csv", sep="")
  out_w_path = paste(dir, out_w_file, sep="/")
  out_t_file = paste(file,"_out_T.csv", sep="")
  out_t_path = paste(dir, out_t_file, sep="/")
  
  X=read.csv(in_path)
  rownames(X)=X[,1]
  X=as.matrix(X[,-1])
  
  result = RetrofitMatch(X, iterations)
  
  write.csv(result["h"], out_h_path)
  write.csv(result["w"], out_w_path)
  write.csv(result["t"], out_t_path)
  
  print("simulation finished")
  paths <- list(in_path, out_w_path, out_h_path, out_t_path)
  names(paths) <- c("in_path", "out_w_path", "out_h_path", "out_t_path")
  return(paths)
}

RetrofitMatchSimulationLocal <- function() {
  print(paste("working directory: ", getwd()))
  setwd("~/Research/retrofit/retrofit/R")
  
  # in_file = c("N=20,M=5_loc_X", "N=10,M=3_loc_X", "extra5_loc_X")
  in_file = "N=10,M=3_loc_X"
  in_dir = "../results"
  K = 10
  
  ref_w_file = "Cerebellum_W_K=10.csv"
  ref_w_path = paste(in_dir, ref_w_file, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  
  decomp_dir = in_dir
  match_dir = in_dir
  decomp_h_file = paste(in_file,"_decomp_H_curved.csv", sep="")
  decomp_h_path = paste(decomp_dir, decomp_h_file, sep="/")
  decomp_w_file = paste(in_file,"_decomp_W_curved.csv", sep="")
  decomp_w_path = paste(decomp_dir, decomp_w_file, sep="/")
  
  match_h_file = paste(in_file,"_match_H_curved.csv", sep="")
  match_h_path = paste(match_dir, match_h_file, sep="/")
  match_w_file = paste(in_file,"_match_W_curved.csv", sep="")
  match_w_path = paste(match_dir, match_w_file, sep="/")
  
  
  w=read.csv(decomp_w_path, row.names = 1)
  h=read.csv(decomp_h_path, row.names = 1)
  
  ret = RetrofitMatchWithRef(ref_w, w, h, K)
  
  write.csv(ret$w, match_w_path)
  write.csv(ret$h, match_h_path)
}

RetrofitMatchSimulationTemp <- function() {
  print(paste("working directory: ", getwd()))
  
  in_r_dir = "~/Research/retrofit/retrofit/results/local/plot"
  in_w_file = "N=20,M=5_loc_W_hat_L=20"
  in_h_file = "N=20,M=5_loc_H_hat_L=20"
  K = 10
  
  w_path = paste(in_r_dir,"/", in_w_file, '.csv', sep="")
  h_path = paste(in_r_dir,"/", in_h_file, '.csv', sep="")
  
  w=read.csv(w_path, row.names = 1)
  h=read.csv(h_path, row.names = 1)
  
  setwd("~/Research/retrofit/retrofit/R")
  in_dir = "../results"
  ref_w_file = "Cerebellum_W_K=10.csv"
  ref_w_path = paste(in_dir, ref_w_file, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1)
  
  ret = RetrofitMatchWithRef(ref_w, w, h, K)
  
  write.csv(ret$w, paste(in_dir, "/", "local", "/", in_w_file, '_match.csv', sep=""))
  write.csv(ret$h, paste(in_dir, "/", "local", "/", in_h_file, '_match.csv', sep=""))
}

RetrofitMatchMarkerTemp <- function() {
  print(paste("working directory: ", getwd()))
  
  w_path = "~/Downloads/Results_Retrofit/iter_10000_seed2/decomposed/extra5_loc_W_decomposed.csv"
  h_path = "~/Downloads/Results_Retrofit/iter_10000_seed2/decomposed/extra5_loc_H_decomposed.csv"
  
  w=read.csv(w_path, row.names = 1)
  h=read.csv(h_path, row.names = 1)
  K=5
  
  ref_w_path = "~/Downloads/Results_Retrofit/iter_10000_seed2/datasets/ref_marker__Cerebellum_pseudomarkers.csv"
  ref_h_path = "~/Downloads/Results_Retrofit/iter_10000_seed2/plot/extra5_H.csv"
  ref_w=read.csv(ref_w_path, row.names = 1)
  ref_marker = list()
  for(r in 1:nrow(ref_w)){
    gene = ref_w[[1]][r]
    cell_type = ref_w[[2]][r]
    if(is.null(ref_marker[[cell_type]])){
      ref_marker[[cell_type]] = c()
    }
    ref_marker[[cell_type]] = c(ref_marker[[cell_type]], gene)
  }
  
  ret = RetrofitMapByMarkers(ref_marker, K, w, h)
  H_mod8=ret$h
  
  H=read.csv(ref_h_path)[,-1]
  H_mod8=H_mod8[1:5,]
  cor_H8=sort(diag(cor(H,H_mod8)),decreasing=T,na.last=T)
  value = round(AUC(x=seq(0,1,length.out = 1000),y=cor_H8),digits=3)
  print(value)
}


RetrofitDecomposeSimulation <- function(dir, file, iterations=2) {
  setwd("~/Research/retrofit/retrofit/R")
  
  in_file = paste(file, ".csv", sep="")
  in_path = paste(dir, in_file, sep="/")
  decomp_h_file = paste(file,"_decomp_H.csv", sep="")
  decomp_h_path = paste(dir, decomp_h_file, sep="/")
  decomp_w_file = paste(file,"_decomp_W.csv", sep="")
  decomp_w_path = paste(dir, decomp_w_file, sep="/")
  decomp_t_file = paste(file,"_decomp_T.csv", sep="")
  decomp_t_path = paste(dir, decomp_t_file, sep="/")
  
  x=read.csv(in_path)
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])
  
  result = RetrofitDecompose(x, iterations=iterations)
  
  write.csv(result["h"], decomp_h_path)
  write.csv(result["w"], decomp_w_path)
  write.csv(result["t"], decomp_t_path)
  
  print("simulation finished")
  paths <- list(in_path, decomp_w_path, decomp_h_path, decomp_t_path)
  names(paths) <- c("decomp_in", "decomp_w", "decomp_h", "decomp_t")
  return(paths)
}

RetrofitDecomposeSimulationLocal <- function() {
  setwd("~/Research/retrofit/retrofit/R")
  # file_source = "N=20,M=5_loc_X"
  file_source = "N=10,M=3_loc_X"
  # file_source = "extra5_loc_X"
  K = 20
  in_dir = "../results"
  decomp_dir = in_dir
  iterations=4000
  in_file = paste(file_source,".csv", sep="")
  in_path = paste(in_dir, in_file, sep="/")
  decomp_h_file = paste(file_source,"_decomp_H.csv", sep="")
  decomp_h_path = paste(decomp_dir, decomp_h_file, sep="/")
  decomp_w_file = paste(file_source,"_decomp_W.csv", sep="")
  decomp_w_path = paste(decomp_dir, decomp_w_file, sep="/")
  decomp_t_file = paste(file_source,"_decomp_T.csv", sep="")
  decomp_t_path = paste(decomp_dir, decomp_t_file, sep="/")
  
  print(in_path)
  x=read.csv(in_path)
  rownames(x)=x[,1]
  x=as.matrix(x[,-1])
  
  result = RetrofitDecompose(x, K=K, iterations=iterations, plot_convergence=TRUE)
  
  write.csv(result["h"], decomp_h_path)
  write.csv(result["w"], decomp_w_path)
  write.csv(result["t"], decomp_t_path)
  
  print(paste("simulation finished saved at: ", decomp_w_path, decomp_h_path, decomp_t_path, sep = " "))
}


RetrofitHippoSimulation <- function(iterations, seed, lambda) {
  setwd("~/storage/home/a/akp6031/retrofit/R")
  
  file_source = "HippoI"
  ref_file_source = "HippoI_W"
  # ref_marker_file_source = "Cerebellum_pseudomarkers"
  L = 16
  K = 8
  in_dir = "~/Downloads/Results_Hippo/"
  out_dir = "~/Downloads/Results_Hippo/"
  iterations=iterations
  lambda=lambda
  seed=seed
  folder = paste("iter_", iterations, "_seed", seed, "_lambda", lambda, sep="")
  
  out_dir = paste(out_dir, folder, sep="/")
  out_decomp_dir = paste(out_dir, "decomposed", sep="/")
  out_map_dir = paste(out_dir, "mapped", sep="/")
  out_dataset_dir = paste(out_dir, "datasets", sep="/")
  out_supplementary_dir = paste(out_dir, "supplementary", sep="/")
  dir.create(out_dir)
  dir.create(out_decomp_dir)
  dir.create(out_map_dir)
  dir.create(out_dataset_dir)
  dir.create(out_supplementary_dir)
  
  ref_w_file = paste(ref_file_source, ".csv", sep="")
  ref_w_path = paste(in_dir, ref_w_file, sep="/")
  ref_w=read.csv(ref_w_path, row.names = 1, check.names = FALSE)
  
  # ref_marker_path = paste(in_dir, paste(ref_marker_file_source, ".csv", sep=""), sep="/")
  # ref_marker_d=read.csv(ref_marker_path, row.names = 1)
  # ref_marker = list()
  # for(r in 1:nrow(ref_marker_d)){
  #   gene = ref_marker_d[[1]][r]
  #   cell_type = ref_marker_d[[2]][r]
  #   if(is.null(ref_marker[[cell_type]])){
  #     ref_marker[[cell_type]] = c()
  #   }
  #   ref_marker[[cell_type]] = c(ref_marker[[cell_type]], gene)
  # }
  
  in_file = paste(file_source, "_X.csv", sep="")
  in_path = paste(in_dir, in_file, sep="/")
  x=read.csv(in_path, row.names = 1, check.names = FALSE)
  
  result = RetrofitMain(x, ref_cor=ref_w, ref_marker=NULL, iterations=iterations, L=L, K=K, lambda=lambda, seed=seed)
  
  dataset_x_path = paste(out_dataset_dir, in_file, sep="/")
  dataset_cor_ref_path = paste(out_dataset_dir, ref_w_file, sep="/")
  decomp_h_path = paste(out_decomp_dir, paste(file_source,"__decomp_h.csv", sep=""), sep="/")
  decomp_w_path = paste(out_decomp_dir, paste(file_source,"__decomp_w.csv", sep=""), sep="/")
  decomp_th_path = paste(out_decomp_dir, paste(file_source,"__decomp_th.csv", sep=""), sep="/")
  match_w_path = paste(out_map_dir, paste(file_source,"__map_cor_w.csv", sep=""), sep="/")
  match_h_path = paste(out_map_dir, paste(file_source,"__map_cor_h.csv", sep=""), sep="/")
  # match_marker_w_path = paste(out_dir, paste(file_source,"__map_marker_w.csv", sep=""), sep="/")
  # match_marker_h_path = paste(out_dir, paste(file_source,"__map_marker_h.csv", sep=""), sep="/")
  cor_w_path = paste(out_supplementary_dir, paste(file_source,"__cor_w.csv", sep=""), sep="/")
  
  write.csv(x, dataset_x_path)
  write.csv(ref_w, dataset_cor_ref_path)
  write.csv(result$h, decomp_h_path)
  write.csv(result$w, decomp_w_path)
  write.csv(result$th, decomp_th_path)
  write.csv(result$h_cor_map, match_h_path)
  write.csv(result$w_cor_map, match_w_path)
  # write.csv(result$h_marker_map, match_marker_h_path)
  # write.csv(result$w_marker_map, match_marker_w_path)
  write.csv(cor(ref_w, result$w_cor_map), cor_w_path)
  
  paths <- list(decomp_w_path, decomp_h_path, decomp_th_path, match_w_path, match_h_path)
  print(paths)
  return(paths)
}

RetrofitSecondaryLambdaSimulation <- function(iterations, seed, lambda, L=10) {
  workdir = getwd()
  setwd(paste("/storage/home/a/akp6031/Results_Secondary/Results_Secondary_L",L,"/",sep=""))
  
  folder = paste('iter_', iterations, '_lamb', lambda, "_seed", seed[[1]], "-",seed[[2]], "-", seed[[3]], sep="")
  dir.create(folder)
  dir.create(paste(folder, "supplementary", sep="/"))
  dir.create(paste(folder, "datasets", sep="/"))
  dir.create(paste(folder, "decomposed", sep="/"))
  dir.create(paste(folder, "mapped", sep="/"))
  
  # common
  output_dir = folder
  ref_cor_file = "Cerebellum_W_K=10.csv"
  ref_marker_file = "Cerebellum_pseudomarkers.csv"
  ref_w=read.csv(ref_cor_file, row.names = 1, check.names = FALSE)
  ref_marker_d=read.csv(ref_marker_file, row.names = 1, check.names = FALSE)
  ref_marker = list()
  for(r in 1:nrow(ref_marker_d)){
    gene = ref_marker_d[[1]][r]
    cell_type = ref_marker_d[[2]][r]
    if(is.null(ref_marker[[cell_type]])){
      ref_marker[[cell_type]] = c()
    }
    ref_marker[[cell_type]] = c(ref_marker[[cell_type]], gene)
  }
  
  # scenario 
  file_source = "N=20,M=5_loc_"
  L = L+10
  K = 10
  
  # run main
  in_path = paste(file_source,"X.csv", sep="")
  x=read.csv(in_path, row.names = 1, check.names = FALSE)
  result = RetrofitMain(x, ref_cor=ref_w, ref_marker=ref_marker, iterations=iterations, L=L, K=K, seed=seed[[1]], lambda=lambda)
  
  # save files
  write.csv(x, paste(output_dir, "/datasets/",file_source,"X.csv", sep=""))
  write.csv(ref_w, paste(output_dir, "/datasets/","Cerebellum_W_K=10.csv", sep="")) 
  write.csv(ref_marker_d, paste(output_dir, "/datasets/","Cerebellum_pseudomarkers.csv", sep=""))
  write.csv(result$w, paste(output_dir, "/decomposed/",file_source,"W_decomposed.csv", sep=""))
  write.csv(result$h, paste(output_dir, "/decomposed/",file_source,"H_decomposed.csv", sep=""))
  write.csv(result$th, paste(output_dir, "/decomposed/",file_source,"TH_decomposed.csv", sep=""))
  write.csv(result$w_cor_map, paste(output_dir, "/mapped/",file_source,"W_mapped_cor.csv", sep=""))
  write.csv(result$h_cor_map, paste(output_dir, "/mapped/",file_source,"H_mapped_cor.csv", sep=""))
  write.csv(result$w_marker_map, paste(output_dir, "/mapped/",file_source,"W_mapped_marker.csv", sep=""))
  write.csv(result$h_marker_map, paste(output_dir, "/mapped/",file_source,"H_mapped_marker.csv", sep=""))
  write.csv(result$w_cor_map_correlation, paste(output_dir, "/supplementary/", file_source,"correlation_performance.csv", sep=""))
  
  # scenario 
  file_source = "N=10,M=3_loc_"
  L = L+10
  K = 10
  
  # run main
  in_path = paste(file_source,"X.csv", sep="")
  x=read.csv(in_path, row.names = 1, check.names = FALSE)
  result = RetrofitMain(x, ref_cor=ref_w, ref_marker=ref_marker, iterations=iterations, L=L, K=K, seed=seed[[2]], lambda=lambda)
  
  # save files
  write.csv(x, paste(output_dir, "/datasets/",file_source,"X.csv", sep=""))
  write.csv(ref_w, paste(output_dir, "/datasets/","Cerebellum_W_K=10.csv", sep="")) 
  write.csv(ref_marker_d, paste(output_dir, "/datasets/","Cerebellum_pseudomarkers.csv", sep=""))
  write.csv(result$w, paste(output_dir, "/decomposed/",file_source,"W_decomposed.csv", sep=""))
  write.csv(result$h, paste(output_dir, "/decomposed/",file_source,"H_decomposed.csv", sep=""))
  write.csv(result$th, paste(output_dir, "/decomposed/",file_source,"TH_decomposed.csv", sep=""))
  write.csv(result$w_cor_map, paste(output_dir, "/mapped/",file_source,"W_mapped_cor.csv", sep=""))
  write.csv(result$h_cor_map, paste(output_dir, "/mapped/",file_source,"H_mapped_cor.csv", sep=""))
  write.csv(result$w_marker_map, paste(output_dir, "/mapped/",file_source,"W_mapped_marker.csv", sep=""))
  write.csv(result$h_marker_map, paste(output_dir, "/mapped/",file_source,"H_mapped_marker.csv", sep=""))
  write.csv(result$w_cor_map_correlation, paste(output_dir, "/supplementary/", file_source,"correlation_performance.csv", sep=""))
  
  # scenario 
  file_source = "extra5_loc_"
  # L = max(L-5, 10)
  # K = 5
  L = L
  K = 10
  
  # run main
  in_path = paste(file_source,"X.csv", sep="")
  x=read.csv(in_path, row.names = 1, check.names = FALSE)
  result = RetrofitMain(x, ref_cor=ref_w, ref_marker=ref_marker, iterations=iterations, L=L, K=K, seed=seed[[3]], lambda=lambda)
  
  # save files
  write.csv(x, paste(output_dir, "/datasets/",file_source,"X.csv", sep=""))
  write.csv(ref_w, paste(output_dir, "/datasets/","Cerebellum_W_K=10.csv", sep="")) 
  write.csv(ref_marker_d, paste(output_dir, "/datasets/","Cerebellum_pseudomarkers.csv", sep=""))
  write.csv(result$w, paste(output_dir, "/decomposed/",file_source,"W_decomposed.csv", sep=""))
  write.csv(result$h, paste(output_dir, "/decomposed/",file_source,"H_decomposed.csv", sep=""))
  write.csv(result$th, paste(output_dir, "/decomposed/",file_source,"TH_decomposed.csv", sep=""))
  write.csv(result$w_cor_map, paste(output_dir, "/mapped/",file_source,"W_mapped_cor.csv", sep=""))
  write.csv(result$h_cor_map, paste(output_dir, "/mapped/",file_source,"H_mapped_cor.csv", sep=""))
  write.csv(result$w_marker_map, paste(output_dir, "/mapped/",file_source,"W_mapped_marker.csv", sep=""))
  write.csv(result$h_marker_map, paste(output_dir, "/mapped/",file_source,"H_mapped_marker.csv", sep=""))
  write.csv(result$w_cor_map_correlation, paste(output_dir, "/supplementary/", file_source,"correlation_performance.csv", sep=""))
  
  setwd(workdir)
}