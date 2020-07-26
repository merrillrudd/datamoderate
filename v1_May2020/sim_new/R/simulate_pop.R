simulate_pop <- function(df, path, itervec, f_list, lh_df, ncores = 1, rewrite){
  registerDoParallel(ncores)
  getDoParWorkers()
  
  for(x in 1:nrow(df)){
    lh <- df[x,"LifeHistory"]
    lh_path <- file.path(path, lh)
    
    f_path <- file.path(lh_path, df[x,"F"])
    
    sig <- df[x,"SigmaR"]
    sig_path <- file.path(f_path, sig)
    if(sig == "Deterministic") dev_mat <- recdev_mat0
    if(sig == "VeryLowSigmaR") dev_mat <- recdev_mat1
    if(sig == "LowSigmaR") dev_mat <- recdev_mat2
    if(sig == "HighSigmaR") dev_mat <- recdev_mat3
    
    
    dfile_path <- file.path(sig_path, "files")
    
    run_iters <- foreach(i = 1:length(itervec), .packages = c('ss3sim', 'r4ss')) %dopar%{
      ipath <- file.path(sig_path, itervec[i])
      dir.create(ipath, showWarnings = FALSE)
      
      om_path <- file.path(ipath, "om")
      dir.create(om_path, showWarnings = FALSE)
      
      if(file.exists(file.path(om_path, "Report.sso")) == FALSE | rewrite == TRUE){
      
      copy_ctl <- file.copy(from = file.path(dfile_path, paste0(lh, "OM.ctl")), to = om_path, overwrite = TRUE)
      copy_dat <- file.copy(from = file.path(dfile_path, paste0(lh, "OM.dat")), to = om_path, overwrite = TRUE)
      copy_f <- file.copy(from = file.path(dfile_path, "forecast.ss"), to = om_path, overwrite = TRUE)
      copy_s <- file.copy(from = file.path(dfile_path, "starter.ss"), to = om_path, overwrite = TRUE)
      
      dat1 <- SS_readdat(file.path(om_path, paste0(lh,"OM.dat")), verbose = FALSE)
      ctl1 <- SS_readctl(file.path(om_path, paste0(lh, "OM.ctl")), use_datlist = TRUE, datlist = dat1, verbose = FALSE)
      for1 <- SS_readforecast(file.path(om_path, "forecast.ss"), verbose = FALSE)
      xyears <- seq(ctl1$recdev_early_start, dat1$endyr+for1$Nforecastyrs)
      devs <- setNames(dev_mat[,i], xyears)
      change_rec_devs(recdevs = devs, ctl_file_in = file.path(om_path, paste0(lh, "OM.ctl")), ctl_file_out = file.path(om_path, paste0(lh,"OM.ctl")))

      fmsy <- as.numeric(lh_info %>% filter(ShortName == lh) %>% select(Fmsy))
      if(df[x,"F"] == "F1") fvec <- f_list[["F1"]]
      if(df[x,"F"] == "F2") fvec <- f_list[["F2"]]
      f1 <- fvec * fmsy
      change_f(years = 1:length(f1), fisheries = 1, fvals = f1, seasons = 1, ses = 0.1, ctl_file_in = file.path(om_path, paste0(lh, "OM.ctl")), ctl_file_out = file.path(om_path, paste0(lh, "OM.ctl")))
      
      #Run model
      navigate <- paste0("cd ", om_path)
      os <- .Platform$OS.type
      ss_bin <- "ss"
      
      bin <- get_bin(ss_bin)
      system(paste0(navigate, ";", bin, " -nohess"), ignore.stdout = TRUE)
    }
  }
  }
}