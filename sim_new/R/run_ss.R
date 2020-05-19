run_ss <- function(df, path, itervec, clean = FALSE, rewrite = TRUE, run_noest = TRUE, ncores, run_hess = FALSE, samp_type = "multinomial"){

  for(x in 1:nrow(df)){
    lh <- df[x,"LifeHistory"]
    lh_path <- file.path(path, lh)
    
    f_path <- file.path(lh_path, df[x,"F"])
    
    sig <- df[x,"SigmaR"]
    sig_path <- file.path(f_path, sig)
    
    dfile_path <- file.path(sig_path, "files")
    
    registerDoParallel(ncores)
    getDoParWorkers()
    
    run_iters <- foreach(i = 1:length(itervec), .packages = c('ss3sim', 'r4ss')) %dopar%{
      ipath <- file.path(sig_path, itervec[i])
      dir.create(ipath, showWarnings = FALSE)
      
      if(df[x,"Rich"] == FALSE) samp_path <- file.path(ipath, df[x,"Samp"])
      if(df[x,"Rich"] == TRUE) samp_path <- file.path(ipath, paste0(df[x,"Samp"],"_Rich"))
      dir.create(samp_path, showWarnings = FALSE)
      
      om_path <- file.path(ipath, "om")
      
      ##################################
      ## setup data file by iteration
      SS_splitdat2(inpath = om_path, outpath = samp_path, inname = "data.ss_new", outpattern = "perfect", number = FALSE, verbose = TRUE, fillblank = TRUE, MLE = TRUE)
      dat <- SS_readdat(file.path(samp_path, "perfect.ss"), version = "3.30", verbose = FALSE)
      if(df[x,"Rich"] == TRUE){
        dat <- sample_index(dat_list = dat, fleets = 1, years = list(1:100), sds_obs = list(0.1))
        dat <- sample_agecomp(dat_list = dat, fleets = 1, Nsamp = 100, years = list(1:100))
      }
      if(df[x,"Rich"] == FALSE){
        dat <- change_data(dat_list = dat, outfile = NULL, years = 1, fleets = 1, types = "index")
        dat$CPUE[1,"se_log"] <- 1
        dat <- sample_agecomp(dat_list = dat, fleets = NULL, Nsamp = NULL, years = NULL, cpar = NULL, ESS = NULL)
      }

      new_comp <- dat$lencomp
      end <- 6
      
      dat <- change_data(dat_list = dat, outfile = NULL, years = 1:100, fleets = 1, types = "len", len_bins = seq(10,82, by = 2), pop_binwidth = 1, pop_minimum_size = 10, pop_maximum_size = 82)
      dat$lencomp[,1:ncol(new_comp)] <- new_comp
      dat$lencomp[,(ncol(new_comp)+1):ncol(dat$lencomp)] <- 0
      
      new_comp <- dat$lencomp
      temp_comp <- dat$lencomp
      
      if(df[x,"Samp"] %in% c("perfect","Ndecline") == FALSE){
        set.seed(456)
        for(y in 1:Nyears){
          year_temp <- new_comp$Yr == y
          if(sum(year_temp)!=0){
            probs <- new_comp[year_temp, -(1:end)]
            probs[which(probs == 0)] <- 0.00001
            probs <- as.numeric(probs/sum(probs))
            Nsamp <- as.numeric(strsplit(as.character(df[x,"Samp"]),"N")[[1]][2])
            if(samp_type == "multinomial") out_comp <- as.matrix(rmultinom(1, size = Nsamp, prob = probs))
            if(samp_type == "dirichlet"){
              out_comp <- matrix(MCMCpack::rdirichlet(1, alpha = Nsamp * probs), ncol = 1)
            }
            temp_comp[year_temp,(end+1):dim(new_comp)[2]] <- out_comp
            temp_comp[,"Nsamp"] <- Nsamp
          }
        }
        test_comp <- cbind(temp_comp[,1:end], temp_comp[,(end+1):dim(new_comp)[2]])
      }
      if(df[x,"Samp"] == "Ndecline"){
        set.seed(456)
        for(y in 1:Nyears){
          year_temp <- new_comp$Yr == y
          if(sum(year_temp)!=0){
            probs <- new_comp[year_temp, -(1:end)]
            Nsamp <- ifelse(y <= 87, 200, 50)
            out_comp <- as.matrix(rmultinom(1, size = Nsamp, prob = probs))
            temp_comp[year_temp,(end+1):dim(new_comp)[2]] <- out_comp
            temp_comp[year_temp,"Nsamp"] <- Nsamp
          }
        }
        test_comp <- cbind(temp_comp[,1:end], temp_comp[,(end+1):dim(new_comp)[2]])
      }
      if(df[x,"Samp"] == "perfect"){
        test_comp <- temp_comp
        test_comp[,"Nsamp"] <- 1000
      }
      dat$lencomp <- test_comp
      if(samp_type == "multinomial") SS_writedat(dat, file.path(samp_path, paste0("ss3.dat")), overwrite = TRUE)
      if(samp_type == "dirichlet") SS_writedat(dat, file.path(samp_path, paste0("ss3_dirichlet.dat")), overwrite = TRUE)
      
      ##################################
      ## setup estimation models
      lyr_path <- file.path(samp_path, df[x,"L"])
      dir.create(lyr_path, showWarnings = FALSE)
      lyears <- as.numeric(strsplit(as.character(df[x,"L"]),"L")[[1]][2])
      lyears_vec <- (Nyears - lyears + 1):Nyears
      
      ########################################
      ## first iter, save unadjusted results
      ########################################
      if(is.na(df[x,"Fix"]) & is.na(df[x,"Adjust"])) rpath1 <- file.path(lyr_path, "R_unadjusted")
      if(is.na(df[x,"Fix"]) == FALSE & is.na(df[x,"Adjust"])) rpath1 <- file.path(lyr_path, paste0("R_unadjusted_fix_", df[x,"Fix"]))
      if(is.na(df[x,"Fix"]) & is.na(df[x,"Adjust"]) == FALSE) rpath1 <- file.path(lyr_path, paste0("R_unadjusted_change_", df[x,"Adjust"]))
      
      if(file.exists(file.path(rpath1, "Report.sso")) == FALSE | rewrite == TRUE){
        unlink(rpath1, TRUE)
        dir.create(rpath1, showWarnings = FALSE)
        copy_ctl <- file.copy(from = file.path(dfile_path, paste0(lh, "EM.ctl")), to = rpath1, overwrite = TRUE)
        if(samp_type == "multinomial") copy_dat <- file.copy(from = file.path(samp_path, 'ss3.dat'), to = rpath1, overwrite = TRUE)
        if(samp_type == "dirichlet") {
          copy_dat <- file.copy(from = file.path(samp_path, 'ss3_dirichlet.dat'), to = rpath1, overwrite = TRUE)
          ignore <- file.rename(from = file.path(rpath1, "ss3_dirichlet.dat"), to = file.path(rpath1, "ss3.dat"))
        }
        copy_f <- file.copy(from = file.path(dfile_path, "forecast.ss"), to = rpath1, overwrite = TRUE)
        copy_s <- file.copy(from = file.path(dfile_path, "starter.ss"), to = rpath1, overwrite = TRUE)
        
        dat <- SS_readdat(file.path(rpath1, "ss3.dat"), verbose = FALSE)
        dat$lencomp <- dat$lencomp %>% filter(Yr %in% lyears_vec)
        SS_writedat(dat, file.path(rpath1, paste0("ss3.dat")), overwrite = TRUE)
        
        ctl <- SS_readctl(file.path(rpath1, paste0(lh,"EM.ctl")), use_datlist = TRUE, datlist = dat, verbose = FALSE)
        
        ctl$Q_setup <- ctl$Q_setup[which(rownames(ctl$Q_setup)=="1 Fishery"),]
        ctl$Q_parms <- ctl$Q_parms[which(rownames(ctl$Q_parms)== "LnQ_base_Fishery(1)"),]
        if(df[x,"Rich"]==FALSE){
          ctl$Q_parms[,"PHASE"] <- -5
        }
        ctl$Q_options <- ctl$Q_options[which(rownames(ctl$Q_options)=="Fishery"),]
        if(df[x,"Rich"] == TRUE){
          ctl$lambdas <- ctl$lambdas[2,]
          ctl$N_lambdas <- 1
        }
        
        ## check settings for other parms
        if(sig == "Deterministic") ctl$SR_parms[which(rownames(ctl$SR_parms) == "SR_sigmaR"),"INIT"] <- 0.01
        if(sig == "VeryLowSigmaR") ctl$SR_parms[which(rownames(ctl$SR_Parms) == "SR_sigmaR"),"INIT"] <- 0.1
        if(sig == "LowSigmaR") ctl$SR_parms[which(rownames(ctl$SR_parms) == "SR_sigmaR"),"INIT"] <- 0.4
        if(sig == "HighSigmaR") ctl$SR_parms[which(rownames(ctl$SR_parms) == "SR_sigmaR"),"INIT"] <- 0.8
        
        ## recruitment
        ## last year lof == age 3 (age associated with length at 5% selectivity)
        ## last year shf == age 2 (age associated with length at 5% selectivity)
        ## last year shs == age 3 (age associated with length at 5% selectivity)
        ## last year los == age 5 (age associated with length at 5% selectivity)
        rmyrs <- ifelse(lh == "short_slow", 3, 
                        ifelse(lh == "short_fast", 2,
                               ifelse(lh == "long_slow", 6,
                                      ifelse(lh == "long_fast", 4, NA))))
        ctl$do_recdev <- 1
        ctl$recdev_early_start <- -29
        ctl$recdev_phase <- 3
        ctl$MainRdevYrFirst <- max(100 - lyears + 1 - dat$Nages,1)
        ctl$MainRdevYrLast <- 100 - rmyrs
        
        ctl$size_selex_parms[1,"LO"] <- 11.5
        ctl$size_selex_parms[1,"HI"] <- 81.5
        ctl$size_selex_parms[7,"LO"] <- 11.5
        ctl$size_selex_parms[7,"HI"] <- 81.5
        
        if(grepl("Selex", df[x,"Fix"])){
          ctl$size_selex_parms[c(1,3),"PHASE"] <- c(-2,-3)
        }
        if(grepl("R0", df[x,"Fix"])){
          ctl$SR_parms[1,"PHASE"] <- -1
          ctl$recdev_phase <- 1
        }
        if(grepl("Linf", df[x,"Adjust"])){
          val <- as.numeric(strsplit(df[x,"Adjust"], "_")[[1]][length(strsplit(df[x,"Adjust"],"_")[[1]])])
          ctl$MG_parms["L_at_Amax_Fem_GP_1","INIT"] <- val
        }
        if(grepl("CV", df[x,"Adjust"])){
          val <- as.numeric(strsplit(df[x,"Adjust"], "_")[[1]][length(strsplit(df[x,"Adjust"],"_")[[1]])])
          ctl$MG_parms[grepl("CV", rownames(ctl$MG_parms)),"INIT"] <- val
        }
        write <- r4ss::SS_writectl(ctllist = ctl, outfile = file.path(rpath1, "ss3.ctl"), overwrite = TRUE, verbose = FALSE)
        file.remove(file.path(rpath1, paste0(lh,"EM.ctl")))
        
        starter <- r4ss::SS_readstarter(file.path(rpath1, 'starter.ss'), verbose = FALSE)
        starter$datfile <- "ss3.dat"
        starter$ctlfile <- "ss3.ctl"
        starter$jitter_fraction <- 0.01
        starter$last_estimation_phase <- 10
        write <- r4ss::SS_writestarter(starter, dir = rpath1, overwrite = TRUE, verbose = FALSE)
        
        #Run model
        navigate <- paste0("cd ", rpath1)
        os <- .Platform$OS.type
        ss_bin <- "ss"
        bin <- get_bin(ss_bin)
        system(paste0(navigate, ";", bin), ignore.stdout = TRUE)

        # r1 <- SS_output(rpath1)
        # d1 <- get_results_derived(r1)
        # t2 <- SS_output(om_path)
        # d2 <- get_results_derived(t2)
        # plot(d2$Value.Recr[which(d2$Yr > 0)], type = "l")
        # lines(d1$Value.Recr[which(d2$Yr > 0)], col = "blue")
        # if(itervec[i] == 1){
        #   out <- r4ss::SS_output(rpath1, verbose = FALSE)
        #   r4ss::SS_plots(dir = rpath1, replist = out)
        # }
      }
      ########################################
      ## second iter, bias adjustment
      ########################################
      if(is.na(df[x,"Fix"]) & is.na(df[x,"Adjust"])) rpath2 <- file.path(lyr_path, "R_biasadjusted")
      if(is.na(df[x,"Fix"]) == FALSE & is.na(df[x,"Adjust"])) rpath2 <- file.path(lyr_path, paste0("R_biasadjusted_fix_", df[x,"Fix"]))
      if(is.na(df[x,"Fix"]) & is.na(df[x,"Adjust"]) == FALSE) rpath2 <- file.path(lyr_path, paste0("R_biasadjusted_change_",df[x,"Adjust"]))
      
      if(file.exists(file.path(rpath2, "Report.sso")) == FALSE | rewrite == TRUE){
        unlink(rpath2, TRUE)
        dir.create(rpath2, showWarnings = FALSE)
        ## copy files 
        copy1 <- file.copy(from = file.path(rpath1, "ss3.dat"), to = rpath2, overwrite = TRUE)
        copy2 <- file.copy(from = file.path(rpath1, "ss3.ctl"), to = rpath2, overwrite = TRUE)
        copy3 <- file.copy(from = file.path(rpath1, "starter.ss"), to = rpath2, overwrite = TRUE)
        copy4 <- file.copy(from = file.path(rpath1, "forecast.ss"), to = rpath2, overwrite = TRUE)
        
        dat <- r4ss::SS_readdat(file.path(rpath2, "ss3.dat"), verbose = FALSE)
        rep_bias <- r4ss::SS_output(rpath1, covar = TRUE, printstats = FALSE, forecast = FALSE)
        new_bias <- tryCatch(r4ss::SS_fitbiasramp(rep_bias), error = function(e) NA)
        if(all(is.na(new_bias))==FALSE){
          ctl <- r4ss::SS_readctl(file.path(rpath2, "ss3.ctl"), use_datlist = TRUE, datlist = dat, verbose = FALSE)
          ctl$do_recdev <- 1
          ctl$last_early_yr_nobias_adj <- new_bias$df[1,1]
          ctl$first_yr_fullbias_adj <- new_bias$df[2,1]
          ctl$last_yr_fullbias_adj <- new_bias$df[3,1]
          ctl$first_recent_yr_nobias_adj <- new_bias$df[4,1]
          ctl$max_bias_adj <- new_bias$df[5,1]
          write <- r4ss::SS_writectl(ctllist = ctl, outfile = file.path(rpath2, "ss3.ctl"), overwrite = TRUE, verbose = FALSE)
          
          #Run model
          navigate <- paste0("cd ", rpath2)
          os <- .Platform$OS.type
          ss_bin <- "ss"
          
          bin <- get_bin(ss_bin)
          if(run_hess == TRUE) system(paste0(navigate, ";", bin), ignore.stdout = TRUE)
          if(run_hess == FALSE) system(paste0(navigate, ";", bin, " -nohess"), ignore.stdout = TRUE)

          # dev.off()
          # r3 <- SS_output(rpath2)
          # d3 <- get_results_derived(r3)
          # t2 <- SS_output(om_path)
          # d2 <- get_results_derived(t2)
          # plot(d2$Value.Recr[which(d2$Yr %in% 1:100)], type = "l")
          # lines(d1$Value.Recr[which(d2$Yr %in% 1:100)], col = "blue")
          # plot(d2$Value.Bratio[which(d2$Yr %in% 1:100)], type = "l")
          # lines(d1$Value.Bratio[which(d2$Yr %in% 1:100)], col = "blue")

          # if(itervec[i] == 1){
          #   out <- r4ss::SS_output(rpath2, verbose = FALSE)
          #   r4ss::SS_plots(dir = rpath2, replist = out)
          # }
        } else {return(NULL)}
      }
      
      ########################################
      ## third iter, no recruitment estimation
      ########################################
      if(run_noest == TRUE){
      rpath3 <- file.path(lyr_path, "R_noest")
      
      if(file.exists(file.path(rpath3, "Report.sso")) == FALSE | rewrite == TRUE){
        unlink(rpath3, TRUE)
        dir.create(rpath3, showWarnings = FALSE)
        ## copy files
        copy1 <- file.copy(from = file.path(rpath2, "ss3.dat"), to = rpath3, overwrite = TRUE)
        copy2 <- file.copy(from = file.path(rpath2, "ss3.ctl"), to = rpath3, overwrite = TRUE)
        copy3 <- file.copy(from = file.path(rpath2, "starter.ss"), to = rpath3, overwrite = TRUE)
        copy4 <- file.copy(from = file.path(rpath2, "forecast.ss"), to = rpath3, overwrite = TRUE)

        dat <- r4ss::SS_readdat(file.path(rpath3, "ss3.dat"), verbose = FALSE)
        ctl <- r4ss::SS_readctl(file.path(rpath3, "ss3.ctl"), use_datlist = TRUE, datlist = dat, verbose = FALSE)
        ctl$do_recdev <- 0
        ctl$recdev_phase <- -3
        write <- r4ss::SS_writectl(ctllist = ctl, outfile = file.path(rpath3, "ss3.ctl"), overwrite = TRUE, verbose = FALSE)

        #Run model
        navigate <- paste0("cd ", rpath3)
        os <- .Platform$OS.type
        ss_bin <- "ss"

        bin <- get_bin(ss_bin)
        if(run_hess == TRUE) system(paste0(navigate, ";", bin), ignore.stdout = TRUE)
        if(run_hess == FALSE) system(paste0(navigate, ";", bin, " -nohess"), ignore.stdout = TRUE)
        
      }
      }
    }
  }
  
}