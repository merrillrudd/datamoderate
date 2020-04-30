#' read either simulation or results files
#' 
#' @param mod_path path to model runs
#' @param df data frame of scenario options with columns LifeHistory, SigmaR, L when reading estimation model
#' @param itervec iterations to read from
#' @param read_truth FALSE to read estimation results, TRUE to read simulation results
#' @param Rscen default NULL, must specify recruitment estimation folder names when read_truth = FALSE, otherwise will look to om folder
#' @param rewrite TRUE will rewrite file if already written
#' @param ncores default = 1
#' @return data frame of results
#' @author Merrill Rudd
get_results <- function(mod_path, df, itervec, read_truth, Rscen = NULL, res_name = NULL, ncores = 1){

  out_dir <- file.path(mod_path, "results")
  dir.create(out_dir, showWarnings = FALSE)
  
  registerDoParallel(ncores)
  getDoParWorkers()
  res <- lapply(1:nrow(df), function(x){
  # res <- for(x in 1:nrow(df)){

    
    lh <- df[x,"LifeHistory"]
    f <- df[x,"F"]
    sig <- df[x,"SigmaR"]
    
    path <- file.path(mod_path, lh, f, sig)
    
      # byIter <- for(y in 1:length(itervec)){
      byIter <-  foreach(y = 1:length(itervec), .packages = c("tidyverse","r4ss")) %dopar%{

        if(read_truth == FALSE){
          lyrs <- df[x,"L"]
          lsamp <- df[x,"Samp"]
          if(df[x,"Rich"] == FALSE) res_dir <- file.path(path, itervec[y], lsamp, lyrs)
          if(df[x,"Rich"] == TRUE) res_dir <- file.path(path, itervec[y], paste0(lsamp,"_Rich"), lyrs)
          if(all(is.null(Rscen)) | any(Rscen %in% c("R_unadjusted","R_biasadjusted","R_noest") == FALSE)) stop("must specify recruitment scenarios to read")
        }
        if(read_truth == TRUE){
          Rscen <- "om"
          res_dir <- file.path(path, itervec[y])
        }
        byR <- lapply(1:length(Rscen), function(z){
        # for(z in 1:length(Rscen)){
          
          if(file.exists(file.path(res_dir, Rscen[z], "Report.sso"))){
            r1 <- r4ss::SS_output(file.path(res_dir, Rscen[z]), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
            
            d1 <- get_results_derived(r1)
            # s1 <- get_results_scalar(r1)
            
            if(any(grepl("StdDev", colnames(d1)))){
              d2 <- d1 %>% 
                rename(year = Yr,
                       SSB = Value.SSB,
                       SSB_sd = StdDev.SSB,
                       Recruit = Value.Recr,
                       Recruit_sd = StdDev.Recr,
                       SPR = Value.SPRratio, 
                       SPR_sd = StdDev.SPRratio,
                       F = Value.F,
                       F_sd = StdDev.F,
                       Depletion = Value.Bratio,
                       Depletion_sd = StdDev.Bratio,
                       OFL = Value.OFLCatch,
                       OFL_sd = StdDev.OFLCatch) %>%
                select(year, SSB, SSB_sd, Recruit, Recruit_sd, SPR, SPR_sd, F, F_sd, Depletion, Depletion_sd, OFL, OFL_sd)
            } else {
              d2 <- d1 %>% 
                rename(year = Yr,
                       SSB = Value.SSB,
                       Recruit = Value.Recr,
                       SPR = Value.SPRratio, 
                       F = Value.F,
                       Depletion = Value.Bratio,
                       OFL = Value.OFLCatch) %>%
                select(year, SSB, Recruit, SPR, F, Depletion, OFL) %>%
                mutate(SSB_sd = NA, Recruit_sd = NA, SPR_sd = NA, F_sd = NA, Depletion_sd = NA, OFL_sd = NA)
            }
            d2 <- d2 %>%
              mutate(SPR = 1-SPR) 
            d2$Depletion[which(d2$year == 1)] <- 1
            d2$Depletion_sd[which(d2$year == 1)] <- 0
            d2$year <- as.numeric(d2$year)
            # d2 <- d2 %>% filter(year <= 100) 
            d2 <- d2 %>% 
              mutate(lifehistory = lh) %>%
              mutate(iteration = itervec[y]) %>%
              mutate(Rest = Rscen[z]) %>%
              mutate(SigmaR = sig) %>%
              mutate(Fscen = f)

            if(read_truth == FALSE){
              d2 <- d2 %>%
              mutate(Lyears = lyrs) %>% 
              mutate(Nsamp = lsamp) %>%
              mutate(Rich = df[x,"Rich"])
            }
              
            s1 <- data.frame('max_grad' = r1$maximum_gradient_component,
                             'LN_R0' = r1$parameters[which(r1$parameters$Label == "SR_LN(R0)"),"Value"],
                             "SSB0" = r1$SBzero)
            s2 <- s1 %>% 
              # select(SSB_MSY, TotYield_MSY, F_MSY, max_grad, SR_LN_R0, Size_DblN_peak_Fishery_1, Size_DblN_ascend_se_Fishery_1) %>%
              # rename(SSBmsy = SSB_MSY,
              #        MSY = TotYield_MSY,
              #        Fmsy = F_MSY,
              #        LN_R0 = SR_LN_R0,
              #        Selex_peak = Size_DblN_peak_Fishery_1,
              #        Selex_shape = Size_DblN_ascend_se_Fishery_1) %>%
              # mutate(SSB0 = r1$SBzero) %>%
              mutate(lifehistory = lh) %>%
              mutate(iteration = itervec[y]) %>%
              mutate(Rest = Rscen[z]) %>%
              mutate(SigmaR = sig) %>%
              mutate(Fscen = f)

            if(read_truth == FALSE){
              s2 <- s2 %>%
              mutate(Lyears = lyrs) %>%
              mutate(Nsamp = lsamp) %>%
              mutate(Rich = df[x,"Rich"])
            }
              
            
            out <- full_join(d2, s2)
          } else {
            out <- NULL
          }
          return(out)
        })
        out2 <- do.call(rbind, byR)
        return(out2)
      }
      byIter <- do.call(rbind, byIter)
    return(byIter)
  })
  out <- do.call(rbind, res)
  if(is.null(res_name) == FALSE) write.csv(out, file.path(out_dir, paste0(res_name,".csv")), row.names = FALSE)
  return(out)
}
