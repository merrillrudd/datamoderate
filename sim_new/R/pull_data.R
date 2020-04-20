pull_data <- function(path, df, itervec){
  byScen <- lapply(1:nrow(df), function(x){
    lh <- df[x,"LifeHistory"]
    lh_path <- file.path(path, lh)
    
    f <- df[x,"F"]
    f_path <- file.path(lh_path, f)
    
    sig <- df[x,"SigmaR"]
    sig_path <- file.path(f_path, sig)
    
    byIter <- lapply(1:length(itervec), function(i){
      ipath <- file.path(sig_path, itervec[i])
      samp_path <- file.path(ipath, df[x,"Samp"])
      lyr_path <- file.path(samp_path, df[x,"L"])
      rpath1 <- file.path(lyr_path, "R_unadjusted")
      dat <- SS_readdat(file.path(rpath1, "ss3.dat"), verbose = FALSE)
      out <- list()
      lencomp <- dat$lencomp[,c(1,7:ncol(dat$lencomp))] %>% rename(year = Yr) %>% filter(year %in% 1:100)
      catch <- dat$catch %>% filter(year >= 1) %>% select(year, catch, catch_se)
      out <- full_join(catch, lencomp, by = "year") %>% mutate(iteration = itervec[i])
      return(out)
    })
    byIter <- do.call(rbind, byIter) 
    out <- cbind.data.frame(byIter, df[x,])
    return(out)
  })
  byScen <- do.call(rbind, byScen)
  return(byScen)
}