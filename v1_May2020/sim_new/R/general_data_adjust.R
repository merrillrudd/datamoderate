#' OM setup for ss3sim: make sure years range from 1-100, adjust age bins, and adjust length bins from 10-72 cm in 2 cm increments
#'
#' @param path path to life history folders
#' @param lh data frame of life history info with columns "ShortName" (3-letter abbreviation for life history) and "Amax" (maximum age)
#' @param sig sigmaR scenario
#' @return list of dat file contents for each scenario
#' @author Merrill Rudd
general_data_adjust <- function(path, lh, sig){
  dat_om <- lapply(1:nrow(lh), function(x){
    lh <- lh[x,"ShortName"]
    lhpath <- file.path(path, lh, sig, "om", paste0(lh,"OM.dat"))
    dat <- r4ss::SS_readdat(lhpath, verbose = FALSE)
    
    dat$styr <- 1
    dat$endyr <- 100
    dat$catch <- dat$catch %>% filter(year %in% dat$styr:dat$endyr)
    
    Amax <- as.numeric(lh_info %>% filter(ShortName == lh) %>% select(Amax))
    dat$Nages <- Amax
    dat <- change_data(dat_list = dat, outfile = NULL, types = "age", age_bins = 1:Amax, fleets = 1, years = 100)
    dat$ageerror <- matrix(rep(c(-1,0.001),length(0:dat$Nages)), nrow = 2)
    colnames(dat$ageerror) <- paste0("age", 0:dat$Nages)
    dat$ageerror <- as.data.frame(dat$ageerror)
    
    ## adjust lengths
    dat <- change_data(dat_list = dat, outfile = NULL, years = 1:dat$endyr, fleets = 1, types = "len", 
                       len_bins = seq(10,72,by=2),
                       pop_binwidth = 1,
                       pop_minimum_size = 10,
                       pop_maximum_size = 72)
    dat$N_lbinspop <- length(seq(10,72,by=1))
    
    # dat <- change_data(dat_list = dat, outfile = NULL, years = 1, fleets = 2, types = "index")
    
    write <- r4ss::SS_writedat(dat, lhpath, overwrite = TRUE, verbose = FALSE)
    return(dat)
  })
}