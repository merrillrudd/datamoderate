rm(list = ls())

####################
## libraries
####################
# devtools::install_github("ss3sim/ss3sim", 
#                          ref = "development", 
#                          dependencies = TRUE, 
#                          build_vignettes = TRUE)

# devtools::install_github("r4ss/r4ss", ref="development")

library(ss3sim)
library(r4ss)
library(tidyverse)
library(foreach)
library(doParallel)

# lime_path <- "C:\\merrill\\LIME"
lime_path <- '~/Packages/LIME'
devtools::load_all(path = lime_path)


####################
## directories
####################

# proj_path <-"C:\\merrill\\datamoderate"
proj_path <- "~/Projects/NOAA/datamoderate"
sim_path <- file.path(proj_path, "sim_new")
dir.create(sim_path, showWarnings = FALSE) 

fig_path <- file.path(sim_path, "figures")
dir.create(fig_path, showWarnings = FALSE)

R_path <- file.path(sim_path, "R")
R_fun <- list.files(R_path)
ignore <- sapply(1:length(R_fun), function(x) source(file.path(R_path, R_fun[x])))

######################################
## read and review life history info
######################################

## read life history info
lh_info <- read.csv(file.path(proj_path, "life_history_info.csv")) %>%
  mutate(ShortName = ifelse(Longevity == "Shorter" & Growth == "Slower", "short_slow",
                            ifelse(Longevity == "Shorter" & Growth == "Faster", "short_fast",
                                   ifelse(Longevity == "Longer" & Growth == "Slower", "long_slow", "long_fast"))))

L_a <- sapply(1:nrow(lh_info), function(x) lh_info[x,"Linf"]*(1 - exp(-lh_info[x,"k"]*(0:lh_info[x,"Amax"] - (-1)))))
lh_info$L_at_Amin <- sapply(1:length(L_a), function(x) L_a[[x]][1])

lh_list <- lapply(1:nrow(lh_info), function(x){
  lh <- create_lh_list(vbk=lh_info[x,"k"], 
            linf=lh_info[x,"Linf"], 
            t0=lh_info[x,"t0"], 
            lwa=6.8e-06, 
            lwb=3.101, 
            M50=0.66 * lh_info[x,"Linf"], 
            S50=20, 
            S95=30, 
            selex_input="length", 
            maturity_input="length", 
            M=lh_info[x,"M"], 
            AgeMax=lh_info[x,"Amax"],
            h = 0.7,  
            binwidth=2, 
            SigmaR=0.4, 
            SigmaF=0.1, 
            SigmaC=0.01, 
            SigmaI=0.01, 
            CVlen=0.1, 
            rho=0, 
            nseasons=1,
            R0 = exp(10))
  return(lh)
})
names(lh_list) <- lh_info$ShortName

#######################################
## pull data
#######################################
data <- read.csv(file.path(sim_path, "results", "data_CL_sample.csv"))
true_wide <- read.csv(file.path(sim_path, "results", "sim_variable_wide.csv"), stringsAsFactors =  FALSE)
true <- true_wide
true <- true %>%
	select(colnames(true)[grepl("_sd", colnames(true)) == FALSE])


#######################################
### ---- LIME -------###

##############################
## Length only
lime_dir <- file.path(sim_path, "LIME_results")
dir.create(lime_dir, showWarnings = FALSE)

data_grid <- unique(data %>% select(SigmaR, LifeHistory, L, Samp, F, Rich))
sig_vec <- unique(data_grid$SigmaR)
lh_vec <- unique(data_grid$LifeHistory)[1]
lyears_vec <- c( "L75")
cyears_vec <- "C10"
samp_vec <- unique(data_grid$Samp)
f_vec <- unique(data_grid$F)


scen_grid <- expand.grid("LifeHistory" = lh_vec,
						"F" = f_vec,
						"SigmaR" = sig_vec,
						"Samp" = samp_vec,
						"C" = cyears_vec,
						"L" = lyears_vec)

itervec <- 1:50
rewrite = FALSE

ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)      

start <- Sys.time()
lime_res <- foreach(x=1:nrow(scen_grid), .packages=c('tidyverse')) %dopar%{
# lime_res <- lapply(1:nrow(scen_grid), function(x){

  lh <- as.character(scen_grid[x,"LifeHistory"])
  dir.create(file.path(lime_dir, lh), showWarnings = FALSE)

  f <- scen_grid[x,"F"]
  dir.create(file.path(lime_dir, lh, f), showWarnings = FALSE)
  
  sig <- scen_grid[x,"SigmaR"]
  dir.create(file.path(lime_dir, lh, f, sig), showWarnings = FALSE)

  samp_dir <- file.path(lime_dir, lh, f, sig, scen_grid[x,"Samp"])
  dir.create(samp_dir, showWarnings = FALSE)

  lh_input <- lh_list[[which(names(lh_list) == lh)]]
  lh_input$SigmaF <- 0.1
  if(sig == "LowSigmaR") lh_input$SigmaR <- 0.4
  if(sig == "HighSigmaR") lh_input$SigmaR <- 0.8
  
  byIter <- lapply(1:length(itervec), function(ii){
  # byIter <- foreach(ii=1:length(itervec), .export = c("itervec","scen_grid","rewrite"), .packages=c('tidyverse')) %dopar%{
  	

	lime_path <- "C:\\merrill\\LIME"
	devtools::load_all(path = lime_path)


  	ipath <- file.path(samp_dir, itervec[ii])
  	dir.create(ipath, showWarnings = FALSE)

  	
  	cyr_dir <- file.path(ipath, scen_grid[x,"C"])
  	dir.create(cyr_dir, showWarnings = FALSE)

 	 lyr_dir <- file.path(cyr_dir, scen_grid[x,"L"])
 	 dir.create(lyr_dir, showWarnings = FALSE)
 	 lyears <- as.numeric(str_split(as.character(scen_grid[x,"L"]),"L")[[1]][2])
 	 cyears <- as.numeric(str_split(as.character(scen_grid[x,"C"]), "C")[[1]][2])
 	 
 	 if(lyears > 0 & cyears > 0) data_type <- "Catch_LC"
 	 if(lyears > 0 & cyears == 0) data_type <- "LC"

 	 	em_path <- file.path(lyr_dir, data_type)
	  	if(rewrite == TRUE) unlink(em_path, TRUE)
 	 	dir.create(em_path, showWarnings = FALSE)

 	 	if(file.exists(file.path(em_path, "LIME_output.rds")) == FALSE){

			idata <- data %>% 
				dplyr::filter(iteration == itervec[ii]) %>%
				dplyr::filter(SigmaR == sig) %>%
				dplyr::filter(LifeHistory == lh) %>%
				dplyr::filter(L == "L75") %>%
				dplyr::filter(Samp == scen_grid[x,"Samp"]) %>%
				dplyr::filter(F == f)

			years <- idata$year
			catch <- matrix(idata$catch,nrow = 1)
			colnames(catch) <- years
			rownames(catch) <- 1
			if(cyears > 0) catch[1:(100 - cyears)] <- -1

			lf <- idata[,4:35]
			bins <- sapply(1:ncol(lf), function(x) as.numeric(strsplit(colnames(lf)[x],"l")[[1]][2]))
			lf_mat <- as.matrix(lf)
			lf_mat[which(is.na(lf_mat))] <- 0
			colnames(lf_mat) <- bins
			rownames(lf_mat) <- years

			## check number of years of length data for scenario
			nodata_vec <- 1:(100 - lyears)
			lf_mat[nodata_vec,] <- 0
			
			first_yr <- 100 - max(c(cyears, lyears)) + 1
			first_lyr <- as.numeric(which(rowSums(lf_mat) > 0)[1])
			agemax <- lh_input$AgeMax
			first_rdev <- max(first_lyr - agemax,1)

			est_rdev_t_inp <- rep(1,length(years))
			if(first_rdev > 1) est_rdev_t_inp[1:(first_rdev-1)] <- 0
			
			years1 <- first_rdev:100
			years <- seq_along(first_rdev:100)
			
			lf_inp <- lf_mat[which(rownames(lf_mat) %in% years1),]
			rownames(lf_inp) <- years
			c_inp <- matrix(catch[1,which(colnames(catch) %in% years1)], nrow = 1)
			colnames(c_inp) <- years
			rownames(c_inp) <- 1
			

			data_list_toUse <- list("years"=years, "LF"=lf_inp, "C_ft"=c_inp)
			inputs <- create_inputs(lh=lh_input, input_data=data_list_toUse)
			
			if(data_type == "LC"){
				data_avail <- "LC"
				C_type <- 0
				SigRpen <- 1
				Fpen <- 0
			}
			if(data_type == "Catch_LC"){
				data_avail <- "Catch_LC"
				C_type <- 2
				SigRpen <- 0
				Fpen <- 0
			}

			# start <- Sys.time()
			out <- run_LIME(modpath=em_path,
							input=inputs,
							data_avail=data_avail,
							newtonsteps=3,
							C_type=C_type,
							LFdist=0,
							rewrite=rewrite,
							derive_quants=TRUE, 
							est_rdev_t = est_rdev_t_inp,
							SigRpen = SigRpen,
							Fpen = Fpen,
							fix_more = "log_sigma_R")

			# end_mod <- Sys.time() - start

			# t1 <- true %>% filter(iteration == 1) %>% filter(lifehistory == as.character(lh))
			# x1 <- t1 %>% filter(variable == "F")

			if(all(is.null(out$df))){
				write("model NA", file.path(em_path,paste0("LIME_modelNA.txt")))
			}
			if(all(is.null(out$df))==FALSE){
				gradient <- out$opt$max_gradient <= 0.001
				pdHess <- out$Sdreport$pdHess
				if(gradient==FALSE){
					write("highgradient", file.path(em_path, paste0( "LIME_highgradient.txt")))
				}
				if(pdHess==FALSE){
					write("Hessian not positive definite", file.path(em_path, paste0("LIME_pdHessFALSE.txt")))
				}
				## save results if converged
				if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(em_path, "LIME_output.rds"))
			}
		}
  })
}#)
end <- Sys.time() - start
stopCluster(cl)

re <- lapply(1:nrow(scen_grid), function(x){
	byIter <- lapply(1:length(itervec), function(y){

		dir <- file.path(lime_dir, scen_grid[x,"LifeHistory"], scen_grid[x,"F"], scen_grid[x,"SigmaR"], scen_grid[x,"Samp"], itervec[y], scen_grid[x,"L"], scen_grid[x,"Data"])
		res <- readRDS(file.path(dir, "LIME_output.rds"))
		truth <- true %>% 
		  filter(lifehistory == scen_grid[x,"LifeHistory"]) %>%
		  filter(iteration == itervec[y]) %>%
		  filter(Fscen == scen_grid[x,"F"]) %>%
		  filter(SigmaR == scen_grid[x,"SigmaR"])
		df <- data.frame("year" = truth$year[1:100], "Depletion_true" = truth$Depletion[1:100], "Depletion_est" = res$Report$D_t)
		plot(df$Depletion_true)
		lines(df$Depletion_est)
		
		plot(truth$Recruit)
		lines(res$Report$R_t)
	
	})
})



plot(t1$F, ylim=c(0,1))

lines(res1$Report$F_ft[1,])
lines(res2$Report$F_ft[1,], col = "blue")

plot(t1$Depletion, ylim = c(0,1))
lines(res1$Report$D_t)
lines(res2$Report$D_t, col = "blue")

plot(t1$Recruit, ylim = c(0, max(t1$Recruit)))
lines(res1$Report$R_t)
lines(res2$Report$R_t, col = "blue")

res$Report$sigma_R