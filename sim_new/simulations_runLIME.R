rm(list = ls())

####################
## libraries
####################
devtools::install_github("ss3sim/ss3sim", 
                         ref = "development", 
                         dependencies = TRUE, 
                         build_vignettes = TRUE)

devtools::install_github("r4ss/r4ss", ref="development")

library(ss3sim)
library(r4ss)
library(tidyverse)
library(foreach)
library(doParallel)

lime_path <- "C:\\merrill\\LIME"
devtools::load_all(path = lime_path)


####################
## directories
####################

proj_path <-"C:\\merrill\\datamoderate"
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
  mutate(ShortName = ifelse(Longevity == "Shorter" & Growth == "Slower", "shs",
                            ifelse(Longevity == "Shorter" & Growth == "Faster", "shf",
                                   ifelse(Longevity == "Longer" & Growth == "Slower", "los", "lof"))))

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
true_wide <- read.csv(file.path(sim_path, "results", "sim_variable_wide.csv"))
true <- true_wide %>%
	select(colnames(true)[grepl("_sd", colnames(true)) == FALSE])


#######################################
### ---- LIME -------###
start <- Sys.time()

lime_dir <- file.path(sim_path, "LIME_results")
dir.create(lime_dir, showWarnings = FALSE)

data_grid <- unique(data %>% select(SigmaR, LifeHistory, L, Samp, F, Rich))
sig_vec <- unique(data_grid$SigmaR)
lh_vec <- unique(data_grid$LifeHistory)[1]
lyears_vec <- unique(data_grid$L)
samp_vec <- unique(data_grid$Samp)
f_vec <- unique(data_grid$F)
data_vec <- c("Catch_LC", "LC")

scen_grid <- expand.grid("LifeHistory" = lh_vec,
						"F" = f_vec,
						"SigmaR" = sig_vec,
						"Samp" = samp_vec,
						"L" = lyears_vec,
						"Data" = data_vec)
scen_grid <- scen_grid[rev(order(scen_grid$LifeHistory)),]

itervec <- 1
rewrite = FALSE

ncores <- 5
registerDoParallel(ncores)
getParWorkers()    
lime_res <- foreach(x=1:nrow(scen_grid), .packages=c('tidyverse')) %dopar%{
 
  devtools::load_all( path = lime_path)

  start <- Sys.time()

  lh <- scen_grid[x,"LifeHistory"]
  f <- scen_grid[x,"F"]
  sig <- scen_grid[x,"SigmaR"]
  modpath <- file.path(sim_path, lh, f, sig, "LIME_EM")
  dir.create(path, showWarnings = FALSE)

  lh_input <- lh_list[[lh]]
  lh_input$SigmaF <- 0.1
  if(sig == "LowSigmaR") lh_input$SigmaR <- 0.4
  if(sig == "HighSigmaR") lh_input$SigmaR <- 0.8
  
  byIter <- lapply(1:length(itervec), function(ii){
  	ipath <- file.path(modpath, itervec[ii])
  	dir.create(ipath, showWarnings = FALSE)
  	if(rewrite==TRUE & file.exists(file.path(ipath,"LIME_output.rds"))) unlink(file.path(ipath, "LIME_output.rds"), TRUE)

			idata <- data %>% 
				filter(iteration == itervec[ii]) %>%
				filter(SigmaR == sig) %>%
				filter(LifeHistory == lh) %>%
				filter(L == scen_grid[x,"L"]) %>%
				filter(Samp == scen_grid[x,"Samp"]) %>%
				filter(F == f)=

			years <- idata$year
			catch <- matrix(idata$catch,nrow = 1)
			colnames(catch) <- years
			rownames(catch) <- 1

			lf <- idata[,4:35]
			bins <- sapply(1:ncol(lf), function(x) as.numeric(strsplit(colnames(lf)[x],"l")[[1]][2]))
			lf_mat <- as.matrix(lf)
			lf_mat[which(is.na(lf_mat))] <- 0
			colnames(lf_mat) <- bins
			rownames(lf_mat) <- years

			data_list_toUse <- list("years"=years, "LF"=lf_mat, "C_ft"=catch)
			inputs <- create_inputs(lh=lh_input, input_data=data_list_toUse)

			first_yr <- as.numeric(which(rowSums(lf_mat) > 0)[1])
			agemat <- lh_input$AgeMax
			first_rdev <- max(first_yr - agemat,1)

			est_rdev_t_inp <- rep(1,length(years))
			if(first_rdev < 1) est_rdev_t_inp[1:(first_rdev-1)] <- 0

			start <- Sys.time()
			out <- run_LIME(modpath=NULL,
							input=inputs,
							data_avail="Catch_LC",
							newtonsteps=3,
							C_type=2,
							LFdist=1,
							rewrite=rewrite,
							derive_quants=TRUE, 
							est_rdev_t = est_rdev_t_inp,
							SigRpen = 0,
							Fpen = 0)
			end_mod <- Sys.time() - start


			t1 <- true %>% filter(iteration == 1) %>% filter(lifehistory == as.character(lh))
			x1 <- t1 %>% filter(variable == "F")

			if(all(is.null(out$df))){
				write("model NA", file.path(modpath, itervec[ii], paste0("LIME_modelNA.txt")))
			}
			if(all(is.null(out$df))==FALSE){
				gradient <- out$opt$max_gradient <= 0.001
				pdHess <- out$Sdreport$pdHess
				if(gradient==FALSE){
					write("highgradient", file.path(modpath, itervec[ii], paste0( "LIME_highgradient.txt")))
				}
				if(pdHess==FALSE){
					write("Hessian not positive definite", file.path(modpath, itervec[ii], paste0("LIME_pdHessFALSE.txt")))
				}
				## save results if converged
				if(gradient == TRUE & pdHess == TRUE) saveRDS(out, file.path(modpath, itervec[ii], "LIME_output.rds"))
			}

  })
  

  return(out)
}
stopCluster(cl)

lime_res <- do.call(rbind, lime_res)



true <- readRDS(file.path(case_path, scen_vec[x], 20, "True.rds"))
res <- readRDS(file.path(case_path, scen_vec[x], 20, "LIME_output.rds"))
