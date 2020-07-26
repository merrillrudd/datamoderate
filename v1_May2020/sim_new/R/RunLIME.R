RunLIME <- function(modpath, itervec, lh, lh_name, LFdist = 1, pool=TRUE, rewrite=TRUE, run=TRUE, data_list=NULL, C_type=0, data_avail="LC", eval_type="SPR", est_selex_f=TRUE, vals_selex_ft=-1, est_rdev_t=TRUE){


		byIter <- lapply(1:length(itervec), function(ii){

			if(rewrite==TRUE & file.exists(file.path(modpath, itervec[ii],"LIME_output.rds"))) unlink(file.path(modpath, itervec[ii], "LIME_output.rds"), TRUE)

			# if(rewrite==FALSE & file.exists(file.path(modpath, itervec[ii],"LIME_output.rds"))) return(NULL)
			if(all(is.null(data_list))){
				true <- readRDS(file.path(modpath, itervec[ii], "True.rds"))
				LFinp <- true$LF[,,1]
				if(C_type %in% c(0,1)) Cinp <- matrix(true$Cn_ft[1,],nrow=1,ncol=ncol(true$Cn_ft))
				if(C_type == 2) Cinp <- matrix(true$Cw_ft[1,],nrow=1,ncol=ncol(true$Cw_ft))
					rownames(Cinp) <- 1
					colnames(Cinp) <- true$years
				Iinp <- matrix(true$I_ft[1,], nrow = 1, ncol = ncol(true$I_ft))
				data_list_toUse <- list("years"=true$years, "LF"=LFinp, "C_ft"=Cinp, I_ft = Iinp)
			}
			if(all(is.null(data_list))==FALSE) data_list_toUse <- data_list
			inputs <- create_inputs(lh=lh, input_data=data_list_toUse)

			if(file.exists(file.path(modpath, itervec[ii], "LF.rds"))==FALSE){
				if(inputs$nfleets==1) LFobs <- inputs$LF
				if(inputs$nfleets>1){
					if(pool==TRUE){
						LFobs <- array(0,dim=c(length(yrs),length(bins),1))
						for(bb in 1:length(bins)){
							for(yy in 1:length(yrs)){
								LFobs[yy,bb,1] <- sum(inputs$LF[yy,bb,])
							}
						}
						rownames(LFobs) <- yrs
						colnames(LFobs) <- bins					
					}
					if(pool==FALSE) LFobs <- inputs$LF
				}
			}
			if(file.exists(file.path(modpath, itervec[ii], "LF.rds"))){
				LFobs <- readRDS(file.path(modpath, itervec[ii], "LF.rds"))
			}
			inputs$LF <- LFobs

			# if(eval_type=="MSY") derive_quants <- TRUE
			# if(eval_type=="SPR") derive_quants <- FALSE

			if(run == TRUE){
				if(rewrite == TRUE | file.exists(file.path(modpath, itervec[ii], "LIME_output.rds"))==FALSE){
				out <- run_LIME(modpath=NULL,
							input=inputs,
							data_avail=data_avail,
							newtonsteps=3,
							C_type=C_type,
							LFdist=LFdist,
							rewrite=rewrite,
							derive_quants=TRUE,
							est_selex_f = est_selex_f,
							vals_selex_ft = vals_selex_ft, 
							est_rdev_t = est_rdev_t)

				## check convergence
				isNA <- all(is.null(out$df))
				if(isNA==TRUE){
					out <- run_LIME(modpath=NULL, input=inputs, data_avail=data_avail, C_type=C_type, LFdist=1, newtonsteps=FALSE, derive_quants=TRUE,
							est_selex_f = est_selex_f,
							vals_selex_ft = vals_selex_ft,
							est_rdev_t = est_rdev_t)
				}
				isNA <- all(is.null(out$df))
				if(isNA==FALSE){
					gradient <- out$opt$max_gradient <= 0.001
					pdHess <- out$Sdreport$pdHess
					if(pdHess==FALSE){
						inputs$theta <- 50
						out <- run_LIME(modpath=NULL, input=inputs, data_avail=data_avail, newtonsteps=3, C_type=C_type, LFdist=1, fix_more="log_theta", derive_quants=TRUE,
							est_selex_f = est_selex_f,
							vals_selex_ft = vals_selex_ft,
							est_rdev_t = est_rdev_t)
					}
					isNA <- all(is.null(out$df))
					if(isNA==TRUE){
						out <- run_LIME(modpath=NULL, input=inputs, data_avail=data_avail, newtonsteps=FALSE, C_type=C_type, LFdist=1,fix_more="log_theta", derive_quants=TRUE,
							est_selex_f = est_selex_f,
							vals_selex_ft = vals_selex_ft,
							est_rdev_t = est_rdev_t)
					}
					isNA <- all(is.null(out$df))
					if(isNA==FALSE){
						gradient <- out$opt$max_gradient <= 0.001
						pdHess <- out$Sdreport$pdHess
					}					
				}
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
			 }
			}

			ests_all <- NULL
			if(data_avail == "Catch_LC" & file.exists(file.path(modpath, itervec[ii], "LIME_output.rds"))){
				output <- readRDS(file.path(modpath, itervec[ii], "LIME_output.rds"))
				sdrep <- summary(output$Sdreport)
				if(eval_type=="SPR"){
	        		est <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "Year"=data_list_toUse$years, "SPR"=output$Report$SPR_t, "FishingMortality"=output$Report$F_ft[1,], "Recruitment"=output$Report$R_t, "RelativeSpawningBiomass"=output$Report$D_t) #"F30"=output$Derived$F30, 
	        		est_long <- est %>% gather(Parameter, Estimate, SPR:RelativeSpawningBiomass)
	        		sd <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "Year"=data_list_toUse$years, "SPR"=sdrep[which(rownames(sdrep)=="SPR_t"),2], "FishingMortality"=sdrep[which(rownames(sdrep)=="lF_ft"),2], "Recruitment"=sdrep[which(rownames(sdrep)=="lR_t"),2], "RelativeSpawningBiomass"=sdrep[which(rownames(sdrep)=="lD_t"),2]) #"F30" = NA, 
	        		sd_long <- sd %>% gather(Parameter, SD, SPR:RelativeSpawningBiomass)					
	        		ests_all <- merge(est_long, sd_long)
				}
				if(eval_type=="MSY"){
	        		ests_all <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "MSY_estimate"=output$Derived$msy, "Fmsy_estimate"=output$Derived$Fmsy, "Bmsy_estimate"=output$Derived$Bmsy)
				}
			} 
			if(data_avail == "LC" & file.exists(file.path(modpath, itervec[ii], "LIME_output.rds"))){
				output <- readRDS(file.path(modpath, itervec[ii], "LIME_output.rds"))
				sdrep <- summary(output$Sdreport)
				if(eval_type=="SPR"){
	        		est <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "Year"=data_list_toUse$years, "SPR"=output$Report$SPR_t, "FishingMortality"=output$Report$F_ft[1,], "Recruitment"=output$Report$R_t, "RelativeSpawningBiomass"=output$Report$D_t) #"F30"=output$Derived$F30, 
	        		est_long <- est %>% gather(Parameter, Estimate, SPR:RelativeSpawningBiomass)
	        		sd <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "Year"=data_list_toUse$years, "SPR"=sdrep[which(rownames(sdrep)=="SPR_t"),2],  "FishingMortality"=sdrep[which(rownames(sdrep)=="lF_ft"),2], "Recruitment"=sdrep[which(rownames(sdrep)=="lR_t"),2], "RelativeSpawningBiomass"=sdrep[which(rownames(sdrep)=="lD_t"),2]) #"F30" = NA,
	        		sd_long <- sd %>% gather(Parameter, SD, SPR:RelativeSpawningBiomass)					
	        		ests_all <- merge(est_long, sd_long)
				}
				if(eval_type=="MSY"){
	        		ests_all <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "MSY_estimate"=output$Derived$msy, "Fmsy_estimate"=output$Derived$Fmsy, "Bmsy_estimate"=output$Derived$Bmsy)
				}
			}			


			if(all(is.null(ests_all))==FALSE) {
        		if(file.exists(file.path(modpath, itervec[ii], "True.rds"))){
        			if(eval_type=="SPR"){
		        		tval <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "Year"=data_list_toUse$years, "SPR"=true$SPR_t, "FishingMortality"=true$F_ft[1,], "Recruitment"=true$R_t, "RelativeSpawningBiomass"=true$D_t) #"F30"=true$F30, 
						tval_long <- tval %>% gather(Parameter, True, SPR:RelativeSpawningBiomass)        				
        			}
        			if(eval_type=="MSY"){
		        		tval_long <- data.frame("LifeHistory"=lh_name, "Iteration"=itervec[ii], "MSY_true"=true$msy)
        			}
        			all_long <- merge(ests_all, tval_long)


					if(eval_type=="SPR"){
						all_long$Year <- as.numeric(all_long$Year)
						all_long_order <- all_long[order(all_long$Year),]
					}
					if(eval_type=="MSY"){
						all_long_order <- all_long
					}
        		} else{
        			all_long_order <- ests_all
        		}

				return(all_long_order)
			} else{
				return(NULL)
			}
		})
		byIter <- do.call(rbind, byIter)
		return(byIter)
}

