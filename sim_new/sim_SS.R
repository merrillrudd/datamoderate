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


####################
## directories
####################

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

###########################################
## Simulation replicates
###########################################
itervec <- 1:100
Nyears <- 100
Fyears <- 10

recdev_mat0 <- matrix(0, nrow = Nyears+Fyears, ncol = length(itervec))

set.seed(123)
recdev_mat1 <- matrix(rnorm((Nyears+Fyears)*length(itervec), mean = 0, 0.1),ncol=length(itervec))

set.seed(123)
recdev_mat2 <- matrix(rnorm((Nyears+Fyears)*length(itervec), mean = 0, 0.4),ncol=length(itervec))

set.seed(123)
recdev_mat3 <- matrix(rnorm((Nyears+Fyears)*length(itervec), mean = 0, 0.8),ncol=length(itervec))


f1 <- c(seq(0.01, 0.05, length.out = 25), seq(0.05, 2, length.out = 30), rep(2, 5), seq(2, 0.6, length.out = 20), rep(0.6,20))
f2 <- c(seq(0.01, 2, length.out = 55), rep(2, 5), seq(2, 0.6, length.out = 40))

###########################################
## Operating model - generate true population
###########################################
########################
### deterministic
########################
## scenarios to simulate
lh_vec <- lh_info$ShortName
SigRscen <- c("Deterministic")
Fscen <- c("F1")
name <- "sim_deterministic"

sim_grid <- expand.grid("SigmaR" = SigRscen, "F"=Fscen, "LifeHistory" = lh_vec)

## simulate scenarios
itervec <- 1
rewrite = FALSE
start <- Sys.time()
sim_det <- simulate_pop(df = sim_grid, 
                         path = sim_path,
                         itervec = itervec, 
                         f_list = list("F1" = f1, "F2" = f2),
                         lh_df = lh_info, 
                         ncores = 4, 
                         rewrite = rewrite)
end_gen <- Sys.time() - start

## read simulations
sim_wide <- get_results(mod_path = sim_path,
                        df = sim_grid,
                        itervec = itervec,
                        read_truth = TRUE,
                        res_name = paste0(name, "_wide"))
sim_long <- sim_wide %>% 
  select(-c(max_grad)) %>%
  pivot_longer(-c(iteration, lifehistory, year, SigmaR, Rest, Fscen), names_to = "variable", values_to = "true") %>%
  filter(grepl("_sd", variable) == FALSE)
write.csv(sim_long, file.path(sim_path, "results", paste0(name,"_long.csv")), row.names = FALSE)


########################
### variable
########################
## scenarios to simulate
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR","HighSigmaR")
Fscen <- c("F1")
name <- "sim_variable"

sim_grid <- expand.grid("SigmaR" = SigRscen, "F"=Fscen, "LifeHistory" = lh_vec)

## simulate scenarios
itervec <- 1:100
rewrite <- FALSE
start <- Sys.time()
sim_var <- simulate_pop(df = sim_grid, 
                        path = sim_path,
                        itervec = itervec, 
                        f_list = list("F1" = f1, "F2" = f2),
                        lh_df = lh_info, 
                        ncores = 14, 
                        rewrite = rewrite)
end_gen <- Sys.time() - start

sim_wide <- get_results(mod_path = sim_path,
                        df = sim_grid,
                        itervec = itervec,
                        read_truth = TRUE,
                        res_name = paste0(name,"_wide"))
sim_long <- sim_wide %>% 
  select(-c(max_grad)) %>%
  pivot_longer(-c(iteration, lifehistory, year, SigmaR, Rest, Fscen), names_to = "variable", values_to = "true") %>%
  filter(grepl("_sd", variable) == FALSE)
write.csv(sim_long, file.path(sim_path, "results", paste0(name,"_long.csv")), row.names = FALSE)

### explore true population
sim_det <- read.csv(file.path(sim_path,"results", "sim_deterministic_long.csv"), stringsAsFactors = FALSE)
sim_var <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv"), stringsAsFactors = FALSE)

sim_long <- full_join(sim_det, sim_var)
write.csv(sim_long, file.path(sim_path, "results", "sim_all_long.csv"), row.names = FALSE)

sim_long <- read.csv(file.path(sim_path, "results", "sim_all_long.csv"))
sim_long$SigmaR <- as.character(sim_long$SigmaR)
sim_long <- sim_long %>%
  mutate(SigmaR_desc = replace(SigmaR, which(SigmaR == "LowSigmaR"), "SigmaR = 0.4"),
         SigmaR_desc = replace(SigmaR_desc, which(SigmaR_desc == "HighSigmaR"), "SigmaR = 0.8")) %>%
  mutate(Longevity = ifelse(grepl("long_", lifehistory), "Long-lived", "Short-lived")) %>%
  mutate(Growth = ifelse(grepl("_slow", lifehistory), "Slow-growing", "Fast-growing"))

check <- sim_long %>% filter(variable == "Recruit") 
p_recom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory, linetype = factor(SigmaR)), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  geom_line(data = check %>% filter(iteration == 1), aes(linetype = factor(SigmaR))) +
  facet_grid(Longevity + Growth ~ SigmaR_desc) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Recruitment") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_path, "Recruit_OM.png"), p_recom, height = 10, width = 12)

check <- sim_long %>% filter(variable == "Depletion") 
p_deplom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory, linetype = factor(SigmaR)), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  geom_line(data = check %>% filter(iteration == 1), aes(linetype = factor(SigmaR))) +
  facet_grid(Longevity + Growth ~ SigmaR_desc) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Depletion") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_path, "Depletion_OM.png"), p_deplom, height = 10, width = 12)


check <- sim_long %>% filter(variable == "F")
p_fom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory, linetype = factor(SigmaR)), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  facet_grid(Longevity + Growth ~ SigmaR_desc) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Fishing mortality rate") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_path, "FishingMortality_OM.png"), p_fom, height = 10, width = 12)


##########################################
## Estimation model
##########################################
#############################
## run perfect deterministic
#############################
lh_vec <- lh_info$ShortName
SigRscen <- c("Deterministic")
Lyears <- c("L100")
Lsamp <- c("perfect")
res_name <- "results_det"
Rich <- FALSE
Fscen <- "F1"

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

itervec <- 1
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = 4)
end_run <- Sys.time() - start

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_unadjusted","R_biasadjusted"),
                        res_name = paste0(res_name, "_wide"))

res_long_all <- res_wide %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)

write.csv(res_long, file.path(sim_path, "results", paste0(res_name, "_long.csv")), row.names = FALSE)

true <- read.csv(file.path(sim_path, "results", "sim_deterministic_long.csv")) %>% select(-Rest)
res <- read.csv(file.path(sim_path, "results", "results_det_long.csv"))
all <- left_join(res, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all_info <- all %>% 
  mutate(re = (estimate - true)/true) #%>%
  # mutate(P25 = estimate - 0.675*sd) %>%
  # mutate(P75 = estimate + 0.675*sd) %>%
  # mutate(cover = ifelse(true >= P25 & true <= P75, 1, 0))

sum_final <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(which(max_grad <= 1))/length(max_grad))#,
            #pcover = sum(cover, na.rm=TRUE)/length(which(is.na(cover) == FALSE)))# %>%
  #filter(variable == "Depletion")
write.csv(sum_final, file.path(fig_path, "RE_summary_deterministic.csv"))

check <- all_info %>% 
  filter(variable %in% c("Recruit","F","SSB","Depletion")) %>% 
  select(-c(sd,re)) %>%
  pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value")
p_det_check <- ggplot(check) +
  geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
  facet_wrap(lifehistory~variable, scales ="free_y", nrow = 4) +
  theme_bw(base_size = 14)

#############################
## run with variation
#############################
## perfect, rich 
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")
Fscen <- c("F1")
res_name <- "results_variable"
Rich <- c(TRUE)

scen_grid_rich <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)
scen_grid <- scen_grid_rich
ncores <- 14
itervec <- 1:14
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       ncores = ncores)
end_run <- Sys.time() - start


## perfect
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")#, "N100")
Fscen <- c("F1")# ,"F2") #
res_name <- "results_variable"
Rich <- c(FALSE)#,FALSE)

scen_grid_perf <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

## sample, 100 years
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR")
Lyears <- c("L100")
Lsamp <- c("N100")
Fscen <- c("F1")# ,"F2") #
res_name <- "results_variable"
Rich <- c(FALSE)

scen_grid_samplong <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

## sample,less than 100 years
lh_vec <- lh_info$ShortName[1]
SigRscen <- c("LowSigmaR")
Lyears <- c("L75","L1")
Lsamp <- c("N100")
Fscen <- c("F1")# ,"F2") #
res_name <- "results_variable"
Rich <- c(FALSE)

scen_grid_samp <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

scen_grid <- scen_grid_perf # rbind.data.frame(scen_grid_perf, scen_grid_samp)
ncores <- 14
itervec <- 1:64
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       ncores = ncores)
end_run <- Sys.time() - start

scen_grid_all <- rbind.data.frame(scen_grid_rich, scen_grid_perf, scen_grid_samplong, scen_grid_samp)
res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid_all, 
                        itervec = 1:100, 
                        read_truth = FALSE, 
                        Rscen = c("R_unadjusted","R_biasadjusted"),
                        rewrite = TRUE, 
                        res_name = paste0(res_name, "_wide"))
res_long_all <- res_wide %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Lyears, Nsamp, Rich, Fscen), names_to = "variable", values_to = "estimate")
res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, paste0(res_name, "_long.csv")), row.names = FALSE)

##########################
true <- read.csv(file.path(sim_path, "sim_variable_long.csv")) %>% select(-Rest)
res <- read.csv(file.path(sim_path, "results_variable_long.csv"))
all <- left_join(res, true, join_by = c(year, lifehistory, iteration, SigmaR, Fscen, variable))
all_info <- all %>% 
  mutate(re = (estimate - true)/true) %>%
  mutate(P25 = estimate - 0.675*sd) %>%
  mutate(P75 = estimate + 0.675*sd) %>%
  mutate(cover = ifelse(true >= P25 & true <= P75, 1, 0)) %>%
  mutate(scenario = paste0(Nsamp,"_",Lyears,"_",Rich)) %>%
  mutate(scenario = replace(scenario, grepl(TRUE, scenario), "data-rich"),
         scenario = replace(scenario, grepl("perfect", scenario), "perfect"))
all_info$scenario <- sapply(1:length(all_info$scenario), function(x) ifelse(grepl("FALSE",all_info$scenario[x]), strsplit(all_info$scenario[x],"_FALSE")[[1]][1], all_info$scenario[x]))

sum_final <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, Rich, Fscen, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(which(max_grad <= 1))/length(max_grad),
            pcover = sum(cover, na.rm=TRUE)/length(which(is.na(cover) == FALSE)))
write.csv(sum_final, file.path(fig_path, "RE_summary.csv"), row.names = FALSE)

sum_final2 <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, scenario, Fscen, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(which(max_grad <= 1))/length(max_grad),
            pcover = sum(cover, na.rm=TRUE)/length(which(is.na(cover) == FALSE)))

check_mre <- lapply(1:length(itervec), function(x){
  sub <- all_info %>% filter(year == 100) %>% filter(variable %in% c("Depletion", "F", "SSB")) %>% filter(iteration <= itervec[x])
  sum <- sub %>%
    group_by(lifehistory, SigmaR, scenario, variable, Rest) %>% 
    summarise(mre = median(re)) %>%
    mutate(iteration = itervec[x])
  return(sum)
})
check_mre <- do.call(rbind, check_mre)
# check_mre$Lyears <- factor(check_mre$Lyears, levels = Lyears)
piter <- ggplot(check_mre) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(x = iteration, y = mre, color = variable, linetype = Rest)) + 
  facet_grid(lifehistory~SigmaR+scenario) +
  coord_cartesian(ylim = c(min(check_mre$mre), quantile(check_mre$mre, 0.99))) +
  theme_bw()
ggsave(file.path(fig_path, "mre_byiter.png"), piter, width = 10, height = 8)

check_pre <- all_info %>% filter(variable %in% c("F","SPR","SSB","Depletion","SSB0")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check_pre$variable <- factor(check_pre$variable, levels = c("F","SPR","SSB","Depletion","SSB0"))
# check_pre$Lyears <- factor(check_pre$Lyears, levels = Lyears) #, "L99","L75","L1"
pre <- ggplot(check_pre) +
  geom_violin( aes(x = scenario, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black') +
  facet_wrap(lifehistory~variable, nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(-1,1)) +
  guides(fill = FALSE) +
  theme_bw()
ggsave(file.path(fig_path, "RE_LowSigmaR.png"), pre, height = 10, width = 18)

check_pret <- all_info %>% filter(variable %in% c("Depletion")) %>% filter(Rest == "R_biasadjusted") %>% filter(year %in% seq(100,1,by=-5))
# check_pret$Lyears <- factor(check_pret$Lyears, levels = Lyears) #, "L99","L75","L1"
pret <- ggplot(check_pret) +
  geom_violin( aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black') +
  facet_wrap(lifehistory~variable+scenario, ncol = 4) +
  xlab("Year") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(-1,1)) + # c(quantile(check_pre$re,0.01,na.rm=TRUE),quantile(check_pre$re,0.99,na.rm=TRUE))) +
  guides(fill = FALSE) +
  theme_bw()
ggsave(file.path(fig_path, "RE_byYear_LowSigmaR.png"), pret, height = 14, width = 14)

check_ic <- sum_final2 %>% filter(variable %in% c("F","SPR","SSB","Depletion")) 
pic <- ggplot(check_ic) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_point(aes(x = scenario, y = pcover, color = Rest), cex = 4, alpha = 0.8) +
  facet_wrap(lifehistory~variable) +
  ylim(c(0,1)) +
  theme_bw()
ggsave(file.path(fig_path, "IC_LowSigmaR.png"), pic, height = 10, width = 18)


r3 <- SS_output(file.path(sim_path, "shs","LowSigmaR","1","om"))
d3 <- get_results_derived(r3) %>% mutate(SigmaR = "LowSigmaR")%>% mutate(type = "om")

r4 <- SS_output(file.path(sim_path, "shs","LowSigmaR",1,"perfect","L100","R_biasadjusted"))
d4 <- get_results_derived(r4) %>% mutate(SigmaR = "LowSigmaR") %>% mutate(type = "catch+length")
d4 <- d4 %>% select(colnames(d4)[which(colnames(d4) %in% colnames(d3))])

r2 <- SS_output(file.path(sim_path, "shs","LowSigmaR",1,"perfect","L100_Rich","R_biasadjusted"))
d2 <- get_results_derived(r2) %>% mutate(SigmaR = "LowSigmaR") %>% mutate(type = "rich")
d2 <- d2 %>% select(colnames(d2)[which(colnames(d2) %in% colnames(d3))])

r5 <- SS_output(file.path(sim_path, "shs","Deterministic","1","om"))
d5 <- get_results_derived(r5) %>% mutate(SigmaR = "Deterministic")%>% mutate(type = "om")

r6 <- SS_output(file.path(sim_path, "shs","Deterministic",1,"perfect","L100","R_biasadjusted"))
d6 <- get_results_derived(r6) %>% mutate(SigmaR = "Deterministic") %>% mutate(type = "catch+length") 
d6 <- d6 %>% select(colnames(d6)[which(colnames(d6) %in% colnames(d5))])

d <- rbind.data.frame(d2,d3,d4,d5,d6)
dlong <- d %>% pivot_longer(-c(Yr,SigmaR,type), names_to = "variable", values_to = "value")
dlong$Yr <- as.numeric(dlong$Yr)
dlong <- dlong %>% filter(Yr <= 100)
dlong$SigmaR <- factor(dlong$SigmaR, levels = c("Deterministic","VeryLowSigmaR","LowSigmaR"))
pi1 <- ggplot(dlong %>% filter(variable %in% c("Value.SSB","Value.Recr","Value.F","Value.Bratio"))) + 
  geom_line(aes(x = Yr, y = value, color = type, linetype = type))+
  facet_wrap(variable~SigmaR, ncol = 2, scales = "free_y") +
  theme_bw()
ggsave(file.path(fig_path, "Example_1iter.png"),  pi1, height = 8, width = 8)
