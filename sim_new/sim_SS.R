rm(list = ls())

####################
## libraries
####################
# devtools::install_github("ss3sim/ss3sim", 
#                          ref = "development", 
#                          dependencies = TRUE, 
#                          build_vignettes = TRUE)
# 
devtools::install_github("r4ss/r4ss", ref="master")

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
Pyears <- 30

recdev_mat0 <- matrix(0, nrow = Pyears+Nyears+Fyears, ncol = length(itervec))

set.seed(123)
recdev_mat1 <- matrix(rnorm((Pyears+Nyears+Fyears)*length(itervec), mean = 0, 0.1),ncol=length(itervec))

set.seed(123)
recdev_mat2 <- matrix(rnorm((Pyears+Nyears+Fyears)*length(itervec), mean = 0, 0.4),ncol=length(itervec))

set.seed(123)
recdev_mat3 <- matrix(rnorm((Pyears+Nyears+Fyears)*length(itervec), mean = 0, 0.8),ncol=length(itervec))


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
                        res_name = paste0(name,"_wide"),
                        ncores = 14)

saveOFL <- sim_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

sim_long <- sim_wide %>% 
  select(-c(max_grad, OFL, OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, SigmaR, Rest, Fscen), names_to = "variable", values_to = "true") %>%
  filter(grepl("_sd", variable) == FALSE)
write.csv(sim_long, file.path(sim_path, "results", paste0(name,"_long.csv")), row.names = FALSE)

### explore true population
# sim_det <- read.csv(file.path(sim_path,"results", "sim_deterministic_long.csv"), stringsAsFactors = FALSE)
sim_var <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv"), stringsAsFactors = FALSE)

sim_long <- sim_var #full_join(sim_det, sim_var)
write.csv(sim_long, file.path(sim_path, "results", "sim_all_long.csv"), row.names = FALSE)

sim_long <- read.csv(file.path(sim_path, "results", "sim_all_long.csv"))
sim_long$SigmaR <- as.character(sim_long$SigmaR)
sim_long <- sim_long %>%
  mutate(SigmaR_desc = replace(SigmaR, which(SigmaR == "LowSigmaR"), "SigmaR = 0.4"),
         SigmaR_desc = replace(SigmaR_desc, which(SigmaR_desc == "HighSigmaR"), "SigmaR = 0.8")) %>%
  mutate(Longevity = ifelse(grepl("long_", lifehistory), "Long-lived", "Short-lived")) %>%
  mutate(Growth = ifelse(grepl("_slow", lifehistory), "Slow-growing", "Fast-growing"))

check <- sim_long %>% filter(variable == "Recruit") %>% filter(year <= 100)
p_recom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  geom_line(data = check %>% filter(iteration == 1)) +
  facet_grid(Longevity + Growth ~ SigmaR_desc) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Recruitment") +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Recruit_OM.png"), p_recom, height = 10, width = 12)

check <- sim_long %>% filter(variable == "Depletion") %>% filter(year <= 100)
p_deplom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  geom_line(data = check %>% filter(iteration == 1)) +
  facet_grid(Longevity + Growth ~ SigmaR_desc) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Depletion") +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Depletion_OM.png"), p_deplom, height = 10, width = 12)


check <- sim_long %>% filter(variable == "F") %>% filter(year <= 100)
p_fom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  facet_grid(Longevity + Growth ~ SigmaR_desc) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Fishing mortality rate") +
  theme_bw(base_size = 24)
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
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

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
  mutate(re = (estimate - true)/true)

sum_converge <- unique(all_info %>% 
                         filter(year == 100) %>%
                         select(iteration, lifehistory, Rest, SigmaR, Lyears, Nsamp, max_grad, LN_R0)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp) %>%
  summarise(pconverge = length(unique(c(which(max_grad <= 1),which(LN_R0<12))))/length(itervec))

sum_final <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE)) %>%
  left_join(sum_converge)
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
########################
## perfect, rich 
########################
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR", "HighSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")
Fscen <- c("F1")
res_name <- "results_variable_rich"
Rich <- c(TRUE)

scen_grid_rich <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)
scen_grid <- scen_grid_rich

ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start


res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_unadjusted","R_biasadjusted"),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)

write.csv(res_long, file.path(sim_path, "results", paste0(res_name, "_long.csv")), row.names = FALSE)

true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
res <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv"))) %>% select(-Rich)
all <- left_join(res, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all_info <- all %>% 
  mutate(re = (estimate - true)/true)

sum_final <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(which(max_grad <= 1))/length(max_grad))
write.csv(sum_final, file.path(fig_path, "RE_summary_rich.csv"))

check <- all_info %>% 
  filter(variable %in% c("Recruit","F","SSB","Depletion")) %>% 
  select(-c(sd,re)) %>%
  pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value") %>%
  filter(iteration == 1) %>%
  filter(year <= 100)
p_rich_check <- ggplot(check) +
  geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
  facet_wrap(lifehistory~variable, scales ="free_y", nrow = 4) +
  theme_bw(base_size = 14)

check <- all_info %>% filter(variable %in% c("F","SPR","SSB","Depletion","SSB0")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("F","SPR","SSB","Depletion","SSB0"))
pre_rich <- ggplot(check) +
  geom_violin( aes(x = lifehistory, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 2) +
  facet_wrap(.~variable, nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(-1,1)) +
  guides(fill = FALSE) +
  theme_bw()

check_mre <- lapply(1:length(itervec), function(x){
  sub <- all_info %>% filter(year == 100) %>% filter(variable %in% c("Depletion", "F", "SSB")) %>% filter(iteration <= itervec[x])
  sum <- sub %>%
    group_by(lifehistory, SigmaR, variable, Rest) %>% 
    summarise(mre = median(re)) %>%
    mutate(iteration = itervec[x])
  return(sum)
})
check_mre <- do.call(rbind, check_mre)
piter <- ggplot(check_mre) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(x = iteration, y = mre, color = variable, linetype = Rest)) + 
  facet_grid(lifehistory~SigmaR) +
  coord_cartesian(ylim = c(min(check_mre$mre), quantile(check_mre$mre, 0.99))) +
  theme_bw()

########################
## perfect, sample 
########################
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR", "HighSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")
Fscen <- c("F1")
res_name <- "results_variable_perfect"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

# for(x in 1:nrow(scen_grid)){
#   for(y in 1:length(itervec)){
#     path <- file.path(sim_path, scen_grid[x,"LifeHistory"], scen_grid[x,"F"], scen_grid[x,"SigmaR"], itervec[y])
#     xx <- file.rename(from = file.path(path, "perfect"), to = file.path(path, "perfect_estLH"))
#   }
# }

ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_unadjusted","R_biasadjusted"),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)

saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)

write.csv(res_long, file.path(sim_path, "results", paste0(res_name, "_long.csv")), row.names = FALSE)

true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
res <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv"))) %>% select(-Rich)
all <- left_join(res, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all_info <- all %>% 
  mutate(re = (estimate - true)/true)

sum_final <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(which(max_grad <= 1))/length(max_grad))
write.csv(sum_final, file.path(fig_path, "RE_summary_perfect.csv"))

check <- all_info %>% 
  filter(variable %in% c("Recruit","F","SSB","Depletion")) %>% 
  select(-c(sd,re)) %>%
  pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value") %>%
  filter(iteration == 1) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(year <= 100)
p_perf_check <- ggplot(check) +
  geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
  facet_wrap(lifehistory~variable, scales ="free_y", nrow = 4) +
  theme_bw(base_size = 14)

check <- all_info %>% filter(SigmaR == "HighSigmaR") %>% filter(variable %in% c("F","SPR","SSB","Depletion","SSB0")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("F","SPR","SSB","Depletion","SSB0"))
pre_perf <- ggplot(check) +
  geom_violin( aes(x = lifehistory, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 2) +
  facet_wrap(.~variable, nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(-1,1)) +
  guides(fill = FALSE) +
  theme_bw()

########################
## sampling
########################
########################
## average sample
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR","HighSigmaR")
Lyears <- c("L75","L20","L10","L1")
Lsamp <- c("N200")
Fscen <- c("F1")
res_name <- "results_variable_avgsample"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_unadjusted","R_biasadjusted"),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)

write.csv(res_long, file.path(sim_path, "results", paste0(res_name, "_long.csv")), row.names = FALSE)

true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
res <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv"))) %>% select(-Rich)
all <- left_join(res, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all_info <- all %>% 
  mutate(re = (estimate - true)/true)

sum_converge <- unique(all_info %>% 
                         filter(year == 100) %>%
                         select(iteration, lifehistory, Rest, SigmaR, Lyears, Nsamp, max_grad, LN_R0)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp) %>%
  summarise(pconverge = length(unique(c(which(max_grad <= 1),which(LN_R0<12))))/100)

sum_final <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(unique(c(which(max_grad <= 1), which(LN_R0<12))))/length(max_grad))
write.csv(sum_final, file.path(fig_path, "RE_summary_avgsampling.csv"))

check <- all_info %>% 
  filter(variable %in% c("Depletion")) %>% 
  select(-c(sd,re)) %>%
  pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value") %>%
  filter(iteration == 1) %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(year <= 100) %>%
  filter(SigmaR == "HighSigmaR")
p_samp_check <- ggplot(check) +
  geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
  facet_wrap(lifehistory~variable+Lyears, scales ="free_y", nrow = 4) +
  theme_bw(base_size = 14)

check <- all_info %>% filter(SigmaR == "HighSigmaR") %>% filter(variable %in% c("F","SPR","SSB","Depletion","SSB0")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("F","SPR","SSB","Depletion","SSB0"))
pre_samp <- ggplot(check) +
  geom_violin( aes(x = Lyears, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 1.2) +
  facet_wrap(lifehistory~variable, nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), 5)) +
  guides(fill = FALSE) +
  theme_bw()


######################
## low sample
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR","HighSigmaR")
Lyears <- c("L75","L20","L10","L1")
Lsamp <- c("N50")
Fscen <- c("F1")
res_name <- "results_variable_lowsample"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)
ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_unadjusted","R_biasadjusted"),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, Lyears, Nsamp, LN_R0) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name, "_long.csv")), row.names = FALSE)

true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
res <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv"))) %>% select(-Rich)
all <- left_join(res, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all_info <- all %>% 
  mutate(re = (estimate - true)/true)

sum_converge <- unique(all_info %>% 
  filter(year == 100) %>%
  select(iteration, lifehistory, Rest, SigmaR, Lyears, Nsamp, max_grad, LN_R0)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp) %>%
  summarise(pconverge = length(unique(c(which(max_grad <= 1),which(LN_R0<12))))/100)

sum_final <- all_info %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(unique(c(which(max_grad <= 1), which(LN_R0<12))))/length(max_grad))
write.csv(sum_final, file.path(fig_path, "RE_summary_lowsampling.csv"))

check <- all_info %>% 
  filter(variable %in% c("Recruit","F","SSB","Depletion")) %>% 
  select(-c(sd,re)) %>%
  pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value") %>%
  filter(iteration == 1) %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(SigmaR == "HighSigmaR")
p_samp_check <- ggplot(check) +
  geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
  facet_wrap(lifehistory~variable+Lyears, scales ="free_y", nrow = 4) +
  theme_bw(base_size = 14)


check <- all_info %>% filter(SigmaR == "HighSigmaR") %>% filter(variable %in% c("F","SPR","SSB","Depletion","SSB0")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("F","SPR","SSB","Depletion","SSB0"))
pre_samp <- ggplot(check) +
  geom_violin( aes(x = Lyears, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 1.2) +
  facet_wrap(lifehistory~variable, nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), 5)) +
  guides(fill = FALSE) +
  theme_bw()


######################
## sample decline
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR","HighSigmaR")
Lyears <- c("L20")
Lsamp <- c("Ndecline")
Fscen <- c("F1")
res_name <- "results_variable_sampledecline"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)
ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_grid, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_unadjusted","R_biasadjusted"),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)

saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, Nsamp, Lyears, LN_R0) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)

write.csv(res_long, file.path(sim_path, "results", paste0(res_name, "_long.csv")), row.names = FALSE)

############################################
## READ ALL
############################################
true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
resr <- read.csv(file.path(sim_path, "results", "results_variable_rich_long.csv"))
resp <- read.csv(file.path(sim_path, "results", "results_variable_perfect_long.csv"))
ress <- read.csv(file.path(sim_path, "results", "results_variable_avgsample_long.csv"))
ress2 <- read.csv(file.path(sim_path, "results", "results_variable_lowsample_long.csv"))
ress3 <- read.csv(file.path(sim_path, "results", "results_variable_sampledecline_long.csv"))
res <- rbind.data.frame(resr, resp, ress, ress2, ress3)
all <- left_join(res, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all$lifehistory <- as.character(all$lifehistory)

all_info <- all %>% 
  mutate(re = (estimate - true)/true) %>%
  mutate(label = ifelse(Rich == TRUE, 'data-rich',
                        ifelse(Rich == FALSE & Nsamp == "perfect", "perfect",
                               paste0(Lyears,"_", Nsamp)))) %>%
  mutate(converge = ifelse(max_grad > 1 | LN_R0 > 12, 0, 1),
         variable = replace(variable, variable == "Depletion", "Fraction of unfished")) %>%
  mutate(longevity = ifelse(grepl("short_",lifehistory), "Short-lived", "Long-lived")) %>%
  mutate(growth = ifelse(grepl("_slow", lifehistory), "Slow-growing", "Fast-growing"))
labels <- unique(all_info$label)
saveRDS(all_info, file.path(sim_path, "results", "results_together.rds"))

sum_converge <- unique(all_info %>% 
                         filter(year == 100) %>%
                         select(iteration, lifehistory, Rest, SigmaR, Lyears, Nsamp, max_grad, LN_R0, label, converge)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp) %>%
  summarise(n_nonconverge = length(which(converge == 0)))

sum_re1 <- unique(all_info %>% 
                    filter(year == 100) %>%
                    filter(variable == "Fraction of unfished") %>%
                    select(iteration, lifehistory, Rest, SigmaR, Lyears, Nsamp, label, re)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, label) %>%
  summarise(n_re1 = length(which(abs(re)>1)))
write.csv(sum_re1, file.path(fig_path, "N_RE_greater_than_1.csv"))

sum_final <- all_info %>%
  filter(converge == 1) %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, variable, year, label) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE)) %>%
  left_join(sum_converge)
write.csv(sum_final, file.path(fig_path, "RE_summary_all.csv"))


# check <- all_info %>% 
#   filter(variable %in% c("Recruit","F","SSB","Depletion")) %>% 
#   select(-c(sd,re)) %>%
#   pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value") %>%
#   filter(iteration == 1) %>%
#   filter(lifehistory == "short_fast") %>%
#   filter(Rest == "R_biasadjusted")
# check$label <- factor(check$label, levels = labels)
# p_samp_check <- ggplot(check) +
#   geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
#   facet_wrap(lifehistory~variable+label, scales ="free_y", nrow = 4) +
#   theme_bw(base_size = 14)
# ggsave(file.path(fig_path, "Check_iter_example.png"), p_samp_check, height = 10, width = 12)
# 

###### Low Sigma R, all scenarios
# vars <- c("SSB0","SSB","Fraction of unfished","OFL")
# check <- all_info %>% filter(SigmaR == "LowSigmaR") %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
# check$variable <- factor(check$variable, levels = vars)
# check$label <- factor(check$label, levels = labels)
# pre_samp <- ggplot(check) +
#   geom_hline(aes(yintercept = 0), lwd = 1) +
#   geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
#   facet_grid(longevity+growth~variable) +
#   xlab("Sampling") +
#   ylab("Relative error") +
#   scale_fill_brewer(palette = "Set1") +
#   coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.995))) +
#   scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
#   guides(fill = FALSE) +
#   theme_bw(base_size = 24)
# ggsave(file.path(fig_path, "RE_LowSigmaR.png"), pre_samp, height = 12, width = 28)
# 
# ###### High Sigma R, all scenarios
# vars <- c("SSB0","SSB","Fraction of unfished","OFL")
# check <- all_info %>% filter(SigmaR == "HighSigmaR") %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
# check$variable <- factor(check$variable, levels = vars)
# check$label <- factor(check$label, levels = labels)
# pre_samp <- ggplot(check) +
#   geom_hline(aes(yintercept = 0), lwd = 1) +
#   geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
#   facet_wrap(lifehistory~variable, nrow = 4) +
#   xlab("Sampling") +
#   ylab("Relative error") +
#   scale_fill_brewer(palette = "Set1") +
#   coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.995))) +
#   scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
#   guides(fill = FALSE) +
#   theme_bw(base_size = 14)
# ggsave(file.path(fig_path, "RE_HighSigmaR.png"), pre_samp, height = 10, width = 20)

## compare low and high sigmaR, Fraction of unfished
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("Fraction of unfished")
check <- all_info %>% filter(label %in% labels) %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = vars)
check$label <- factor(check$label, levels = labels)
check$SigmaR <- as.character(check$SigmaR)
check <- check %>% mutate(SigmaR = replace(SigmaR, SigmaR == "LowSigmaR", 0.4),
                          SigmaR = replace(SigmaR, SigmaR == "HighSigmaR", 0.8))
check$Lyears <- factor(check$Lyears, levels = c("L100","L75","L20","L10","L1"))
pre_samp <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  # geom_hline(aes(yintercept = 1), lty = 2) +
  # geom_hline(aes(yintercept = 2), lty = 2) +
  geom_violin( aes(x = factor(SigmaR), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(Nsamp+Lyears~longevity+growth) +
  xlab("SigmaR") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  # scale_x_discrete(name = "SigmaR", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Depl_SigmaR_REcompare.png"), height = 15, width =10)

## compare low and high sigmaR, OFL
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("OFL")
check <- all_info %>% filter(label %in% labels) %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = vars)
check$label <- factor(check$label, levels = labels)
check$SigmaR <- as.character(check$SigmaR)
check <- check %>% mutate(SigmaR = replace(SigmaR, SigmaR == "LowSigmaR", 0.4),
                          SigmaR = replace(SigmaR, SigmaR == "HighSigmaR", 0.8))
check$Lyears <- factor(check$Lyears, levels = c("L100","L75","L20","L10","L1"))
pre_samp <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  # geom_hline(aes(yintercept = 1), lty = 2) +
  # geom_hline(aes(yintercept = 2), lty = 2) +
  geom_violin( aes(x = factor(SigmaR), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(Nsamp+Lyears~longevity+growth) +
  xlab("SigmaR") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  # scale_x_discrete(name = "SigmaR", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "OFL_SigmaR_REcompare.png"), height = 15, width =10)

###### Low SigmaR , only N50
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("SSB0", "SSB","Fraction of unfished","OFL")
check <- all_info %>% filter(label %in% labels) %>% filter(SigmaR == "LowSigmaR") %>% filter(Nsamp %in% c("perfect","N50")) %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = vars)
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
pre_samp <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(longevity+growth~variable) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.995))) +
  scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "RE_LowSigmaR_N50.png"), pre_samp, height = 10, width = 18)

### low SigmaR, only N200
vars <- c("SSB0", "SSB","Fraction of unfished","OFL")
check <- all_info %>% filter(SigmaR == "LowSigmaR") %>% filter(Nsamp %in% c("perfect","N200")) %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = vars)
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
pre_samp <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(longevity+growth~variable) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "RE_LowSigmaR_N200.png"), pre_samp, height = 10, width = 18)

###### High SigmaR , only N50
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("SSB0", "SSB","Fraction of unfished","OFL")
check <- all_info %>% filter(label %in% labels) %>% filter(SigmaR == "HighSigmaR") %>% filter(Nsamp %in% c("perfect","N50")) %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = vars)
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
pre_samp <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(longevity+growth~variable) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.995))) +
  scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "RE_HighSigmaR_N50.png"), pre_samp, height = 10, width = 18)

### high SigmaR, only N200
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("SSB0", "SSB","Fraction of unfished","OFL")
check <- all_info %>% filter(label %in% labels) %>% filter(SigmaR == "HighSigmaR") %>% filter(Nsamp %in% c("perfect","N200")) %>% filter(variable %in% vars) %>% filter(converge == 1) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = vars)
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
pre_samp <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(longevity+growth~variable) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "RE_HighSigmaR_N200.png"), pre_samp, height = 10, width = 18)

### 20 years
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("SSB0","SSB","Fraction of unfished","OFL")
check <- all_info %>%
  filter(converge == 1) %>%
  filter(Lyears %in% c("L100", "L20")) %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(year == 100) %>%
  filter(variable %in% vars) %>%
  filter(label %in% labels) %>%
  filter(SigmaR == "LowSigmaR")
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
check$variable <- factor(check$variable, levels = vars)
pre20 <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(longevity+growth~variable) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.0001), quantile(check$re,0.9999))) +
  scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24) 
ggsave(file.path(fig_path, "RE_LowSigmaR_L20.png"), height = 10, width = 20)

if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("SSB0","SSB","Fraction of unfished","OFL")
check <- all_info %>%
  filter(converge == 1) %>%
  filter(Lyears %in% c("L100", "L20")) %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(year == 100) %>%
  filter(variable %in% vars) %>%
  filter(label %in% labels) %>%
  filter(SigmaR == "HighSigmaR")
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
check$variable <- factor(check$variable, levels = vars)
pre20 <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_grid(longevity+growth~variable) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.0001), quantile(check$re,0.9999))) +
  scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24) 
ggsave(file.path(fig_path, "RE_HighSigmaR_L20.png"), height = 10, width = 20)

ggsave(file.path(fig_path, "RE_HighSigmaR_L20.png"), height = 10, width = 12)

#############
## by year
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
check <- all_info %>% filter(Nsamp == "N50") %>% 
  filter(variable == "Fraction of unfished") %>% 
  filter(converge == 1) %>% 
  filter(label %in% labels) %>%
  filter(Rest == "R_biasadjusted") %>% 
  filter(year %in% seq(100,1,by = -5))
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
pret <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin( aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_wrap(lifehistory~label, nrow = 4) +
  xlab("Year") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re, 0.999))) +
  guides(fill = FALSE) +
  theme_bw()
ggsave(file.path(fig_path, "RE_byYear_Depletion_LowSigmaR_N50.png"), pret, height = 10, width = 18)

if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
check <- all_info %>% filter(Nsamp == "N200") %>% 
  filter(variable == "Fraction of unfished") %>% 
  filter(converge == 1) %>% 
  filter(Rest == "R_biasadjusted") %>% 
  filter(label %in% labels) %>%
  filter(year %in% seq(100,1,by = -5)) 
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
pret <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin( aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_wrap(lifehistory~label, nrow = 4) +
  xlab("Year") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re, 0.999))) +
  guides(fill = FALSE) +
  theme_bw()
ggsave(file.path(fig_path, "RE_byYear_Depletion_LowSigmaR_N200.png"), pret, height = 10, width = 18)

### 20 years
if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
vars <- c("SSB0","SSB","Fraction of unfished","OFL")
check <- all_info %>%
  filter(converge == 1) %>%
  filter(Lyears %in% c("L100", "L20")) %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(label %in% labels) %>%
  filter(year == 100) %>%
  filter(variable %in% vars)
labels <- as.character(unique(check$label))
check$label <- factor(check$label, levels = labels)
check$variable <- factor(check$variable, levels = vars)
pre20 <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin( aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_wrap(lifehistory~variable, nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.0001), quantile(check$re,0.9999))) +
  scale_x_discrete(name = "Sampling", labels = str_replace(labels,"_","\n")) +
  guides(fill = FALSE) +
  theme_bw() 
ggsave(file.path(fig_path, "RE_HighSigmaR_L20.png"), height = 10, width = 12)

if(any(labels == "data-rich")) labels <- labels[-which(labels == "data-rich")]
check <- all_info %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(converge == 1) %>%
  filter(variable %in% c("Fraction of unfished", "Recruitment")) %>%
  filter(Nsamp %in% c("perfect","N200"))
check$Rest <- as.character(check$Rest)
check <- check %>%
  mutate(Rest = replace(Rest, Rest == "R_unadjusted", "Unadjusted"),
         Rest = replace(Rest, Rest == "R_biasadjusted", "Bias-adjusted"))
labels <- as.character(unique(check$label)) 
check$label <- factor(check$label, levels = labels)
prest <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin(aes(x = Rest, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25, 0.5,0.75)) +
  facet_wrap(lifehistory~label, nrow = 4) +
  xlab("Recruitment estimation model") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(quantile(check$re,0.001), quantile(check$re,0.999))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 14)
ggsave(file.path(fig_path, "RE_byRec_Depletion_LowSigmaR_N200.png"), prest, height = 10, width = 14)

check <- all_info %>% filter(variable %in% c("F","SPR","SSB","Depletion","SSB0")) %>% 
  filter(converge == 1) %>% 
  filter(Rest == "R_biasadjusted") %>% 
  filter(year %in% seq(100,1,by = -5)) %>%
  filter(variable == "F")
check$variable <- factor(check$variable, levels = c("F","SPR","SSB","Depletion","SSB0"))
check$label <- factor(check$label, levels = c("data-rich", "perfect", "L75_N200", "L1_N200"))
pret <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin( aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_wrap(lifehistory~label, nrow = 4) +
  xlab("Year") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re, 0.999))) +
  guides(fill = FALSE) +
  theme_bw()
ggsave(file.path(fig_path, "RE_byYear_F_LowSigmaR.png"), pret, height = 10, width = 12)



###########################
### Pull data
###########################
lh_vec <- lh_info$ShortName
SigRscen <- c("LowSigmaR")
Lyears <- c("L75")
Lsamp <- c("N200")
Fscen <- c("F1")
res_name <- "results_variable_sample"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

data <- pull_data(path = sim_path, df = scen_grid, itervec = 1:100)
write.csv(data, file.path(sim_path, "results", "data_CL_sample.csv"), row.names = FALSE)



check_ic <- sum_final2 %>% filter(variable %in% c("F","SPR","SSB","Depletion")) 
pic <- ggplot(check_ic) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_point(aes(x = scenario, y = pcover, color = Rest), cex = 4, alpha = 0.8) +
  facet_wrap(lifehistory~variable) +
  ylim(c(0,1)) +
  theme_bw()
ggsave(file.path(fig_path, "IC_LowSigmaR.png"), pic, height = 10, width = 18)


check_mre <- lapply(1:length(itervec), function(x){
  sub <- all_info %>% filter(year == 100) %>% filter(variable %in% c("Depletion", "F", "SSB")) %>% filter(iteration <= itervec[x])
  sum <- sub %>%
    group_by(lifehistory, SigmaR, label, variable, Rest) %>% 
    summarise(mre = median(re)) %>%
    mutate(iteration = itervec[x])
  return(sum)
})
check_mre <- do.call(rbind, check_mre)
piter <- ggplot(check_mre) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(x = iteration, y = mre, color = variable, linetype = Rest)) + 
  facet_grid(lifehistory~SigmaR+label) +
  coord_cartesian(ylim = c(min(check_mre$mre), quantile(check_mre$mre, 0.99))) +
  theme_bw()