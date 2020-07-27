rm(list = ls())

####################
## libraries
####################
devtools::install_github("ss3sim/ss3sim", 
                         ref = "v1.0.0", #"v1.0.3", 
                         dependencies = TRUE)
# 
# devtools::install_github("r4ss/r4ss", ref="master")

library(ss3sim)
library(r4ss)
library(tidyverse)
library(foreach)
library(doParallel)


####################
## directories
####################

proj_path <- "~/NOAA/datamoderate"
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
# zero <- matrix(0, nrow = Pyears, ncol = length(itervec))
recdev_mat3 <- matrix(rnorm((Pyears+Nyears+Fyears)*length(itervec), mean = 0, 0.8),ncol=length(itervec))

f1 <- c(seq(0.01, 0.05, length.out = 25), seq(0.05, 2, length.out = 30), rep(2, 5), seq(2, 0.6, length.out = 20), rep(0.6,20))
# f2 <- c(seq(0.01, 0.05, length.out = 25), seq(0.05, 1, length.out = 30), rep(1, 5), seq(1, 0.3, length.out = 20), rep(0.3,20))
# f3 <- c(seq(0.01, 0.05, length.out = 25), seq(0.05, 2, length.out = 30), rep(2, 5), seq(2, 0.6, length.out = 40))

###########################################
## Operating model - generate true population
###########################################
########################
### deterministic
########################
## scenarios to simulate
lh_vec <- lh_info$ShortName[4]
SigRscen <- c("Deterministic")
Fscen <- c("F1")
name <- "sim_deterministic"

sim_grid <- expand.grid("SigmaR" = SigRscen, "F"=Fscen, "LifeHistory" = lh_vec)

## simulate scenarios
itervec <- 1
rewrite = TRUE
start <- Sys.time()
sim_det <- simulate_pop(df = sim_grid, 
                         path = sim_path,
                         itervec = itervec, 
                         f_list = list("F1" = f1),
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

p <- ggplot(sim_long) +
geom_line(aes(x = year, y = true)) +
facet_wrap(~variable, scales = "free_y") +
theme_bw()
 

########################
### variable
########################
## scenarios to simulate
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR")#,"LowSigmaR")
Fscen <- c("F1")#,"F2","F3")
name <- "sim_variable"

sim_grid <- expand.grid("SigmaR" = SigRscen, "F"=Fscen, "LifeHistory" = lh_vec)

## simulate scenarios
itervec <- 1:100
rewrite <- FALSE
start <- Sys.time()
sim_var <- simulate_pop(df = sim_grid, 
                        path = sim_path,
                        itervec = itervec, 
                        f_list = list("F1" = f1),
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
  dplyr::select(-c(max_grad, OFL, OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, SigmaR, Rest, Fscen), names_to = "variable", values_to = "true") %>%
  filter(grepl("_sd", variable) == FALSE)
write.csv(sim_long, file.path(sim_path, "results", paste0(name,"_long.csv")), row.names = FALSE)

### explore true population
# sim_det <- read.csv(file.path(sim_path,"results", "sim_deterministic_long.csv"), stringsAsFactors = FALSE)
# sim_var <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv"), stringsAsFactors = FALSE)

# sim_long <- sim_var #full_join(sim_det, sim_var)
write.csv(sim_long, file.path(sim_path, "results", "sim_all_long.csv"), row.names = FALSE)

sim_long <- read.csv(file.path(sim_path, "results", "sim_all_long.csv"))
sim_long$SigmaR <- as.character(sim_long$SigmaR)
sim_long <- sim_long %>%
  mutate(SigmaR_desc = replace(SigmaR, which(SigmaR == "LowSigmaR"), "SigmaR = 0.4"),
         SigmaR_desc = replace(SigmaR_desc, which(SigmaR_desc == "HighSigmaR"), "SigmaR = 0.8")) %>%
  mutate(Longevity = ifelse(grepl("long_", lifehistory), "Longer-lived", "Shorter-lived")) %>%
  mutate(Growth = ifelse(grepl("_slow", lifehistory), "Slower-to-Linf", "Faster-to-Linf"))

check <- sim_long %>% filter(variable == "Recruit") %>% filter(year <= 100) %>% filter(Fscen == "F1")
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

check <- sim_long %>% filter(variable == "Depletion") %>% filter(year <= 100) %>% filter(Fscen == "F1")
p_deplom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  geom_line(data = check %>% filter(iteration == 1)) +
  geom_line(data = check %>% filter(iteration == 2)) +
  geom_line(data = check %>% filter(iteration == 3)) +
  geom_line(data = check %>% filter(iteration == 4)) +
  geom_line(data = check %>% filter(iteration == 5)) +
  facet_grid(Longevity + Growth ~ SigmaR_desc + Fscen) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Depletion") +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Depletion_OM.png"), p_deplom, height = 10, width = 12)

# check <- sim_long %>% filter(variable == "Depletion") %>% filter(year <= 100) %>% filter(Fscen == "F3")
# p_deplom <- ggplot(check, aes(x = year, y = true)) +
#   stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
#   stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
#   stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
#   geom_line(data = check %>% filter(iteration == 1)) +
#   facet_grid(Longevity + Growth ~ SigmaR_desc + Fscen) +
#   scale_fill_brewer(palette = "Set1") +
#   scale_color_brewer(palette = "Set1") +
#   guides(fill = FALSE, color = FALSE, linetype = FALSE) +
#   expand_limits(y = 0) +
#   xlab("Year") + ylab("Depletion") +
#   theme_bw(base_size = 24)
# ggsave(file.path(fig_path, "Depletion_OM_F3.png"), p_deplom, height = 10, width = 12)


# check <- sim_long %>% filter(variable == "Depletion") %>% filter(year <= 100)
# p_deplom_v2 <- ggplot(check, aes(x = year, y = true)) +
#   stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
#   stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
#   stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
#   geom_line(data = check %>% filter(iteration == 1)) +
#   facet_grid(Longevity + Growth ~ SigmaR_desc + Fscen) +
#   scale_fill_brewer(palette = "Set1") +
#   scale_color_brewer(palette = "Set1") +
#   guides(fill = FALSE, color = FALSE, linetype = FALSE) +
#   expand_limits(y = 0) +
#   xlab("Year") + ylab("Depletion") +
#   theme_bw(base_size = 24)
# ggsave(file.path(fig_path, "Depletion_OM_F2.png"), p_deplom_v2, height = 10, width = 18)

sub <- check%>%
  filter(year == 100) %>%
  group_by(lifehistory, SigmaR, Fscen) %>%
  summarise(med = median(true, na.rm = TRUE))

check <- sim_long %>% filter(variable == "F") %>% filter(year <= 100) %>% filter(Fscen == "F1")
p_fom <- ggplot(check, aes(x = year, y = true)) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
  stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
  facet_grid(Longevity + Growth ~ SigmaR_desc + Fscen) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE) +
  expand_limits(y = 0) +
  xlab("Year") + ylab("Fishing mortality rate") +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "FishingMortality_OM.png"), p_fom, height = 10, width = 12)

# check <- sim_long %>% filter(variable == "F") %>% filter(year <= 100) 
# p_fom_v2 <- ggplot(check, aes(x = year, y = true)) +
#   stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.05), fun.max = function(x) stats::quantile(x, 0.95), geom = "ribbon", alpha = 0.125) +
#   stat_summary(aes(fill = lifehistory), fun.min = function(x) stats::quantile(x, 0.25), fun.max = function(x) stats::quantile(x, 0.75), geom = "ribbon", alpha = 0.25) +
#   stat_summary(aes(colour = lifehistory), fun = function(x) stats::quantile(x, 0.5), geom = "line", lwd = 2) +
#   facet_grid(Longevity + Growth ~ SigmaR_desc + Fscen) +
#   scale_fill_brewer(palette = "Set1") +
#   scale_color_brewer(palette = "Set1") +
#   guides(fill = FALSE, color = FALSE, linetype = FALSE) +
#   expand_limits(y = 0) +
#   xlab("Year") + ylab("Fishing mortality rate") +
#   theme_bw(base_size = 24)
# ggsave(file.path(fig_path, "FishingMortality_OM_F2.png"), p_fom_v2, height = 10, width = 18)


##########################################
## Estimation model
##########################################
#############################
## run perfect deterministic
#############################
lh_vec <- lh_info$ShortName[4]
SigRscen <- c("Deterministic")
Lyears <- c("L100")
Lsamp <- c("perfect")
res_name <- "results_det"
Rich <- FALSE
Fscen <- "F1"

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich, "Fix" = NA, "Adjust" = NA)

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


res_long <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv")))
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
lh_vec <- lh_info$ShortName[4]
SigRscen <- c("HighSigmaR")#,"LowSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")
Fscen <- c("F1")#,"F2")
res_name <- "results_variable_rich"
Rich <- c(TRUE)

scen_grid_rich <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich, "Fix" = NA, "Adjust" = NA)
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

ncores <- 1
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

res_long <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv")))
true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
all <- left_join(res_long, true, join_by = c(year, lifehistory, iteration, SigmaR, Fscen, variable))
all$variable <- as.character(all$variable)
all_info <- all %>% 
  mutate(re = (estimate - true)/true) %>%
  mutate(label = ifelse(Rich == TRUE, 'data-rich',
                        ifelse(Rich == FALSE & Nsamp == "perfect", "perfect",
                               paste0(Lyears,"_", Nsamp)))) %>%
  mutate(converge = ifelse(max_grad > 1 | LN_R0 > 12, 0, 1),
         variable = replace(variable, variable == "Depletion", "Fraction of unfished")) %>%
  mutate(longevity = ifelse(grepl("short_",lifehistory), "Shorter-lived", "Longer-lived")) %>%
  mutate(growth = ifelse(grepl("_slow", lifehistory), "Slower-to-Linf", "Faster-to-Linf"))

sum_final <- all_info %>%
  filter(converge == 1) %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, Fscen, variable, year, label) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE))

sum_selex <- sum_final %>% filter(Rest == "R_biasadjusted") %>% filter(grepl("Selex", variable))

check <- all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  filter(year %in% seq(100,1,by=-5)) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(Fscen == "F1")
p_rec <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin(aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = FALSE) +
  ylab("Year") + xlab("RE") +
  facet_grid(longevity~growth) +
  theme_bw(base_size = 20)

check <- all_info %>%
  filter(variable == "Recruit") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(iteration == min(iteration)) %>%
  filter(year <= 100)
p_rec2 <- ggplot(check) +
  geom_line(aes(x = year, y = true), lwd = 2) +
  geom_line(aes(x = year, y = estimate, col = "red"), lwd = 2) +
  theme_bw()

check <-  all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  # filter(year %in% 75:100) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(iteration %in% 1) %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(Fscen == "F1")
p_iter <- ggplot(check) +
  geom_line(aes(x = year, y = true), alpha = 0.8, lwd = 1) +
  geom_line(aes(x = year, y = estimate), alpha = 0.8, lwd = 2, color = "red") +
  facet_wrap(iteration~longevity+growth) +
  theme_bw(base_size = 20)

check <- all_info %>% 
  filter(variable %in% c("Recruit","F","Fraction of unfished", "SSB")) %>% 
  select(-c(sd,re)) %>%
  pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value") %>%
  filter(iteration == 1) %>%
  filter(year <= 100) %>%
  filter(Fscen == "F1") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(SigmaR == "HighSigmaR")
p_rich_check <- ggplot(check) +
  geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
  facet_wrap(lifehistory~variable, scales ="free_y", nrow = 4) +
  theme_bw(base_size = 14)

check <- all_info %>% 
      filter(variable %in% c("Fraction of unfished","OFL","Recruit", "SSB0","Selex_peak")) %>% 
      # filter(year == 100) %>% 
      filter(Rest == "R_biasadjusted") #%>% filter(abs(max_grad) <= 1e-4)
check$variable <- factor(check$variable, levels = c("Fraction of unfished","OFL","Recruit", "SSB0","Selex_peak"))
pre_rich <- ggplot(check) +
  geom_violin( aes(x = lifehistory, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 2) +
  facet_wrap(Fscen~variable, scales = "free_y") +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(-1,1)) +
  guides(fill = FALSE) +
  theme_bw()

check_mre <- lapply(1:length(itervec), function(x){
  sub <- all_info %>% filter(year == 100) %>% filter(variable %in% c("OFL","Fraction of unfished","SSB0")) %>% filter(iteration <= itervec[x])
  sum <- sub %>%
    group_by(lifehistory, SigmaR, Fscen, variable, Rest) %>% 
    summarise(mre = median(re)) %>%
    mutate(iteration = itervec[x])
  return(sum)
})
check_mre <- do.call(rbind, check_mre)
piter <- ggplot(check_mre %>% filter(Fscen == "F1") %>% filter(SigmaR == "HighSigmaR")) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(x = iteration, y = mre, color = variable, linetype = Rest)) + 
  facet_grid(lifehistory~SigmaR) +
  coord_cartesian(ylim = c(min(check_mre$mre), quantile(check_mre$mre, 0.99))) +
  theme_bw()

########################
## perfect
########################
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR")#,"LowSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")
Fscen <- c("F1")#,"F2")
res_name <- "results_variable_perfect"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich, "Fix" = NA, "Adjust" = NA)


# r1 <- SS_output(file.path(sim_path, "v2", "long_fast", "F1", "HighSigmaR", 2, "om"))
# r2 <- SS_output(file.path(sim_path, "v2", "long_fast", "F1", "HighSigmaR", 2, "perfect", "L100", "R_biasadjusted"))

# plot(r1$recruit$dev[which(r1$recruit$Yr %in% 1:100)], type = "l")
# lines(r2$recruit$dev[which(r2$recruit$Yr %in% 1:100)], col = "blue")
# itervec <- 1:100
# for(i in 1:nrow(scen_grid)){
#   for(j in itervec){
#     path <- file.path(sim_path, scen_grid[i,"LifeHistory"], scen_grid[i,"F"], scen_grid[i,"SigmaR"], j)
#     # file.rename(from = file.path(path, "perfect"), to = file.path(path, "perfect_known"))
#     unlink(file.path(path, "perfect"), TRUE)
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


ncores <- 1
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
                        Rscen = c("R_biasadjusted"),
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
  
  res_long <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv")))
  true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
  all <- left_join(res_long, true, join_by = c(year, lifehistory, iteration, SigmaR, Fscen, variable))
  all_info <- all %>% 
    mutate(re = (estimate - true)/true) %>%
    mutate(label = ifelse(Rich == TRUE, 'data-rich',
                          ifelse(Rich == FALSE & Nsamp == "perfect", "perfect",
                                 paste0(Lyears,"_", Nsamp)))) %>%
    mutate(converge = ifelse(max_grad > 1 | LN_R0 > 12, 0, 1)) %>%
    # mutate(variable = replace(variable, variable == "Depletion", "Fraction of unfished")) %>%
    mutate(longevity = ifelse(grepl("short_",lifehistory), "Shorter-lived", "Longer-lived")) %>%
    mutate(growth = ifelse(grepl("_slow", lifehistory), "Slower-to-Linf", "Faster-to-Linf"))

sum_final <- all_info %>%
  filter(converge == 1) %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, Fscen, variable, year, label) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE))

sum_selex <- sum_final %>% filter(Rest == "R_biasadjusted") %>% filter(grepl("Selex", variable))

check <- all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  filter(year %in% seq(100,1,by=-5)) %>%
  filter(SigmaR == "HighSigmaR") %>%
  # filter(Rest == "R_biasadjusted") %>%
  filter(Fscen == "F1")
p_rec <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin(aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = FALSE) +
  ylab("Year") + xlab("RE") +
  facet_grid(longevity~growth) +
  theme_bw(base_size = 20)

check <-  all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  filter(year %in% 75:100) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(Fscen == "F1") %>%
  filter(iteration == 2)
p_iter <- ggplot(check) +
  geom_line(aes(x = year, y = true), alpha = 0.8, lwd = 1) +
  geom_line(aes(x = year, y = estimate), alpha = 0.8, lwd = 2, color = "red") +
  facet_wrap(iteration~longevity+growth) +
  theme_bw(base_size = 20)

check <- all_info %>% 
  filter(variable %in% c("Depletion", "OFL")) %>% 
  select(-c(sd,re)) %>%
  pivot_longer(cols=c(estimate,true), names_to = "model", values_to = "value") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(year <= 100) %>%
  filter(Fscen == "F1") %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(iteration == min(iteration))
p_check <- ggplot(check) +
  geom_line(aes(x = year, y = value, color = model, linetype = model), lwd = 2) +
  facet_wrap(lifehistory~variable, scales ="free_y", nrow = 4) +
  theme_bw(base_size = 14)

check <- all_info %>% 
  filter(variable %in% c("OFL","Depletion","SSB0", "Selex_peak")) %>% 
  filter(year == 100) %>% 
  filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("OFL","Depletion","SSB0", "Selex_peak"))
pre <- ggplot(check) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 2) +
  geom_violin( aes(x = lifehistory, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_wrap(variable~., nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(-1,1)) +
  guides(fill = FALSE) +
  theme_bw()

check <- all_info %>% filter(variable %in% c("OFL", "Depletion")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("OFL", "Depletion"))
pre <- ggplot(check) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 2) +
  geom_violin( aes(x = lifehistory, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  facet_wrap(Fscen+SigmaR~variable, nrow = 4) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(-1,1)) +
  guides(fill = FALSE) +
  theme_bw()

check_mre <- lapply(1:length(itervec), function(x){
  sub <- all_info %>% filter(year %in% 100) %>% filter(variable %in% c("OFL", "Fraction of unfished")) %>% filter(iteration <= itervec[x])
  sum <- sub %>%
    group_by(lifehistory, SigmaR, Fscen, variable, Rest, year) %>% 
    summarise(mre = median(re)) %>%
    mutate(iteration = itervec[x])
  return(sum)
})
check_mre <- do.call(rbind, check_mre)
piter <- ggplot(check_mre %>% filter(Fscen == "F1") %>% filter(SigmaR == "HighSigmaR") %>% filter(Rest == "R_biasadjusted")) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(x = iteration, y = mre, color = factor(year), linetype = variable)) + 
  facet_grid(lifehistory~SigmaR) +
  coord_cartesian(ylim = c(min(check_mre$mre), quantile(check_mre$mre, 0.99))) +
  theme_bw()

########################
## sampling
########################
########################
## average sample
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR")#,"LowSigmaR")
Lyears <- c("L75","L20","L10","L1")
Lsamp <- c("N200")
Fscen <- c("F1")#,"F2")
res_name <- "results_variable_avgsample"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich, "Fix" = NA, "Adjust" = NA)

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

ncores <- 1
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
                        Rscen = c("R_biasadjusted"),
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
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, Fscen, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(unique(c(which(max_grad <= 1), which(LN_R0<12))))/length(max_grad))
write.csv(sum_final, file.path(fig_path, "RE_summary_avgsampling.csv"))

sub <- sum_final %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(variable == "Depletion")

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

check <- all_info %>% filter(SigmaR == "HighSigmaR") %>% filter(variable %in% c("OFL","Depletion")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("OFL","Depletion"))
check$Lyears <- factor(check$Lyears, levels = c("L75","L20","L10","L1"))
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
SigRscen <- c("HighSigmaR")#,"LowSigmaR")
Lyears <- c("L75","L20","L10","L1")
Lsamp <- c("N50")
Fscen <- c("F1")#,"F2")
res_name <- "results_variable_lowsample"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich, "Fix" = NA, "Adjust" = NA)

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


ncores <- 1
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

res_long <- read.csv(file.path(sim_path, "results", paste0(res_name, "_long.csv")))
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
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp,Fscen, variable, year) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE),
            pconverge = length(unique(c(which(max_grad <= 1), which(LN_R0<12))))/length(max_grad))
write.csv(sum_final, file.path(fig_path, "RE_summary_lowsampling.csv"))


check <- all_info %>% filter(SigmaR == "HighSigmaR") %>% filter(variable %in% c("OFL","Depletion", "SSB0","Selex_peak")) %>% filter(year == 100) %>% filter(Rest == "R_biasadjusted")
check$variable <- factor(check$variable, levels = c("OFL","Depletion", "SSB0","Selex_peak"))
check$Lyears <- factor(check$Lyears, levels = c("L75","L20","L10","L5","L2","L1"))
pre_samp <- ggplot(check) +
  geom_violin( aes(x = Lyears, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd = 1.2) +
  facet_grid(lifehistory~variable) +
  xlab("Sampling") +
  ylab("Relative error") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), 5)) +
  guides(fill = FALSE) +
  theme_bw()


######################
## sample decline
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR")#,"LowSigmaR")
Lyears <- c("L20")
Lsamp <- c("Ndecline")
Fscen <- c("F1")#,"F2")
res_name <- "results_variable_sampledecline"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich, "Fix" = NA, "Adjust" = NA)

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

ncores <- 1
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


######################
## sensitivities to Linf and CV
## linf low
lh_vec <- lh_info$ShortName[4]

sens1 <- expand.grid("SigmaR" = "HighSigmaR",
                    "LifeHistory" = lh_vec,
                    "L" = "L100", 
                    "Samp" = "perfect",
                    "F" = "F1",
                    "Rich" = FALSE, 
                    "Fix" = NA,
                    "Adjust" = c("CV_0.125", "CV_0.075", "Linf_49.5","Linf_60.5")) #c()) #, ))

sens2 <- expand.grid("SigmaR" = "HighSigmaR",
                    "LifeHistory" = lh_vec,
                    "L" = c("L75","L20"), 
                    "Samp" = "N200",
                    "F" = "F1",
                    "Rich" = FALSE, 
                    "Fix" = NA,
                    "Adjust" = c("CV_0.125", "CV_0.075", "Linf_49.5", "Linf_60.5")) #c()) #, "Linf_60.5", "CV_0.125", "CV_0.075")) 

scen_grid <- rbind.data.frame(sens1, sens2)

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

ncores <- 1
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


res_name <- "CV_0.075"
scen <- scen_grid %>% filter(grepl(res_name, Adjust))
res_wide1 <- get_results(mod_path = sim_path, 
                        df = scen, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c(paste0("R_unadjusted_change_", res_name)),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)

res_name <- "CV_0.125"
scen <- scen_grid %>% filter(grepl(res_name, Adjust))
res_wide2 <- get_results(mod_path = sim_path, 
                        df = scen, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c(paste0("R_unadjusted_change_", res_name)),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)

res_name <- "Linf_49.5"
scen <- scen_grid %>% filter(grepl(res_name, Adjust))
res_wide3 <- get_results(mod_path = sim_path, 
                        df = scen, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c(paste0("R_unadjusted_change_", res_name)),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)

res_name <- "Linf_60.5"
scen <- scen_grid %>% filter(grepl(res_name, Adjust))
res_wide4 <- get_results(mod_path = sim_path, 
                        df = scen, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c(paste0("R_unadjusted_change_", res_name)),
                        res_name = paste0(res_name, "_wide"),
                        ncores = 14)

res_wide <- rbind.data.frame(res_wide1, res_wide2, res_wide3, res_wide4)

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

write.csv(res_long, file.path(sim_path, "results", "Sensitivities_long.csv"), row.names = FALSE)

############################################
## READ ALL
############################################
true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% select(-Rest)
# resr <- read.csv(file.path(sim_path, "results", "results_variable_rich_long.csv"))
resp <- read.csv(file.path(sim_path, "results", "results_variable_perfect_long.csv"))
ress <- read.csv(file.path(sim_path, "results", "results_variable_avgsample_long.csv"))
ress2 <- read.csv(file.path(sim_path, "results", "results_variable_lowsample_long.csv"))
ress3 <- read.csv(file.path(sim_path, "results", "results_variable_sampledecline_long.csv"))
ressens <- read.csv(file.path(sim_path, "results", "Sensitivities_long.csv"))
res <- rbind.data.frame(resp, ress, ress2, ress3, ressens) #resr, 
res1 <- res %>% select(-c(variable, estimate)) %>% rename(estimate = LN_R0) %>% mutate(variable = "LN_R0") %>% mutate(LN_R0 = estimate)
res2 <- full_join(res, res1)

all <- left_join(res2, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all$lifehistory <- as.character(all$lifehistory)
all$variable <- as.character(all$variable)
all_info <- all %>% 
  mutate(re = (estimate - true)/true) %>%
  mutate(label = ifelse(Rich == TRUE, 'data-rich',
                        ifelse(Rich == FALSE & Nsamp == "perfect", "perfect",
                               paste0(Lyears,"_", Nsamp)))) %>%
  mutate(converge = ifelse(max_grad > 1 | LN_R0 > 12, 0, 1),
        variable = replace(variable, variable == "Depletion", "Fraction of unfished")) %>%
  mutate(longevity = ifelse(grepl("short_",lifehistory), "Shorter-lived", "Longer-lived")) %>%
  mutate(growth = ifelse(grepl("_slow", lifehistory), "Slower-to-Linf", "Faster-to-Linf"))
labels <- unique(all_info$label)
saveRDS(all_info, file.path(sim_path, "results", "results_together.rds"))

############################################
## RESULTS
############################################
all_info <- readRDS(file.path(sim_path, "results", "results_together.rds"))

sum_converge <- unique(all_info %>% 
                         filter(year == 100) %>%
                         select(iteration, lifehistory, Rest, SigmaR, Lyears, Nsamp, Fscen, max_grad, LN_R0, label, converge)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp) %>%
  summarise(n_nonconverge = length(which(converge == 0)))

sum_re1 <- unique(all_info %>% 
                    filter(year == 100) %>%
                    filter(variable == "Fraction of unfished") %>%
                    select(iteration, lifehistory, Rest, SigmaR, Lyears, Nsamp, Fscen, label, re)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, Fscen, label) %>%
  summarise(n_re1 = length(which(abs(re)>1)))
write.csv(sum_re1, file.path(fig_path, "N_RE_greater_than_1.csv"))

sum_final <- all_info %>%
  filter(converge == 1) %>%
  filter(year %in% c(100)) %>%
  group_by(lifehistory, Rest, SigmaR, Lyears, Nsamp, Fscen, variable, year, label) %>%
  summarise(mre = median(re, na.rm = TRUE),
            mare = median(abs(re), na.rm = TRUE)) %>%
  left_join(sum_converge)
write.csv(sum_final, file.path(fig_path, "RE_summary_all.csv"))

sub <- sum_final %>%
  filter(lifehistory == "long_fast") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(variable == "Fraction of unfished") %>%
  filter(Fscen == "F1")


##############################
## figures for presentation
##############################
## Depl, OFL
labels <- unique(all_info$label)
labels_inp <- labels#[-which(labels == "data-rich")]
vars <- c("Fraction of unfished", "OFL")
check <- all_info %>% 
  filter(Fscen == "F1") %>% 
  filter(label %in% labels_inp) %>% 
  filter(variable %in% vars) %>% 
  filter(converge == 1) %>% 
  filter(year == 100) %>%  
  filter(Rest == "R_biasadjusted") %>%
  mutate(label = str_replace(label, "_", "\n"))
check$variable <- factor(check$variable, levels = vars)
labels_inp2 <- unique(check$label)
check$label <- factor(check$label, levels = labels_inp2)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(longevity+growth~variable) +
  ylab("Relative error") + xlab("Data scenario") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.995))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Depl_OFL_finalyear.png"), p, height = 10, width = 18)

## Recruitment, final year
labels <- unique(all_info$label)
# labels_inp <- labels#[-which(labels == "data-rich")]
vars <- c("Recruit")
check <- all_info %>% 
  filter(Fscen == "F1") %>% 
  filter(label %in% labels_inp) %>% 
  filter(variable %in% vars) %>% 
  filter(converge == 1) %>% 
  filter(year %in% c(seq(5,100,by=5))) %>%  
  filter(Rest == "R_biasadjusted") %>%
  filter(label == "perfect")
  # mutate(label = str_replace(label, "_", "\n")) 
# check$variable <- factor(check$variable, levels = vars)
# labels_inp2 <- unique(check$label)
# check$label <- factor(check$label, levels = labels_inp2)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(longevity+growth~variable) +
  ylab("Relative error") + xlab("Year") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.995))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Recruit_perfect.png"), p, height = 10, width = 18)

labels <- unique(all_info$label)
# labels_inp <- labels#[-which(labels == "data-rich")]
vars <- c("Recruit")
check <- all_info %>% 
  filter(Fscen == "F1") %>% 
  filter(label %in% labels_inp) %>% 
  filter(variable %in% vars) %>% 
  filter(converge == 1) %>% 
  filter(year %in% c(seq(5,100,by=5))) %>%  
  filter(Rest == "R_biasadjusted") %>%
  filter(label == "L20_N200")
  # mutate(label = str_replace(label, "_", "\n")) 
# check$variable <- factor(check$variable, levels = vars)
# labels_inp2 <- unique(check$label)
# check$label <- factor(check$label, levels = labels_inp2)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(longevity+growth~variable) +
  ylab("Relative error") + xlab("Year") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.995))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Recruit_L20_N200.png"), p, height = 10, width = 18)


## Rec, Selex
labels <- unique(all_info$label)
labels_inp <- labels#[-which(labels == "data-rich")]
vars <- c("LN_R0")#, "Selex_peak", "Selex_ascend")
check <- all_info %>%
  filter(Fscen == "F1") %>% 
  filter(label %in% labels_inp) %>% 
  filter(variable %in% vars) %>% 
  filter(converge == 1) %>% 
  filter(year == 100) %>%  
  filter(Rest == "R_biasadjusted") %>%
  mutate(label = str_replace(label, "_", "\n"))
check$variable <- factor(check$variable, levels = vars)
labels_inp2 <- unique(check$label)
check$label <- factor(check$label, levels = labels_inp2)
check$lifehistory2 <- paste0(check$longevity, ", ", check$growth)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(longevity+growth~variable, scales = "free_y") +
  ylab("Relative error") + xlab("Data scenario") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,1))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "EstPar_R0.png"), p, height = 10, width = 10)

labels <- unique(all_info$label)
labels_inp <- labels#[-which(labels == "data-rich")]
vars <- c("Selex_peak", "Selex_ascend")
check <- all_info %>%
  filter(Fscen == "F1") %>% 
  filter(label %in% labels_inp) %>% 
  filter(variable %in% vars) %>% 
  filter(converge == 1) %>% 
  filter(year == 100) %>%  
  filter(Rest == "R_biasadjusted") %>%
  mutate(label = str_replace(label, "_", "\n"))
check$variable <- factor(check$variable, levels = vars)
labels_inp2 <- unique(check$label)
check$label <- factor(check$label, levels = labels_inp2)
check$lifehistory2 <- paste0(check$longevity, ", ", check$growth)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(longevity+growth~variable, scales = "free_y") +
  ylab("Relative error") + xlab("Data scenario") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,1))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "EstPar_Selex.png"), p, height = 10, width = 18)

##############################
### sensitivities
##############################
## Depl, OFL
vars <- c("Fraction of unfished", "OFL")
perf <- all_info %>% filter(label == "perfect")
sens <- all_info %>% filter(grepl("change", Rest))
check <- rbind.data.frame(perf, sens)
check$Rest <- as.character(check$Rest)
check <- check %>%
  filter(lifehistory == "long_fast") %>%
  filter(variable %in% vars) %>%
  filter(converge == 1) %>%
  filter(year == 100) %>%
  mutate(Rest = replace(Rest, Rest == 'R_biasadjusted', "perfect")) %>%
  mutate(Rest = replace(Rest, grepl("CV_0.075", Rest), "0.75 * CV")) %>%
  mutate(Rest = replace(Rest, grepl("CV_0.125", Rest), "1.25 * CV")) %>%
  mutate(Rest = replace(Rest, grepl("Linf_49.5", Rest), "0.9 * Linf")) %>%
  mutate(Rest = replace(Rest, grepl("Linf_60.5", Rest), "1.1 * Linf"))

check1 <- check %>% filter(Rest %in% c("perfect", '0.75 * CV','1.25 * CV'))
Rest <- unique(check1$Rest)
check1$Rest <- factor(check1$Rest, levels = Rest)
p <- ggplot(check1) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = variable, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(~Rest) +
  ylab("Relative error") + xlab("Variable") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Depl_OFL_finalyear_CV.png"), p, height = 10, width = 15)

check1 <- check %>% filter(Rest %in% c("perfect", '0.9 * Linf','1.1 * Linf'))
Rest <- unique(check1$Rest)
check1$Rest <- factor(check1$Rest, levels = Rest)
p <- ggplot(check1) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = variable, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(~Rest) +
  ylab("Relative error") + xlab("Variable") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Depl_OFL_finalyear_Linf.png"), p, height = 10, width = 15)


vars <- c("LN_R0", "Selex_peak", "Selex_ascend")
check <- rbind.data.frame(perf, sens)
check$Rest <- as.character(check$Rest)
check <- check %>%
  filter(lifehistory == "long_fast") %>%
  filter(variable %in% vars) %>%
  filter(converge == 1) %>%
  filter(year == 100) %>%
  mutate(Rest = replace(Rest, Rest == 'R_biasadjusted', "perfect")) %>%
  mutate(Rest = replace(Rest, grepl("CV_0.075", Rest), "0.75 * CV")) %>%
  mutate(Rest = replace(Rest, grepl("CV_0.125", Rest), "1.25 * CV")) %>%
  mutate(Rest = replace(Rest, grepl("Linf_49.5", Rest), "0.9 * Linf")) %>%
  mutate(Rest = replace(Rest, grepl("Linf_60.5", Rest), "1.1 * Linf"))

check1 <- check %>% filter(Rest %in% c("perfect", '0.75 * CV','1.25 * CV'))
Rest <- unique(check1$Rest)
check1$Rest <- factor(check1$Rest, levels = Rest)
p <- ggplot(check1) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = variable, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(~Rest) +
  ylab("Relative error") + xlab("Variable") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Est_pars_CV.png"), p, height = 10, width = 18)

check1 <- check %>% filter(Rest %in% c("perfect", '0.9 * Linf','1.1 * Linf'))
Rest <- unique(check1$Rest)
check1$Rest <- factor(check1$Rest, levels = Rest)
p <- ggplot(check1) +
  geom_hline(aes(yintercept = 0), lwd = 1) +
  geom_violin(aes(x = variable, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.05, 0.5, 0.95)) +
  facet_grid(~Rest) +
  ylab("Relative error") + xlab("Variable") +
  scale_fill_brewer(palette = "Set1") +
  # coord_cartesian(ylim=c(quantile(check$re, 0.001), quantile(check$re,0.999))) +
  guides(fill = FALSE) +
  theme_bw(base_size = 24)
ggsave(file.path(fig_path, "Est_pars_Linf.png"), p, height = 10, width = 18)
