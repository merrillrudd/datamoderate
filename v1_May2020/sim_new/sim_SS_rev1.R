rm(list = ls())

####################
## libraries
####################

library(ss3sim)
library(r4ss)
library(tidyverse)
library(foreach)
library(doParallel)
library(MCMCpack)

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

itervec <- 1:100
Nyears <- 100

###############################
## Request 1: Additional outputs
###############################
########################
## perfect
########################
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR")#,"LowSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")
Fscen <- c("F1")#,"F2")
res_name_inp <- "results_variable_rich_all"
Rich <- c(TRUE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)

saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)

write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)


########################
## perfect
########################
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR","LowSigmaR")
Lyears <- c("L100")
Lsamp <- c("perfect")
Fscen <- c("F1")#,"F2")
res_name_inp <- "results_variable_perfect_all"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted"),
                        res_name = paste0(res_name_inp, "_wide"),
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

write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)

########################
## sampling
########################
########################
## average sample
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR","LowSigmaR")
Lyears <- c("L75","L20","L10","L1") #,"L5","L2"
Lsamp <- c("N200")
Fscen <- c("F1")#,"F2")
res_name_inp <- "results_variable_avgsample_all"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted"),
                        res_name = paste0(res_name_inp, "_wide"),
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

write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)

######################
## low sample
lh_vec <- lh_info$ShortName
SigRscen <- c("HighSigmaR","LowSigmaR")
Lyears <- c("L75","L20","L10","L1") #"L5","L2",
Lsamp <- c("N50")
Fscen <- c("F1")#,"F2")
res_name_inp <- "results_variable_lowsample_all"
Rich <- c(FALSE)

scen_grid <- expand.grid("SigmaR" = SigRscen, "LifeHistory" = lh_vec, "L" = Lyears, "Samp" = Lsamp, "F" = Fscen, "Rich" = Rich)

res_wide <- get_results(mod_path = sim_path, 
                        df = scen_grid, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted"),
                        res_name = paste0(res_name_inp, "_wide"),
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

write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)

true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% dplyr::select(-Rest)
resp <- read.csv(file.path(sim_path, "results", "results_variable_perfect_all_long.csv"))
ress <- read.csv(file.path(sim_path, "results", "results_variable_avgsample_all_long.csv"))
ress2 <- read.csv(file.path(sim_path, "results", "results_variable_lowsample_all_long.csv"))
res <- rbind.data.frame(resp, ress, ress2)
res2 <- res %>%
  dplyr::select(-c(variable, estimate, sd)) %>%
  pivot_longer(cols = LN_R0, names_to = "variable", values_to = "estimate") %>%
  mutate(sd = NA) %>%
  mutate(LN_R0 = estimate)
res_all <- rbind.data.frame(res, res2)
all <- left_join(res_all, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all$lifehistory <- as.character(all$lifehistory)
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

sum_selex <- sum_final %>% filter(grepl("Selex", variable))

labels <- unique(all_info$label)
labels_inp2 <- str_replace(labels, "_", "\n")
check <- all_info %>% 
  filter(converge == 1) %>% 
  filter(year == 100) %>% 
  filter(grepl("Selex",variable)) %>%
  filter(Fscen == "F1") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(SigmaR == "HighSigmaR") %>%
  mutate(label = str_replace(label, "_", "\n"))
check$label <- factor(check$label, levels = labels_inp2)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin(aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(longevity+growth~variable) +
  guides(fill = FALSE) +
  xlab("Scenario") + ylab("RE") +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Selex.png"), p, height = 10, width = 18)

sum_R0 <- sum_final %>% filter(grepl("LN_R0", variable))
labels <- unique(all_info$label)
labels_inp2 <- str_replace(labels, "_", "\n")
check <- all_info %>% 
  filter(converge == 1) %>% 
  filter(year == 100) %>% 
  filter(grepl("R0",variable)) %>%
  filter(Fscen == "F1") %>%
  filter(Rest == "R_biasadjusted") %>%
  filter(SigmaR == "HighSigmaR") %>%
  mutate(label = str_replace(label, "_", "\n"))
check$label <- factor(check$label, levels = labels_inp2)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin(aes(x = label, y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(longevity~growth) +
  guides(fill = FALSE) +
  xlab("Scenario") + ylab("RE") +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "R0.png"), p, height = 10, width = 18)


check <- all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  filter(year %in% seq(60,100,by=10)) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(label == "perfect")
p_rec <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin(aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = FALSE) +
  coord_cartesian(ylim = c(min(check$re),quantile(check$re,0.999))) +
  ylab("RE") + xlab("Year") +
  facet_grid(longevity~growth) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "RYr_perfect.png"), p_rec, height = 10, width = 18)

check <- all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  filter(year %in% seq(60,100,by=10)) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(label == "L75_N200")
p_rec <- ggplot(check) +
  geom_hline(aes(yintercept = 0)) +
  geom_violin(aes(x = factor(year), y = re, fill = lifehistory), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = FALSE) +
  coord_cartesian(ylim = c(min(check$re),quantile(check$re,0.999))) +
  ylab("RE") + xlab("Year") +
  facet_grid(longevity~growth) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "RYr_L75_N200.png"), p_rec, height = 10, width = 18)

check <-  all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  filter(year %in% 50:100) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(iteration %in% sample(1:100,1)) %>%
  filter(label == "perfect")
p_iter <- ggplot(check) +
  geom_line(aes(x = year, y = true), alpha = 0.8, lwd = 1) +
  geom_line(aes(x = year, y = estimate), alpha = 0.8, lwd = 2, color = "red") +
  facet_wrap(iteration~longevity+growth) +
  xlab("Year") + ylab("Recruitment") +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Riter1_perfect.png"), p_iter, height = 10, width = 12)

check <-  all_info %>%
  filter(converge == 1) %>%
  filter(variable == "Recruit") %>%
  filter(year %in% 50:100) %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(iteration %in% 1) %>%
  filter(label == "L75_N200")
p_iter <- ggplot(check) +
  geom_line(aes(x = year, y = true), alpha = 0.8, lwd = 1) +
  geom_line(aes(x = year, y = estimate), alpha = 0.8, lwd = 2, color = "red") +
  facet_wrap(iteration~longevity+growth) +
  xlab("Year") + ylab("Recruitment") +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Riter1_L75_N200.png"), p_iter, height = 10, width = 12)


check_mre <- lapply(1:length(itervec), function(x){
  sub <- all_info %>% filter(year == 100) %>% filter(variable %in% c("LN_R0","OFL")) %>% filter(iteration <= itervec[x])
  sum <- sub %>%
    group_by(lifehistory, SigmaR, Fscen, variable, Rest, label) %>% 
    summarise(mre = median(re)) %>%
    mutate(iteration = itervec[x])
  return(sum)
})
check_mre <- do.call(rbind, check_mre)
piter <- ggplot(check_mre %>% filter(Fscen == "F1") %>% filter(SigmaR == "HighSigmaR")) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(x = iteration, y = mre, color = variable, linetype = Rest)) + 
  facet_grid(lifehistory~label) +
  coord_cartesian(ylim = c(min(check_mre$mre), quantile(check_mre$mre, 0.99))) +
  theme_bw()



###############################
## Requests
###############################
lh_vec <- lh_info$ShortName[4]

# ## request 1: rerun 1 iter with adjusted estimation of end recruitment
s1 <- expand.grid("SigmaR" = "HighSigmaR",
                  "LifeHistory" = lh_vec,
                  "L" = "L100",
                  "Samp" = "perfect",
                  "F" = "F1",
                  "Rich" = FALSE,
                  "Fix" = NA,
                  "Adjust" = NA)

set.seed(123)
Pyears <- 30
Nyears <- 100
Fyears <- 1
recdev_mat3 <- matrix(rnorm((Pyears+Nyears+Fyears)*length(itervec), mean = 0, 0.8),ncol=length(itervec))

#Run model
rpath1 <- "~/Projects/NOAA/datamoderate/sim_new/long_fast/F1/HighSigmaR/1/perfect/L100/R_input_recdevs"
navigate <- paste0("cd ", rpath1)
os <- .Platform$OS.type
ss_bin <- "ss"
bin <- get_bin(ss_bin)
system(paste0(navigate, ";", bin), ignore.stdout = TRUE)


ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = s1,
       path = sim_path,
       itervec = itervec,
       clean = FALSE,
       rewrite = rewrite,
       run_noest = FALSE,
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start


## data rich selectivity R0 check
s0a <- expand.grid("SigmaR" = "HighSigmaR", 
                  "LifeHistory" = lh_vec, 
                  "L" = "L100", 
                  "Samp" = "perfect", 
                  "F" = "F1", 
                  "Rich" = TRUE, 
                  "Fix" = c("Selex"),
                  "Adjust" = NA)

s0b <- expand.grid("SigmaR" = "HighSigmaR", 
                   "LifeHistory" = lh_vec, 
                   "L" = "L100", 
                   "Samp" = "perfect", 
                   "F" = "F1", 
                   "Rich" = TRUE, 
                   "Fix" = c("Selex_R0"),
                   "Adjust" = NA)

s0 <- rbind.data.frame(s0a, s0b)

ncores <- 14
itervec <- 1:42
rewrite = FALSE
start <- Sys.time()
run_ss(df = s0,
       path = sim_path,
       itervec = itervec,
       clean = FALSE,
       rewrite = rewrite,
       run_noest = FALSE,
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

## fixing selectivity
s1 <- expand.grid("SigmaR" = "HighSigmaR", 
                  "LifeHistory" = lh_vec, 
                  "L" = "L100", 
                  "Samp" = "perfect", 
                  "F" = "F1", 
                  "Rich" = FALSE, 
                  "Fix" = c("Selex"),
                  "Adjust" = NA)

s2 <- expand.grid("SigmaR" = "HighSigmaR", 
                  "LifeHistory" = lh_vec, 
                  "L" = "L75", 
                  "Samp" = "N200", 
                  "F" = "F1", 
                  "Rich" = FALSE, 
                  "Fix" = c("Selex"),
                  "Adjust" = NA)

s_selex <- rbind.data.frame(s1,s2,s0a)

ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = s_selex,
       path = sim_path,
       itervec = itervec,
       clean = FALSE,
       rewrite = rewrite,
       run_noest = FALSE,
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

res_name_inp <- "Fix_selex"
res_wide <- get_results(mod_path = sim_path, 
                        df = s_selex, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_fix_Selex"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)

## fixing selectivity and R0
s3 <- expand.grid("SigmaR" = "HighSigmaR",
                  "LifeHistory" = lh_vec,
                  "L" = "L100",
                  "Samp" = "perfect",
                  "F" = "F1",
                  "Rich" = c(FALSE), 
                  "Fix" = c("Selex_R0"),
                  "Adjust" = NA)

s4 <- expand.grid("SigmaR" = "HighSigmaR",
                  "LifeHistory" = lh_vec,
                  "L" = "L75",
                  "Samp" = "N200",
                  "F" = "F1",
                  "Rich" = c(FALSE), 
                  "Fix" = c("Selex_R0"),
                  "Adjust" = NA)

s_selexR0 <- rbind.data.frame(s3, s4, s0b)

ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = s_selexR0,
       path = sim_path,
       itervec = itervec,
       clean = FALSE,
       rewrite = rewrite,
       run_noest = FALSE,
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

res_name_inp <- "Fix_selexR0"
res_wide <- get_results(mod_path = sim_path, 
                        df = s_selexR0, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_fix_Selex_R0"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp, Rich) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)


## linf low
s3b <- expand.grid("SigmaR" = "HighSigmaR",
                  "LifeHistory" = lh_vec,
                  "L" = "L100",
                  "Samp" = "perfect",
                  "F" = "F1",
                  "Rich" = c(TRUE,FALSE), 
                  "Fix" = NA,
                  "Adjust" = c("Linf_49.5"))

s4 <- expand.grid("SigmaR" = "HighSigmaR",
                  "LifeHistory" = lh_vec,
                  "L" = c("L75","L20"),
                  "Samp" = "N200",
                  "F" = "F1",
                  "Rich" = FALSE, 
                  "Fix" = NA,
                  "Adjust" = c("Linf_49.5"))

slinflow <- rbind.data.frame(s3b, s4)

ncores <- 14
itervec <- 1:100
rewrite = FALSE
start <- Sys.time()
run_ss(df = slinflow, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

sub <- slinflow
res_name_inp <- "Linf_49.5"
itervec <- 1:100
res_wide <- get_results(mod_path = sim_path, 
                        df = sub, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_change_Linf_49.5"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)


## linf low
s3b <- expand.grid("SigmaR" = "HighSigmaR",
                  "LifeHistory" = lh_vec,
                  "L" = "L100",
                  "Samp" = "perfect",
                  "F" = "F1",
                  "Rich" = c(TRUE,FALSE), 
                  "Fix" = NA,
                  "Adjust" = c("Linf_60.5"))


slinfhigh <- rbind.data.frame(s3b)

ncores <- 14
itervec <- 1:50
rewrite = FALSE
start <- Sys.time()
run_ss(df = slinflow, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

sub <- slinfhigh
res_name_inp <- "Linf_60.5"
itervec <- 1:50
res_wide <- get_results(mod_path = sim_path, 
                        df = sub, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_change_Linf_60.5"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
gc()
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)


## adjusted cv
s4a <-  expand.grid("SigmaR" = "HighSigmaR",
                    "LifeHistory" = lh_vec,
                    "L" = c("L75"),
                    "Samp" = "N200",
                    "F" = "F1",
                    "Rich" = FALSE, 
                    "Fix" = NA,
                    "Adjust" = c("CV_0.125"))

s4b <- expand.grid("SigmaR" = "HighSigmaR",
                   "LifeHistory" = lh_vec,
                   "L" = c("L75"),
                   "Samp" = "N200",
                   "F" = "F1",
                   "Rich" = FALSE, 
                   "Fix" = NA,
                   "Adjust" = c("CV_0.075"))

scen_cv <- rbind.data.frame(s4a, s4b)
scen_cv[,"Fix"] <- as.character(scen_cv[,"Fix"])
scen_cv[,"Adjust"] <- as.character(scen_cv[,"Adjust"])

ncores <- 14
itervec <- 1:56
rewrite = FALSE
start <- Sys.time()
run_ss(df = scen_cv, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores)
end_run <- Sys.time() - start

sub <- scen_cv %>% filter(Adjust == "CV_0.075")
res_name_inp <- "CV_0.075"
res_wide <- get_results(mod_path = sim_path, 
                        df = sub, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_change_CV_0.075"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)

sub <- scen_cv %>% filter(Adjust == "CV_0.125")
res_name_inp <- "CV_0.125"
res_wide <- get_results(mod_path = sim_path, 
                        df = sub, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_change_CV_0.125"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)


### dirichlet
sd <- expand.grid("SigmaR" = "HighSigmaR",
                  "LifeHistory" = "long_fast",
                  "L" = "L75",
                  "Samp" = "N50",
                  "F" = "F1",
                  "Rich" = FALSE, 
                  "Fix" = NA,
                  "Adjust" = "Dirichlet")

ncores <- 14
itervec <- 1:84
rewrite = FALSE
start <- Sys.time()
run_ss(df = sd, 
       path = sim_path, 
       itervec = itervec,
       clean = FALSE, 
       rewrite = rewrite, 
       run_noest = FALSE, 
       run_hess = FALSE,
       ncores = ncores,
       samp_type = "dirichlet")
end_run <- Sys.time() - start

sub <- scen_cv %>% filter(Adjust == "CV_0.075")
res_name_inp <- "CV_0.075"
res_wide <- get_results(mod_path = sim_path, 
                        df = sub, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_change_CV_0.075"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)

itervec <- 1:28
res_name_inp <- "Dirichlet"
res_wide <- get_results(mod_path = sim_path, 
                        df = sd, 
                        itervec = itervec, 
                        read_truth = FALSE, 
                        Rscen = c("R_biasadjusted_change_Dirichlet"),
                        res_name = paste0(res_name_inp, "_wide"),
                        ncores = 14)
res_wide <- unique(res_wide)
saveOFL <- res_wide %>% 
  group_by(lifehistory, iteration, Rest, SigmaR, Fscen, LN_R0, Lyears, Nsamp) %>%
  summarise(OFL = OFL[which(year == Nyears + 1)])

res_long_all <- res_wide %>%
  dplyr::select(-c(OFL,OFL_sd)) %>%
  left_join(saveOFL) %>%
  pivot_longer(-c(iteration, lifehistory, year, max_grad, SigmaR, Rest, Rich, Lyears, Nsamp, Fscen, LN_R0), names_to = "variable", values_to = "estimate")

res_long_mle <- res_long_all %>% filter(grepl("_sd", variable) == FALSE)
res_long_sd <- res_long_all %>% filter(grepl("_sd", variable)) %>% 
  rename(sd = estimate) %>%
  separate(col = variable, into = c("variable", "rm"), sep = "_") %>%
  dplyr::select(-rm)
res_long <- full_join(res_long_mle, res_long_sd)
write.csv(res_long, file.path(sim_path, "results", paste0(res_name_inp, "_long.csv")), row.names = FALSE)



####################
true <- read.csv(file.path(sim_path, "results", "sim_variable_long.csv")) %>% dplyr::select(-Rest)
res1 <- read.csv(file.path(sim_path, "results", "Fix_selex_long.csv"))
res2 <- read.csv(file.path(sim_path, "results", "Fix_selexR0_long.csv"))
res3 <- read.csv(file.path(sim_path, "results", "Linf_49.5_long.csv"))
res3b <- read.csv(file.path(sim_path, "results", "Linf_60.5_long.csv"))
res4 <- read.csv(file.path(sim_path, "results", "CV_0.075_long.csv"))
res5 <- read.csv(file.path(sim_path, "results", "CV_0.125_long.csv"))
res6 <- read.csv(file.path(sim_path, "results", "Dirichlet_long.csv"))
resr <- read.csv(file.path(sim_path, "results", "results_variable_rich_all_long.csv"))
resp <- read.csv(file.path(sim_path, "results", "results_variable_perfect_all_long.csv"))
ress <- read.csv(file.path(sim_path, "results", "results_variable_avgsample_all_long.csv"))
resl <- read.csv(file.path(sim_path, "results", "results_variable_lowsample_all_long.csv"))


res <- rbind.data.frame(res1, res2, res3, res3b, res4, res5, res6, resr, resp, ress, resl)
resx <- res %>%
  dplyr::select(-c(variable, estimate, sd)) %>%
  pivot_longer(cols = LN_R0, names_to = "variable", values_to = "estimate") %>%
  mutate(sd = NA) %>%
  mutate(LN_R0 = estimate)
res_all <- rbind.data.frame(res, resx)
all <- left_join(res_all, true, join_by = c(year, lifehistory, iteration, SigmaR, variable))
all$lifehistory <- as.character(all$lifehistory)
all_info <- all %>% 
  mutate(re = (estimate - true)/true) %>%
  mutate(label = ifelse(Rich == TRUE, 'data-rich',
                        ifelse(Rich == FALSE & Nsamp == "perfect", "perfect",
                               paste0(Lyears,"_", Nsamp)))) %>%
  mutate(converge = ifelse(max_grad > 1 | LN_R0 > 12, 0, 1),
         variable = replace(variable, variable == "Depletion", "Fraction of unfished")) %>%
  mutate(longevity = ifelse(grepl("short_",lifehistory), "Shorter-lived", "Longer-lived")) %>%
  mutate(growth = ifelse(grepl("_slow", lifehistory), "Slower-to-Linf", "Faster-to-Linf"))

all_info$Rest <- as.character(all_info$Rest)
all_info2 <- all_info %>%
  mutate(Rest = replace(Rest, Rest == "R_biasadjusted", "base")) %>%
  mutate(Rest = str_remove(Rest, "R_biasadjusted_")) %>%
  mutate(Rest = str_replace(Rest, "_", "\n")) %>%
  mutate(Rest = str_remove(Rest, "change\n"))

find_selex <- all_info2 %>% filter(variable == "Selex_peak") %>%
  dplyr::select(-c(variable,sd,true,re,converge)) %>%
  rename(Selex_peak = estimate)

all_info3 <- left_join(all_info2, find_selex) %>%
  mutate(converge = replace(converge, Selex_peak > 80, 0))


labels <- c("data-rich", "perfect", "L75_N200")
check <- all_info3 %>%
  filter(converge == 1) %>%
  filter(Rest %in% c("base", "fix\nSelex", "fix\nSelex_R0")) %>%
  filter(variable %in% c("OFL","Fraction of unfished", "Recruit", "LN_R0")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR")
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = Rest, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(variable~label, scales = "free_y", ncol = 3) +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Fix_selectivity.png"), p, height = 14, width = 18)

tab <- check %>% 
  filter(variable %in% c("Recruit","OFL","LN_R0")) %>%
  group_by(Rest, variable, lifehistory, SigmaR, label) %>%
  summarise(mre = median(re))

labels <- c("data-rich","perfect")
check <- all_info3 %>%
  filter(converge == 1) %>%
  filter(Rest %in% c("base", "Linf_49.5", "Linf_60.5")) %>%
  filter(variable %in% c("OFL","Fraction of unfished")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR")
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = label, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable+Rest, ncol = 3) +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Adjust_Linf_RE_rich.png"), p, height = 8, width = 10)


check_mre <- lapply(1:length(itervec), function(x){
  sub <- check %>% filter(year == 100) %>% filter(variable %in% c("LN_R0","Selex_peak","OFL")) %>% filter(iteration <= itervec[x])
  sum <- sub %>%
    group_by(lifehistory, SigmaR, Fscen, variable, Rest, label) %>% 
    summarise(mre = median(re)) %>%
    mutate(iteration = itervec[x])
  return(sum)
})
check_mre <- do.call(rbind, check_mre)
piter <- ggplot(check_mre %>% filter(Fscen == "F1") %>% filter(SigmaR == "HighSigmaR")) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(aes(x = iteration, y = mre, color = variable, linetype = Rest)) + 
  facet_grid(lifehistory~label) +
  coord_cartesian(ylim = c(min(check_mre$mre), quantile(check_mre$mre, 0.99))) +
  theme_bw()
  
  
labels <- c("data-rich","perfect")
check <- all_info3 %>%
  filter(converge == 1) %>%
  filter(Rest %in% c("base", "Linf_49.5", "Linf_60.5")) %>%
  filter(variable %in% c("Selex_peak","LN_R0")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR")
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = label, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable+Rest, ncol = 3, scales = "free_y") +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Adjust_Linf_RE_rich_selex.png"), p, height = 8, width = 10)

tab <- check %>% 
  group_by(Rest, variable, lifehistory, SigmaR, label) %>%
  summarise(iter = max(iteration))




labels <- c("data-rich","perfect")
check <- all_info3 %>%
  # filter(converge == 1) %>%
  filter(Rest %in% c("base", "Linf_49.5", "Linf_60.5")) %>%
  filter(variable %in% c("Selex_peak","LN_R0")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR")
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = label, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable+Rest, ncol = 3, scales = "free_y") +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Adjust_Linf_RE_rich_selex_nonconverge.png"), p, height = 8, width = 10)


# labels <- c("data-rich","L75_N200")
# check <- all_info3 %>%
#   filter(converge == 1) %>%
#   filter(Rest %in% c("base", "Linf_49.5")) %>%
#   filter(variable %in% c("OFL","Fraction of unfished","Selex_peak","Selex_ascend")) %>%
#   filter(year == 100) %>%
#   filter(label %in% labels) %>%
#   filter(lifehistory == "long_fast") %>%
#   filter(SigmaR == "HighSigmaR")
# check$label <- factor(check$label, levels = labels)
# p <- ggplot(check) +
#   geom_hline(aes(yintercept = 0), lwd = 1.2) +
#   geom_violin(aes(x = label, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
#   scale_fill_brewer(palette = "Set1") +
#   facet_wrap(~variable+Rest, scales = "free_y", ncol = 2) +
#   xlab("Scenario") + ylab("RE") +
#   guides(fill = FALSE) +
#   theme_bw(base_size = 20)
# ggsave(file.path(fig_path, "Adjust_Linf_RE.png"), p, height = 12, width = 10)


labels <- c("L75_N200")
check <- all_info3 %>%
  filter(converge == 1) %>%
  filter(Rest %in% c("base", "Linf_49.5","CV_0.075", "CV_0.125")) %>%
  filter(variable %in% c("OFL","Fraction of unfished")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR")# %>%
  # filter(iteration %in% 1:56)
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = Rest, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Adjust_Linf_CV_RE.png"), p, height = 8, width = 12)

labels <- c("L75_N200")
check <- all_info3 %>%
  filter(converge == 1) %>%
  filter(Rest %in% c("base", "Linf_49.5","CV_0.075", "CV_0.125")) %>%
  filter(variable %in% c("Selex_peak","LN_R0")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR")# %>%
  # filter(iteration %in% 1:56)
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = Rest, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Adjust_Linf_CV_RE_selex.png"), p, height = 8, width = 12)

labels <- c("L75_N200")
check <- all_info3 %>%
  # filter(converge == 1) %>%
  filter(Rest %in% c("base", "Linf_49.5","CV_0.075", "CV_0.125")) %>%
  filter(variable %in% c("Selex_peak","LN_R0")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR")# %>%
  # filter(iteration %in% 1:56)
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = Rest, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)
ggsave(file.path(fig_path, "Adjust_Linf_CV_RE_selex_nonconverge.png"), p, height = 8, width = 12)




labels <- c("L75_N50")
check <- all_info3 %>%
  filter(converge == 1) %>%
  filter(Rest %in% c("base", "Dirichlet")) %>%
  filter(variable %in% c("OFL","Fraction of unfished")) %>%
  filter(year == 100) %>%
  filter(label %in% labels) %>%
  filter(lifehistory == "long_fast") %>%
  filter(SigmaR == "HighSigmaR") %>%
  filter(iteration %in% 1:28)
check$label <- factor(check$label, levels = labels)
p <- ggplot(check) +
  geom_hline(aes(yintercept = 0), lwd = 1.2) +
  geom_violin(aes(x = Rest, y = re, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  xlab("Scenario") + ylab("RE") +
  guides(fill = FALSE) +
  theme_bw(base_size = 20)


# p <- ggplot(check) +
#   # geom_hline(aes(yintercept = 0), lwd = 1.2) +
#   geom_violin(aes(x = label, y = estimate, fill = Rest), scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
#   scale_fill_brewer(palette = "Set1") +
#   facet_wrap(~variable+Rest, scales = "free_y", ncol = 2) +
#   xlab("Scenario") + ylab("RE") +
#   guides(fill = FALSE) +
#   theme_bw(base_size = 20)
# ggsave(file.path(fig_path, "Adjust_Linf_Value.png"), p, height = 8, width = 18)

### request 5 -- plot china rockfish
path <- file.path(sim_path, "compare", "china_run")
rc <- SS_output(path, verbose = FALSE)
SS_plots(dir = path, replist = rc)


path <- file.path(sim_path, "short_slow", "F1", "HighSigmaR", 1, "N50", "L10", "R_biasadjusted")
r1 <- SS_output(path, verbose = FALSE)
SS_plots(dir = path, replist = r1)

path <- file.path(sim_path, "short_fast", "F1", "HighSigmaR", 1, "N50", "L10", "R_biasadjusted")
r1 <- SS_output(path, verbose = FALSE)
SS_plots(dir = path, replist = r1)


path <- file.path(sim_path, "long_fast", "F1", "HighSigmaR", 1, "N50", "L10", "R_biasadjusted")
r1 <- SS_output(path, verbose = FALSE)
SS_plots(dir = path, replist = r1)

path <- file.path(sim_path, "long_slow", "F1", "HighSigmaR", 1, "N50", "L10", "R_biasadjusted")
r1 <- SS_output(path, verbose = FALSE)
SS_plots(dir = path, replist = r1)






### testing dirichlet
om_path <- file.path(sim_path, "long_fast", "F1", "HighSigmaR", 1, "om")
samp_path <- file.path(sim_path, "long_fast", "F1", "HighSigmaR", 1, "N200")
      
      ##################################
      ## setup data file by iteration
      SS_splitdat2(inpath = om_path, outpath = samp_path, inname = "data.ss_new", outpattern = "perfect", number = FALSE, verbose = TRUE, fillblank = TRUE, MLE = TRUE)
      dat <- SS_readdat(file.path(samp_path, "perfect.ss"), version = "3.30", verbose = FALSE)

        dat <- change_data(dat_list = dat, outfile = NULL, years = 1, fleets = 1, types = "index")
        dat$CPUE[1,"se_log"] <- 1
        dat <- sample_agecomp(dat_list = dat, fleets = NULL, Nsamp = NULL, years = NULL, cpar = NULL, ESS = NULL)


      new_comp <- dat$lencomp
      end <- 6
      
      dat <- change_data(dat_list = dat, outfile = NULL, years = 1:100, fleets = 1, types = "len", len_bins = seq(10,82, by = 2), pop_binwidth = 1, pop_minimum_size = 10, pop_maximum_size = 82)
      dat$lencomp[,1:ncol(new_comp)] <- new_comp
      dat$lencomp[,(ncol(new_comp)+1):ncol(dat$lencomp)] <- 0
      
      new_comp <- dat$lencomp
      temp_comp <- dat$lencomp
      

      library(MCMCpack)
        set.seed(456)
        for(y in 1:Nyears){
          year_temp <- new_comp$Yr == y
          if(sum(year_temp)!=0){
            probs <- new_comp[year_temp, -(1:end)]
            probs[which(probs == 0)] <- 0.00001
            probs <- as.numeric(probs / sum(probs))
            Nsamp <- 200
            out_comp <- rdirichlet(1, alpha =  Nsamp * probs)
            temp_comp[year_temp,(end+1):dim(new_comp)[2]] <- out_comp
            temp_comp[,"Nsamp"] <- Nsamp
          }
        }
        test_comp1 <- cbind(temp_comp[,1:end], temp_comp[,(end+1):dim(new_comp)[2]])


      temp_comp <- dat$lencomp

        set.seed(456)
        for(y in 1:Nyears){
          year_temp <- new_comp$Yr == y
          if(sum(year_temp)!=0){
            probs <- new_comp[year_temp, -(1:end)]
            probs[which(probs == 0)] <- 0.00001
            probs <- as.numeric(probs / sum(probs))
            Nsamp <- 200

            out_comp <- t(rmultinom(1, size = Nsamp, prob = probs))
            out_comp <- out_comp / sum(out_comp)
            temp_comp[year_temp,(end+1):dim(new_comp)[2]] <- out_comp
            temp_comp[,"Nsamp"] <- Nsamp
          }
        }
        test_comp2 <- cbind(temp_comp[,1:end], temp_comp[,(end+1):dim(new_comp)[2]])

p1 <- test_comp1 %>%
  pivot_longer(cols = l10:l82, names_to = "bin", values_to = "value") %>%
  mutate(Distribution = "Dirichlet")

p2 <- test_comp2 %>%
 pivot_longer(cols = l10:l82, names_to = "bin", values_to = "value") %>%
 mutate(Distribution = "Multinomial")

pdf <- full_join(p1, p2)

p <- ggplot(pdf %>% filter(Yr %in% 96:100)) +
geom_line(aes(x = bin, y = value)) +
facet_wrap(Distribution~Yr) +
theme_bw()