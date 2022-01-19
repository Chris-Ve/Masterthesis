library(tidyverse)
library(tidymodels)
library(lmtest)
library(sandwich)
library(stargazer)

df <- readxl::read_excel("uni_data_master.xlsx")

### Transform data for simplicity

uni_data <- recipe(~ ., data = df) %>%
  step_log(avg_fac_sal, applicants, admissions, enrolled, tuition, instr_exp_per_stud) %>% 
  step_mutate_at(applicants, admissions, enrolled, admrate,yield,tuition,
                 fn = ~ .*100) %>% 
  step_rename(grad = grad_rate, ret=ret_rate, sal=avg_fac_sal, ft=fac_fulltime,
              instr = instr_exp_per_stud, sfr=sf_ratio) %>% 
  step_lag(grad, ret, sal, ft, instr, sfr, 
           lag = 2*118) %>%
  step_rename_at(starts_with("lag"), fn = ~ gsub("lag_236_", "l", ., fixed = TRUE)) %>% 
  step_naomit(rank) %>% 
  prep(df) %>%
  bake(df)


# Since the rank is published one year ahead, it is already lagged by one year by default.
# Therefore, taking an additional lag would be redundant.

##### Baseline Estimation Results #######################

# for correct ordering of results
myorder = c("sat" ,"applicants", "admissions", "enrolled", "admrate", "yield", "tuition")

# estimation
mods.baseline <- uni_data %>%
  gather(variable, value, sat, applicants, admissions, enrolled, admrate, yield, tuition) %>%
  group_by(variable) %>%
  do(mod1 = lm(value ~ rank+lgrad+lret+linstr, .),
     mod2 = lm(value ~ rank+lgrad+lret+linstr+lsal+lft+lsfr, .),
     mod3 = lm(value ~ rank+lgrad+lret+linstr+year, .),
     mod4 = lm(value ~ rank+lgrad+lret+linstr+lsal+lft+lsfr+year, .),
     mod5 = lm(value ~ rank+name+year, .),
     mod6 = lm(value ~ rank+lgrad+lret+linstr+name+year, .),
     mod7 = lm(value ~ rank+lgrad+lret+linstr+lsal+lft+lsfr+name+year, .)) %>% 
  mutate(variable =  factor(variable, levels = myorder)) %>%
  arrange(variable)

# robust standard errors
mods.baseline.robust <- lapply(t(mods.baseline)[-1,],
                                  function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# this reproduces the results for the SAT score in table 6 for example
stargazer(mods.baseline.robust[1:7], keep = "rank", type="text")

# for the yield rate in table 6
stargazer(mods.baseline.robust[36:42], keep = "rank", type="text")

# R^2 can be obtained by using the model with uncorrected SEs, e.g.
stargazer(as.matrix(mods.baseline)[1,-1], keep = "rank", keep.stat = c("n", "adj.rsq"),
          type="text")


##### Flexible Timing Results #######################

mods.timing <- uni_data %>%
  gather(variable, value, sat, applicants, admissions, enrolled, admrate, yield, tuition) %>%
  group_by(variable) %>%
  do(mod1 = lm(value ~ rank+lag(rank,118)+lgrad+lret+linstr+name+year, .),
     mod2 = lm(value ~ rank+lag(rank,118)+lag(rank,2*118)+lgrad+lret+linstr+name+year, .)) %>% 
  mutate(variable =  factor(variable, levels = myorder)) %>%
  arrange(variable)

# robust standard errors
mods.timing.robust <- lapply(t(mods.timing)[-1,],
                               function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 7
stargazer(mods.timing.robust[1:8], keep = "rank", type="text")

# table 8
stargazer(mods.timing.robust[9:14], keep = "rank", type="text")



##### Heterogeneity Results #######################

### by control
mods.control <- uni_data %>%
  gather(variable, value, sat, applicants, admissions, enrolled, admrate, yield, tuition) %>%
  group_by(variable) %>%
  do(mod1 = lm(value ~ rank+rank:private+lgrad+lret+linstr+name+year, .)) %>% 
  mutate(variable =  factor(variable, levels = myorder)) %>%
  arrange(variable)

# robust standard errors
mods.control.robust <- lapply(t(mods.control)[-1,],
                             function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 9, panel A
stargazer(mods.control.robust, keep = "rank", type="text")


### by size

uni_data$size_ug_bin <- fct_relevel(uni_data$size_ug_bin, "small", "med", "large")

# estimation
mods.size <- uni_data %>%
  gather(variable, value, sat, applicants, admissions, enrolled, admrate, yield, tuition) %>%
  group_by(variable) %>%
  do(mod1 = lm(value ~ rank+rank:size_ug_bin+lgrad+lret+linstr+name+year, .)) %>% 
  mutate(variable =  factor(variable, levels = myorder)) %>%
  arrange(variable)

# robust standard errors
mods.size.robust <- lapply(t(mods.size)[-1,],
                              function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 9, panel B
stargazer(mods.size.robust, keep = "rank", type="text")



### by relative rank position

uni_data$rank_bin <- fct_relevel(uni_data$rank_bin, "25", "50", "75", "100", "125", "150")

# estimation
mods.rank <- uni_data %>%
  gather(variable, value, sat, applicants, admissions, enrolled, admrate, yield, tuition) %>%
  group_by(variable) %>%
  do(mod1 = lm(value ~ rank+rank:rank_bin+lgrad+lret+linstr+name+year, .)) %>% 
  mutate(variable =  factor(variable, levels = myorder)) %>%
  arrange(variable)

# robust standard errors
mods.rank.robust <- lapply(t(mods.rank)[-1,],
                           function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 9, panel C
stargazer(mods.rank.robust, keep = "rank", type="text")


##### Sensitivity Results ###################

# Impuation for missing year observations of the student faculty ratio
uni_data.imp <- uni_data %>%
  mutate(lsfr.imp = lsfr) %>% 
  group_by(name) %>%
  fill(lsfr.imp, .direction = "downup") %>%
  dplyr::ungroup()

# estimation
mods.sens <- uni_data.imp %>%
  gather(variable, value, sat, applicants, admissions, enrolled, admrate, yield, tuition) %>%
  group_by(variable) %>%
  do(mod1 = lm(value ~ rank+lgrad+lret+linstr+lsal+lft+lsfr.imp+name+year, .)) %>% 
  mutate(variable =  factor(variable, levels = myorder)) %>%
  arrange(variable)

# robust standard errors
mods.sens.robust <- lapply(t(mods.sens)[-1,],
                           function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table A1
stargazer(mods.sens.robust, keep = "rank", type="text")


##### Control Set Selection Results ##############################

# results are not exactly reproducible due to the random split of the data (!)

uni_data2 <- recipe(~ ., data = df) %>%
  step_rename(grad = grad_rate, ret=ret_rate, sal=avg_fac_sal, ft=fac_fulltime,
              instr = instr_exp_per_stud, sfr=sf_ratio) %>% 
  step_lag(grad, ret, sal, ft, instr, sfr, 
           lag = 2*118) %>%
  step_rename_at(starts_with("lag"), fn = ~ gsub("lag_236_", "l", ., fixed = TRUE)) %>% 
  step_mutate_at(lsal, linstr, fn = ~ ./1000) %>% 
  step_naomit(rank) %>% 
  prep(df) %>%
  bake(df)


data_split <- initial_split(uni_data2, prop = 0.75)
train_data <- training(data_split)
test_data  <- testing(data_split)

# estimation 
mods.css <- train_data %>%
  do(mod1 = lm(rank~lgrad, .),
     mod2 = lm(rank~lgrad+lret, .),
     mod3 = lm(rank~lgrad+lret+linstr, .),
     mod4 = lm(rank~lgrad+lret+linstr+lsal, .),
     mod5 = lm(rank~lgrad+lret+linstr+lsal+lsfr, .),
     mod6 = lm(rank~lgrad+lret+linstr+lsal+lsfr+lft, .))

# robust standard errors
mods.css.robust <- lapply(t(mods.css),
                           function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 5
stargazer(mods.css.robust, type="text", omit = "Constant")


##### RDD Results ##################################

# data prep
uni_data_rdd <- recipe(~ ., data = df) %>%
  step_rename(grad = grad_rate, ret=ret_rate, sal=avg_fac_sal, ft=fac_fulltime,
              instr = instr_exp_per_stud, sfr=sf_ratio) %>% 
  step_mutate(top50=ifelse(rank<=50,TRUE,FALSE), rankc = rank-50) %>% 
  step_log(instr) %>% 
  step_mutate_at(admrate,yield,
                 fn = ~ .*100) %>% 
  step_mutate_at(applicants, admissions, enrolled, tuition,
                 fn = ~ ./1000) %>% 
  step_lag(grad,sfr,instr,ret,ft,sal, lag = 2*118) %>%
  step_rename_at(starts_with("lag"), fn = ~ gsub("lag_236_", "l", ., fixed = TRUE)) %>% 
  step_naomit(rank) %>% 
  prep(df) %>%
  bake(df)


# estimation
mods.rdd <- uni_data_rdd %>%
  gather(variable, value, sat, applicants, admissions, enrolled, admrate, yield) %>%
  group_by(variable) %>%
  do(mod1 = lm(value ~ top50*rankc+top50*I(rankc^2)+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -50 & rankc <50),
     mod2 = lm(value ~ top50*rankc+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -50 & rankc <50),
     mod3 = lm(value ~ top50*rankc+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -20 & rankc <20),
     mod4 = lm(value ~ top50*rankc+lgrad+lret+linstr+name+year, .,
               subset = rankc > -15 & rankc <15),
     mod5 = lm(value ~ top50*rankc+lgrad+lret+linstr+name+year, .,
               subset = rankc > -10 & rankc <10)) %>% 
  mutate(variable =  factor(variable, levels = myorder)) %>%
  arrange(variable)


# robust standard errors
mods.rdd.robust <- lapply(t(mods.rdd)[-1,],
                          function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 10, SAT score
stargazer(mods.rdd.robust[1:5], type="text", keep = "top50")

# table 10, Applications
stargazer(mods.rdd.robust[6:10], type="text", keep = "top50")


### Private & Public Tuition ########################

#### public tuition
mods.rdd.public <- uni_data_rdd %>%
  filter(private==0) %>% 
  do(mod1 = lm(tuition ~ top50*rankc+top50*I(rankc^2)+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -50 & rankc <50),
     mod2 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -50 & rankc <50),
     mod3 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -20 & rankc <20),
     mod4 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, .,
               subset = rankc > -15 & rankc <15),
     mod5 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, .,
               subset = rankc > -10 & rankc <10))

# robust standard errors
mods.rdd.public.robust <- lapply(t(mods.rdd.public),
                                  function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 10, panel G
stargazer(mods.rdd.public.robust, type="text", keep = "top50")



### private tuition
mods.rdd.private <- uni_data_rdd %>%
  filter(private==1) %>% 
  do(mod1 = lm(tuition ~ top50*rankc+top50*I(rankc^2)+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -50 & rankc <50),
     mod2 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -50 & rankc <50),
     mod3 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, ., 
               subset = rankc > -20 & rankc <20),
     mod4 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, .,
               subset = rankc > -15 & rankc <15),
     mod5 = lm(tuition ~ top50*rankc+lgrad+lret+linstr+name+year, .,
               subset = rankc > -10 & rankc <10))

# robust standard errors
mods.rdd.private.robust <- lapply(t(mods.rdd.private),
                          function(x) coeftest(x, vcov=vcovCL(x, type = "HC1", cluster = ~name)))

# table 10, panel H
stargazer(mods.rdd.private.robust, type="text", keep = "top50")








