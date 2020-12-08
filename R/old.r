###### INITIAL STUFF ==========
setwd("E:/Bolnick lab/Thesis/chapter 2 copepod local adapt")
copepods <- read.csv("chapter2.copepods.cleaned.csv")
#copepodsurvival <- read.csv("chapter2.copepodsurvival.cleaned.csv")

#making a new column to join copepods to its lake (i.e. tapeworm came from same copepod lake true or false?)
copepods$native = as.character(copepods$cop.lake) == 
  as.character(copepods$worm.lake)
# View(copepods)
library("plyr") #if not loaded before dplyr, it'll cause porblems 
# to dplyr. I need this for the paakage Rmisc to get means and CI 
# in a dataset accoding to different groupings
library(tidyverse) # has dplyr, ggplot2
library(cowplot)
library("ggpubr")
library("brms")
library(loo)
# Use this data for subsets
prediction_data1 = copepods %>% select(cop.lake, worm.lake, native, plate, worm.fam) %>% distinct


# Hurdle Poisson model
# There is some hurdle that the data-generating process needs to overcome 
# Before you get non-zero data
# in this case, that hurdle is successful infection

# Here's the first model to try
# model_formula1 = bf(numb.worm ~ cop.lake * worm.lake + native,
# hu ~ cop.lake * worm.lake + native)
# Note that hu is the probability of infection, 
# which is part of a logistic regression sub-model

# Make a list of formulae to be run later and compared using LOOIC
#list(
#bf(numb.worm ~ cop.lake * worm.lake + (1|plate) + native,
#hu ~ cop.lake * worm.lake + (1|plate) + native),
#bf(numb.worm ~ cop.lake + worm.lake + (1|plate) + native,
#hu ~ cop.lake * worm.lake + (1|plate) + native),
#bf(numb.worm ~ cop.lake * worm.lake + (1|plate) + native,
#hu ~ cop.lake + worm.lake + (1|plate) + native))


# Prior distributions
# First, use the model specification to figure out which priors need to be set
# View(get_prior(model_formula1, family = hurdle_poisson(), data = copepods))

# Run this function to define your hurdle models =======
run_hurdle_model = function(formula, priors, name, data = copepods, cores = 4,
                            save_dir = "E:/Bolnick lab/Thesis/chapter 2 copepod local adapt/hurdle_models", 
                            pred_dat = prediction_data1 , ...) {
  browser()
  model_fit = brm(formula,
                  family = hurdle_poisson(), data = data, cores = cores,
                  prior = priors, sample_prior = "yes", ...)
  
  saveRDS(model_fit, file = file.path(save_dir, paste0(name, "_model.rds")))
  
  
  hurdle_pred = posterior_linpred(model_fit, newdata = pred_dat, transform = TRUE,dpar = "hu")
  pp_hurdle = 
    pred_dat %>% mutate(key = paste0("V", 1:n())) %>% 
    left_join(hurdle_pred %>% as_tibble %>% gather(), by = "key") %>% 
    group_by(cop.lake, worm.lake) %>% mutate(value = 1-value) %>% # value was non-infection rate
    summarize(mean = mean(value), sd = sd(value),
              L_95 = quantile(value, .025), L_50 = quantile(value, .25), 
              U_95 = quantile(value, .975), U_50 = quantile(value, .75)) %>% 
    mutate(type = "hurdle")
  
  pois_pred = posterior_linpred(model_fit, newdata = pred_dat, transform = TRUE, dpar = "mu")
  pp_pois = 
    pred_dat %>% mutate(key = paste0("V", 1:n())) %>% 
    left_join(pois_pred %>% as_tibble %>% gather(), by = "key") %>% 
    group_by(cop.lake, worm.lake) %>% mutate(value = value) %>% # value was non-infection rate
    summarize(mean = mean(value), sd = sd(value),
              L_95 = quantile(value, .025), L_50 = quantile(value, .25), 
              U_95 = quantile(value, .975), U_50 = quantile(value, .75)) %>% 
    mutate(type = "pois")
  bind_rows(pp_pois, pp_hurdle) %>% mutate(name = name) %>% 
    write_csv(file.path(save_dir, paste0(name, "smry.csv")))
  message(name)
  invisible()
}

### Don't run from here until next heading ============
prior_strings1 = 
  c(prior_string("normal(0, 6)", class = "Intercept"),
    prior_string("normal(0, 1.5)", class = "b") #try SD priors from 1.5 to 10
    ,prior_string("normal(0, 1)", class = "b", dpar = "hu")#,
    #prior_string("student_t(7, 0, 3)", class = "sd"), #this and the one code below are for the random variable (1|plate)
    #prior_string("student_t(7, 0, 3)", class = "sd", dpar = "hu")
  )



#####start running several models here 


# make sure that the save_dir directory actually exists

formula_list = list(
  name1 = bf(numb.worm ~ cop.lake * worm.lake + native,
             hu ~ cop.lake * worm.lake + native), 
  name2 = bf(numb.worm ~ cop.lake * worm.lake + native,
             hu ~ cop.lake * worm.lake + native),
  name3 = bf(numb.worm ~ cop.lake * worm.lake + native,
             hu ~ cop.lake * worm.lake + native),
  name4 = bf(numb.worm ~ cop.lake * worm.lake + native,
             hu ~ cop.lake * worm.lake + native),
  name5 = bf(numb.worm ~ cop.lake * worm.lake + native,
             hu ~ cop.lake * worm.lake + native))

prior_list = list(
  prior_strings1 = c(prior_string("normal(0, 6)", class = "Intercept") #the prior in the intercept always has to be bigger than any other priors below 
                     ,prior_string("normal(0, 1.5)", class = "b") #try SD priors from 1.5 to 10
                     ,prior_string("normal(0, 1)", class = "b", dpar = "hu")),
  prior_strings2 = c(prior_string("normal(0, 6)", class = "Intercept") #the prior in the intercept always has to be bigger than any other priors below
                     ,prior_string("normal(0, 3)", class = "b") #try SD priors from 1.5 to 10
                     ,prior_string("normal(0, 1)", class = "b", dpar = "hu")),
  prior_strings3 = c(prior_string("normal(0, 6)", class = "Intercept") #the prior in the intercept always has to be bigger than any other priors below
                     ,prior_string("normal(0, 5)", class = "b") #try SD priors from 1.5 to 10
                     ,prior_string("normal(0, 1)", class = "b", dpar = "hu")),
  prior_strings4 = c(prior_string("normal(0, 8)", class = "Intercept") #the prior in the intercept always has to be bigger than any other priors below
                     ,prior_string("normal(0, 7.5)", class = "b") #try SD priors from 1.5 to 10
                     ,prior_string("normal(0, 1)", class = "b", dpar = "hu")),
  prior_strings5 = c(prior_string("normal(0, 11)", class = "Intercept") #the prior in the intercept always has to be bigger than any other priors below
                     ,prior_string("normal(0, 10)", class = "b") #try SD priors from 1.5 to 10
                     ,prior_string("normal(0, 1)", class = "b", dpar = "hu"))
)

length()

library(purrr)
list(formula = formula_list, priors = prior_list, name = names(formula_list)) %>% 
  pmap(run_hurdle_model )
# readRDS
#loaded_file = readRDS("something.rds") #change names here according to the file names recorded after running functions

#after running:

out_file = "/Volumes/SP UFD U3/Bolnick lab/Thesis/chapter 2 copepod local adapt/hurdle_models/testing for priors/"

dir.exists(out_file)
brms_fits = map(dir(out_file, pattern = "rds"), ~readRDS(file.path(out_file, .x)))

brms_summary = map_dfr(dir(out_file, pattern = "csv"), ~read_csv(file.path(out_file, .x)) %>% mutate(name = .x))

library(cowplot)
brms_summary %>% filter(type == "hurdle") %>% 
  ggplot(aes(x = name, y = mean, color = name)) + 
  facet_grid(cop.lake ~ worm.lake, switch = "both") + 
  geom_linerange(aes(ymin = L_95, ymax = U_95), size = .9) + 
  geom_linerange(aes(ymin = L_50, ymax = U_50), size = 2) + 
  geom_point(aes(y = mean), size = 5, shape = 4, color = "black")  + 
  theme(axis.text.x = element_blank()) 


#### Explore all variable combinations ================

# Define model terms
terms = c("cop.lake", "worm.lake", "worm.lake:cop.lake", "native", "(1|plate)", "(1|worm.fam)")
# make sure you don't include an interaction w/o both main effects
fix_interaction = function(arg) {
  if("worm.lake:cop.lake" %in% arg) {
    if(! (("cop.lake" %in% arg) && ("worm.lake" %in% arg)))
      arg = ""
  }
  arg
}
remove_null = function(.x) .x [.x != "~ "] # discard removed parts
# Create a combination of all model terms on rhs
rhs_combos = map(1:5, ~combn(seq_along(terms), m = .x, simplify = FALSE) %>% map(~terms[.x]) ) %>% 
  flatten %>% map(fix_interaction) %>% map_chr(~paste("~", paste(.x, collapse = " + "))) %>% 
  remove_null
# Since this is a two part model, we need all 2-part combos of the above
model_formulas_part = 
  list(pois = paste("numb.worm", rhs_combos),
       hu = paste("hu", rhs_combos)) %>% 
  cross() 
# Convert them into brms formulae
model_formulas = model_formulas_part %>% 
  map(~bf(as.formula(.x$pois), as.formula(.x$hu)))

# Determine if randome effect of plate is present; if so, add the prior to the appropriate formula section
model_priors = model_formulas_part %>% map(function(.x){
  re_hurdle = grepl("plate", x = .x$hu, fixed = TRUE) | grepl("worm.fam", x = .x$hu, fixed = TRUE)
  re_pois = grepl("plate", x = .x$pois, fixed = TRUE) | grepl("worm.fam", x = .x$pois, fixed = TRUE)
  b_hurdle = grepl("lake", x = .x$hu, fixed = TRUE)|grepl("native", x = .x$hu, fixed = TRUE)
  b_pois = grepl("lake", x = .x$pois, fixed = TRUE)|grepl("native", x = .x$pois, fixed = TRUE)
  
  prior = prior_string("normal(0, 6)", class = "Intercept") #the prior in the intercept always has to be bigger than any other priors below 
  if(isTRUE(re_hurdle)) prior = prior + 
    prior_string("student_t(7, 0, 1)", class = "sd", dpar = "hu")
  if(isTRUE(re_pois)) prior = prior +
    prior_string("student_t(7, 0, 1.5)", class = "sd")
  if(isTRUE(b_hurdle)) prior = prior+ 
    prior_string("normal(0, 1)", class = "b", dpar = "hu")
  if(isTRUE(b_pois)) prior = prior + 
    prior_string("normal(0, 1.5)", class = "b")
  prior
})

# Run here ==================
if(!exists("run_hurdle_model")) stop ("DUDE, Run the code above defining run_hurdle model")
# Change to where you want to save it
save_direc = "E:/Bolnick lab/Thesis/chapter 2 copepod local adapt/hurdle_models/combination_run"
dir.create(save_direc)

# Format for the run_hurdle_model function
parm_list = 
  list(formula = model_formulas, priors = model_priors, 
       name = paste0("combo_run_",seq_along(model_formulas)))
# Run it!
parm_list %>%  as_tibble %>% #.[-(1:934),] %>%  # Comment in and change numbers to run/exclude specific cases   
  pmap(run_hurdle_model, save_dir = save_direc, control = list(adapt_delta = .95))



# Check which files had divergences ###############
rds_files = dir(save_direc, pattern = "model.rds", full.names = TRUE)
(n_files = length(rds_files))

divergences = imap_dbl(rds_files, ~readRDS(.x) %>%
                         nuts_params(pars = "divergent__") %>% pull(Value) %>% sum)
repeat_runs = 
  tibble(rds_files = basename(rds_files), divergences) %>% filter(divergences>0) %>% 
  mutate(model_num = str_remove(rds_files, fixed("combo_run_")) %>% 
           str_remove(fixed("_model.rds")) %>% as.integer)
# re-run w/ higher adapt_delta
parm_list %>%  as_tibble %>% .[repeat_runs$model_num,] %>%  # Comment in and change numbers to run/exclude specific cases   
  pmap(run_hurdle_model, save_dir = save_direc, control = list(adapt_delta = .99))

## Calculate LOO-IC for each model #################
library(loo)
N_CORES = 4 # double check
looic_smry = map_dfr(rds_files, function(.x) {
  fit = readRDS(.x) 
  looic = brms::loo(fit, cores = N_CORES)
  # Look to see if there's a K problem?
  k_problems = length(pareto_k_ids(looic, threshold = .07))
  if(isTRUE(k_problems > 0)) {
    # re-fit the model to handle the bad k values
    looic = brms::loo(fit, k_threshold = 0.7, cores = N_CORES)
    k_problems = length(pareto_k_ids(looic, threshold = .07))
  }
  out_file = str_replace(.x, fixed("_model"), fixed("_loo"))
  saveRDS(looic, file = out_file)
  
  tibble(file = out_file, k_problems = k_problems)
})

# Check to see if any of them have k-problems
# If not...
loo_weights = map(looic_smry$out_file, readRDS) %>% 
  loo_model_weights(cores = 4)
saveRDS(loo_weights, "loo_weights.rds")

# maybe need to format this as a tibble or something
loo_weights # Rope in the formulas for the largest ones, etc.

done_loos = dir("hurdle_models/combination_run/", pattern = "loo.rds", full.names = TRUE)

loo_weights = map(done_loos, readRDS) %>% loo_model_weights(cores=4)
saveRDS(loo_weights, "loo_weights.rds")
loo_weights_pretty = loo_weights %>% as.matrix %>% as_tibble(rownames = "model_ID") %>% rename(LOO_Weight = V1) %>% arrange(desc(LOO_Weight)) 
loo_weights_pretty
loo_weights_pretty %>% ggplot(aes(x = seq_along(model_ID), y = LOO_Weight))+geom_line() + xlim(1, 5)

# So it seems like model 1444 was by far the best at prediction (Note: make sure that's 1444 by the original count, not alphabetically)
# What's the looic value of it compared with the next closest, model 1108?
done_loos[1108]
model1444_loo = readRDS("hurdle_models/combination_run/combo_run_1444_loo.rds")
model999_loo = readRDS("hurdle_models/combination_run/combo_run_999_loo.rds") # This was the last one by position
model1108_loo = readRDS("hurdle_models/combination_run/combo_run_1108_loo.rds") # This was the last one by position
model696_loo = readRDS("hurdle_models/combination_run/combo_run_696_loo.rds") # This was the last one by position


# Full model seems to be best 

## Posterior Predictive Checks ##########
good_fit = read_rds(file.path(save_direc, "combo_run_1444_model.rds"))

library(tidybayes)
# Looking at the posterior predictive distribution, creating posterior intervals around it, and comparing that w/ the observed values
preds = predicted_draws(good_fit, newdata = copepods %>% select(numb.worm, cop.lake, worm.fam, worm.lake, plate, native)) 
hdis = point_interval(preds, .prediction, .width = c(.95, .5), .interval = hdci)

hdis %>% arrange(numb.worm) %>% mutate(jitter = runif(n(), -.4, .4)) %>% 
  ggplot(aes(x = numb.worm + jitter, y = .prediction)) + 
  geom_interval(size = .05) + geom_point(size = 1) + scale_color_grey(start = .8, end = .5)

library(moments)
# Let's do a posterior predictive check
pp_dist = 
  preds %>% ungroup %>% group_by(.draw) %>% 
  summarize(mean = mean(.prediction), sd = sd(.prediction),
            skew = skewness(.prediction), kurt = kurtosis(.prediction)) %>% 
  select(-.draw) %>% gather(key = 'Summary', value = 'value')

pp_real = copepods %>% rename(.prediction = numb.worm) %>% 
  summarize(mean = mean(.prediction), sd = sd(.prediction),
            skew = skewness(.prediction), kurt = kurtosis(.prediction)) %>% 
  gather(key = 'Summary', value = 'value')
# PP check does decent job of predicting the moments of the data
# Slightly over-skewed, slightly under-kurtotic
pp_dist %>% ggplot(aes(x = value)) + facet_wrap(~Summary, scale = "free") +
  geom_density() + geom_vline(aes(xintercept = value), data = pp_real, color = "red")

# Let's visualize expected values ##############

betas = get_variables(good_fit)[1:16] %>% # betas are the first set
  map(as.symbol)
beta_draws = gather_draws(good_fit, !!!betas)

bd2 = beta_draws %>% ungroup %>% 
  mutate(hurdle = str_detect(.variable, "_hu_"), 
         .variable = str_remove(.variable, "_hu") %>% str_remove("b_")) %>% group_by(draw, hurdle)
# New plan: create a list of all possible combos of lake, stream, then get preds

intercept_and_native = bd2 %>% filter((.variable %in% c("Intercept", "nativeTRUE"))) %>%
  spread(key = .variable, value = .value)

lakes = bd2 %>% filter(!(.variable %in% c("Intercept", "nativeTRUE"))) %>% 
  separate(.variable, into = c("type", "lake"), by = ".lake") %>% 
  mutate(lake = str_remove(lake, "lake"))

semi_spread = left_join(intercept_and_native, lakes)

semi_spread %>% ungroup %>% distinct(lake, type)
# boot is the reference lake
# add it back in
boot_vals = semi_spread %>% select(-type, -lake, -.value) %>% distinct %>% mutate(lake = "boo", .value = 0)
full_spread = bind_rows(semi_spread, boot_vals %>% mutate(type = "cop"), boot_vals %>% mutate(type = "worm")) %>% 
  group_by(.draw) %>% nest

mutate(value = if_else(hurdle, 1- inv.logit(value), trunc_pois_expect(value))
trunc_pois_expect = function(log_lambda) {
  lambda = exp(log_lambda)
  lambda / (1 - exp(-lambda)) # see wikipeida page on zero-truncated poisson
}

library(gtools)
combos = copepods %>% as_tibble %>% distinct(worm.lake, cop.lake) %>% mutate_all(as.character) %>%  mutate(native = worm.lake == cop.lake)
expected_values = full_spread %>%
  mutate(etas = map(data, ~pmap_df(combos, function(worm.lake, cop.lake, ...) {
    .x %>% filter((type == "worm"& lake == worm.lake) | (type == "cop"& lake == cop.lake)) %>% 
      select(-lake) %>% spread(key = type, value = .value) %>% 
      mutate(value = if_else(rep(worm.lake == cop.lake, 2), nativeTRUE, 0) + Intercept + cop + worm) %>% 
      select(hurdle, value) %>% mutate(value = if_else(hurdle, 1- inv.logit(value), trunc_pois_expect(value)), 
                                       cop = cop.lake, worm = worm.lake)
  })))

hdis = expected_values$etas %>% bind_rows %>%  
  mutate(value = if_else(hurdle, 1-value, 
                         trunc_pois_expect(logit(value)))) %>% # do this only until you re-run the above expected_value code
  group_by(cop, worm, hurdle) %>% 
  point_interval(value, .width = c(.95, .5), .interval = hdci) %>% 
  mutate(Uncertainty = if_else(.width == .5, "50% HDI", "95% HDI") %>% factor(., levels = unique(.))) %>% 
  ungroup %>% 
  mutate(worm = paste(worm, "worms")) %>% 
  mutate(hurdle = if_else(hurdle, "Probability of infection", "Expected Infection intensity")) %>% 
  mutate(cop_pos = as.integer(as.factor(cop)))

raw_data_for_plot = 
  copepods %>% as_tibble %>% select(cop = cop.lake, worm = worm.lake, numb.worm) %>% 
  mutate(infected = as.integer(numb.worm > 0)) %>% 
  gather(numb.worm, infected, key = "hurdle", value = "value") %>% 
  filter(value > 0 | hurdle == "infected") %>% # drop uninfected from the intensity column
  mutate(worm = paste(worm, "worms"),
         hurdle = if_else(hurdle == "infected",  "Probability of infection", "Expected Infection intensity"))

# These all seem low for the infection intensity. I'm betting that they should be higher by 1.  Need to think on it

median_bar_width = .1
infection_plot = 
  hdis %>% filter(hurdle ==  "Probability of infection") %>% 
  ggplot(aes(x = cop, y = value)) +
  geom_jitter(data = raw_data_for_plot %>% filter(hurdle ==  "Probability of infection"), 
              width = .3, height = .1, size = .5, shape = 21, alpha = .15, fill = "red", color = "red")+
  geom_interval(aes(size = Uncertainty, color = Uncertainty)) +
  # geom_point(size = 3, shape = 4) + 
  geom_intervalh(aes(xmin = cop_pos - median_bar_width, xmax = cop_pos + median_bar_width), color = "black", size = .8)+
  facet_wrap( ~worm, ncol = 1) +
  scale_color_grey(start = .8, end = .6) + 
  theme(strip.placement = "outside", strip.background.y = element_blank())+
  scale_size_manual(values = c(1, 1.5)) + ylab("Probability of infection") + xlab("Copepod Lake")
intensity_plot = 
  hdis %>% filter(hurdle == "Expected Infection intensity") %>% 
  ggplot(aes(x = cop, y = value)) +
  geom_jitter(data = raw_data_for_plot %>% filter(hurdle == "Expected Infection intensity"),
              width = .3, height = .3, size = .5, shape = 21, alpha = .3, fill = "blue", color = "blue")+
  geom_interval(aes(size = Uncertainty, color = Uncertainty)) +
  # geom_point(size = 3, shape = 4) + 
  geom_intervalh(aes(xmin = cop_pos - median_bar_width, xmax = cop_pos + median_bar_width), color = "black", size = .8)+
  facet_wrap( ~worm, ncol = 1) +
  scale_color_grey(start = .8, end = .6) + 
  theme(strip.placement = "outside", strip.background.y = element_blank())+
  scale_size_manual(values = c(1, 1.5)) + ylab( "Expected Infection intensity") + xlab("Copepod Lake")
ci_legend = get_legend(intensity_plot)
no_legend = theme(legend.position = "none")
joint_plot = 
  plot_grid(plot_grid(infection_plot+no_legend, intensity_plot+no_legend, align = "hv", axis = "tb", ncol = 2),
            ci_legend, ncol = 2, rel_widths = c(1, .2))
# cairo_png = partial(png, type = "cairo") # let's use a nicer png renderer # not working atm
ggsave(filename = "hurdle_models/figures/bayes_figure.pdf", joint_plot, dpi = 500, width = 8, height = 6)

write_csv(hdis, "hurdle_model/run_1444_transformed_hdis.csv")

## What about well number? ###########

as_tibble(copepods) %>% pull(well.numb) %>% unique

model_formulas[[1444]]
model_formulas[[1108]]   # OH, right, this is the best one
fit_1108 = readRDS(file.path(save_direc,"combo_run_1108_model.rds")) 

add_well_form = 
  bf(numb.worm ~ cop.lake + worm.lake ,
     hu ~ cop.lake + worm.lake + native )

run_hurdle_model(add_well_form, priors = model_priors[[1444]], name = "combo_run_1108_plus_well",  
                 save_dir = save_direc, control = list(adapt_delta = .99), data = copepods)
extra_fit2 = readRDS(file.path(save_direc,"combo_run_1108_plus_well_model.rds")) 
extra_loo2 = brms::loo(extra_fit, cores = 4)

loo_model_weights(lst(model1444_loo,model999_loo,
                      model1108_loo, extra_loo, extra_loo2, model696_loo))
# It looks like the 

