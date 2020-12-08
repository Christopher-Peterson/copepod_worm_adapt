# Run a version of the hurdle model on just the M. albidus 
# Boot lake worms vs Boot & Law copepods
library(dplyr)
library(readr)
library(brms)
library(loo)  

cop_dat = read_csv("data/chapter2.copepods.cleaned.csv")

# Q1: Do boot worms do better with boot copepods vs lau copepods (both M. albidus)? ####

q1_data = cop_dat %>% 
  filter(cop.lake %in% c("lau", "boo"), worm.lake == "boo") %>% 
  mutate(infected = as.integer(infected.yes.no == "y")) %>% 
  select(cop.lake, worm.fam, infected, n_worm = numb.worm, plate) %>% 
  mutate(plate = if_else(is.na(plate), "na_plate", plate))

q1_form_full = brmsformula(
  n_worm ~ cop.lake + (1 | worm.fam) + (1 | plate),
  hu ~ cop.lake + (1 | worm.fam) + (1 | plate),
  family = hurdle_poisson())
q1_form_simp = brmsformula(
  n_worm ~ cop.lake, hu ~ cop.lake,
  family = hurdle_poisson())
q1_form_null = brmsformula(
  n_worm ~ 1, hu ~ 1,
  family = hurdle_poisson())
q1_form_re = brmsformula(
  n_worm ~ (1 | worm.fam) + (1 | plate),
  hu ~ (1 | worm.fam) + (1 | plate),
  family = hurdle_poisson())

q1_priors_simp = q1_priors_null + 
  prior_string("normal(0, .5)", class = "b", dpar = "hu") +
  prior_string("normal(0, .25)", class = "b")
q1_priors_re = 
  prior_string("student_t(7, 0, 0.25)", class = "sd", dpar = "hu") +
  prior_string("student_t(7, 0, 0.125)", class = "sd")
q1_priors_re_norm = 
  prior_string("normal(0, 0.25)", class = "sd", dpar = "hu") +
  prior_string("normal(0, 0.125)", class = "sd")

q1_model_full = 
  brm(q1_form_full, data = q1_data, 
      prior = q1_priors_simp + q1_priors_re, 
      control = list(adapt_delta = .99))
q1_model_full_norm = 
  brm(q1_form_full, data = q1_data, 
      prior = q1_priors_simp + q1_priors_re_norm, 
      control = list(adapt_delta = .99))
q1_model_simp = 
  brm(q1_form_simp, data = q1_data, 
      prior = q1_priors_simp, 
      control = list(adapt_delta = .99))
q1_model_re = 
  brm(q1_form_re, data = q1_data, 
      prior = q1_priors_null + q1_priors_re, 
      control = list(adapt_delta = .99))

# Get stacking weights:
q1_stacked_models = brms::loo_model_weights(
  q1_model_full,
  q1_model_full_norm,
  q1_model_simp,
  q1_model_re,
  model_names = c("Full", "Compromise",
                  "Fixed Effect", "Random Effects")
)


# Q2: Do boot copepods have a preference for M. albidus over other copepods? ####
q2_data = cop_dat %>% 
  filter(worm.lake == "boo") %>% 
  mutate(infected = as.integer(infected.yes.no == "y")) %>% 
  select(cop.lake, worm.fam, infected, n_worm = numb.worm, plate) %>% 
  mutate(different_spp = !cop.lake %in% c("lau", "boo"),
         different_genus = cop.lake %in% c("ech", "gos"))
  
q2_form_full = brmsformula(
  n_worm ~ different_genus + different_spp + (1 | cop.lake),
  hu ~ different_genus + different_spp + (1 | cop.lake),
  family = hurdle_poisson())
q2_form_spp_re = brmsformula(
  n_worm ~ different_spp + (1 | cop.lake),
  hu ~ different_spp + (1 | cop.lake),
  family = hurdle_poisson())
q2_form_gen_re = brmsformula(
  n_worm ~ different_genus + (1 | cop.lake),
  hu ~ different_genus + (1 | cop.lake),
  family = hurdle_poisson())
q2_form_spp = brmsformula(
  n_worm ~ different_spp ,
  hu ~ different_spp ,
  family = hurdle_poisson())
q2_form_gen = brmsformula(
  n_worm ~ different_genus ,
  hu     ~ different_genus ,
  family = hurdle_poisson())
q2_form_spp_gen = brmsformula(
  n_worm ~ different_genus + different_species,
  hu     ~ different_genus + different_species,
  family = hurdle_poisson())
q2_form_re = brmsformula(
  n_worm ~ (1 | cop.lake),
  hu ~ (1 | cop.lake),
  family = hurdle_poisson())

q2_prior_int = 
  prior_string("normal(0,1)", class = "Intercept") +
  prior_string("normal(0,1)", class = "Intercept", dpar = "hu") 
q2_prior_re_nml = 
  prior_string("normal(0, 0.25)", class = "sd", dpar = "hu") +
  prior_string("normal(0, 0.125)", class = "sd") 
q2_prior_re_t = 
  prior_string("student_t(7, 0, 0.25)", class = "sd", dpar = "hu") +
  prior_string("student_t(7, 0, 0.125)", class = "sd") 
q2_prior_cont =  
  prior_string("normal(0, .5)", class = "b", dpar = "hu") +
  prior_string("normal(0, .25)", class = "b")

q2_model_full = 
  brm(q2_form_full,  data = q2_data, 
      prior = q2_prior_int + q2_prior_re_t + q2_prior_cont,
      control = list(adapt_delta = .99))
q2_model_full_nml = 
  brm(q2_form_full,  data = q2_data, 
      prior = q2_prior_int + q2_prior_re_nml + q2_prior_cont,
      control = list(adapt_delta = .99))
q2_model_re_spp = 
  brm(q2_form_spp_re,  data = q2_data, 
      prior = q2_prior_int + q2_prior_re_t + q2_prior_cont,
      control = list(adapt_delta = .99))
q2_model_re_genus = 
  brm(q2_form_genus_re,  data = q2_data, 
      prior = q2_prior_int + q2_prior_re_t + q2_prior_cont,
      control = list(adapt_delta = .99))
q2_model_full_nml = 
  brm(q2_form_full,  data = q2_data, 
      prior = q2_prior_int + q2_prior_re_nml + q2_prior_cont,
      control = list(adapt_delta = .99))


q2_model_nocont = 
  brm(q2_form_nocont,  data = q2_data, 
      prior = q2_prior_int + q2_prior_re_t,
      control = list(adapt_delta = .99))
q2_model_nore = 
  brm(q2_form_nore,  data = q2_data, 
      prior = q2_prior_int + qt_prior_cont,
      control = list(adapt_delta = .99))

# Get stacking weights:
q2_stacked_models = brms::loo_model_weights(
  q2_model_nore,
  q2_model_full_t,
  q2_model_full_nml,
  q2_model_nocont,
  model_names = c("Fixed Effect", "Full", "Compromise", "Random Effect")
)



### Q2: Investigate & Plot #####
# Investigate q2_model_full_t
set.seed(9842892)
q2_scaled_stacks = q2_stacked_models/max(q2_stacked_models) 
q2_newdat = q2_data %>% select(cop.lake, different_spp) %>% distinct()
q2_re_draws = sample(1:4000, size = 4000 * q2_scaled_stacks[4])
# Helper functions
ev_freq = function(.x, ...) {
  1 - gtools::inv.logit(.x)
}
ev_num = function(.x, ...) {
  # x = .x[draws, par_ids , drop = FALSE]
  # Truncated poisson expected value (see Wikipedia entry for formula)
  # .x is log_lambda
  lambda = exp(.x)
  lambda / (1 - exp(-lambda))
}
get_ci = function(.x, probs = c(.025, .05, .25, .5, .75, .95, .975)){
  quants = quantile(.x, probs, na.rm = TRUE)
  names(quants) = paste0("Q", substring(probs, 2))
  tibble(!!!quants)
}
order_cop_lake = function(cop_lakes) {
  # Convert cop.lake into a factor
  lake_order = list(
    boo = "Boot",  lau = "Lawier",
    rob = "Roberts",
    ech = "Echo",
    gos = "Gossling") 
  rlang::exec(recode, .x = cop_lakes, !!!lake_order) %>%
    factor(levels = lake_order)
}
extract_linpred = function(model, ..., newdata = q2_newdat) {
  brms::posterior_linpred(model, scale = "linear", ...,
                          summary = FALSE, newdata = newdata) 
}


# Get model-stacked expected values
num_worms = rbind(
    extract_linpred(q2_model_full_t),
    extract_linpred(q2_model_nocont, subset = q2_re_draws) ) %>% ev_num()
worm_freq = rbind(
  extract_linpred(q2_model_full_t, dpar = "hu"),
  extract_linpred(q2_model_nocont, dpar = "hu", 
                  subset = q2_re_draws) ) %>% ev_freq()
# Connect EVs to their predictors & get CIs
ev_ci = q2_newdat %>% 
  mutate(Intensity  = map(1:n(), ~num_worms[, .x]),
         Prevalence = map(1:n(), ~worm_freq[, .x])) %>% 
  pivot_longer(c(Intensity, Prevalence), 
               names_to = "component", values_to = ".value") %>% 
  mutate(ci = map(.value, get_ci)) %>% 
  select(-.value) %>% 
  unnest(ci) %>% 
  # select(median = Q.5, lo = Q.025, hi = Q.975) %>% 
  ungroup() %>% 
  mutate(species = if_else(different_spp, "Other", "M. albidus")) %>% 
  mutate(cop.lake = order_cop_lake(cop.lake))

q2_infection_emp = q2_data %>% 
  mutate(species = if_else(different_spp, "Other", "M. albidus")) %>% 
  group_by(cop.lake, species) %>% 
  summarize(prevalence = mean(infected)) %>% 
  mutate(cop.lake = order_cop_lake(cop.lake)) %>% 
  mutate(component = "Prevalence")

# Visualize it
ev_ci %>% ggplot() + 
  aes(x = cop.lake, y = Q.5, color = species) + 
  facet_grid(component~species, scale = "free", space = "free_x", switch = "both") + 
  geom_col(aes(y = prevalence), color = "black", fill = "#ffffff00",
           data = q2_infection_emp)+
  geom_linerange(aes(ymin = Q.025, ymax = Q.975), size = .8, alpha = .5) +
  geom_linerange(aes(ymin = Q.25, ymax = Q.75), size = 1.5) +
  geom_point(size = 5, shape = 3) + 
  scale_color_viridis_d(begin = .2, end = .7) + 
  ylab("Effect Size (Boot Tapeworms)") + 
  xlab("Copepod lake & Species") +
  cowplot::theme_cowplot() + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(face = "italic"),
        strip.placement = "outside",
        legend.position = "none")
ggsave("figures/boo_worm_comparison.png", dpi = 300, width = 6, height = 5)


