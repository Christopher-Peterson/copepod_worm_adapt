suppressPackageStartupMessages({
  library(dplyr);library(tidyr); library(purrr);
  library(tibble); library(readr); library(rlang);
  library(brms)
})

argv = commandArgs(TRUE)
job_file = argv[1] %|% "jobs/run_hurdle.job"

# This sets up the hurdle model run

set.seed(93093) # to set seeds for all bayes runs
dir.create("out")
dir.create("out/model_runs")
dir.create("out/loo")
## Create a list of all predictor combinations ========
terms = c("cop.lake", "worm.lake", "worm.lake:cop.lake", 
          "genus", "genus:worm.lake", # We're explicitly identifying the interactions here so they have the possibility of being included
          "native", "(1|plate)", "(1|worm.fam)")

## Removes combinations where interaction is present 
## without corresponding main effect
fix_interaction = function(arg) {
  if("worm.lake:cop.lake" %in% arg) {
    if(! (("cop.lake" %in% arg) && ("worm.lake" %in% arg)))
      arg = "" # empty string will later be filtered out
  } # copy this nested if statement & make it work for genus:worm.lake
  if("worm.lake:genus" %in% arg) {
    if(! (("genus" %in% arg) && ("worm.lake" %in% arg)))
      arg = ""
  }
  arg
}
remove_null = function(.x) .x [.x != "~ "] # discard removed parts

# Create a combination of all predictors
pick_n_predictors = function(n, preds = terms) {
  combn(seq_along(preds), m = n, simplify = FALSE) %>% 
    map(~preds[.x])
}

rhs_combos = 
  map(seq_along(terms), pick_n_predictors) %>%  # 1 - 8 possible predictors in a model
  flatten %>% map(fix_interaction) %>%
  # Convert to character version of a formula
  map_chr(~paste("~", paste(.x, collapse = " + "))) %>% 
  remove_null()
# Since this is a two part model, we need all 2-part combos of the above
model_formulas_txt = 
  list(pois = paste("numb.worm", rhs_combos),
       hu = paste("hu", rhs_combos)) %>% cross() 
# Convert them into brms formulae
model_formulas = model_formulas_txt %>% 
  map(~bf(as.formula(.x$pois), as.formula(.x$hu)))

### Set Priors ===========

# One challenge is determining what effects are present
detect_rand_eff = function(val) {
  grepl("plate", x = val, fixed = TRUE) |
    grepl("worm.fam", x = val, fixed = TRUE)
}
detect_fixed_eff = function(val) {
  grepl("lake", x = val, fixed = TRUE) | 
  grepl("native", x = val, fixed = TRUE) |
  grepl("genus", x = val, fixed = TRUE)
}
determine_prior = function(.x){
  re_hurdle = detect_rand_eff(.x$hu)
  re_pois   = detect_rand_eff(.x$pois)
  b_hurdle  = detect_fixed_eff(.x$hu)
  b_pois    = detect_fixed_eff(.x$pois)
  
  prior = prior_string("normal(0, 6)", class = "Intercept") #the prior in the intercept always has to be bigger than any other priors below 
  if(isTRUE(re_hurdle)) prior = prior + prior_string("student_t(7, 0, 1)", class = "sd", dpar = "hu")
  if(isTRUE(re_pois)) prior = prior + prior_string("student_t(7, 0, 1.5)", class = "sd")
  if(isTRUE(b_hurdle)) prior = prior+ prior_string("normal(0, 1)", class = "b", dpar = "hu")
  if(isTRUE(b_pois)) prior = prior + prior_string("normal(0, 1.5)", class = "b")
  prior
}
# has_re = function(model_txt) {
#   browser()
#   any(grepl(pattern = "(1|", x = unlist(model_txt), fixed = TRUE))
# }

run_table = tibble(
    model_txt = model_formulas_txt,
    formula = model_formulas) %>% 
  mutate(priors = map(model_txt, determine_prior),
         id = 1:n(), name = paste0("combo_run_", id),
         adapt_delta = .95, 
         # Create rng seeds for model
         rng_seed = purrr::rdunif(n(), as.integer(2L ^ 31 - 1)))
write_rds(run_table, file = "model_settings.rds")

# Write job file
exec_cmd = "docker_stan" # "Rscript"
run_script = "R/run_hurdle.r"
paste(exec_cmd, run_script, 1:nrow(run_table)) %>% 
  write_lines(job_file)
