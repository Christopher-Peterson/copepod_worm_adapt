suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(rlang)
} )
if(!exists('argv')) argv = commandArgs(TRUE)
job_file = argv[1] %|% 'jobs/continue_run_hurdle.job'

model_settings = read_rds("/home/peterson/Downloads/model_settings.rds")
completed = dir("out/loo") %>% str_remove(".rds")

models_expanded = model_settings %>% 
  mutate(model_txt2 = map(model_txt, as_tibble)) %>% 
  select(id, name, model_txt2) %>% 
  unnest(model_txt2)


still_to_do = models_expanded %>% 
  filter(!(name %in% completed) )
      # If this fails, replace it w/ the subsequent filter line
      # THere was an issue with plate encoding that should have been resolved
      # but may not have been. 
  # filter(!(name %in% completed) | 
  #     str_detect(hu, "plate") | str_detect(hu, "plate")) 

cmd    = "docker_stan" # WIll, change this to whatever you called you version of it
script = "run_hurdle.r" # same
new_job = glue("{cmd} {script} {still_to_do$id}")

write_lines(new_job,  job_file)