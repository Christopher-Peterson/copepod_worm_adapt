library(tidyverse)

model_settings = read_rds("/home/peterson/Downloads/model_settings.rds")
completed = dir("out/loo") %>% str_remove(".rds")

models_expanded = model_settings %>% 
  mutate(model_txt2 = map(model_txt, as_tibble)) %>% 
  select(id, name, model_txt2) %>% 
  unnest(model_txt2)


still_to_do = models_expanded %>% 
  filter(!(name %in% completed) | 
      # This is a one-time fix to deal w/ the plate issue
      # Remove the | from the previous line and next line and
      # close the parentheses on the filter command once you've run this once.
      str_detect(hu, "plate") | str_detect(hu, "plate")) 

cmd    = "docker_stan" # WIll, change this to whatever you called you version of it
script = "run_hurdle.r" # same
new_job = glue("{cmd} {script} {still_to_do$id}")

write_lines(new_job, "long_run_hurdle.job")