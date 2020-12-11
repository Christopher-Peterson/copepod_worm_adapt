library(tidyverse)
cleaned_data = read_csv("data/chapter2.copepod.cleaned.Aug2020.csv") %>% 
  select(number, cop.lake, worm.fam, worm.lake, plate, numb.worm) %>% 
  replace_na(list(plate = "p100")) %>% 
  mutate(native = cop.lake == worm.lake,
         genus = if_else(cop.lake %in% c("gos", "ech"), "A", "M" )) # This is where you'd want to define "genus"

write_csv(cleaned_data, "data/chapter_2_copepod_for_bayes.csv")
