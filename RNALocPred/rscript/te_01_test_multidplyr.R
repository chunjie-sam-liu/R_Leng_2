library(tidyverse)
library(ggthemes)
library(stringr)
library(rvest)
library(multidplyr)


root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")

hsa_mirna_clean <- read_rds(path = file.path(human_path_mirna, "rds_02_hsa_mirna_clean.rds.gz"))

# ----------------------------------------get sequence from RNALocate-------------------------
RNALocate_seq <- function(oid, return_format = "tibble", ...){
  url <- "http://www.rna-society.org/rnalocate/sequence.php?Gene_ID"
  url <- str_c(url, oid,sep = "=")
  
  content <- read_html(url) %>% 
    html_nodes("pre") %>% 
    html_text() 
  
  if (length(content) == 0){
    tibble(names = NA, seq = NA)
  } else{
    content %>%
      as_tibble() %>% 
      mutate(type = c("names", "seq")) %>% 
      spread(type, value) %>% 
      mutate(
        names = str_replace(names, ">", ""),
        seq = str_replace(seq, "\\r","")
      )
  }
}

# test multidplyr
cl <- parallel::detectCores()
cluster <- create_cluster(floor(cl * 5 / 6))

hsa_mirna_clean %>% 
  head(40) %>% 
  partition(cluster = cluster) -> te

te %>% 
  cluster_library("tidyverse") %>% 
  cluster_library("stringr") %>% 
  cluster_library("rvest") %>% 
  cluster_assign_value("RNALocate_seq", RNALocate_seq)

te
cluster_eval(te, search())[[1]]
cluster_get(te, "RNALocate_seq")[[1]]

start <- proc.time()
te %>% 
  mutate(
    out = map(.x = oid, ~ RNALocate_seq(oid = .x, return_format = "tibble"))
  ) %>% 
  collect() %>% 
  as_tibble() %>% 
  unnest()
time_elapsed_series <- proc.time() - start
time_elapsed_series
parallel::stopCluster(cluster)

# in seria
start <- proc.time()
hsa_mirna_clean %>% 
  head(40) %>% 
  mutate(
    out = map(.x = oid, ~ RNALocate_seq(oid = .x, return_format = "tibble"))
  )
time_elapsed_series <- proc.time() - start
time_elapsed_series

