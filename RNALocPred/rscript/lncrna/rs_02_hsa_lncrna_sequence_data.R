library(tidyverse)
library(ggthemes)
library(stringr)
library(rvest)
library(multidplyr)


root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")

human_path_mirna <- file.path(human_path, "lncrna")

hsa_mirna_clean <- read_rds(path = file.path(human_path_mirna, "rds_02_hsa_lncrna_clean.rds.gz"))


RNALocate_seq <- function(oid, return_format = "tibble", ...){
  url <- "http://www.rna-society.org/rnalocate/sequence.php?Gene_ID"
  url <- str_c(url, oid,sep = "=")
  print(oid)
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

cl <- parallel::detectCores()
cluster <- create_cluster(floor(cl * 5 / 6))
hsa_mirna_clean %>% 
  # register cores
  partition(cluster = cluster) %>% 
  # load library to cluster
  cluster_library("tidyverse") %>% 
  cluster_library("stringr") %>% 
  cluster_library("rvest") %>% 
  cluster_assign_value("RNALocate_seq", RNALocate_seq) ->
  hsa_mirna_clean_shards


hsa_mirna_clean_seq <-
  hsa_mirna_clean_shards %>% 
  # get names and seq
  mutate(out = map(.x = oid, ~ RNALocate_seq(oid = .x, return_format = "tibble"))) %>% 
  collect() %>% 
  as_tibble() %>% 
  unnest()
parallel::stopCluster(cluster) #


hsa_mirna_clean_seq_ready <-
  hsa_mirna_clean_seq %>% 
  ungroup() %>% 
  dplyr::select(-PARTITION_ID) %>% 
  mutate(seq_len = str_length(seq)) 





























