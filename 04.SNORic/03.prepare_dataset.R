#!/usr/bin/Rscript
#
library(methods)
library(tidyverse)
library(stringr)

# without chol clinical information
dataset_path = '/home/cliu18/liucj/projects/4.SNORic/database/dataset/dataset'

dataset <- 
  read_tsv(
    dataset_path,
    col_names = c(
      "dataset_description",
      "cancer_type",
      "total",
      "tumor_n",
      "normal_n",
      "average_mappable_reads",
      "snorna_n",
      "snorna_rpkm_n"
      )
           )  %>%
  mutate(
    dataset_id = paste("TCGA", cancer_type, sep = "-"),
    source = "TCGA",
    build = "hg19"
  ) %>%
  select(
    dataset_id,
    source,
    cancer_type,
    dataset_description,
    normal_n,
    tumor_n,
    build,
    average_mappable_reads,
    snorna_n,
    snorna_rpkm_n
  ) 

dataset %>%
  write_tsv(paste(dataset_path, 'tsv', sep = '.'))
