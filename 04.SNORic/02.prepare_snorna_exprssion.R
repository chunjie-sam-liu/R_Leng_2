#!/usr/bin/Rscript
#
library(methods)
library(tidyverse)
library(stringr)

# without chol clinical information
expression_path = '/home/cliu18/liucj/projects/4.SNORic/database/snorna_exprssion'

expression_files = list.files(expression_path, pattern = 'matrix$')


load_expression <- function(x){
  path = file.path(expression_path, x)
  label = str_split(x,'\\.', simplify = T)[1,2]
  dataset_id <- paste('TCGA',label,sep='-')
  print(dataset_id)
  
  d <- 
    read_tsv(path) %>%
    gather(sample_id, snorna_expression, -snoRNA) %>%
    mutate(dataset_id = dataset_id) %>%
    select(
      snorna = snoRNA,
      dataset_id,
      sample_id,
      snorna_expression
    )
  
  d %>%
    mutate(
      code = str_match(sample_id, pattern = "-(\\d\\d)$")[,2],
      barcode = str_replace(sample_id, pattern='-\\d\\d$',''),
      code = ifelse(code=="01", "Tumor", "Normal"),
      sample_id = paste(label, code, barcode, sep = '-')
  ) %>%
    select(
      snorna,
      dataset_id,
      sample_id,
      snorna_expression
    )
}

expression_files %>% 
  lapply(load_expression) %>%
  bind_rows() %>%
  write_tsv(file.path(expression_path, 'snorna_expression.tsv'))