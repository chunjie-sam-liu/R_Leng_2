#!/usr/bin/Rscript
#
library(methods)
library(tidyverse)
library(stringr)


clinical_path = '/home/cliu18/liucj/projects/4.SNORic/database/clinical'

clinical_files = list.files(path = clinical_path, pattern = 'txt$')

load_data <- 
  function(x){
    
    path = file.path(clinical_path, x)
    dataset_id <- paste('TCGA',str_split(x,'_', simplify = T)[1,1],sep='-')
    print(dataset_id)
    label = paste(str_split(x,'_', simplify = T)[1,1],'Tumor',sep='-')
    print(label)
    
    d <- 
      read_tsv(file = path)
    
    d %>%
      mutate(
        dataset_id = dataset_id,
        sample_id = paste(label, barcode, sep = "-")
      ) %>%
      dplyr::select(
        matches('Subtype|stage|grade|tobacco'),
        dataset_id,
        sample_id,
        time = os_days,
        status = os_status
      ) ->
      d_select

    d_select %>%
      gather(subtype,stage, 
             -c(dataset_id, sample_id, time, status)) %>%
      mutate(
        status = plyr::revalue(status, c("Alive" = 0, "Dead" = 1))
      )
  }

clinical_files %>%
  lapply(load_data) %>%
  bind_rows() %>%
  write_tsv(file.path(clinical_path, 'clinical_dataset.tsv'))


read_tsv(file.path(clinical_path, 'clinical_dataset.tsv')) -> clinical_data

clinical_data %>%
  group_by(dataset_id, subtype,stage) %>% 
  count() %>%
  drop_na() %>%
  View()
