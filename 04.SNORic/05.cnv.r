data_path <- "/home/jgong/gongj_snoRNA/database_data/CNV_data"

file_name <- list.files(data_path, "_snoRNA_TCGA_CNV_20161216.match")
library(tidyverse)
library(stringr)

get_cnv<- function(x){
  name <- str_replace(x, "_.*", "")
  dataset_id <- paste("TCGA", name, sep='-')
  
  data <- read_tsv(file.path(data_path, x)) 
  type <- str_sub(colnames(data),14,14)
  col_name <- str_sub(colnames(data),1,12)
  col_name[type=="1"] = paste(name,"Normal", col_name[type=="1"], sep="-")
  col_name[type=="0"] = paste(name,"Tumor", col_name[type=="0"], sep="-")
  
  colnames(data) <- col_name
  
  data %>% 
    subset(select = which(!duplicated(colnames(.)))) %>% 
    mutate(dataset_id, dataset_id) %>%
    gather(sample_id, copy, -c(cnv, dataset_id))
    
}

file_name %>%
  lapply(get_cnv) %>%
  bind_rows() %>%
  write_tsv(path = "/home/cliu18/liucj/projects/4.SNORic/database/cna/cna.tsv", col_names = F)
