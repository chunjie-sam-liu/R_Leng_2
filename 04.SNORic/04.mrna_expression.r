data_path <- '/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp'

file_name <- list.files(data_path, "_mRNA_each_exp_20160513")
file_name <- file_name[!grepl(file_name, pattern = "GBM|LAML")]
library(tidyverse)
library(stringr)

database_path <- "/home/jgong/gongj_snoRNA/database_data/mRNA"

host_data <- read_tsv("/home/jgong/gongj_snoRNA/database_data/snoRNA_overlap_gene.position.type", col_names = c("sno", "gene", "h", "tyep")) %>% 
  mutate(host = "h") %>%
  select(sno, gene, host)

get_mrna_expr <- function(x){
  print(x)
  name <- str_replace(x, "_.*", "")
  snor_gene <- read_tsv(file.path(database_path, list.files(path = database_path, pattern = name)))
  dataset_id <- paste("TCGA", name, sep='-')
  
  snor_gene %>% filter(!grepl("X\\.\\.", gene), abs(r) > 0.3) %>%
    separate(gene, into= c("gene", "entrez"), sep = "\\.") %>%
    mutate(dataset_id = dataset_id) %>% 
    left_join(host_data, by = c("sno", "gene")) %>%
    select(snorna = sno, dataset_id, gene, host, r, p, fdr) ->
    snor_gene_sep
  
  
  data <- read_tsv(file.path(data_path, x)) %>% select(-length(.))
  type <- str_sub(colnames(data),14,14)
  col_name <- str_sub(colnames(data),1,12)
  col_name[type=="1"] = paste(name,"Normal", col_name[type=="1"], sep="-")
  col_name[type=="0"] = paste(name,"Tumor", col_name[type=="0"], sep="-")
  colnames(data) <- col_name
  
  data %>% 
    subset(select = which(!duplicated(colnames(.)))) %>%
    filter(!grepl("\\?", gene)) %>%
    separate(gene, into=c("gene","entrez"), sep="\\|") %>% 
    select(-entrez) ->
    data_sep
  
  snor_gene_sep %>% 
    left_join(data_sep, by = "gene") %>%
    gather(sample_id, rna_expression, -c(snorna, dataset_id, gene, host, r, p, fdr)) %>%
    replace_na(list(host = "\\N")) -> result_data
  return(result_data)
}

file_name %>% 
  lapply(get_mrna_expr) %>%
  bind_rows() %>%
  write_tsv(path = "/home/cliu18/liucj/projects/4.SNORic/database/mrna_expression/mrna_expression.tsv", col_names = F)



