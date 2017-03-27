data_path <- "/extraspace/TCGA/TCGA_methylation450k"

file_name <- list.files(data_path, "_methy_probe_20160801")
file_name <- file_name[!grepl(file_name, pattern = "GBM|LAML")]

library(tidyverse)
library(stringr)

whole_meth <- read_tsv(
  file = "/extraspace/jgong/gongj_snoRNA/database_data/methylation_data/whole_snoRNA_UCSC_GENECODE.up2k.methyprobes", 
  col_names = c("snorna", "chrom", "strand", "start", "end", "source", "distance", "chip", "fuck",'pos' )) %>%
  select(-fuck)


get_meth <- function(x) {
  print(x)
  name <- str_replace(x, "_.*", "")
  dataset_id <- paste("TCGA", name, sep='-')
  
  data <- read_tsv(file.path(data_path, x)) %>% select(-length(.))
  type <- str_sub(colnames(data),14,14)
  col_name <- str_sub(colnames(data),1,12)
  col_name[type=="1"] = paste(name,"Normal", col_name[type=="1"], sep="-")
  col_name[type=="0"] = paste(name,"Tumor", col_name[type=="0"], sep="-")
  colnames(data) <- col_name
  
  data %>% 
    subset(select = c(1,which(grepl(name, colnames(.))))) %>%
    subset(select = which(!duplicated(colnames(.)))) %>%
    drop_na() %>%
    separate(gene, into = c('chip', 'gene'), sep = '_') %>%
    select(-gene) ->
    data_sep
  
  whole_meth %>%
    left_join(data_sep, by = "chip") %>%
    gather(sample_id, methy, -c(snorna, chrom, strand, 
                                start, end, source, distance, chip, pos)) %>%
    mutate(dataset_id = dataset_id) %>%
    select(snorna, dataset_id, sample_id, chrom, strand, start, end, source, distance, meth_id = chip, pos, level = methy)
  
  
}

file_name %>%
  lapply(get_meth) %>%
  bind_rows() %>%
  write_tsv(path = "/home/cliu18/liucj/projects/4.SNORic/database/methylation/methylation.tsv", col_names = F)
