data_path <- "/extraspace/TCGA/TCGA_protein"

file_name <- list.files(data_path, "_protein.re_20160627")
file_name <- file_name[!grepl(file_name, pattern = "GBM|LAML")]

library(tidyverse)
library(stringr)


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
    drop_na() 
  
  snor_pro <- read_tsv(
    file.path(
      "/home/jgong/gongj_snoRNA/database_data/protein",
      list.files(
        path = "/home/jgong/gongj_snoRNA/database_data/protein", 
        pattern = paste(name,".*_id", sep = "")
        )
      )
    ) %>%
    mutate(dataset_id = dataset_id) %>%
    select(snorna = sno, dataset_id, protein = gene, spearman_corr = r, p_value=p, fdr) 
  
  snor_pro %>%
    left_join(data, by = "protein") %>%
    gather(sample_id, expression, -c(snorna,dataset_id, protein, spearman_corr, p_value, fdr)) %>%
    drop_na() ->
    result_data
  result_data %>% 
    write_tsv(path=paste("/home/cliu18/liucj/projects/4.SNORic/database/protein_expresion/protein_expression", name,"tsv", sep="."),na="", col_names=F)
}

file_name %>%
  lapply(get_meth) 

# %>%
  # write_tsv(path = "/home/cliu18/liucj/projects/4.SNORic/database/protein_expresion/protein_expression.tsv", col_names = F)
