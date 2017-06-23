library(tidyverse)
library(stringr)
library(Biostrings)
library(seqLogo)


root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "snorna")
hsa_mirna_clean <- read_rds(path = file.path(human_path_mirna, "rds_03_hsa_snorna_clean_seq_ready.rds.gz"))
hsa_mirna_clean %>% 
  mutate(subloc = ifelse(subloc %in% c("Nucleolus", "Nucleus"), "Nucleus", "non-Nucleus")) %>% 
  drop_na() -> hsa_mirna_clean_1

hsa_mirna_clean_1 %>% filter(subloc == "Nucleus") %>% mutate(seq = str_replace_all(seq, "T", "U"))-> nucleus
hsa_mirna_clean_1 %>% filter(subloc != "Nucleus") %>% mutate(seq = str_replace_all(seq, "T", "U"))-> non_nucleus

nucleus_rnastringset <- RNAStringSet(x = nucleus$seq)
names(nucleus_rnastringset) <- nucleus$oid
writeXStringSet(x = nucleus_rnastringset, filepath = file.path(human_path_mirna, "nucleus.fa"))
non_nucleus_rnastringset <- RNAStringSet(x = non_nucleus$seq)
names(non_nucleus_rnastringset) <- non_nucleus$oid
writeXStringSet(x = non_nucleus_rnastringset, filepath = file.path(human_path_mirna, "non_nucleus.fa"))

nucleus_motif <- read_tsv(file = file.path(human_path_mirna, "nucleus_motif.tsv"), col_names = F) %>% dplyr::rename(rowname = X1)
nucleus_rnastringset %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  inner_join(nucleus_motif, by = "rowname") %>% 
  mutate(subx = str_sub(x, start = X3, end = X3 + 5)) -> nucleus_ready
DNAStringSet(RNAStringSet(x = nucleus_ready$subx)) -> nucleus_ready_set
names(nucleus_ready_set) <- nucleus_ready$rowname
pwm <- PWM(nucleus_ready_set, type = "prob")
p <- makePWM(consensusMatrix(nucleus_ready_set, as.prob = T) %>% head(4))
p
seqLogo(p)
write_rds(pwm, path = file.path(human_path_mirna, "rds_04_nucleus_pwm.rds.gz"), compress = "gz")



