library(rGADEM)
library(tidyverse)
library(stringr)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")
hsa_mirna_clean <- read_rds(path = file.path(human_path_mirna, "rds_03_hsa_mirna_clean_seq_ready.rds.gz"))

hsa_mirna_clean %>% 
  filter(subloc == "Exosome", type == "Mature") %>% 
  drop_na() %>% 
  mutate(seq = str_replace_all(seq, "T", "U")) %>%
  distinct(seq, .keep_all = T) -> exosome

exosome_rnastringset <- RNAStringSet(x = exosome$seq) 
names(exosome_rnastringset) <- exosome$osymbol
writeXStringSet(x = exosome_rnastringset, filepath = file.path(human_path_mirna, "exosome.fa"))


hsa_mirna_clean %>% 
  filter(subloc == "Mitochondrion") %>%
  drop_na() %>% 
  mutate(seq = str_replace_all(seq, "T", "U")) %>%
  distinct(seq, .keep_all = T) -> mitochondrion

mitochondrion_rnastringset <- RNAStringSet(x = mitochondrion$seq) 
names(mitochondrion_rnastringset) <- mitochondrion$osymbol
writeXStringSet(x = mitochondrion_rnastringset, filepath = file.path(human_path_mirna, "mitochondrion.fa"))

gadem <- GADEM(Sequences = exosome_rnastringset, verbose = T)
#

# miRBase mature mirna
mirbase_rnastringset <- readRNAStringSet(filepath = "/home/cliu18/liucj/reference/miRBase/mature.fa")

mirbase_rnastringset %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  separate(rowname, into = c("symbol", "acc", "species1", "species2", "name"), sep = " ") -> mirbase

hnRNPQ <- unique(c("mmu-let-7b-5p", "mmu-miR-127-3p", "mmu-miR-127-5p", "mmu-miR-15b-5p", "mmu-miR-1934-3p", "mmu-miR-1941-5p", "mmu-miR-194-2-3p", "mmu-miR-195a-3p", "mmu-miR-208b-3p", "mmu-miR-3099-3p", "mmu-miR-3110-5p", "mmu-miR-328-5p", "mmu-miR-3470a", "mmu-miR-3470b", "mmu-miR-3473a", "mmu-miR-3473b", "mmu-miR-3473g", "mmu-miR-3620-5p", "mmu-miR-369-3p", "mmu-miR-433-3p", "mmu-miR-493-5p", "mmu-miR-500-3p", "mmu-miR-504-3p", "mmu-miR-5108", "mmu-miR-5119", "mmu-miR-5128", "mmu-miR-5128", "mmu-miR-6337", "mmu-miR-6538", "mmu-miR-671-5p", "mmu-miR-690", "mmu-miR-690", "mmu-miR-6909-5p", "mmu-miR-6930-5p", "mmu-miR-6935-5p", "mmu-miR-6959-5p", "mmu-miR-6981-5p", "mmu-miR-7015-5p", "mmu-miR-7015-5p", "mmu-miR-7056-5p", "mmu-miR-7221-3p", "mmu-miR-7687-5p", "mmu-miR-133a-3p", "mmu-miR-150-5p", "mmu-miR-193b-5p", "mmu-miR-23a-5p", "mmu-miR-3099-3p", "mmu-miR-3102-5p.2-5p", "mmu-miR-3473b", "mmu-miR-3620-5p", "mmu-miR-409-3p", "mmu-miR-411-3p", "mmu-miR-451a", "mmu-miR-485-5p", "mmu-miR-5124b", "mmu-miR-6337", "mmu-miR-6715-3p", "mmu-miR-671-5p", "mmu-miR-6896-5p", "mmu-miR-690", "mmu-miR-6959-5p", "mmu-miR-6986-5p", "mmu-miR-6994-5p", "mmu-miR-6995-5p", "mmu-miR-706", "mmu-miR-7068-5p", "mmu-miR-7578", "mmu-miR-760-3p", "mmu-miR-7648-3p", "mmu-miR-7667-5p", "mmu-miR-7667-5p", "mmu-miR-7684-5p", "mmu-miR-7687-5p"))

mirbase %>% filter(symbol %in% hnRNPQ) %>% mutate(x = str_replace_all(x, "U", "T"))   -> hnRNPQ_seq

DNAStringSet(x = hnRNPQ_seq$x) -> hnRNPQ_seq_set
names(hnRNPQ_seq_set) <- hnRNPQ_seq$symbol
# gadem <- GADEM(Sequences = hnRNPQ_seq_set, verbose = T)
writeXStringSet(x = hnRNPQ_seq_set, filepath = file.path(human_path_mirna, "hnrnpq.fa"))
# can't find motif using rGADEM
hnRNPQ_motif <- read_tsv(file.path(human_path_mirna, "hnRNPQ_motif.tsv"),col_names = F) %>% rename(rowname = X1)
hnRNPQ_seq_set %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% inner_join(hnRNPQ_motif, by = "rowname") %>% mutate(subx = str_sub(x, start = X3, end = X3 + 5)) ->hnRNPQ_ready

DNAStringSet(x = hnRNPQ_ready$subx) -> hnRNPQ_ready_set
names(hnRNPQ_ready_set) <- hnRNPQ_ready$rowname
##
# get pwm Cell Reports paper hnRNPQ
# #
pwm <- PWM(hnRNPQ_ready_set, type = "prob")
write_rds(pwm, path = file.path(human_path_mirna, "rds_05_hnrnpq_pwm.rds.gz"), compress = "gz")
p <- makePWM(consensusMatrix(DNAStringSet(hnRNPQ_ready_set), as.prob = T) %>% head(4))
seqLogo(p)





# hnrnpa2b1
hnrnpa2b1 <- c("hsa-miR-654-5p", "hsa-miR-339-3p", "hsa-miR-769-3p", "hsa-miR-135a-3p", "hsa-miR-638", "hsa-miR-1225-5p", "hsa-miR-877-5p", "hsa-miR-134-5p", "hsa-miR-765", "hsa-miR-575", "hsa-miR-887-3p", "hsa-miR-188-5p", "hsa-miR-513a-5p", "hsa-miR-513b-5p", "hsa-miR-513c-5p", "hsa-miR-630", "hsa-miR-583", "hsa-miR-483-5p", "hsa-miR-125a-3p", "hsa-miR-125b-1-3p", "hsa-miR-1224-5p", "hsa-miR-601", "hsa-miR-198", "hsa-miR-193b-5p", "hsa-miR-422a", "hsa-miR-1226-5p", "hsa-miR-671-5p", "hsa-miR-451a", "hsa-miR-520b", "hsa-miR-520e")
mirbase  %>% filter(symbol %in% hnrnpa2b1) %>% mutate(x = str_replace_all(x, "U", "T"))  -> hnrnpa2b1_seq
  
DNAStringSet(x = hnrnpa2b1_seq$x)  -> hnrnpa2b1_seq_set
names(hnrnpa2b1_seq_set) <- hnrnpa2b1_seq$symbol
writeXStringSet(x = hnrnpa2b1_seq_set, filepath = file.path(human_path_mirna, "hnrnpa2b1.fa"))

hnrnpa2b1_motif <- read_tsv(file.path(human_path_mirna, "hnRNPA2B1_motif.tsv"),col_names = F) %>% rename(rowname = X1)
hnrnpa2b1_seq_set %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% inner_join(hnrnpa2b1_motif, by = "rowname") %>% mutate(subx = str_sub(x, start = X3, end = X3 + 5)) ->hnrnpa2b1_ready

DNAStringSet(x = hnrnpa2b1_ready$subx) -> hnrnpa2b1_ready_set
names(hnrnpa2b1_ready_set) <- hnrnpa2b1_ready$rowname
##
# get pwm Cell Reports paper hnRNPQ
# #
pwm <- PWM(hnrnpa2b1_ready_set, type = "prob")
write_rds(pwm, path = file.path(human_path_mirna, "rds_05_hnrnpa2b1_pwm.rds.gz"), compress = "gz")
p <- makePWM(consensusMatrix(DNAStringSet(hnrnpa2b1_ready_set), as.prob = T) %>% head(4))
seqLogo(p)



mito_motif <- read_tsv(file.path(human_path_mirna, "mito_motif.tsv"),col_names = F) %>% rename(rowname = X1)
DNAStringSet(mitochondrion_rnastringset) %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% inner_join(mito_motif, by = "rowname") %>% mutate(subx = str_sub(x, start = X3, end = X3 + 5)) ->mito_ready
  
DNAStringSet(x = mito_ready$subx) -> mito_ready_set
names(mito_ready_set) <- mito_ready$rowname
##
# get pwm Cell Reports paper hnRNPQ
# #
pwm <- PWM(mito_ready_set, type = "prob")
write_rds(pwm, path = file.path(human_path_mirna, "rds_05_mito_pwm.rds.gz"), compress = "gz")
p <- makePWM(consensusMatrix(DNAStringSet(mito_ready_set), as.prob = T) %>% head(4))
seqLogo(p)


for(i in seq_along(mitochondrion_rnastringset)) {
  matchPWM(pwm,
           DNAStringSet(mitochondrion_rnastringset)[[i]],
           with.score = T) -> hits
  s <-
    mcols(hits)$score

  # print(paste(i, max(s)))
  if (length(s) > 0) {
    print(max(s))
  } else{
    print(0)
  }
}
  