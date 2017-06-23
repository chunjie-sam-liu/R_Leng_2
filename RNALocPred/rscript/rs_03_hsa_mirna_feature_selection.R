library(tidyverse)
library(ggthemes)
library(stringr)
library(multidplyr)
library(Biostrings)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")

hsa_mirna <- 
  read_rds(path = file.path(human_path_mirna, "rds_03_hsa_mirna_clean_seq_ready.rds.gz")) %>% 
  drop_na() %>% 
  mutate(seq = str_replace_all(str_to_upper(seq), "T", "U")) %>% 
  drop_na() %>% 
  distinct(oid, .keep_all = T)

# GC content
hsa_mirna_seq_stat <-
  hsa_mirna %>% 
  mutate(
    A = str_count(seq, "A") ,
    G = str_count(seq, "G") ,
    C = str_count(seq, "C") ,
    U = str_count(seq, "U") 
  ) %>% 
  gather(nt,pt, A:U) %>% 
  mutate(
    subloc = forcats::as_factor(subloc), 
    type = forcats::as_factor(type),
    nt = forcats::as_factor(nt)
    ) 

hsa_mirna_seq_stat %>% 
  ggplot(aes(x = subloc, y = pt, fill = nt)) +
  geom_bar(stat = "identity") +
  theme_gdocs()

# encoding based on the sequence
# main feature is the total length to encoding and 
encodeKMerSeq <- function (k, rnaStringSet) 
{
  lookupTable <- diag(4^k)
  row.names(lookupTable) <- mkAllStrings(c("A", "C", "G", "U"), k)
  l <- nchar(toString(rnaStringSet[1]))
  n <- (l - k + 1) * (4^k)
  m <- length(rnaStringSet)
  featureVector <- matrix(0L, m, n)
  for (i in 1:length(rnaStringSet)) {
    features <- c()
    seq <- toupper(toString(rnaStringSet[i]))
    for (j in 1:(nchar(seq) - k + 1)) {
      if (is.na(match(substr(seq, j, j + k - 1), row.names(lookupTable)))) {
      }
      else {
        featureVector[i, ((j - 1) * (4^k) + 1):(j * (4^k))] <- lookupTable[substr(seq, j, j + k - 1), ]
      }
    }
  }
  row.names(featureVector) <- names(rnaStringSet)
  return(featureVector)
}

#
# test for seed region 2-8
# 
hsa_mirna_set <- RNAStringSet(x = hsa_mirna$seq, start = 2, end = 8)
names(hsa_mirna_set) <- hsa_mirna$names

# make features vector
hsa_mirna_feature <- 
  as_tibble(rownames_to_column(as.data.frame(encodeKMerSeq(1, hsa_mirna_set)), var = "names")) %>% 
  inner_join(hsa_mirna, by = "names") %>% 
  mutate(gc_content = str_count(seq, "[CG]") / seq_len) %>% 
  dplyr::select(-c(names, oid, osymbol, type, seq, seq_len))
write_rds(hsa_mirna_feature, path = file.path(human_path_mirna, "rds_04_hsa_mirna_feature.rds.gz"), compress = "gz")













