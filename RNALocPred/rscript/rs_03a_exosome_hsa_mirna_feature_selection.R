library(tidyverse)
library(ggthemes)
library(stringr)
library(multidplyr)
library(Biostrings)
library(broom)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")

hsa_mirna <-
  read_rds(path = file.path(human_path_mirna, "rds_03_hsa_mirna_clean_seq_ready.rds.gz")) %>%
  drop_na() %>%
  mutate(seq = str_replace_all(str_to_upper(seq), "T", "U")) %>%
  drop_na() %>%
  distinct(oid, .keep_all = T) %>%
  mutate(subloc = ifelse(
    subloc %in% c("Exosome", "Extracellular vesicle", "Microvesicle"),
    "EV",
    "non_EV"
  ))


#----------------------nucleotide composition significance-------------------
seq_composition <- function(seq) {
  seq %>% RNAString() %>% DNAString() -> d_string
  c(
    oligonucleotideFrequency(
      x = d_string,
      as.prob = T,
      width = 1
    ),
    oligonucleotideFrequency(
      x = d_string,
      as.prob = T,
      width = 2
    )
  ) -> comp
  comp %>% enframe() %>% spread(name, value)
}

hsa_mirna %>%
  select(osymbol, subloc, seq, names) %>%
  mutate(comp = map(.x = seq, .f = seq_composition)) %>%
  unnest() -> hsa_mirna_comp

colnames(hsa_mirna_comp)[-c(1, 2, 3, 4)] %>%
  map( ~ t.test(get(.x) ~ subloc, data = hsa_mirna_comp) %>% glance()) %>%
  enframe() %>%
  unnest() %>%
  mutate(
    name = colnames(hsa_mirna_comp)[-c(1, 2, 3, 4)],
    p.value = signif(-log10(p.value), digits = 3),
    mean_ev = signif(estimate1, 3),
    mean_non_ev = signif(estimate2, 3),
    ratio = log2(mean_ev / mean_non_ev)
  ) %>%
  select(name, p.value, mean_ev, mean_non_ev, ratio) %>%
  arrange(-p.value) -> significant_nt_count
# GC,A,AA,GT filterout

# significance of count
significant_nt_count %>%
  ggplot(aes(x = ratio, y = p.value)) +
  geom_point() +
  geom_point(data = significant_nt_count %>% filter(abs(ratio) > 0.2),
             color = "red") +
  geom_text(
    aes(label = name),
    hjust = 0,
    nudge_x = 0.03,
    data = significant_nt_count %>% filter(abs(ratio) > 0.2),
    color = "red"
  ) +
  labs(x = "Ratio of EV vs non-EV",
       y = "P-value (log10)",
       title = "Significant nucleotide counts. EV vs non-EV") +
  theme_gdocs()

ggsave(filename = "fig_12_hsa_feature_selection_significant_nucleotide.png", device = "png", path = human_path_mirna)

hsa_mirna %>% select(osymbol, subloc, seq) -> hsa_mirna_pos

pos_test <- list()
for (i in 1:16) {
  chisq.test(
    hsa_mirna_pos$subloc,
    str_sub(hsa_mirna_pos$seq, start = i, end = i),
    simulate.p.value = T,
    B = 10000
  ) -> r_t
  pos_test <- c(pos_test, list(r_t))
}

pos_test %>% 
  map(tidy) %>% 
  enframe() %>% 
  unnest() %>%
  mutate(p.value = signif(-log10(p.value)), digits = 3) %>% 
  ggplot(aes(x = name, y = p.value)) +
  geom_rect(aes(xmin = 2, xmax = 8, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.05) + 
  geom_line() + 
  geom_point() + 
  scale_x_discrete(limits = 1:16) +
  theme_gdocs() +
  labs(
    x = "miRNA Position",
    y = "P-value for base composition (log10)"
  ) 
ggsave(filename = "fig_12_hsa_feature_selection_position_difference.png", device = "png", path = human_path_mirna)



# sequence encoding------------------------------------
encodeKMerSeq <- function (rnaStringSet, k)
{
  lookupTable <- diag(4 ^ k)
  row.names(lookupTable) <- mkAllStrings(c("A", "C", "G", "U"), k)
  l <- nchar(toString(rnaStringSet[1]))
  n <- (l - k + 1) * (4 ^ k)
  m <- length(rnaStringSet)
  featureVector <- matrix(0L, m, n)
  for (i in 1:length(rnaStringSet)) {
    features <- c()
    seq <- toupper(toString(rnaStringSet[i]))
    for (j in 1:(nchar(seq) - k + 1)) {
      if (is.na(match(substr(seq, j, j + k - 1), row.names(lookupTable)))) {
        
      }
      else {
        featureVector[i, ((j - 1) * (4 ^ k) + 1):(j * (4 ^ k))] <-
          lookupTable[substr(seq, j, j + k - 1),]
      }
    }
  }
  row.names(featureVector) <- names(rnaStringSet)
  return(featureVector)
}

hsa_mirna_set <- RNAStringSet(x = hsa_mirna$seq, start = 2, end = 8)
names(hsa_mirna_set) <- hsa_mirna$names

# create seed sequence as feature
seed_feature <-
  hsa_mirna_set %>%
  encodeKMerSeq(1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "names") %>%
  as_tibble()

# ---------------------motif pwm------------------------------------
# Add pwm to features
pwm_score <- function(seq, pwm) {
  seq %>% RNAString() %>% DNAString() -> d_string
  matchPWM(pwm = pwm,
           subject = d_string,
           with.score = T) -> hits
  s <- mcols(hits)$score
  if (length(s)) {
    max(s)
  } else{
    0
  }
}

# the pwm was calculate from rs_03a.R
# hnrnpq pwm
hnrnpq_pwm <- read_rds(path = file.path(human_path_mirna, "rds_05_hnrnpq_pwm.rds.gz"))
# hnrnpa2b1 pwm
hnrnpa2b1_pwm <- read_rds(path = file.path(human_path_mirna, "rds_05_hnrnpa2b1_pwm.rds.gz"))

hsa_mirna %>%
  select(subloc, names, seq, seq_len) %>%
  mutate(hnrnpq_pwm_score = map_dbl(.x = seq, .f = pwm_score, pwm = hnrnpq_pwm)) %>%
  mutate(hnrnpa2b1_pwm_score = map_dbl(.x = seq, .f = pwm_score, pwm = hnrnpa2b1_pwm)) %>%
  mutate(gc_content = str_count(seq, "[CG]") / seq_len) -> pwm_feature

# mirna composition features
comp_feature <- hsa_mirna_comp %>% select(-osymbol, - subloc, -seq)


features <-
  seed_feature %>%
  inner_join(pwm_feature, by = "names") %>%
  inner_join(comp_feature, by = "names") %>% 
  select(-seq, -names, -GC, -A, -AA, -GT)

write_rds(
  features,
  path = file.path(
    human_path_mirna,
    "rds_06_hsa_mirna_features_exosome.rds.gz"
  ),
  compress = "gz"
)
