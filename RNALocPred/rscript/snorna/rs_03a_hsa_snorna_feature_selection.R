library(tidyverse)
library(ggthemes)
library(stringr)
library(multidplyr)
library(Biostrings)
library(broom)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "snorna")


hsa_mirna <-
  read_rds(path = file.path(human_path_mirna, "rds_03_hsa_snorna_clean_seq_ready.rds.gz")) %>%
  mutate(seq = str_replace_all(str_to_upper(seq), "T", "U")) %>%
  distinct(oid, .keep_all = T) %>%
  mutate(subloc = ifelse(subloc %in% c("Nucleolus", "Nucleus"), "Nucleus", "non-Nucleus")) %>% 
  drop_na()

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
  select(oid, osymbol, subloc, seq, names) %>%
  mutate(comp = map(.x = seq, .f = seq_composition)) %>%
  unnest() -> hsa_mirna_comp

colnames(hsa_mirna_comp)[-c(1, 2, 3, 4,5)] %>%
  map( ~ t.test(get(.x) ~ subloc, data = hsa_mirna_comp) %>% glance()) %>%
  enframe() %>%
  unnest() %>%
  mutate(
    name = colnames(hsa_mirna_comp)[-c(1, 2, 3, 4, 5)],
    p.value = signif(-log10(p.value), digits = 3),
    mean_ev = signif(estimate1, 3),
    mean_non_ev = signif(estimate2, 3),
    ratio = log2(mean_ev / mean_non_ev)
  ) %>%
  select(name, p.value, mean_ev, mean_non_ev, ratio) %>%
  arrange(-p.value) -> significant_nt_count


significant_nt_count %>%
  ggplot(aes(x = ratio, y = p.value)) +
  geom_point() +
  geom_point(data = significant_nt_count %>% filter(abs(ratio) > 0.2, p.value > -log10(0.05)),
             color = "red") +
  geom_text(
    aes(label = name),
    hjust = 0,
    nudge_x = 0.03,
    data = significant_nt_count %>% filter(abs(ratio) > 0.2, p.value > -log10(0.05)),
    color = "red"
  ) +
  labs(x = "Ratio of EV vs non-EV",
       y = "P-value (log10)",
       title = "Significant nucleotide counts. EV vs non-EV") +
  theme_gdocs()
ggsave(filename = "fig_12_hsa_snorna_feature_selection_significant_nucleotide.png", device = "png", path = human_path_mirna)


hsa_mirna %>% select(osymbol, subloc, seq) -> hsa_mirna_pos














