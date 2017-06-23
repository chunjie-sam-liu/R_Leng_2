library(tidyverse)
library(ggthemes)
library(stringr)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "lncrna")

hsa_lncrna_raw <- read_rds(path = file.path(human_path_mirna, "rds_01_hsa_lncrna_raw.rds.gz")) %>% 
  select(
    rlid = RLID,
    oid = `Official ID`, 
    osymbol = `Official Symbol`, 
    subloc = `Subcellular Localization`, 
    tissue = Tissue
  )

hsa_lncrna_raw %>% 
  ggplot(aes(x = reorder(subloc, subloc, function(x)-length(x)))) +
  geom_bar(aes(fill = subloc)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  scale_fill_gdocs() +
  labs(
    x = "Subcellular Localization",
    y = "Count",
    title = "miRNA total data"
  ) +
  theme_gdocs() +
  theme(
    axis.text.x = element_blank()
  )

ggsave(filename = "fig_01_hsa_all_lncrna_subcell_distribution.png", path = human_path_mirna, device = "png")

hsa_lncrna_raw %>% 
  ggplot(aes(x = reorder(tissue, tissue, function(x) -length(x)))) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  theme_gdocs() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Tissue")

ggsave(filename = "fig_02_hsa_all_lrnrna_tissue_distribution.png", path = human_path_mirna, device = "png")


hsa_lncrna_raw %>% 
  group_by(oid, osymbol, subloc) %>% 
  count() %>% 
  ungroup() %>% 
  ggplot(aes(x = as.factor(n))) +
  geom_bar(aes(fill = as.factor(n))) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  theme_gdocs() +
  labs(x = "Same lncRNA and subcell in different tissue") +
  scale_fill_gdocs(name = "NO.")
ggsave(filename = "fig_03_hsa_same_lrnrna_subcell_in_different_tissue.png", path = human_path_mirna, device = "png")


hsa_lncrna_raw %>% 
  distinct(oid, osymbol, subloc) ->
  hsa_lncrna_dedup_tissue


hsa_lncrna_dedup_tissue %>% 
  group_by(oid, osymbol)  %>% 
  mutate(subloc = ifelse(subloc %in% c("Cytoplasm", "Cytosol", "Endoplasmic reticulum", "Exosome", "Mitochondrion", "Ribosome", "Circulating"), "Cytosol", subloc)) %>% 
  mutate(subloc = ifelse(subloc %in% c("Nucleus", "Nucleolus"), "Nucleus", subloc)) %>% 
  ungroup() %>% 
  distinct() ->
  hsa_lncrna_dedup_tissue5

hsa_lncrna_dedup_tissue5 %>% ggplot(aes(x = subloc)) + geom_bar(aes(fill = subloc))  + geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + theme_gdocs() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hsa_lncrna_dedup_tissue5 %>% 
  group_by(oid, osymbol) %>% 
  count() %>% 
  ungroup() %>% 
  ggplot(aes(x = as.factor(n))) +
  geom_bar(aes(fill = as.factor(n))) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  scale_fill_gdocs(name = "NO.") +
  theme_gdocs() +
  labs(x = "The same lncRNA in different subcellular compartment")

ggsave(filename = "fig_04_hsa_same_lrnrna_different_subcell.png", path = human_path_mirna, device = "png")

hsa_lrnrna_dedup_tissue_monoplace <-
  hsa_lncrna_dedup_tissue5 %>% 
  group_by(oid, osymbol) %>% 
  filter(n() == 1) %>% 
  ungroup()

hsa_lrnrna_dedup_tissue_monoplace %>% 
  ggplot(aes(x = as.factor(subloc))) +
  geom_bar(aes(fill = subloc)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  scale_fill_gdocs() +
  theme_gdocs()

ggsave(filename = "fig_07_hsa_sequence_length_dist.png", device = "png", path = human_path_mirna)
write_rds(hsa_lrnrna_dedup_tissue_monoplace, path = file.path(human_path_mirna, "rds_02_hsa_lncrna_clean.rds.gz"), compress = "gz")


















