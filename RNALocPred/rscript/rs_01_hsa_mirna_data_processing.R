library(tidyverse)
library(ggthemes)
library(stringr)

root_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"
human_path <- file.path(root_path, "human")
human_path_mirna <- file.path(human_path, "mirna")

# ---------------------Load raw data-----------------------------
hsa_mirna_raw <- 
  read_rds(path = file.path(human_path_mirna, "rds_01_hsa_mirna_raw.rds.gz")) %>% 
  select(
    rlid = RLID,
    oid = `Official ID`, 
    osymbol = `Official Symbol`, 
    subloc = `Subcellular Localization`, 
    tissue = Tissue
  ) 

# ----------------------mirna distribution in subcellular localization-------------------
hsa_mirna_raw %>% 
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

ggsave(filename = "fig_01_hsa_all_mirna_subcell_distribution.png", path = human_path_mirna, device = "png")

# -----------------tissue distribution------------------
hsa_mirna_raw %>% 
  ggplot(aes(x = reorder(tissue, tissue, function(x) -length(x)))) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  theme_gdocs() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Tissue")

ggsave(filename = "fig_02_hsa_all_mirna_tissue_distribution.png", path = human_path_mirna, device = "png")


# ----------------same mirna in the same subcell in different tissue --------------
hsa_mirna_raw %>% 
  group_by(oid, osymbol, subloc) %>% 
  count() %>% 
  ungroup() %>% 
  ggplot(aes(x = as.factor(n))) +
  geom_bar(aes(fill = as.factor(n))) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  theme_gdocs() +
  labs(x = "Same miRNA and subcell in different tissue") +
  scale_fill_gdocs(name = "NO.")
ggsave(filename = "fig_03_hsa_same_mirna_subcell_in_different_tissue.png", path = human_path_mirna, device = "png")


# -------------dedup the mirna_subcellular in different tissue--------------------
hsa_mirna_raw %>% 
  distinct(oid, osymbol, subloc) ->
  hsa_mirna_dedup_tissue

hsa_mirna_dedup_tissue %>% 
  group_by(oid, osymbol) %>% 
  count() %>% 
  ungroup() %>% 
  ggplot(aes(x = as.factor(n))) +
  geom_bar(aes(fill = as.factor(n))) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  scale_fill_gdocs(name = "NO.") +
  theme_gdocs() +
  labs(x = "The same miRNA in different subcellular compartment")
ggsave(filename = "fig_04_hsa_same_mirna_different_subcell.png", path = human_path_mirna, device = "png")


# ----------------------mirna in multiple place are not suitable for trainning-------------
hsa_mirna_dedup_tissue_monoplace <-
  hsa_mirna_dedup_tissue %>% 
  group_by(oid, osymbol) %>% 
  filter(n() == 1) %>% 
  ungroup()

hsa_mirna_dedup_tissue_monoplace %>% 
  mutate(id_abbr = str_sub(oid, start = 1, end = 4)) %>% 
  ggplot(aes(x = as.factor(id_abbr))) +
  geom_bar(aes(fill = subloc)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  scale_fill_gdocs() +
  theme_gdocs() +
  labs(x = "Sequence Type") +
  scale_x_discrete(labels = c("Precursor", "Mature", "Novel"))

ggsave(filename = "fig_05_hsa_sequence_type.png", path = human_path_mirna, device = "png")

#------------------------exclude precursor----------------------------------------
hsa_mirna_dedup_tissue_monoplace_precursor <- 
  hsa_mirna_dedup_tissue_monoplace %>% 
  mutate(type = str_sub(oid, start = 1, end = 4)) %>% 
  mutate(type = plyr::revalue(type, c(MI00 = "Precursor", MIMA = "Mature", RLGI = "Novel"))) %>% 
  filter(type != "Precursor")

hsa_mirna_dedup_tissue_monoplace_precursor %>% 
  ggplot(aes(x = subloc)) +
  geom_bar(aes(fill = type)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) +
  theme_gdocs() +
  scale_fill_gdocs(name = "Type") +
  labs(x = "Subcellular Localization", y = "Count") +
  theme(axis.text.x = element_text(color = rep(c("green", "black"), times = c(4,2))))
ggsave(filename = "fig_06_hsa_mature_novel_stat.png", path = human_path_mirna, device = "png")
write_rds(hsa_mirna_dedup_tissue_monoplace_precursor, path = file.path(human_path_mirna, "rds_02_hsa_mirna_clean.rds.gz"), compress = "gz")

# 1. I don't believe novel mirna is real.and most of novel mirna are in mitochondrion. Key question is wether we use this novel mirna data?
# 2. Based on the fig 07, I can classify mirna into in_cell/out_cell, in_exosome/out_exosome, or in_mitochondrion/out_mitchondrion.
save(list = ls(), file = file.path(human_path_mirna, "rda_01_mirna_raw_data_scale.rda.gz"), compress = "gzip")



