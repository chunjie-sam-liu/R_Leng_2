library(tidyverse)
library(ggthemes)

source_path <- "/home/cliu18/liucj/reference/RNALocate"

data_file <- file.path(source_path, "RNALocate_all_data.txt")

script_path <- "/extraspace/liucj/github/RstudioWithGit/RNALocPred"


data <-
  read_tsv(data_file, col_types = cols(`Official ID` = col_character()))

data %>% 
  # group_by(Organism) %>% 
  # filter(n() > 300) %>% 
  # ungroup() %>% 
  ggplot(aes(x = Organism)) +
  geom_bar(aes(fill = `RNA Category`)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.3) + 
  theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  scale_fill_gdocs() -> 
  p

# figure_1
ggsave(filename = "01.fig_all_rnalocate_data_distribution.png", plot = p, device = "png", path = script_path)


# --------------------------for human--------------------
human <- 
  data %>% filter(
    Organism == "Homo sapiens"
  )

write_rds(human, path = file.path(script_path, "human/01.data_hsa_rna.rds.gz"), compress = "gz")

human %>%
  ggplot(aes(x = `RNA Category`)) +
  geom_bar(aes(fill = `Subcellular Localization`)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.3) +
  scale_fill_gdocs() +
  theme_gdocs() ->
  plot

ggsave(filename = "02.fig_hsa_rna_category_distribution.pdf", plot = plot, device = "pdf", path = file.path(script_path, "human"))
ggsave(filename = "02.fig_hsa_rna_category_distribution.png", plot = plot, device = "png", path = file.path(script_path, "human"))

# --------------------focus on the lncrna, mirna and mrna-----------------
lncrna_mirna_mrna_human <-
  human %>% filter(`RNA Category` %in% c("lncRNA", "miRNA", "mRNA"))

lncrna_mirna_mrna_human %>%
  ggplot(aes(x = `RNA Category`)) +
  geom_bar(aes(fill = `Subcellular Localization`)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.3)
# ---------------------------snorna---------------------
human %>% filter(`RNA Category` %in% c("snoRNA")) %>% write_rds(
  path  = file.path(script_path, "human/snorna/rds_01_hsa_snorna_raw.rds.gz"),
  compress = "gz"
)


# ------------------save human lncrna mirna and mrna for further analysis--------------------
write_rds(lncrna_mirna_mrna_human, path = file.path(script_path, 'human/03.data_hsa_lncrna_mirna_mrna.rds.gz'), compress = 'gz')


# 1. considering mirna
# rm mirna residing more than one subcellular

lncrna_mirna_mrna_human %>% 
  filter(`RNA Category` == "miRNA") %>% 
  write_rds(path = file.path(script_path, "human/mirna/01.data_hsa_mirna_raw.rds.gz"), compress = "gz")

lncrna_mirna_mrna_human %>% 
  filter(`RNA Category` == "lncRNA") %>% 
  write_rds(path = file.path(script_path, "human/lncrna/01.data_hsa_lncrna_raw.rds.gz"), compress = "gz")

lncrna_mirna_mrna_human %>% 
  filter(`RNA Category` == "mRNA") %>% 
  write_rds(path = file.path(script_path, "human/mrna/01.data_hsa_mrna_raw.rds.gz"), compress = "gz")
  

lncrna_mirna_mrna_human %>%
  filter(`RNA Category` == "miRNA") %>%
  group_by(`Official ID`) %>%
  filter(n() < 2) %>%
  ungroup() ->
  mirna


mirna %>% 
  ggplot(aes(x = reorder(`Subcellular Localization`,`Subcellular Localization`, function(x)-length(x)))) +
  geom_bar(aes(fill = `Subcellular Localization`)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3) + 
  scale_fill_gdocs() +
  labs(
    x = "Subcellular Localization",
    y = "Count",
    title = "miRNA filtered data"
  ) +
  theme_gdocs() 
write_rds(mirna, path = file.path(script_path, 'mirna_rnalocate_human.rds.gz'), compress = 'gz')


# -----------------mouse ---------------------------------
mouse_rna <- 
  data %>% 
  filter(Organism == "Mus musculus")
write_rds(mouse_rna, path = file.path(script_path, "mouse/01.data_mmu_rna.rds.gz"), compress = "gz")

mouse_rna %>%
  ggplot(aes(x = `RNA Category`)) +
  geom_bar(aes(fill = `Subcellular Localization`)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.3) +
  theme_gdocs() ->
  plot
ggsave(filename = "02.fig_mmu_rna_category_distribution.pdf", plot = plot, device = "pdf", path = file.path(script_path, "mouse"))
ggsave(filename = "02.fig_mmu_rna_category_distribution.png", plot = plot, device = "png", path = file.path(script_path, "mouse"))


mouse_lncrna_mirna_mrna <-
  mouse_rna %>% filter(`RNA Category` %in% c("lncRNA", "miRNA", "mRNA"))

mouse_lncrna_mirna_mrna %>%
  ggplot(aes(x = `RNA Category`)) +
  geom_bar(aes(fill = `Subcellular Localization`)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.3)

# ------------------save human lncrna mirna and mrna for further analysis--------------------
write_rds(mouse_lncrna_mirna_mrna, path = file.path(script_path, 'mouse/03.data_mmu_lncrna_mirna_mrna.rds.gz'), compress = 'gz')


# 1. considering mirna
# rm mirna residing more than one subcellular

mouse_lncrna_mirna_mrna %>% 
  filter(`RNA Category` == "miRNA") %>% 
  write_rds(path = file.path(script_path, "mouse/mirna/01.data_mmu_mirna_raw.rds.gz"), compress = "gz")

mouse_lncrna_mirna_mrna %>% 
  filter(`RNA Category` == "lncRNA") %>% 
  write_rds(path = file.path(script_path, "mouse/lncrna/01.data_mmu_lncrna_raw.rds.gz"), compress = "gz")

mouse_lncrna_mirna_mrna %>% 
  filter(`RNA Category` == "mRNA") %>% 
  write_rds(path = file.path(script_path, "mouse/mrna/01.data_mmu_mrna_raw.rds.gz"), compress = "gz")


# 
save(list = ls(), file = file.path(script_path, 'rnalocate_data.rda'))



