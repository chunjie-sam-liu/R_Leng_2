#!/usr/bin/R
library(RMySQL)
library(tidyverse)
library(broom)
library(stringr)
library(jsonlite)
library(latex2exp)

#"snorna_clinical"
#"snorna_dataset"
#"snorna_snorna_expression"
#"snorna_clinical_subtype"
#"snorna_genomic_analysis"
#"snorna_snorna_info"
args <- commandArgs(TRUE)
#args <- 'VTQ2O0VOU0cwMDAwMDIwMDkxMztTTk9SRDQ2I0tJUlA='
arg_encode <- args[1]
arg_decode <- str_split(rawToChar(base64_dec(arg_encode)),"#", simplify=T)

# output dir
root <- "/home/liucj/web/snorna_data_portal"
resource_jsons = file.path(root, 'snorna/resource/jsons')
resource_pngs = file.path(root, 'snorna/resource/pngs')

q_name <- arg_decode[1,1]
dataset_id <- arg_decode[1,2]
database <- paste("snorna_snorna_expression", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")

con <- dbConnect(RMySQL::MySQL(), username = "username", password="password", dbname="db", host="localhost")

clinical <- dbReadTable(con, "snorna_clinical") %>% as_tibble()

res <- dbSendQuery(con, statement = paste(
                        "SELECT * FROM", database,
                        "WHERE snorna = ", paste('"',q_name,'"', sep=''),
                        "AND dataset_id = ",
                        paste('"',dataset_id,'"', sep = ""), sep=" "))

expr <- dbFetch(res, n = -1) %>% as_tibble()
dbClearResult(res)

expr %>%
    dplyr::select(snorna, sample_id, expression = snorna_expression) %>%
    left_join(clinical, by="sample_id") %>%
    dplyr::select(snorna, dataset_id,subtype,stage, expression) %>%
    mutate(expression = log2(expression + 1)) %>%
    drop_na() ->
    diff_subtype

diff_subtype %>%
    group_by(snorna, dataset_id, subtype) %>%
    do(tidy(anova(lm(expression ~ stage, data = .)))) %>%
    filter(term == "stage") %>%
    ungroup() %>%
    mutate(encode = arg_encode)->
    diff_subtype_pvalue
write_json(diff_subtype_pvalue, path=file.path(resource_jsons, paste(arg_encode,'diff_subtype','json',sep='.')))

diff_subtype_distinct <- unique(diff_subtype_pvalue$subtype)
plot_survival <- function(x){
    diff_subtype %>%
    filter(subtype == x) %>%
    ggplot(aes(x = stage, y = expression, color = stage)) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    geom_boxplot( outlier.shape = NA, width = 0.6) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle=45, hjust = 1)
        ) +
    labs(
        y = TeX("snoRNA expression ($log_2(RPKM)$)"),
        x = str_replace_all(x,"_"," ")
        ) -> plot
    ggsave(plot = plot, filename = file.path(resource_pngs, paste(arg_encode,x,"diff_subtype","png",sep = ".")), device = "png")
}
diff_subtype_distinct %>% sapply(plot_survival)

# disconnect
dbDisconnect(con)