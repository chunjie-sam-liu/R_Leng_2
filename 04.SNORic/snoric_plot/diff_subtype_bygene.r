#!/usr/bin/R
library(RMySQL)
library(stringr)
library(jsonlite)
library(tidyverse)
library(broom)
library(latex2exp)

#snorna_gene_snorna_pair

args <- commandArgs(TRUE)
#args <- 'R0FETDEjVENHQS1BQ0M='
arg_encode <- args[1]
arg_decode <- str_split(rawToChar(base64_dec(arg_encode)),"#", simplify=T)

# output dir
root <- "/home/liucj/web/snorna_data_portal"
resource_jsons = file.path(root, 'snorna/resource/jsons')
#resource_jsons = "./"
resource_pngs = file.path(root, 'snorna/resource/pngs')
#resource_pngs = "./"

q_name <- arg_decode[1,1]
dataset_id <- arg_decode[1,2]
table_pair <- paste("snorna_gene_snorna_pair", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")
table_snorna_expression <- paste("snorna_snorna_expression", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")

con <- dbConnect(RMySQL::MySQL(), username = "username", password="password", dbname="db", host="localhost")

res_pair <- dbSendQuery(con, statement = paste(
                       "SELECT snorna FROM ", table_pair,
                       "WHERE gene_symbol = ", paste('"',q_name,'"', sep=''),
                       "AND dataset_id = ",
                       paste('"',dataset_id,'"', sep = ""), sep=" "))
pair_info <- dbFetch(res_pair, n = -1) %>% as_tibble()
dbClearResult(res_pair)


if(nrow(pair_info) < 1){
    write_json(
        NULL,
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'diff_subtype_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}


res_snorna_expr <- dbSendQuery(con, statement = paste(
                       "SELECT * FROM ", table_snorna_expression,
                       "WHERE snorna IN (", str_c("'", str_c(pair_info$snorna, collapse="', '"), "'"),
                       ")", sep=" "))
snorna_expr <- dbFetch(res_snorna_expr, n = -1) %>% as_tibble()
dbClearResult(res_snorna_expr)

clinical <- dbReadTable(con, "snorna_clinical") %>% as_tibble()

snorna_expr %>%
    dplyr::select(snorna, sample_id, expression = snorna_expression) %>%
    left_join(clinical, by="sample_id") %>%
    dplyr::select(snorna, dataset_id,subtype,stage, expression) %>%
    mutate(
        expression = log2(expression + 1),
        stage = as.factor(stage)) %>%
    drop_na() %>%
    group_by(snorna) %>%
    filter(n() >=10)->
    diff_subtype


if(nrow(diff_subtype) < 1){
    write_json(
        NULL,
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'diff_subtype_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}


diff_subtype %>%
    group_by(snorna, dataset_id, subtype) %>%
    do(tidy(anova(lm(expression ~ stage, data = .)))) %>%
    filter(term == "stage") %>%
    ungroup() %>%
    mutate(encode = arg_encode)->
    diff_subtype_pvalue
write_json(diff_subtype_pvalue, path=file.path(resource_jsons, paste(arg_encode,'diff_subtype_bygene','json',sep='.')))

plot_diff_subtype <- function(df){
    ggplot(data = df, aes(x = stage, y = expression, color = stage)) +
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
        y = TeX("snoRNA expression ($log_2(RPKM)$)")
        )
}

diff_subtype %>%
    group_by(snorna, subtype) %>%
    by_slice(plot_diff_subtype) ->
    plots

paths <- str_c(
            arg_encode,
            plots$snorna,
            plots$subtype,
            "diff_subtype_bygene",
            "png",
            sep = "."
            )

pwalk(list(paths, plots$.out), ggsave, device = "png", path= resource_pngs)
dbDisconnect(con)