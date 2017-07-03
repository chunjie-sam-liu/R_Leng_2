#!/usr/bin/R
library(stringr)
library(jsonlite)
library(tidyverse)

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
table_mrna_expression <- paste("snorna_mrna_expression", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")

db <- src_mysql(dbname="db", host="localhost", user = "username", password="password")

pair_info <- tbl(db, table_pair) %>%
    filter(gene_symbol == q_name) %>%
    collect()

if(nrow(pair_info) < 1){
    write_json(
        NULL,
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'corr_geneexpr_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}


snorna_expr <- tbl(db, table_snorna_expression) %>%
    collect(n = Inf) %>%
    semi_join(pair_info, by = "snorna")


mrna_expr <- tbl(db, table_mrna_expression) %>%
    filter(gene_symbol == q_name) %>%
    collect(n = Inf)

corr_data <-
    mrna_expr%>%
    inner_join(snorna_expr, by = c("snorna", "dataset_id", "sample_id")) %>%
    group_by(snorna) %>%
    filter(n() >=10) %>%
    ungroup()

if(nrow(corr_data) < 1){
    write_json(
        NULL,
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'corr_geneexpr_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}

corr_data %>%
    select(dataset_id, snorna, gene_symbol, spearman_corr, p_value) %>%
    mutate(p_value = signif(p_value, digits = 3), encode = arg_encode) %>%
    distinct(snorna, .keep_all = T) %>%
    write_json(
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'corr_geneexpr_bygene', 'json', sep="."))
                    )

plot_point <- function(df){
    df %>%
    ggplot(aes(x = rank(log2(snorna_expression + 1)), y = rank(log2(rna_expression + 1)))) +
        geom_point(color = "#ca2015") +
        geom_smooth(method=lm, se = FALSE, fullrange=TRUE, color = "#1abc9c") +
        theme_bw() +
        theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black")
            ) +
        labs(
            x = "snoRNA expression (Rank)",
            y = "mRNA expression (Rank)"
            )->p
    return(p)
}

corr_data %>%
    group_by(snorna) %>%
    by_slice(plot_point) %>%
    transmute(
            filename = str_c(
            arg_encode,
            snorna,
            "corr_geneexpr_bygene",
            "png",
            sep = "."
            ),
            plot = .out
    ) %>%
    pwalk(ggsave, path = resource_pngs, device = "png")





