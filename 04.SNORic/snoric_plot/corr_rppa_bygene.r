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

# prepare table
table_pair <- paste("snorna_gene_snorna_pair", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")
table_snorna_expression <- paste("snorna_snorna_expression", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")
table_protein_expression <- paste("snorna_protein_expression", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")

#connect database
db <- src_mysql(dbname="db", host="localhost", user = "username", password="password")

#query data
pair_info <- tbl(db, table_pair) %>%
    filter(gene_symbol == q_name) %>%
    collect()

if(nrow(pair_info) < 1){
    write_json(
        NULL,
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'corr_rppa_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}
snorna_expr <- tbl(db, table_snorna_expression) %>%
    collect(n = Inf) %>%
    semi_join(pair_info, by = "snorna")
rppa_expr <- tbl(db, table_protein_expression) %>%
    collect(n = Inf)

# data wrangling
corr_data <-
    rppa_expr%>%
    inner_join(snorna_expr, by = c("snorna", "dataset_id", "sample_id")) %>%
    rename(protein_expression = rna_expression)%>%
    group_by(snorna) %>%
    filter(n() >=10) %>%
    ungroup()

if(nrow(corr_data) < 1){
    write_json(
        NULL,
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'corr_rppa_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}


# keep json
corr_data %>%
    select(dataset_id, snorna, protein, spearman_corr, p_value) %>%
    mutate(p_value = signif(p_value, digits = 3), encode = arg_encode) %>%
    distinct(snorna, .keep_all = T) %>%
    write_json(
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'corr_rppa_bygene', 'json', sep="."))
                    )

# keep picture
plot_point <- function(df){
    df %>%
    ggplot(aes(x = rank(log2(snorna_expression + 1)), y = rank(log2(protein_expression + 1)))) +
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
            y = "protein expression (Rank)"
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
            "corr_rppa_bygene",
            "png",
            sep = "."
            ),
            plot = .out
    ) %>%
    pwalk(ggsave, path = resource_pngs, device = "png")
