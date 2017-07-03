#!/usr/bin/R
library(RMySQL)
library(stringr)
library(jsonlite)
library(tidyverse)

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
                    str_c(arg_encode, 'snorna_expr_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}

res_snorna_expr <- dbSendQuery(con, statement = paste(
                       "SELECT * FROM ", table_snorna_expression,
                       "WHERE snorna IN (", str_c("'", str_c(pair_info$snorna, collapse="', '"), "'"),
                       ")", sep=" "))
snorna_expr <- dbFetch(res_snorna_expr, n = -1) %>% as_tibble()
dbClearResult(res_snorna_expr)



snorna_expr %>%
    mutate(snorna_expression = log2(snorna_expression + 1)) %>%
    group_by(snorna, dataset_id) %>%
    summarise(mean_expression = mean(snorna_expression)) %>%
    filter(mean_expression >= 1) %>%
    write_json( path = file.path(resource_jsons, paste(arg_encode, 'snorna_expr_bygene', 'json', sep=".")))
dbDisconnect(con)