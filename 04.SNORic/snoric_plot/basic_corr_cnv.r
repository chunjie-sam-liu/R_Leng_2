#!/usr/bin/R
library(RMySQL)
library(tidyverse)
library(broom)
library(stringr)
library(jsonlite)
library(broom)

args <- commandArgs(TRUE)
# all
# args <- 'VTg7RU5TRzAwMDAwMjAwNDYzO1NOT1JEMTE4I1RDR0EtQUNDI2FsbA=='
# 0
# args <- 'VTg7RU5TRzAwMDAwMjAwNDYzO1NOT1JEMTE4I1RDR0EtQ0hPTCMw'
# Stage I
# args <- 'VTg7RU5TRzAwMDAwMjAwNDYzO1NOT1JEMTE4I1RDR0EtQUNDI1N0YWdlIEk='
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
subtype_id <- arg_decode[1,3]


database_cnv <- paste("snorna_cnv", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")
database_snorna_expression <- paste("snorna_snorna_expression", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")

# connect db
con <- dbConnect(RMySQL::MySQL(), username = "username", password="password", dbname="db", host="localhost")

# mrna_expression
res_cnv <- dbSendQuery(con, statement = paste(
                       "SELECT * FROM ", database_cnv,
                       "WHERE snorna = ", paste('"',q_name,'"', sep=''),
                       "AND dataset_id = ",
                       paste('"',dataset_id,'"', sep = ""), sep=" "))
expr_cnv <- dbFetch(res_cnv, n = -1) %>% as_tibble()
dbClearResult(res_cnv)

# snorna_expression
res_snorna_expression <- dbSendQuery(con, statement = paste(
                       "SELECT * FROM ", database_snorna_expression,
                       "WHERE snorna = ", paste('"',q_name,'"', sep=''),
                       "AND dataset_id = ",
                       paste('"',dataset_id,'"', sep = ""), sep=" "))
expr_snorna_expression <- dbFetch(res_snorna_expression, n = -1) %>% as_tibble()
dbClearResult(res_snorna_expression)

if(! subtype_id %in% c("all", "0")) {
    res_clinical <- dbSendQuery(con, statement = str_c(
                    "SELECT * FROM snorna_clinical WHERE stage = ",
                    str_c('"', subtype_id, '"', sep = ''),
                    "AND dataset_id = ",
                    str_c('"', dataset_id, '"', sep = ''),
                    sep = " "))
    clinical <- dbFetch(res_clinical, n = -1) %>% as_tibble()
    dbClearResult(res_clinical)
}else{
    res_clinical <- dbSendQuery(con, statement = str_c(
                    "SELECT * FROM snorna_clinical WHERE ",
                    "dataset_id = ",
                    str_c('"', dataset_id, '"', sep = ''),
                    sep = " "))
    clinical <- dbFetch(res_clinical, n = -1) %>% as_tibble()
    dbClearResult(res_clinical)
}


# disconnect
dbDisconnect(con)

if (nrow(clinical)>1) {
    expr_snorna_expression %>%
        semi_join(clinical, by="sample_id") ->
        expr_snorna_expression
}

if(nrow(expr_snorna_expression) < 10 || mean(expr_snorna_expression$snorna_expression) < 1){
    write_json(NULL, path = file.path(resource_jsons, paste(arg_encode, 'corr_cnv', 'json', sep=".")))
    print("Not enough samples OR mean expression less than 1!")
    quit("no", status = 0)
}


corr_data <-
    expr_snorna_expression%>%
    inner_join(expr_cnv, by = c("snorna", "dataset_id", "sample_id")) %>%
    drop_na(snorna_expression, copy)


corr_data %>%
    group_by(snorna, dataset_id ) %>%
    do(tidy(cor.test(.$snorna_expression,.$copy, method = 'spearman'))) %>%
    ungroup() %>%
    mutate(encode = arg_encode)->
    corr_data_json
write_json(corr_data_json, path = file.path(resource_jsons, paste(arg_encode, 'corr_cnv', 'json', sep=".")))


corr_data %>%
    ggplot(aes(x = rank(log2(snorna_expression + 1)), y = rank(copy))) +
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
        y = "Copy (Rank)"
        ) -> p
ggsave(filename = file.path(resource_pngs, paste(arg_encode,"corr_cnv","png",sep = ".")), device = "png", plot = p )


