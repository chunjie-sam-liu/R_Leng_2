#!/usr/bin/R
library(RMySQL)
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(stringr)
library(jsonlite)


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

database <- paste("snorna_snorna_expression", str_to_lower(str_replace(dataset_id, "TCGA-", "")), sep = "_")

# connect mysql
con <- dbConnect(RMySQL::MySQL(), username = "username", password="password", dbname="db", host="localhost")
res <- dbSendQuery(con, statement = paste(
                        "SELECT * FROM", database,
                        "WHERE snorna = ", paste('"',q_name,'"', sep=''),
                        "AND dataset_id = ",
                        paste('"',dataset_id,'"', sep = ""), sep=" "))
expr <- dbFetch(res, n = -1) %>% as_tibble()
dbClearResult(res)

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

# calculate coxp
expr %>%
    dplyr::select(snorna, sample_id, expression = snorna_expression) %>%
    left_join(clinical, by="sample_id") %>%
    filter(!is.na(time), time > 0, !is.na(status)) %>%
    group_by(snorna) %>%
    mutate(expression = log2(expression + 1)) %>%
    mutate(group = as.factor(ifelse(expression <= median(expression),"Low", "High"))) %>%
    mutate(months = ifelse(time / 30 > 60, 60, time / 30)) %>%
    mutate(status = ifelse(months >= 60 & status == 1, 0, status)) %>%
    ungroup() ->
    for_survival


if(nrow(for_survival) < 20 || nlevels(for_survival$group) != 2 || mean(for_survival$expression) < 1){
    write_json(NULL, path=file.path(resource_jsons, paste(arg_encode,'survival','json',sep='.')))
    print("Not enough samples OR can't cut group OR mean expression less than 1!")
    quit("no", status = 0)
}



survival_coxph_model <-
    for_survival %>%
    group_by(snorna) %>%
    do(tidy(coxph(Surv(months, status) ~ expression, data = ., na.action = na.exclude))) %>%
    ungroup() %>%
    mutate(dataset_id = dataset_id,
            encode = arg_encode)
write_json(survival_coxph_model, path=file.path(resource_jsons, paste(arg_encode,'survival','json',sep='.')))


for_survival_distinct <- unique(for_survival$snorna)
plot_survival <- function(x){
   fit_x <- survfit(Surv(months,status) ~ group, data = for_survival %>%
        filter(snorna == x), na.action = na.exclude)
    plot_x <- ggsurvplot(fit_x,
                data= for_survival %>% filter(snorna == x),
                pval=T,
                pval.method = T,
                title = paste("5-year Survival, Coxph =",
                        ifelse(survival_coxph_model$p.value < 1e-04, 0.0001,
                        signif(survival_coxph_model$p.value, 3))),
                xlab = "5-year survival (months)",
                ylab = 'Probability of survival')
    ggsave(filename = file.path(resource_pngs, paste(arg_encode,"survival","png",sep = ".")), device = "png")
}
for_survival_distinct %>% sapply(plot_survival)
