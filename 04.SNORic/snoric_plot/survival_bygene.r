#!/usr/bin/R
library(RMySQL)
library(stringr)
library(jsonlite)
library(tidyverse)
library(broom)
library(latex2exp)
library(survival)
library(survminer)
#snorna_gene_snorna_pair

args <- commandArgs(TRUE)
# args <- 'UFRFTiNUQ0dBLUFDQw=='
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
                    str_c(arg_encode, 'survival_bygene', 'json', sep="."))
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
    filter(!is.na(time), time > 0, !is.na(status)) %>%
    group_by(snorna) %>%
    mutate(expression = log2(expression + 1)) %>%
    mutate(group = ifelse(expression <= median(expression),"Low", "High")) %>%
    mutate(months = ifelse(time / 30 > 60, 60, time / 30)) %>%
    mutate(status = ifelse(months >= 60 & status == 1, 0, status)) %>%
    filter(n() >=20)->
    for_survival


if(nrow(for_survival) < 1){
    write_json(
        NULL,
        path = file.path(
                    resource_jsons,
                    str_c(arg_encode, 'survival_bygene', 'json', sep="."))
                    )
    quit("no", status = 0)
}

survival_coxph_model <-
    for_survival %>%
    group_by(snorna) %>%
    do(tidy(coxph(Surv(months, status) ~ expression, data = ., na.action = na.exclude))) %>%
    ungroup() %>%
    mutate(dataset_id = dataset_id,
            encode = arg_encode,
            p.value = ifelse(p.value < 1e-04, 0.0001, signif(p.value, 3)))
write_json(survival_coxph_model, path=file.path(resource_jsons, paste(arg_encode,'survival_bygene','json',sep='.')))

for_survival_distinct <- unique(for_survival$snorna)
plot_survival <- function(x){
   fit_x <- survfit(Surv(months,status) ~ group,
                    data = filter(for_survival, snorna == x),
                    na.action = na.exclude)
   coxp <- coxph(Surv(months, status) ~ expression,
                    data = filter(for_survival, snorna == x),
                    na.action = na.exclude) %>%
                tidy() %>%
                .$p.value
    plot_x <- ggsurvplot(fit_x,
                data= for_survival %>% filter(snorna == x),
                pval=T,
                pval.method = T,
                title = paste("5-year Survival, Coxph =",
                        ifelse(coxp < 1e-04, 0.0001,
                        signif(coxp, 3))),
                xlab = "5-year survival (months)",
                ylab = 'Probability of survival')
    ggsave(filename = file.path(resource_pngs, paste(arg_encode,x, "survival_bygene","png",sep = ".")), device = "png")
}
for_survival_distinct %>% sapply(plot_survival)


#plot_survival <- function(df){
#    fit <- survfit(Surv(months,status) ~ group, data = df, na.action = na.exclude)
#    ggsurvplot(fit,
#        data= df,
#        pval=T,
#        pval.method = T,
#        title = paste("5-year Survival, Coxph =",
#                ifelse(survival_coxph_model$p.value < 1e-04, 0.0001,
#                signif(survival_coxph_model$p.value, 3))),
#        xlab = "5-year survival (months)",
#        ylab = 'Probability of survival') ->
#    plot
#    res <- arrange_ggsurvplots(list(plot), ncol = 1, nrow = 1, print = FALSE)
#    return(res)
#}

#for_survival %>%
#    group_by(snorna) %>%
#    by_slice(plot_survival, .labels = TRUE, .to="plot") %>%
#    transmute(
#                filename = str_c(
#                arg_encode,
#                snorna,
#                "survival_bygene",
#                "png",
#                sep = "."
#                ),
#                plot = plot
#            ) ->
#    survplots
#survplots %>% pwalk(ggsave, path = resource_pngs, device = "png", width = 7, height = 7)
dbDisconnect(con)
