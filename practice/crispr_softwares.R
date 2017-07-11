library(tidyverse)
library(ggthemes)
library(stringr)
library(multidplyr)
# CRISPR chronological distribution
library(RISmed)
path <- "/home/cliu18/liucj/projects/5.CRISPR/softwares/omictools"
# res <- EUtilsSummary("crispr", type = 'esearch', db = 'pubmed')
# QueryCount(res)
# 
# tally <- array()
# x <- 1
# for(i in 1986:2017){
#   Sys.sleep(1)
#   r <- EUtilsSummary('crispr', type = 'esearch', db = 'pubmed', mindate = i, maxdate = i)
#   tally[x] <- QueryCount(r)
#   x <- x+ 1
# }

query_count <- function(year, query = "enhancer"){
  res <- EUtilsSummary(query, type = 'esearch', db = 'pubmed', mindate = year, maxdate = year)
  QueryCount(res)
}

pub <- tibble(year = 1956:2017)

# pub %>% head(2) %>%  mutate( num = map_dbl(.x = year, .f = query_count))

cl <- parallel::detectCores()
cluster <- create_cluster(floor(cl * 5 / 6))
pub %>% 
  # register cores
  partition(cluster = cluster) %>% 
  # load library to cluster
  cluster_library("RISmed") %>% 
  cluster_library("tidyverse") %>% 
  cluster_assign_value("query_count", query_count) ->
  pub_shards

pub_num <-
  pub_shards %>% 
  # get names and seq
  mutate(num = map_dbl(.x = year, .f = query_count)) %>% 
  collect() %>% 
  as_tibble() %>% 
  ungroup()
parallel::stopCluster(cluster) 


pub_num %>% 
  arrange(num) %>% filter(year > 1977 ) %>%
  ggplot(mapping = aes(x = year, y = num)) +
  geom_point() +
  geom_line() +
  theme_gdocs() +
  labs(x = "Year", y = "Number")  -> p
ggsave(plot = p, filename = "02_crispr_publication_journal.pdf", device = "pdf", path = path)

read_csv(file = file.path(path, "crispr_tool_list.csv")) %>% 
  select(time) %>% 
  group_by(time) %>% 
  count() -> cal

ggplot(data = cal, mapping = aes(x = time, y = n))+
  geom_point(color = "red") + 
  geom_line(color = "red") + 
  theme_few() +
  theme(axis.title = element_blank()) -> p_cal
ggsave(plot = p_cal, filename = "01_cal_publication_journal.pdf", device = "pdf", path = path)

# tally %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   as_tibble() %>% 
#   rename(., number = `.`) %>% mutate(rownames = as.numeric(rowname)) %>% 
#   ggplot(mapping = aes(x = rowname, y = number)) +
#   geom_text(aes(label = number), size = 3, vjust = -0.3) + 
#   geom_bar(stat = 'identity') + 
#   theme_gdocs() +
#   theme(axis.text.x = element_text(angle = 45)) +
#   labs(x = "Year", y = "Number")
# 
# ggsave(filename = "CRISPR_publications.png", device = "png", path = path)

# get tools list from omictools
library(rvest)
url <- "https://omictools.com"
url_crispr <- str_c(str_c(url, "/search?q=crispr"), c('', str_c("&page", 2:5, sep = "=")))

get_description <- function(x){
  x %>% html_node('article p.tool-section-grid-item-content-description-long') %>% 
    html_attr(name = "title")
}

get_tool_name <- function(x){
  x %>% html_node('header a.js-tool-card-link') %>% html_attrs() %>% unlist()
}

crispr_list <- function(x){
  desc <- get_description(x)
  tool_name <- get_tool_name(x)
  c(tool_name['title'], tool_name['href'], 'desc' = desc) %>% enframe() %>% spread(name, value) 
}

get_basic_crispr_info <- function(x){
  x %>% 
    read_html() %>% 
    html_nodes('div.tool-section-grid-item') %>% 
    map(crispr_list) %>% 
    bind_rows()
}

url_crispr %>% 
  map(get_basic_crispr_info) %>% 
  bind_rows() -> crispr_tools_info

get_doc <- function(x){
  # href
  x %>% html_node("a.tool-heading-title-link") %>% html_attr(name = "href") -> href
  # interfeace
  x %>% html_nodes("div.section-row-column dt") %>% html_text() %>% 
    str_replace_all(pattern = "\\r\\n|  |\\:", replacement = "") -> dt
  x %>% html_nodes("div.section-row-column dd") %>% html_text() %>% 
    str_replace_all(pattern = "\\r\\n|  ", replacement = "")-> dd
  interface <- tibble(dt = dt, dd = dd)
  
  x %>% html_nodes("li.tool-publication-list-item.is-primary") %>% 
    html_text() %>% 
    str_split(pattern = "\\n", simplify = T) %>% 
    str_replace_all("\\  ", "") -> pub
  pub <- pub[pub != ""]
  journal <- pub[3]
  pmid <- str_split(pub[4], pattern = "\\:", simplify = T)[2] %>% str_trim()
  author <- pub[1]
  time <- str_extract(string = pub[1], pattern = "\\d+") %>% as.numeric()
  
  publi <- tibble(time = time, journal = journal, pmid = pmid)
  
  x %>% html_nodes("li.literature-link-list-item") %>% 
    html_text() %>% 
    str_replace_all(pattern = "\\r\\n|  ", replacement = "")  -> link_lit
  
  tibble(url = href, time = time, journal = journal, pmid = pmid, interface = list(interface), link_lit = list(link_lit))
}

crispr_detail <- function(x){
  #href, Interface, stability, publication, link literature
  str_c(url, x) %>% read_html() %>% get_doc()
}

# crispr_tools_info %>% 
#   mutate(details = map(href, crispr_detail)) %>% 
#   unnest() -> crispr_tools_all

cl <- parallel::detectCores()
cluster <- create_cluster(floor(cl * 5 / 6))
crispr_tools_info %>% 
  # register cores
  partition(cluster = cluster) %>% 
  # load library to cluster
  cluster_library("tidyverse") %>% 
  cluster_library("stringr") %>% 
  cluster_library("rvest") %>% 
  cluster_assign_value("crispr_detail", crispr_detail) %>% 
  cluster_assign_value("get_doc", get_doc) %>% 
  cluster_assign_value("url", url)->
  crispr_tools_info_shards

crispr_tools_all <-
  crispr_tools_info_shards %>% 
  # get names and seq
  mutate(details = map(href, crispr_detail)) %>% 
  collect() %>% 
  as_tibble() %>% 
  unnest() %>% 
  ungroup()
parallel::stopCluster(cluster) 

crispr_tools_all %>% 
  mutate(journal = ifelse(str_detect(journal, "DOI"), "BMC Bioinformatics.", journal)) ->
  crispr_tools_all

crispr_tools_all %>% 
  ggplot(aes(x = as.factor(time))) +
  geom_bar() +
  geom_text(stat='count', aes(label = ..count..), size = 3, vjust = -0.3) +
  theme_gdocs() +
  labs(x = "Year", y = "Number")
ggsave(filename = "publication_time.png", device = "png", path = path)

crispr_tools_all %>% 
  ggplot(aes(x = reorder(journal, journal, function(x) -length(x)))) + 
  geom_bar(aes(fill = as.factor(time))) +
  geom_text(stat='count', aes(label = ..count..), size = 3, vjust = -0.3) + 
  theme_gdocs() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Journal", y = "Number") +
  scale_fill_gdocs(name = "Year")
ggsave(filename = "publication_journal.png", device = "png", path = path)

crispr_tools_all  %>% filter(time >= 2013) %>% select(title, interface, link_lit) %>%  filter(link_lit %>% map_dbl(length) > 0) %>% print(n = Inf)

crispr_tools_all  %>% filter(time >= 2013) %>% mutate(reviewed = link_lit %>% map_dbl(length)) %>% select(-c(PARTITION_ID, interface, link_lit)) -> xls 

write_excel_csv(x = xls, path = file.path(path, "crispr_tool_list.csv"))
 